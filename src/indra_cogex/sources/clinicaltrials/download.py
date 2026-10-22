"""
Download and parse the ClinicalTrials.gov data using Trialsynth.
"""
import itertools
import json
import logging
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set, Tuple

import pystow
import pandas as pd
from tqdm import tqdm

from indra.databases import drugbank_client
from indra.ontology.bio import bio_ontology
from indra_cogex.client import process_identifier
from indra_cogex.representation import dump_norm_id
from indra_cogex.sources.clinicaltrials.grounders import ClinicalTrialsDrugGrounder
from trialsynth.ctgov import config, process
from trialsynth.base.extract.paths import (
    RESULTS_GROUNDED_CTGOV_DIR,
    RESULTS_GROUNDED_DIR,
)


def _clean(text: str) -> str:
    """Strip control characters from free-text fields before TSV serialization."""
    if not text:
        return text
    return text.replace("\x00", "").replace("\t", " ").replace("\n", " ").replace("\r", " ")

__all__ = [
    "ensure_clinical_trials_df",
    "process_trialsynth_edges",
    "process_trialsynth_trial_publication_edges",
    "process_trialsynth_bioentity_nodes",
    "process_trialsynth_trial_nodes",
    "load_all",
]

CLINICAL_TRIALS_MODULE = pystow.module(
    "indra",
    "cogex",
    "clinicaltrials",
)
TRIALSYNTH_PATH = CLINICAL_TRIALS_MODULE.module("trialsynth")
os.environ["DATA_DIR"] = str(TRIALSYNTH_PATH.base.absolute())

# Initializes the configuration for the clinical trials module in trialsynth
ctconfig = config.CTConfig()


#: Labels to keep for interventions in the bioentity nodes.
INTERVENTION_LABELS = {
    "DRUG",
    "DIETARY_SUPPLEMENT",
}


logger = logging.getLogger(__name__)

#: Corpus name -> (grounded JSON dir, record ID field, namespace of the node
#: owning the TrialResult)
RESULT_CORPORA = {
    "pubmed": (RESULTS_GROUNDED_DIR.base, "pmid", "PUBMED"),
    "ctgov": (RESULTS_GROUNDED_CTGOV_DIR.base, "nct_id", "CLINICALTRIALS"),
}


def ensure_clinical_trials_df(
    *,
    redownload: bool = False,
    reprocess: bool = False,
    max_pages: Optional[int] = None
):
    """Download and parse the ClinicalTrials.gov data using Trialsynth.

    Parameters
    ----------
    redownload :
        If True, redownload the raw data, even if it already exists. This will
        also set `reprocess` to True, since redownloading the data means that it
        needs to be reprocessed. Default: False.
    reprocess :
        If True, reprocess the data, even if it already exists. This will
        reprocess the raw data, but not redownload it if it already exists.
        Default: False.
    max_pages :
        The maximum number of pages to download from the ClinicalTrials.gov API.
        If None, all pages will be downloaded. Default: None.
    """
    if redownload:
        reprocess = True

    # Check the processed trialsynth data directory
    if not reprocess and all(
        p.exists() for p in (
            ctconfig.edges_path,
            ctconfig.bio_entities_path,
            ctconfig.trials_path,
            ctconfig.trial_publication_edges_path,
        )
    ):
        logger.info("ClinicalTrials.gov data already processed, skipping download.")
        return

    drugbank_grounder = ClinicalTrialsDrugGrounder()
    ctp = process.CTProcessor(
        reload_api_data=redownload,
        store_samples=True,
        validate=False,
        intervention_grounder=drugbank_grounder,
    )
    ctp.run(max_pages=max_pages)


def _mesh_to_chebi(row) -> str:
    """Convert a MeSH CURIE to a DrugBank or ChEBI CURIE if possible for interventions"""
    # Some interventions in the trialsynth data are directly from mesh
    # annotations and don't go through grounding, where drugbank and chebi
    # are prioritized over mesh
    curie = row["bioentity"]
    if curie is None:
        raise ValueError(
            "The row does not have a 'bioentity' column, cannot convert to CHEBI."
        )
    if row["rel_type:string"] == "condition" or not curie.lower().startswith("mesh:"):
        # If it's not mesh or it's not an intervention row just return the
        # original CURIE
        return curie

    # At this point we know that the CURIE is a MeSH CURIE and we want to
    # upgrade it to a DrugBank CURIE if possible, falling back to ChEBI.
    # The bio_ontology has chebi nodes stored as ("CHEBI", "CHEBI:12345")
    mesh_id = curie.split(":")[1]
    chebi_ns, chebi_id = bio_ontology.map_to(
        ns1="MESH", id1=mesh_id, ns2="CHEBI"
    ) or (None, None)

    # There is no direct MeSH-to-DrugBank cross-reference in the bio_ontology,
    # so this has to go through CHEBI as an intermediate step: MESH -> CHEBI
    # -> DRUGBANK, using the same DrugBank-to-CHEBI mapping data the drug
    # grounder itself is built from.
    if chebi_id:
        drugbank_id = drugbank_client.get_drugbank_id_from_chebi_id(chebi_id)
        if drugbank_id:
            return f"DRUGBANK:{drugbank_id}"

    return chebi_id if chebi_id else curie


def _normalize_bioentity_curie(curie: str) -> str:
    if pd.isna(curie):
        return curie
    db_ns, db_id = process_identifier(curie)
    return dump_norm_id(db_ns, db_id)


def _merge_grounding_sources(sources: pd.Series) -> str:
    merged: set[str] = set()
    for value in sources:
        if pd.isna(value) or value == "":
            continue
        merged.update(
            item.strip() for item in str(value).split(";") if item.strip()
        )
    return ";".join(sorted(merged))


def process_trialsynth_edges() -> pd.DataFrame:
    """Convert the edge file from trialsynth to CoGEx format

    Returns
    -------
    :
        The converted edges DataFrame with CoGEx format headers and intervention
         values converted to Chebi identifiers.
    """
    headers_translation = {
        "from:CURIE": "trial",
        "to:CURIE": "bioentity",
    }

    # Read the edges file from trialsynth
    edges_df = pd.read_csv(ctconfig.edges_path, sep="\t", compression="gzip")

    # Rename the columns to match CoGEx format
    edges_df.rename(columns=headers_translation, inplace=True)

    # Only translate has_intervention edges from mesh to chebi, but use a new column
    # for the resulting values and fill in with the untranslated values for the rows
    # that were not translated
    edges_df["bioentity_mapped"] = edges_df.apply(_mesh_to_chebi, axis=1)

    # Drop the "source_registry:string" column
    if "source_registry:string" in edges_df.columns:
        edges_df.drop(columns=["source_registry:string"], inplace=True)

    # Normalize bioentity CURIEs so equivalent identifiers group together
    edges_df["bioentity_mapped"] = edges_df["bioentity_mapped"].map(
        _normalize_bioentity_curie
    )

    # Merge duplicate rows on trial, mapped bioentity, and relation type while
    # unioning semicolon-separated grounding source lists
    merge_keys = ["trial", "bioentity_mapped", "rel_type:string"]
    agg_columns = {
        col: "first"
        for col in edges_df.columns
        if col not in merge_keys + ["grounding_sources:string[]"]
    }
    agg_columns["grounding_sources:string[]"] = _merge_grounding_sources
    edges_df = edges_df.groupby(merge_keys, as_index=False).agg(agg_columns)

    return edges_df


def process_trialsynth_trial_publication_edges() -> pd.DataFrame:
    """Load the publication-trial edge file from trialsynth

    Returns
    -------
    :
        A dataframe with CoGEx formatted edges
    """
    return pd.read_csv(
        ctconfig.trial_publication_edges_path, sep="\t", compression="gzip"
    )


def process_trialsynth_bioentity_nodes(mesh_chebi_map: Dict[str, str]) -> pd.DataFrame:
    """Convert the bioentity nodes file from the trialsynth to CoGEx format

    Returns
    -------
    :
        writeme
    """
    headers_translation = {
        "curie:CURIE": "bioentity",
        "term:string": "name",
    }

    # Read the bioentity nodes file from trialsynth
    bioentity_nodes_df = pd.read_csv(
        ctconfig.bio_entities_path, sep="\t", compression="gzip"
    )

    # Rename the columns to match CoGEx format
    bioentity_nodes_df.rename(columns=headers_translation, inplace=True)

    def _map_to_chebi(row) -> str:
        if "intervention" in row["labels:LABEL[]"].lower():
            return mesh_chebi_map.get(
                row["bioentity"], row["bioentity"]
            )
        return row["bioentity"]

    # Translate the same rows that were translated in the edges file:
    # Any mesh id that is an intervention should be attempted to be converted to chebi
    bioentity_nodes_df["bioentity_mapped"] = bioentity_nodes_df.apply(_map_to_chebi, axis=1)

    # Map names for the chebi mapped bioentities
    def _nsid_to_name(row) -> str:
        if pd.isna(row["bioentity_mapped"]) or row["bioentity_mapped"] is None:
            return row["name"]
        if "intervention" in row["labels:LABEL[]"].lower() and \
            row["bioentity_mapped"].lower().startswith("chebi:"):
            # If it's an intervention, get the name from the chebi id, use the
            # existing name as default
            mapped_name = bio_ontology.get_name(
                ns="CHEBI", id=row["bioentity_mapped"].upper()
            )
            return mapped_name or row["name"]
        return row["name"]

    bioentity_nodes_df["name"] = bioentity_nodes_df.apply(_nsid_to_name, axis=1)

    # Filter to interventions that are labeled as drugs or dietary supplements
    def _intervention_filter(row) -> bool:
        """Check if the row is an intervention that should be kept."""
        labels = row["labels:LABEL[]"].lower()
        if "condition" in labels:
            # We keep conditions as they are not interventions
            return True
        # Check if the labels contain any of the allowed intervention labels
        if any(label.lower() in labels for label in INTERVENTION_LABELS):
            return True
        # Check if the curie is from chebi
        if row["bioentity"].lower().startswith("chebi:"):
            return True
        return False

    # Filter the interventions to the ones that are labeled as drugs or dietary
    # supplements or where the curie is from chebi
    bioentity_nodes_df = bioentity_nodes_df[
        bioentity_nodes_df.apply(_intervention_filter, axis=1)
    ]

    return bioentity_nodes_df[["bioentity_mapped", "name"]]


def process_trialsynth_trial_nodes() -> pd.DataFrame:
    """Convert the trial nodes file from the trialsynth to CoGEx format

    Returns
    -------
    :
        writeme
    """
    # Cogex headers:
    # id:ID
    # :LABEL
    # phase:int  <-- -1 for unknown phase
    # randomized:boolean
    # start_year:int
    # start_year_anticipated:boolean
    # completion_year:int
    # completion_year_anticipated:boolean
    # last_update_submit_year:int
    # status  # Completed, terminated etc...
    # study_type  # Observational, interventional etc...
    # why_stopped

    # Trialsynth has the following headers:
    # curie:CURIE
    # title:string <- Corresponds to BriefTitle
    # official_title:string <- Corresponds to OfficialTitle
    # brief_summary:string,
    # detailed_description:string,
    # labels:LABEL[] - scsv of clinical_trials, intervention, observational, expanded_access
    # design:DESIGN - e.g. 'Purpose: PREVENTION; Allocation: RANDOMIZED;Masking: DOUBLE; Assignment: '
    # conditions:CURIE[] - scsv of mesh conditions, e.g. 'mesh:D000001;mesh:D000002'
    # interventions:CURIE[]
    # primary_outcome:OUTCOME[]
    # secondary_outcome:OUTCOME[]
    # secondary_ids:CURIE[]
    # source_registry:string
    # phases:PHASE[]
    # start_year:integer
    # start_year_anticipated:boolean
    # primary_completion_year:integer
    # primary_completion_year_type:string
    # completion_year:integer
    # completion_year_type:string
    # last_update_submit_year:integer
    # status:string
    # why_stopped:string
    # references:string[]

    # Translate the headers to CoGEx format, these map 1:1
    headers_translation = {
        "curie:CURIE": "id:ID",
        "why_stopped:string": "why_stopped",
        "status:string": "status",
    }

    def _get_phase(phase_string: str) -> int:
        if pd.notna(phase_string):
            # The phase string is a list of phases, e.g. "PHASE1|PHASE2"
            # Get the highest phase number
            max_phase = max(
                int(p[-1]) if p[-1].isdigit() else -1
                for p in phase_string.split("|")
            )
            return max_phase
        return -1

    # Read the trial nodes file from trialsynth
    trials_nodes_df = pd.read_csv(ctconfig.trials_path, sep="\t", compression="gzip")

    # Rename the columns to match CoGEx format
    trials_nodes_df.rename(columns=headers_translation, inplace=True)

    # Add the :LABEL column
    trials_nodes_df[":LABEL"] = "ClinicalTrial"

    # Add the phase column, defaulting to -1 (unknown), pick the max phase
    trials_nodes_df["phase:int"] = trials_nodes_df["phases:PHASE[]"].apply(_get_phase)

    # Set the randomized boolean column based on the labels containing 'Allocation: RANDOMIZED'
    trials_nodes_df["randomized:boolean"] = trials_nodes_df["labels:LABEL[]"].str.contains(
        "Allocation: RANDOMIZED")

    # Copy the year columns as nullable ints
    for col in ("start_year", "completion_year", "last_update_submit_year"):
        trials_nodes_df[col] = trials_nodes_df[f"{col}:integer"].astype("Int64")

    # fixme: need a better way to get the study type out
    def _get_study_type(labels: str) -> str:
        if pd.isna(labels):
            return "unknown"
        # The study type is in the labels, e.g. "Study Type: Interventional"
        for label in labels.split(";"):
            for study_type in ["interventional", "observational", "expanded access"]:
                if study_type in label.lower():
                    return study_type
        return "unknown"
    # Add study type column
    trials_nodes_df["study_type"] = trials_nodes_df["labels:LABEL[]"].apply(_get_study_type)

    # Clean up
    # Remove newlines and 2+ whitespaces from the why_stopped column
    trials_nodes_df["why_stopped"] = trials_nodes_df["why_stopped"].str.replace(
        r"\s+", " ", regex=True
    ).str.strip()

    return trials_nodes_df


#: Node label -> name of its ID counter. The counters are shared by all
#: corpora so node IDs stay unique across them.
ID_COUNTERS = {
    "TrialResult": "result_id",
    "TrialArm": "arm_id",
    "TrialAsRunArm": "as_run_arm_id",
    "TrialAdverseEvent": "ae_id",
    "TrialCriterion": "criterion_id",
    "TrialMetric": "metric_id",
    "TrialOutcome": "outcome_id",
    "TrialStatisticalComparison": "stat_id",
}

#: Keys of the DataFrames returned by load_all
ROW_KEYS = (
    "result_nodes",
    "arms",
    "as_run_arms",
    "metrics",
    "adverse_events",
    "criteria",
    "outcomes",
    "stat_comparisons",
    "genetic_edges",
    "ae_bioentity_edges",
    "criterion_bioentity_edges",
    "source_edges",
)

#: Edge DataFrames that are deduplicated on these columns
DEDUPLICATE_ON = {
    "genetic_edges": ["result_id", "hgnc_id"],
    "ae_bioentity_edges": ["adverseevent_id", "db", "id"],
    "criterion_bioentity_edges": ["criterion_id", "db", "id"],
}


def _load_jsons(
    json_dir: Path,
    id_field: str,
    allowed_ids: Optional[Set[str]] = None,
) -> Iterator[Tuple[str, dict]]:
    """Yield (record_id, data) tuples from the grounded JSON files.

    Records are yielded one at a time so the caller never holds the whole
    corpus in memory. Files are read by a thread pool, which hides the
    per-file I/O latency, a chunk at a time so read-ahead stays bounded.

    Parameters
    ----------
    json_dir :
        Path to directory containing grounded JSON files.
    id_field :
        Key holding the record ID (e.g. pmid or nct_id) in each JSON.
    allowed_ids :
        If given, records whose ID is not in this set are skipped.

    Yields
    ------
    :
        A (record_id, data) tuple per JSON file.
    """
    paths = sorted(json_dir.glob("*.json"))
    n_loaded = 0
    n_skipped = 0
    chunk_size = 1000
    with ThreadPoolExecutor(max_workers=8) as executor, tqdm(
        total=len(paths), unit="file", desc=f"Loading {json_dir.name}"
    ) as progress:
        for start in range(0, len(paths), chunk_size):
            chunk = paths[start:start + chunk_size]
            progress.update(len(chunk))
            for path, raw in zip(chunk, executor.map(Path.read_bytes, chunk)):
                try:
                    data = json.loads(raw)
                except Exception as e:
                    logger.warning("Failed to load %s: %s", path.name, e)
                    continue
                record_id = str(data.get(id_field, path.stem))
                if allowed_ids is not None and record_id not in allowed_ids:
                    n_skipped += 1
                    continue
                n_loaded += 1
                yield record_id, data
    logger.info("Loaded %d JSON files from %s", n_loaded, json_dir)
    if n_skipped:
        logger.info("Skipped %d records with an ID not in allowed_ids", n_skipped)


def _get_trial_nct_ids() -> Set[str]:
    """Return the NCT IDs of the ClinicalTrial nodes built from trialsynth."""
    curies = pd.read_csv(
        ctconfig.trials_path, sep="\t", compression="gzip", usecols=["curie:CURIE"]
    )["curie:CURIE"]
    return {curie.split(":", 1)[1] for curie in curies.dropna()}


def _get_id_starts(look_up_current_max_ids: bool) -> Dict[str, int]:
    """Return the first ID to assign for each ID counter."""
    if not look_up_current_max_ids:
        return {name: 1 for name in ID_COUNTERS.values()}

    from indra_cogex.client import Neo4jClient
    client = Neo4jClient()
    starts = {}
    for label, name in ID_COUNTERS.items():
        rows = client.query_tx(
            f"MATCH (n:{label}) "
            "RETURN max(toInteger(split(n.id, ':')[-1])) AS max_id"
        )
        starts[name] = (int(rows[0][0] or 0) if rows else 0) + 1
    return starts


def _metric_row(metric: dict, parent_ns: str, parent_id: int, metric_id: int) -> dict:
    return {
        "parent_ns": parent_ns,
        "parent_id": parent_id,
        "metric_id": metric_id,
        "name": _clean(metric.get("name", "")),
        "value_numeric": metric.get("value_numeric"),
        "unit": _clean(metric.get("unit", "")),
        "value_text": _clean(metric.get("value_text", "")),
        "source_sentence": _clean(metric.get("source_sentence", "")),
    }


def _add_grounding_edge(edges: List[dict], id_key: str, node_id: int, item: dict):
    """Append an edge to the item's grounding, if it has one."""
    grounding = item.get("grounding") or {}
    if grounding.get("db") and grounding.get("id"):
        edges.append({id_key: node_id, "db": grounding["db"], "id": grounding["id"]})


def _extract_arms(data: dict, result_id: int, ids: dict, rows: dict):
    for arm in data.get("arms", []):
        arm_id = next(ids["arm_id"])
        rows["arms"].append({
            "result_id": result_id,
            "arm_id": arm_id,
            "arm_name": _clean(arm.get("arm_name", "")),
            "arm_type": _clean(arm.get("arm_type") or ""),
            "n": arm.get("n"),
            "dosage": _clean(arm.get("dosage") or ""),
            "source_sentence": _clean(arm.get("source_sentence", "")),
        })
        for metric in arm.get("metrics", []):
            rows["metrics"].append(
                _metric_row(metric, "arm", arm_id, next(ids["metric_id"]))
            )
        for ae in arm.get("adverse_events", []):
            ae_id = next(ids["ae_id"])
            rows["adverse_events"].append({
                "arm_id": arm_id,
                "adverseevent_id": ae_id,
                "event_name": _clean(ae.get("event_name", "")),
                "incidence_numeric": ae.get("incidence_numeric"),
                "unit": _clean(ae.get("unit", "")),
                "value_text": _clean(ae.get("value_text", "")),
                "source_sentence": _clean(ae.get("source_sentence", "")),
            })
            _add_grounding_edge(
                rows["ae_bioentity_edges"], "adverseevent_id", ae_id, ae
            )


def _extract_as_run_arms(data: dict, result_id: int, ids: dict, rows: dict):
    for arm in data.get("as_run_arms", []):
        rows["as_run_arms"].append({
            "result_id": result_id,
            "as_run_arm_id": next(ids["as_run_arm_id"]),
            "group_label": _clean(arm.get("group_label") or ""),
            "intervention": _clean(arm.get("drug") or ""),
            "dosage": _clean(arm.get("dose") or ""),
            "route": _clean(arm.get("route") or ""),
            "schedule": _clean(arm.get("schedule") or ""),
            "evidence_text": _clean(arm.get("evidence_text") or ""),
        })


def _extract_criteria(data: dict, result_id: int, ids: dict, rows: dict):
    for criterion_type in ("inclusion", "exclusion"):
        # Prefer the grounded lists, which carry the same text plus a
        # grounding, and fall back to the plain ones
        items = data.get(f"grounded_{criterion_type}_criteria") or data.get(
            f"{criterion_type}_criteria", []
        )
        for item in items:
            criterion_id = next(ids["criterion_id"])
            rows["criteria"].append({
                "result_id": result_id,
                "criterion_id": criterion_id,
                "text": _clean(item.get("text", "")),
                "criterion_type": criterion_type,
                "evidence_text": _clean(item.get("evidence_text", "")),
            })
            _add_grounding_edge(
                rows["criterion_bioentity_edges"], "criterion_id", criterion_id, item
            )


def _extract_outcomes(data: dict, result_id: int, ids: dict, rows: dict):
    for item in data.get("results", []):
        rows["outcomes"].append({
            "result_id": result_id,
            "outcome_id": next(ids["outcome_id"]),
            "text": _clean(item.get("text", "")),
            "evidence_text": _clean(item.get("evidence_text", "")),
        })


def _extract_stat_comparisons(data: dict, result_id: int, ids: dict, rows: dict):
    for comp in data.get("statistical_comparisons", []):
        stat_id = next(ids["stat_id"])
        rows["stat_comparisons"].append({
            "result_id": result_id,
            "statcomparison_id": stat_id,
            "comparison_name": _clean(comp.get("comparison_name", "")),
        })
        for metric in comp.get("metrics", []):
            rows["metrics"].append(
                _metric_row(metric, "statcomparison", stat_id, next(ids["metric_id"]))
            )


def _extract_genetic_edges(data: dict, result_id: int, rows: dict):
    for entry in data.get("genetic", {}).get("grounded_inclusion", []):
        for grounding in entry.get("groundings", []):
            info = grounding.get("info", {})
            if info.get("db") == "HGNC" and info.get("id"):
                rows["genetic_edges"].append({
                    "result_id": result_id,
                    "hgnc_id": info["id"],
                    "symbol": info.get("entry_name", grounding.get("symbol", "")),
                    "variant": entry.get("variant"),
                    "evidence_text": _clean(entry.get("evidence_text", "")),
                })


def _extract_record(
    data: dict, record_id: str, source_ns: str, ids: dict, rows: dict
):
    """Append the rows for one grounded JSON record."""
    result_id = next(ids["result_id"])
    randomization = data.get("randomization") or {}
    recruitment = data.get("recruitment") or {}
    rows["result_nodes"].append({
        "result_id": result_id,
        "study_info": _clean(data.get("study_info", "")),
        "trial_ids:string[]": ";".join(data.get("trial_ids", [])),
        "phase": data.get("phase") or "",
        "locations:string[]": ";".join(data.get("locations", [])),
        "randomization_ratio": _clean(randomization.get("ratio") or ""),
        "randomization_stratification_factors:string[]": ";".join(
            _clean(f).replace(";", ",")
            for f in randomization.get("stratification_factors") or []
        ),
        "recruitment_channel": _clean(recruitment.get("channel") or ""),
        "recruitment_window": _clean(recruitment.get("window") or ""),
    })
    rows["source_edges"].append({
        "source_ns": source_ns,
        "source_id": record_id,
        "result_id": result_id,
    })
    _extract_arms(data, result_id, ids, rows)
    _extract_criteria(data, result_id, ids, rows)
    _extract_as_run_arms(data, result_id, ids, rows)
    _extract_outcomes(data, result_id, ids, rows)
    _extract_stat_comparisons(data, result_id, ids, rows)
    _extract_genetic_edges(data, result_id, rows)


def load_all(
    corpora: Tuple[str, ...] = ("pubmed", "ctgov"),
    look_up_current_max_ids: bool = False,
) -> Dict[str, pd.DataFrame]:
    """Load all grounded JSONs in a single pass and return all DataFrames.

    All corpora share one set of ID counters so node IDs stay unique across
    them.

    Parameters
    ----------
    corpora :
        Names of the corpora in RESULT_CORPORA to load.
    look_up_current_max_ids :
        If True, look up the current max IDs in the database and use them to
        set the IDs for the nodes. This is useful for loading data into a
        database that already has some nodes. Default: False.

    Returns
    -------
    :
        Dictionary of DataFrames keyed by the names in ROW_KEYS.
    """
    ids = {
        name: itertools.count(start)
        for name, start in _get_id_starts(look_up_current_max_ids).items()
    }
    rows = {key: [] for key in ROW_KEYS}

    for corpus in corpora:
        json_dir, id_field, source_ns = RESULT_CORPORA[corpus]
        # Registry results attach directly to a ClinicalTrial node, so only
        # keep those whose trial exists in the trialsynth trial nodes
        allowed_ids = _get_trial_nct_ids() if source_ns == "CLINICALTRIALS" else None
        for record_id, data in _load_jsons(json_dir, id_field, allowed_ids):
            _extract_record(data, record_id, source_ns, ids, rows)

    frames = {key: pd.DataFrame(key_rows) for key, key_rows in rows.items()}
    for key, subset in DEDUPLICATE_ON.items():
        if not frames[key].empty:
            frames[key] = frames[key].drop_duplicates(subset=subset)
    for key, df in frames.items():
        logger.info("Extracted %d %s rows", len(df), key)
    return frames
