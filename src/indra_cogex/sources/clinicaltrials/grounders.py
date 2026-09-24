from functools import lru_cache
from typing import Optional, Literal

import gilda

from indra.databases import drugbank_client as db
from gilda import make_grounder, Term, Annotation
from gilda.process import normalize
from gilda.ner import annotate as gilda_annotate
from indra.ontology.bio import bio_ontology
from trialsynth.base.ground import InterventionGrounder, Annotator
from trialsynth.base.models import Intervention, BioEntity

ANNOTATE_MIN_LEN = 4
ANNOTATE_MIN_SCORE = 0.7
ANNOTATE_STOPLIST = {
    "NT", "CON", "GCA", "TAB", "CPAP", "COPE", "CHCPE", "TPF", "PF",
    "JIA", "OR", "IV", "WT", "HR", "CI", "RR", "OS", "PFS", "CR", "PR",
    "SD", "PD", "CT", "MRI", "PCR", "IHC", "AE", "SAE", "PS", "ECOG",
}
MESH_PREFIX = "MESH"
DRUG_NAMESPACES = ["DRUGBANK", "CHEBI", "MESH"]  # Todo: expand?
SKIP_INTERVENTIONS = set()  # Todo: fill out after first try

_drugbank_grounder = None


def build_drugbank_terms():
    """Parse INDRA's DrugBank ID to name mappings into Gilda Terms."""
    terms = []
    for drugbank_id, name in db.drugbank_names.items():
        terms.append(Term(
            norm_text=normalize(name),
            text=name,
            db="DRUGBANK",
            id=drugbank_id,
            entry_name=name,
            status="name",
            source="drugbank",
        ))
    return terms


def get_drugbank_grounder():
    """Build once and return a DrugBank-only grounder, without touching Gilda's global grounder."""
    global _drugbank_grounder
    if _drugbank_grounder is None:
        _drugbank_grounder = make_grounder(build_drugbank_terms())
    return _drugbank_grounder


def _term_result(scored_match):
    """Return db, id, name, and score for a scored match."""
    top = scored_match.term
    return {"db": top.db, "id": top.id, "entry_name": top.entry_name, "score": scored_match.score}


def _first_valid_annotation(text, grounder, namespaces):
    """Fallback: annotate text with the given grounder, return the first hit above the length, stoplist, and score filters."""
    annotations = gilda_annotate(text, grounder=grounder, namespaces=namespaces)
    for annotation in annotations:
        matched_text = annotation.text.strip()
        if len(matched_text) < ANNOTATE_MIN_LEN:
            continue
        if matched_text.upper() in ANNOTATE_STOPLIST:
            continue
        if annotation.matches[0].score < ANNOTATE_MIN_SCORE:
            continue
        return annotation.matches[0]
    return None


def drugbank_ground(text):
    """Ground text with DrugBank first, falling back to the default Gilda grounder."""
    drugbank_results = get_drugbank_grounder().ground(text)
    if drugbank_results:
        return _term_result(drugbank_results[0])
    # namespaces=None here: get_drugbank_grounder() only ever holds DRUGBANK terms, so
    # there is nothing else to restrict against.
    drugbank_match = _first_valid_annotation(text, get_drugbank_grounder(), namespaces=None)
    if drugbank_match:
        return _term_result(drugbank_match)
    fallback_results = gilda.get_grounder().ground(text, namespaces=DRUG_NAMESPACES)
    if fallback_results:
        return _term_result(fallback_results[0])
    fallback_match = _first_valid_annotation(text, gilda.get_grounder(), namespaces=DRUG_NAMESPACES)
    if fallback_match:
        return _term_result(fallback_match)
    return None


@lru_cache(1)
def get_drug_grounder():
    terms = build_drugbank_terms()
    grounder = make_grounder(terms)
    return grounder


class ClinicalTrialsDrugAnnotator(Annotator):
    """Annotator for drug interventions in clinical trials."""

    def __init__(self, namespaces):
        super().__init__(namespaces=namespaces, mesh_prefix=MESH_PREFIX)

    def annotate(self, text: str, *, context: str = None) -> list[Annotation]:
        drugbank_grounder = get_drug_grounder()
        annotations = gilda_annotate(
            text,
            grounder=drugbank_grounder,
            namespaces=self.namespaces,
            context_text=context,
        )

        # Filter out annotations
        filtered_annotations = []
        for annotation in annotations:
            matched_text = annotation.text.strip()
            if len(matched_text) < ANNOTATE_MIN_LEN:
                continue
            if matched_text.upper() in ANNOTATE_STOPLIST:
                continue
            if annotation.matches[0].score < ANNOTATE_MIN_SCORE:
                continue
            filtered_annotations.append(annotation)

        return filtered_annotations


class ClinicalTrialsDrugGrounder(InterventionGrounder):
    """Grounder for drug interventions in clinical trials."""

    def __init__(self):
        self.curations = {}
        self.namespaces = DRUG_NAMESPACES
        annotator = ClinicalTrialsDrugAnnotator(namespaces=self.namespaces)
        self.drugbank_grounder = get_drug_grounder()
        super().__init__(
            namespaces=self.namespaces,
            annotator=annotator,
            grounder_func=self.drug_grounder,
            mesh_prefix=MESH_PREFIX
        )

    def drug_grounder(
        self, text: str, *, namespaces: Optional[list[str]], context: str = None
    ) -> list[dict]:
        """Ground the text to a drug term."""
        if namespaces is None:
            namespaces = self.namespaces
        matches = self.drugbank_grounder.ground(text, context=context, namespaces=namespaces)

        # Filter matches
        if matches and matches[0].term.get_curie() in SKIP_INTERVENTIONS:
            # If it is, return an empty list
            return []

        # Trialsynth does filtering, so return all matches here
        return matches

    def preprocess(self, intervention: Intervention) -> Intervention:
        # Remove grounding for interventions part of SKIP_INTERVENTIONS
        if (
            intervention.ns
            and intervention.ns_id
            and intervention.curie in SKIP_INTERVENTIONS
        ):
            intervention.ns = None
            intervention.ns_id = None
            intervention.grounded_term = None
            return intervention

        if intervention.text:
            # Clean the text of the entity from ®, ™ and ©, which is known to cause
            # issues when grounding brand names
            clean_text = _remove_symbols(intervention.text)
            intervention.text = clean_text

            if clean_text in self.curations:
                curie = self.curations[clean_text]
                # CHEBI:CHEBI:1234 -> CHEBI, CHEBI:1234
                db_ns, db_id = curie.split(":", maxsplit=1)
                intervention.ns = db_ns
                intervention.ns_id = db_id
                intervention.grounded_term = bio_ontology.get_name(intervention.curie.split())
                return intervention

        # If not curated, just return the entity
        return intervention


def _remove_symbols(s: str) -> str:
    # Clean a string from special characters

    # Remove special characters
    special = ["®", "™", "©"]
    for char in special:
        s = s.replace(char, "")

    return s
