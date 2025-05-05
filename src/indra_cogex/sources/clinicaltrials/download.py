"""
Download and parse the ClinicalTrials.gov data using the v2 API.
See https://clinicaltrials.gov/data-api/api
"""

from typing import Optional, List

import pystow
import requests
from tqdm.auto import tqdm, trange
import pandas as pd
import io

__all__ = [
    "CLINICAL_TRIALS_PATH",
    "ensure_clinical_trials_df",
]

API_URL = "https://clinicaltrials.gov/api/v2/studies"
CLINICAL_TRIALS_MODULE = pystow.join(
    "indra",
    "cogex",
    "clinicaltrials"
)
CLINICAL_TRIALS_PATH = CLINICAL_TRIALS_MODULE.joinpath(
    "clinical_trials.tsv.gz"
)

#: The fields that are used by default. A full list can be found
#: here: https://classic.clinicaltrials.gov/api/info/study_fields_list
DEFAULT_FIELDS = [
    "NCTId",
    "BriefTitle",
    "Condition",
    "ConditionMeshTerm",
    "ConditionMeshId",
    "InterventionName",
    "InterventionType",
    "InterventionMeshTerm",
    "InterventionMeshId",
    "StudyType",
    "DesignAllocation",
    "OverallStatus",
    "Phase",
    "WhyStopped",
    "SecondaryIdType",
    "SecondaryId",
    "StartDate",  # Month [day], year: "November 1, 2023", "May 1984" or NaN
    "StartDateType",  # "Actual" or "Anticipated" (or NaN)
    "ReferencePMID",  # these are tagged as relevant by the author, but not necessarily about the trial
]

def ensure_clinical_trials_df(
    *,
    refresh: bool = False,
    page_size: int = 1000,
    max_pages: Optional[int] = None
) -> pd.DataFrame:
    """Download and parse data from ClinicalTrials.gov using the v2 API

    Parameters
    ----------
    refresh :
        If True, will re-download the data from ClinicalTrials.gov. If False,
        will use the cached data if available.
    page_size :
        The number of trials to fetch per page. Defaults to 1000, which is also
        the maximum allowed by the rest API.
    max_pages :
        The maximum number of pages to fetch. If None, will fetch all pages.

    Returns
    -------
    :
        A pandas DataFrame containing the clinical trials data.
    """
    if CLINICAL_TRIALS_PATH.exists() and not refresh:
        return pd.read_csv(CLINICAL_TRIALS_PATH, compression="gzip")

    params = {"format": "csv", "pageSize": page_size}

    # Get the total number of trials available from the API's json response
    resp_json = requests.get(API_URL, params={"pageSize": 1, "countTotal": "true"})
    resp_json.raise_for_status()
    resp_json = resp_json.json()
    total_trials = resp_json.get("totalCount", 0)

    pbar = tqdm(
        desc="Downloading ClinicalTrials.gov data",
        total=total_trials if max_pages is None else max_pages * page_size,
        unit="trial",
        unit_scale=True,
    )

    # First request: get header + first page of data
    resp = requests.get(API_URL, params=params)
    resp.raise_for_status()
    text = resp.text
    lines = text.splitlines()

    # Extract header from the very first line
    header = lines[0].split(",")

    # read first page normally (header=0)
    df_pages = [pd.read_csv(io.StringIO(text))]
    pbar.update(df_pages[0].shape[0])

    # Page through the rest
    token = resp.headers.get("X-Next-Page-Token")
    while token:
        params["pageToken"] = token
        resp = requests.get(API_URL, params=params)
        resp.raise_for_status()

        # read with our saved header, and tell pandas there is no header row in this chunk
        df = pd.read_csv(io.StringIO(resp.text), names=header, header=None)
        df_pages.append(df)

        # update the progress bar with the number of rows in this chunk
        pbar.update(df.shape[0])

        if max_pages is not None and len(df_pages) >= max_pages:
            tqdm.write(f"Reached max pages: {max_pages}. Stopping download.")
            break
        token = resp.headers.get("X-Next-Page-Token")

    pbar.close()

    # Concatenate
    df = pd.concat(df_pages, ignore_index=True)

    # Save the dataframe
    df.to_csv(CLINICAL_TRIALS_PATH, index=False, compression="gzip")

    return df


if __name__ == "__main__":
    ensure_clinical_trials_df(refresh=True)
