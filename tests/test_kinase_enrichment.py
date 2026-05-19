from indra_cogex.client.enrichment import discrete


class FakeClient:
    pass


def test_kinase_ora_uses_annotation_overlap_and_lazy_statements(monkeypatch):
    kinase_to_phosphosites = {
        ("hgnc:1", "KINASE1"): {("MAPK1", "T202"), ("AKT1", "S473")},
        ("hgnc:2", "KINASE2"): {("MAPK1", "Y187")},
    }

    def fail_count(*args, **kwargs):
        raise AssertionError("count_phosphosites should not be called")

    monkeypatch.setattr(discrete, "count_phosphosites", fail_count)
    monkeypatch.setattr(
        discrete,
        "get_kinase_phosphosites",
        lambda **kwargs: kinase_to_phosphosites,
    )

    df = discrete.kinase_ora(
        client=FakeClient(),
        phosphosite_ids=[("MAPK1", "T202"), ("NOPE", "S1")],
        alpha=0.05,
        keep_insignificant=True,
    )

    diagnostics = df.attrs["diagnostics"]
    assert diagnostics["enrichment_universe_source"] == "kinase_annotated_phosphosites"
    assert diagnostics["enrichment_universe_count"] == 3
    assert diagnostics["known_kinase_annotated_overlap"] == 1
    assert diagnostics["used_in_enrichment_universe_count"] == 1
    assert diagnostics["tested_kinases"] == 2
    assert df["statements"].tolist() == [[], []]
    assert dict(zip(df["curie"], df["matched_phosphosites"])) == {
        "hgnc:1": ["MAPK1-T202"],
        "hgnc:2": [],
    }
