from indra_cogex.client.enrichment.utils import extend_by_ontology


def test_extend_by_ontology():
    # Create a mapping with some of the hgnc IDs
    gene_set_mapping = {
        (
            "go:0000978",
            "RNA polymerase II cis-regulatory region sequence-specific DNA binding",
        ): {
            "123",
            "456",
            "789",
            "120",
            "654",
        },
        ("go:0035497", "cAMP response element binding"): {
            "711",
            "852",
            "963",
            "159",
            "357",
        },
    }

    # Extend the mapping by the GO ontology
    len_before = sum(map(lambda x: len(x), gene_set_mapping.values()))
    extend_by_ontology(gene_set_mapping)
    assert sum(map(lambda x: len(x), gene_set_mapping.values())) > len_before
