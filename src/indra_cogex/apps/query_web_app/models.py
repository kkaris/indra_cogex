from typing import Tuple, List, Any, Dict
from pydantic import BaseModel
from indra_cogex.representation import standardize, Node, norm_id


################
# Query Models #
################

class DefaultQuery(BaseModel):
    term: Tuple[str, str]

    def grounding(self) -> Tuple[str, str, str]:
        return standardize(*self.term)

    def normalize(self) -> str:
        _, ns, id_ = self.grounding()
        return norm_id(ns, id_)


class CheckConnection(BaseModel):
    source_term: Tuple[str, str]
    target_term: Tuple[str, str]

    def grounding(self) -> Tuple[Tuple[str, str, str], Tuple[str, str, str]]:
        return standardize(*self.source_term), standardize(*self.target_term)


class GoForGene(DefaultQuery):
    include_indirect: bool = False


class GeneForGo(DefaultQuery):
    include_indirect: bool = False


###################
# Response Models #
###################


# A json representation of indra_cogex.representation.Node
class NodeJson(BaseModel):
    labels: List[str]
    data: Dict[str, Any]


query_model_map = {
    "default": DefaultQuery,
    "get_go_terms_for_gene": GoForGene,
    "get_genes_for_go_term": GeneForGo,
    "is_query": CheckConnection,

}
