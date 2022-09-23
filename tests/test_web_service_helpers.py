"""
Tests functionalities related to the CoGEx web service serving INDRA Discovery
"""
import json
from typing import Iterable

from indra.statements import Evidence, Agent, Activation
from indra_cogex.apps.utils import unicode_escape, _stmt_to_row
from indra_cogex.apps.queries_web.helpers import get_docstring


def test_unicode_double_escape():
    """Test unicode_double_escape function"""
    true_beta = "β"
    single_escaped_beta = r"\u03b2"
    double_escaped_beta = r"\\u03b2"

    true_alpha = "α"
    quadruple_escaped = r"\\\\u03b1"

    unequal_escaped = r"\\\\u03b1 and \\u03b2"
    true_alpha_and_beta = r"α and β"

    # Test with unicode
    assert unicode_escape(single_escaped_beta) == true_beta
    assert unicode_escape(double_escaped_beta) == true_beta
    assert unicode_escape(quadruple_escaped) == true_alpha
    assert unicode_escape(unequal_escaped) == true_alpha_and_beta

    # Test with non-unicode
    assert unicode_escape("a") == "a"
    assert unicode_escape("no unicode in here") == "no unicode in here"


def test__stmt_to_row():
    db_ev = Evidence._from_json(
        {
            "source_api": "biopax",
            "pmid": "12917261",
            "source_id": "http://pathwaycommons.org/pc12/Catalysis_8049495032c7bba740de082d7bf6c3da",
            "annotations": {"source_sub_id": "pid"},
            "epistemics": {"direct": True},
            "text_refs": {"PMID": "12917261"},
            "source_hash": 7478359958559662154,
        }
    )
    a = Agent("a")
    b = Agent("b")
    db_stmt = Activation(a, b, evidence=[db_ev])
    source_counts = {"biopax": 1}
    ev_array, english, stmt_hash, sources, total_evidence, badges, = _stmt_to_row(
        stmt=db_stmt,
        cur_dict={},
        cur_counts={},
        source_counts=source_counts,
        include_belief_badge=True,
    )
    assert int(total_evidence) == 1
    assert "biopax" in ev_array
    assert "7478359958559662154" in ev_array
    assert 'null' in ev_array
    assert sources == json.dumps(source_counts)
    assert english == '"<b>A</b> activates <b>b</b>."'


def test_rest_api_type_annotation():
    def my_func(mystr: str, my_iter: Iterable[str]) -> str:
        """Top level description

        Parameters
        ----------
        mystr :
            A string.
        my_iter :
            An Iterable over strings.

        Returns
        -------
        :
            A longer string.
        """
        return mystr + ",".join(my_iter)

    # The extra space between parameter descriptions is
    parsed_docstr = """Top level description

Parameters
----------
mystr : str
    A string.

my_iter : List[str]
    An Iterable over strings.

Returns
-------
str
    A longer string.
"""
    short, full_docstr = get_docstring(my_func)
    assert short == parsed_docstr.split("\n")[0], short
    assert full_docstr == parsed_docstr, full_docstr
