from flask import Blueprint, render_template, request
from flask_jwt_extended import jwt_required
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField
from wtforms.validators import DataRequired

from indra_cogex.apps.utils import render_statements, format_stmts
from indra_cogex.client.queries import *
__all__ = ["search_blueprint"]

search_blueprint = Blueprint("search", __name__, url_prefix="/search")

class SearchForm(FlaskForm):
    agent_name = StringField("Agent Name", validators=[DataRequired()])
    agent_role = StringField("Agent Role")
    source_type = StringField("Source Type")
    rel_type = StringField("Relationship Type")
    submit = SubmitField("Search")

@search_blueprint.route("/", methods=['GET','POST'])
@jwt_required(optional=True)
def search():
    form = SearchForm()
    if form.validate_on_submit():
        agent_name = form.agent_name.data
        agent_role = form.agent_role.data
        source_type = form.source_type.data
        rel_type = form.rel_type.data 
        statements, evidence_count = get_statements(
            agent=agent_name,
            agent_role=agent_role,
            stmt_sources=source_type,
            rel_types=rel_type,
            limit=50,
            evidence_limit=2000,
            return_evidence_counts=True
        )
        return render_statements(stmts=statements, evidence_count=evidence_count)

    return render_template("search/search_page.html", form=form)
