"""Reusable fields for the INDRA CoGEx analysis application."""

from wtforms import BooleanField, FileField, FloatField, IntegerField, RadioField
from wtforms.validators import DataRequired

indra_path_analysis_field = BooleanField("Include INDRA path-based analysis (slow)")
keep_insignificant_field = BooleanField(
    "Keep insignificant results (leads to long results lists)"
)
alpha_field = FloatField(
    "Alpha",
    default=0.05,
    validators=[DataRequired()],
    description="The alpha is the threshold for significance in the"
    " Fisher's exact test with which multiple hypothesis"
    " testing correction will be executed.",
)
correction_field = RadioField(
    "Multiple Hypothesis Test Correction",
    choices=[
        ("fdr_bh", "Family-wise Correction with Benjamini/Hochberg"),
        ("bonferroni", "Bonferroni (one-step correction)"),
        ("sidak", "Sidak (one-step correction)"),
        ("holm-sidak", "Holm-Sidak (step down method using Sidak adjustments)"),
        ("holm", "Holm (step-down method using Bonferroni adjustments)"),
        ("fdr_tsbh", "Two step Benjamini and Hochberg procedure"),
        (
            "fdr_tsbky",
            "Two step estimation method of Benjamini, Krieger, and Yekutieli",
        ),
    ],
    default="fdr_bh",
)
file_field = FileField("File", validators=[DataRequired()])
species_field = RadioField(
    "Species",
    choices=[
        ("human", "Human"),
        ("rat", "Rat"),
        ("mouse", "Mouse"),
    ],
    default="human",
)
permutations_field = IntegerField(
    "Permutations",
    default=100,
    validators=[DataRequired()],
    description="The number of permutations used with GSEA",
)
