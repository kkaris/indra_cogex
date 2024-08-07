# -*- coding: utf-8 -*-

"""Processors for generating nodes and relations for upload to Neo4j."""

from class_resolver import Resolver

from .bgee import BgeeProcessor
from .cbioportal import (
    CcleCnaProcessor,
    CcleDrugResponseProcessor,
    CcleMutationsProcessor,
)
from .cellmarker import CellMarkerProcessor
from .chembl import ChemblIndicationsProcessor
from .clinicaltrials import ClinicaltrialsProcessor
from .disgenet import DisgenetProcessor
from .ec import HGNCEnzymeProcessor
from .goa import GoaProcessor
from .hpoa import HpDiseasePhenotypeProcessor, HpPhenotypeGeneProcessor
from .indra_db import DbProcessor, EvidenceProcessor
from .indra_ontology import OntologyProcessor
from .interpro import InterproProcessor
from .nih_reporter import NihReporterProcessor
from .pathways import ReactomeProcessor, WikipathwaysProcessor
from .processor import Processor
from .pubmed import PublicationProcessor, JournalProcessor
from .sider import SIDERSideEffectProcessor
from .wikidata import JournalPublisherProcessor
from .gwas import GWASProcessor
from .depmap import DepmapProcessor

__all__ = [
    "processor_resolver",
    "Processor",
    "BgeeProcessor",
    "ReactomeProcessor",
    "WikipathwaysProcessor",
    "GoaProcessor",
    "DbProcessor",
    "OntologyProcessor",
    "CcleCnaProcessor",
    "CcleMutationsProcessor",
    "CcleDrugResponseProcessor",
    "ClinicaltrialsProcessor",
    "ChemblIndicationsProcessor",
    "SIDERSideEffectProcessor",
    "EvidenceProcessor",
    "PublicationProcessor",
    "JournalProcessor",
    "HpDiseasePhenotypeProcessor",
    "HpPhenotypeGeneProcessor",
    "NihReporterProcessor",
    "InterproProcessor",
    "CellMarkerProcessor",
    "JournalPublisherProcessor",
    "DisgenetProcessor",
    "GWASProcessor",
    "HGNCEnzymeProcessor",
    "DepmapProcessor",
]

processor_resolver: Resolver[Processor] = Resolver.from_subclasses(Processor)
