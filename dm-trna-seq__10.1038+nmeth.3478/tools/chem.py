"""Source-owned structured-RNA library declaration and reaction graph."""
from __future__ import annotations
from types import SimpleNamespace
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from chemdraw import Construct, Scene, Segment, Workflow, feature
from circular import CircularizationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from protocol_families.structured_rna import StructuredRNAInput, exact, _finish_rt_pcr, _finish_pcr

TITLE = 'DM-tRNA-seq (DM-TGIRT-seq)'
DOI = '10.1038/nmeth.3478'
CITATION = 'Zheng et al., Nature Methods (2015)'
SUMMARY = 'Wild-type and D135S AlkB remove major Watson–Crick-face methylations before TGIRT template switching, cDNA circularization and indexed PCR.'
NOTES = '01_dm_trna_seq.html'
SEQUENCING_ENDING = 'The source gives both library PCR primers and states Illumina HiSeq sequencing, but not a separate sequencing-primer oligo.'
CAVEAT = ''
SOURCE = 'Zheng et al., Nature Methods (2015); <a href="https://doi.org/10.1038/nmeth.3478">doi:10.1038/nmeth.3478</a>.'
FINAL_CAPTION = 'Role-level final library. Exact source-published oligos are listed above; undisclosed kit junctions remain role tokens.'
SEQUENCING_INTRO = ''
SEQUENCING_UNAVAILABLE = 'The defining source does not publish a distinct sequencing-primer oligo or enough kit sequence to compute its exact landing site.'
SEQ_PRIMERS = ()
READ_LENGTHS = None
oligos = lambda: ('TGIRT RNA template — 5′-AGAUCGGAAGAGCACACGUCUAGUUCUACAGUCCGACGAUC/3SpC3/-3′',
 'multiplex PCR primer — 5′-AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC-3′',
 'indexed PCR primer — '
 '5′-CAAGCAGAAGACGGCATACGAGAT[BARCODE]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT-3′')
_definition = SimpleNamespace(title=TITLE, doi=DOI, citation=CITATION, endpoint=SEQUENCING_ENDING)

DM_TEMPLATE = exact("DM-TGIRT RNA template", "AGATCGGAAGAGCACACGTCTAGTTCTACAGTCCGACGATC", "r1")
DM_FORWARD = exact("DM-TGIRT multiplex PCR primer", "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC", "p5")
DM_REVERSE = exact("DM-TGIRT indexed PCR primer", "CAAGCAGAAGACGGCATACGAGAT[BARCODE]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT", "p7", placeholder=True)

FINAL_LIBRARY = Construct((role("P5 flow-cell end"), *(role("DM-TGIRT starter-derived entry"), role("tRNA cDNA")),
    role("sample index", tag="cbc", length_bp=6,
         molecular_feature=feature("sample_index", "sample_index", "unknown")),
    role("P7 flow-cell end")), name=TITLE)

def _dm() -> Workflow:
    source = StructuredRNAInput("deacylated mature tRNA", "free 3-prime CCA",
                                "wild-type AlkB plus AlkB-D135S", False)
    return chain(source.scene(), (
        ("Demethylate Watson–Crick-face m1A, m3C and m1G",
         molecule("demethylated tRNA", (role("structured tRNA with relieved RT blocks"),
                                         role("free 3-prime CCA"))),
         "Wild-type AlkB targets m1A/m3C and D135S improves m1G removal."),
        ("Template-switch TGIRT from a blocked RNA/DNA starter duplex",
         molecule("full-length tRNA cDNA", (DM_TEMPLATE, role("tRNA-derived cDNA"))),
         "The DNA primer carries a one-base 3′ T overhang; TGIRT switches onto the tRNA CCA end."),
        ("Circularize the single-stranded cDNA",
         molecule("circular DM-TGIRT cDNA", (role("CircLigase-II-closed cDNA"),)),
         "Circularization places the two PCR sites around one continuous cDNA template."),
        ("Amplify with the multiplex and indexed primers",
         molecule("DM-tRNA-seq library", (DM_FORWARD, role("tRNA cDNA"), DM_REVERSE)),
         "The defining paper sequenced the indexed products on an Illumina HiSeq."),
    ))

workflow = _dm

def _dm_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    trna = graph.input(
        "trna", "deacylated mature tRNA",
        MaterialKind.BIOLOGICAL,
        state=StructuredRNAInput(
            "deacylated mature tRNA", "free 3-prime CCA",
            "wild-type AlkB plus AlkB-D135S", False).scene(), primary=True)
    alkb = graph.input("alkb", "wild-type AlkB plus AlkB-D135S",
                       MaterialKind.CATALYST)
    demethylated = graph.product(
        "demethylated_trna", "demethylated structured tRNA",
        MaterialKind.MOLECULAR,
        state=molecule("demethylated tRNA", (
            role("structured tRNA with relieved RT blocks"), role("free 3-prime CCA"))))
    graph.react("demethylate", "Demethylate m1A, m3C and m1G",
                StepKind.CHEMICAL_CONVERSION,
                ((trna, InputRole.TEMPLATE), (alkb, InputRole.CATALYST)),
                (demethylated,))
    starter = graph.input("starter_duplex", "blocked DM-TGIRT RNA/DNA starter duplex",
                          MaterialKind.OLIGO, state=molecule("TGIRT starter duplex", (DM_TEMPLATE,)))
    tgirt = graph.input("tgirt", "thermostable group-II-intron reverse transcriptase",
                        MaterialKind.CATALYST)
    cdna = graph.product(
        "full_length_cdna", "full-length tRNA cDNA",
        MaterialKind.MOLECULAR,
        state=molecule("full-length tRNA cDNA", (DM_TEMPLATE, role("tRNA-derived cDNA"))))
    graph.react("template_switch", "Template-switch TGIRT onto the tRNA CCA end",
                StepKind.TEMPLATE_SWITCHING,
                ((demethylated, InputRole.TEMPLATE), (starter, InputRole.STRUCTURAL),
                 (tgirt, InputRole.CATALYST)), (cdna,))
    circligase = graph.input("circligase", "CircLigase II", MaterialKind.CATALYST)
    circle = graph.product(
        "circular_cdna", "circular DM-TGIRT cDNA", MaterialKind.MOLECULAR,
        state=molecule("circular DM-TGIRT cDNA", (role("CircLigase-II-closed cDNA"),)))
    graph.react("circularize", "Circularize the single-stranded cDNA",
                StepKind.CIRCULARIZATION,
                ((cdna, InputRole.TEMPLATE), (circligase, InputRole.CATALYST)),
                (circle,), geometry=CircularizationGeometry.known(
                    "single-strand intramolecular ligation",
                    statement="The protocol circularizes single-stranded cDNA with CircLigase II."))
    reviewed = [
        "trna", "alkb", "demethylate", "starter_duplex", "tgirt",
        "template_switch", "circligase", "circularize",
    ]
    _finish_pcr(graph, "library", circle, final, d.endpoint, reviewed)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    return _dm_graph(_definition, FINAL_LIBRARY)
