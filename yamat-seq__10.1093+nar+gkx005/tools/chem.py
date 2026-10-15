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

TITLE = 'YAMAT-seq — Y-shaped adapter mature-tRNA sequencing'
DOI = '10.1093/nar/gkx005'
CITATION = 'Shigematsu et al., Nucleic Acids Research (2017)'
SUMMARY = 'A four-way Y adapter recognizes the conserved mature-tRNA NCCA overhang and T4 RNA ligase 2 seals both acceptor-stem nicks in one reaction.'
NOTES = '01_yamat_seq.html'
SEQUENCING_ENDING = 'A custom small-RNA sequencing primer reads through the adapter into mature-tRNA cDNA; the experiment used a 100-nt single read.'
CAVEAT = ''
SOURCE = 'Shigematsu et al., Nucleic Acids Research (2017); <a href="https://doi.org/10.1093/nar/gkx005">doi:10.1093/nar/gkx005</a>.'
FINAL_CAPTION = 'Role-level final library. Exact source-published oligos are listed above; undisclosed kit junctions remain role tokens.'
SEQUENCING_INTRO = ''
SEQUENCING_UNAVAILABLE = 'The defining source does not publish a distinct sequencing-primer oligo or enough kit sequence to compute its exact landing site.'
SEQ_PRIMERS = ()
READ_LENGTHS = None
oligos = lambda: ('Y-3′-AD — 5′-/5Phos/GTATCCAGTTGGAATTCTCGGGTGCCAAGG/3ddC/-3′',
 'Y-5′-AD-A — 5′-GTTCAGAGTTCTACAGTCCGACGATCACTGGATACTGga-3′',
 'Y-5′-AD-G — 5′-GTTCAGAGTTCTACAGTCCGACGATCACTGGATACTGgg-3′',
 'Y-5′-AD-C — 5′-GTTCAGAGTTCTACAGTCCGACGATCACTGGATACTGgc-3′',
 'Y-5′-AD-U — 5′-GTTCAGAGTTCTACAGTCCGACGATCACTGGATACTGgu-3′')
_definition = SimpleNamespace(title=TITLE, doi=DOI, citation=CITATION, endpoint=SEQUENCING_ENDING)

Y3 = exact("Y-3-prime-AD", "GTATCCAGTTGGAATTCTCGGGTGCCAAGG", "r2")
Y5_A = exact("Y-5-prime-AD-A", "GTTCAGAGTTCTACAGTCCGACGATCACTGGATACTGGA", "r1")
Y5_G = exact("Y-5-prime-AD-G", "GTTCAGAGTTCTACAGTCCGACGATCACTGGATACTGGG", "r1")
Y5_C = exact("Y-5-prime-AD-C", "GTTCAGAGTTCTACAGTCCGACGATCACTGGATACTGGC", "r1")
Y5_U = exact("Y-5-prime-AD-U", "GTTCAGAGTTCTACAGTCCGACGATCACTGGATACTGGT", "r1")

FINAL_LIBRARY = Construct((role("P5 flow-cell end"), *(role("Y-5-prime adapter-derived read entry"), role("mature tRNA cDNA"), role("Y-3-prime adapter-derived end")),
    role("sample index", tag="cbc", length_bp=6,
         molecular_feature=feature("sample_index", "sample_index", "unknown")),
    role("P7 flow-cell end")), name=TITLE)

def _yamat() -> Workflow:
    source = StructuredRNAInput("mature aminoacylated tRNA", "3-prime NCCA end")
    return chain(source.scene(), (
        ("Deacylate mature tRNA at alkaline pH",
         molecule("deacylated mature tRNA", (role("paired acceptor stem"), role("free 3-prime NCCA"))),
         "The free 3′-OH and conserved NCCA overhang define the ligation substrate."),
        ("Anneal the four Y-5′ adapters with Y-3′-AD",
         molecule("Y-adapter–tRNA complex", (Y5_A, role("mature tRNA"), Y3)),
         "The four terminal adapter variants complement the four possible discriminator bases."),
        ("Seal both nicks with T4 RNA ligase 2",
         molecule("covalently Y-adapted tRNA", (Y5_A, role("mature tRNA"), Y3),
                  junctions=(("Y-5-prime-AD-A", "mature tRNA", "Rnl2 nick ligation"),
                             ("mature tRNA", "Y-3-prime-AD", "Rnl2 nick ligation"))),
         "Rnl2 seals the RNA/RNA and RNA/DNA nicks in the annealed acceptor-stem complex."),
        ("Reverse-transcribe, PCR and gel-select the tRNA product",
         molecule("YAMAT-seq library", (role("P5 and small-RNA read arm"),
                  role("mature tRNA cDNA"), role("indexed P7 arm"))),
         "The defining experiment used TruSeq Small RNA reagents and 100-nt single-read HiSeq sequencing."),
    ))

workflow = _yamat

def _yamat_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    trna = graph.input(
        "aminoacylated_trna", "mature aminoacylated tRNA",
        MaterialKind.BIOLOGICAL,
        state=StructuredRNAInput("mature aminoacylated tRNA", "3-prime NCCA end").scene(),
        primary=True)
    alkaline = graph.input("alkaline_condition", "alkaline deacylation condition",
                           MaterialKind.CONDITION)
    deacylated = graph.product(
        "deacylated_trna", "deacylated mature tRNA with free NCCA",
        MaterialKind.MOLECULAR,
        state=molecule("deacylated mature tRNA", (
            role("paired acceptor stem"), role("free 3-prime NCCA"))))
    graph.react("deacylate", "Deacylate mature tRNA at alkaline pH",
                StepKind.CHEMICAL_CONVERSION,
                ((trna, InputRole.TEMPLATE), (alkaline, InputRole.CONDITION)),
                (deacylated,))
    y_adapter = graph.input("y_adapter", "Y-3-prime-AD plus four Y-5-prime adapter variants",
                            MaterialKind.ADAPTER, state=molecule("Y adapter set", (
                                Y5_A, Y5_G, Y5_C, Y5_U, Y3)))
    complex_ = graph.product(
        "y_complex", "annealed Y-adapter–tRNA complex",
        MaterialKind.MOLECULAR,
        state=molecule("Y-adapter–tRNA complex", (Y5_A, role("mature tRNA"), Y3)))
    graph.react("anneal", "Anneal the discriminator-matched Y adapter",
                StepKind.HYBRIDIZATION,
                ((deacylated, InputRole.TEMPLATE), (y_adapter, InputRole.STRUCTURAL)),
                (complex_,))
    rnl2 = graph.input("rnl2", "T4 RNA ligase 2", MaterialKind.CATALYST)
    adapted = graph.product(
        "adapted_trna", "covalently Y-adapted tRNA", MaterialKind.MOLECULAR,
        state=molecule("covalently Y-adapted tRNA", (Y5_A, role("mature tRNA"), Y3),
            junctions=(("Y-5-prime-AD-A", "mature tRNA", "Rnl2 nick ligation"),
                       ("mature tRNA", "Y-3-prime-AD", "Rnl2 nick ligation"))))
    graph.react("seal_nicks", "Seal both Y-adapter nicks", StepKind.NICK_SEALING,
                ((complex_, InputRole.TEMPLATE), (y_adapter, InputRole.STRUCTURAL),
                 (rnl2, InputRole.CATALYST)), (adapted,))
    reviewed = [
        "aminoacylated_trna", "alkaline_condition", "deacylate", "y_adapter",
        "anneal", "rnl2", "seal_nicks",
    ]
    _finish_rt_pcr(graph, "library", adapted, final, d.endpoint, reviewed,
                   rt_primer_label="Y-adapter reverse-transcription primer")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    return _yamat_graph(_definition, FINAL_LIBRARY)
