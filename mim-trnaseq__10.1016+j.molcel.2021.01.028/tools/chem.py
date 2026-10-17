"""Source-owned structured-RNA library declaration and reaction graph."""
from __future__ import annotations
from dataclasses import replace
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

TITLE = 'mim-tRNAseq — modification-induced misincorporation tRNA sequencing'
DOI = '10.1016/j.molcel.2021.01.028'
CITATION = 'Behrens et al., Molecular Cell (2021)'
SUMMARY = 'Barcoded 3′-adapter ligation and prolonged low-salt TGIRT reverse transcription preserve modification-dependent mismatches while producing full-length tRNA cDNA.'
NOTES = '01_mim_trnaseq.html'
SEQUENCING_ENDING = 'The indexed library was sequenced for 150 cycles on an Illumina NextSeq; the source does not list a separate custom sequencing primer.'
CAVEAT = ''
SOURCE = 'Behrens et al., Molecular Cell (2021); <a href="https://doi.org/10.1016/j.molcel.2021.01.028">doi:10.1016/j.molcel.2021.01.028</a>.'
FINAL_CAPTION = 'Role-level final library. Exact source-published oligos are listed above; undisclosed kit junctions remain role tokens.'
SEQUENCING_INTRO = ''
SEQUENCING_UNAVAILABLE = 'The defining source does not publish a distinct sequencing-primer oligo or enough kit sequence to compute its exact landing site.'
SEQ_PRIMERS = ()
READ_LENGTHS = None
oligos = lambda: ('3′ adapter I1 — 5′-/5Phos/GATATCGTCAAGATCGGAAGAGCACACGTCTGAA/ddC/-3′',
 'primer-dependent RT oligo — '
 '5′-/5Phos/RNAGATCGGAAGAGCGTCGTGTAGGGAAAGAG/iSp18/GTGACTGGAGTTCAGACGTGTGCTC-3′',
 'library forward — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCT*C-3′',
 'indexed reverse — 5′-CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGT*G-3′')
_definition = SimpleNamespace(title=TITLE, doi=DOI, citation=CITATION, endpoint=SEQUENCING_ENDING)

MIM_ADAPTER_I1 = exact("mim 3-prime adapter I1", "GATATCGTCAAGATCGGAAGAGCACACGTCTGAA", "r2")
MIM_RT = exact("mim primer-dependent RT oligo", "[5-prime-pRN]AGATCGGAAGAGCGTCGTGTAGGGAAAGAG[iSp18]GTGACTGGAGTTCAGACGTGTGCTC", "r2", placeholder=True)
MIM_FORWARD = exact("mim library forward primer", "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC", "p5")
MIM_REVERSE = exact("mim indexed reverse primer", "CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTG", "p7", placeholder=True,
                    feature=feature("sample_index", "sample_index", "random"))

FINAL_LIBRARY = Construct((role("P5 flow-cell end"), *(role("mim adapter-derived entry"), role("native-modification tRNA cDNA")),
    role("sample index", tag="cbc", length_bp=6,
         molecular_feature=feature("sample_index", "sample_index", "unknown")),
    role("P7 flow-cell end")), name=TITLE)

def _mim() -> Workflow:
    source = StructuredRNAInput("60–100-nt native mature tRNA", "native 3-prime CCA")
    return chain(source.scene(), (
        ("Ligate one of four preadenylated barcoded 3′ adapters",
         molecule("adapter-ligated tRNA", (role("native modified tRNA"), MIM_ADAPTER_I1),
                  junctions=(("native modified tRNA", "mim 3-prime adapter I1", "Rnl2tr KQ ligation"),)),
         "The adapter barcode permits samples to be pooled before reverse transcription."),
        ("Reverse-transcribe for 16 h with TGIRT under low-salt conditions",
         molecule("full-length misincorporation-bearing cDNA", (MIM_RT,
                  role("tRNA cDNA retaining modification signatures"))),
         "Conditions favour readthrough rather than enzymatically removing the native modifications."),
        ("Hydrolyse RNA, size-select cDNA and circularize with CircLigase",
         molecule("circular mim-tRNAseq cDNA", (role("circular full-length tRNA cDNA"),)),
         "The RT primer's 5′ RN dinucleotide supports circularization."),
        ("PCR with common P5 and indexed P7 primers",
         molecule("mim-tRNAseq library", (MIM_FORWARD,
                  role("tRNA cDNA with misincorporations"), MIM_REVERSE)),
         "Five to six cycles precede 150-cycle Illumina NextSeq sequencing."),
    ))

workflow = _mim

def _mim_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    trna = graph.input(
        "native_trna", "60–100-nt native mature tRNA",
        MaterialKind.BIOLOGICAL,
        state=StructuredRNAInput("60–100-nt native mature tRNA",
                                 "native 3-prime CCA").scene(), primary=True)
    adapter = graph.input("barcoded_adapter", "preadenylated barcoded 3-prime adapter",
                          MaterialKind.ADAPTER, state=molecule("mim 3-prime adapter", (
                              MIM_ADAPTER_I1,)))
    ligase = graph.input("adapter_ligase", "truncated T4 RNA ligase 2 KQ",
                         MaterialKind.CATALYST)
    adapted = graph.product(
        "adapted_trna", "barcoded 3-prime-adapter-ligated tRNA",
        MaterialKind.MOLECULAR,
        state=molecule("adapter-ligated tRNA", (
            role("native modified tRNA"), MIM_ADAPTER_I1),
            junctions=(("native modified tRNA", "mim 3-prime adapter I1",
                        "Rnl2tr KQ ligation"),)))
    known_end = reported_adapter_end(
        "mim preadenylated 3-prime adapter", kind="single-strand",
        polarity="5-prime", terminal_sequence="G",
        statement="The source publishes the preadenylated, 3-prime-blocked adapter sequence and its RNA-ligation role.")
    graph.react("adapter_ligation", "Ligate a preadenylated barcoded 3-prime adapter",
                StepKind.LIGATION,
                ((trna, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,), geometry=known_end)
    rt_oligo = graph.input("rt_oligo", "published mim primer-dependent RT oligo",
                           MaterialKind.OLIGO, state=molecule("mim RT oligo", (MIM_RT,)))
    tgirt = graph.input("tgirt", "TGIRT reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product(
        "full_length_cdna", "full-length modification-signature cDNA",
        MaterialKind.MOLECULAR,
        state=molecule("full-length misincorporation-bearing cDNA", (
            MIM_RT, role("tRNA cDNA retaining modification signatures"))))
    graph.react("reverse_transcription", "Reverse-transcribe for 16 hours at low salt",
                StepKind.REVERSE_TRANSCRIPTION,
                ((adapted, InputRole.TEMPLATE), (rt_oligo, InputRole.STRUCTURAL),
                 (tgirt, InputRole.CATALYST)), (cdna,))
    hydrolysis = graph.input("rna_hydrolysis", "RNA hydrolysis condition",
                             MaterialKind.CONDITION)
    clean_cdna = graph.product("clean_cdna", "RNA-free full-length tRNA cDNA",
                               MaterialKind.MOLECULAR, state=cdna.state)
    degraded_rna = graph.product("degraded_rna", "hydrolysed RNA fragments",
                                 MaterialKind.WASTE)
    graph.react("hydrolyse_rna", "Hydrolyse RNA with NaOH (95°C, 3 min)",
                StepKind.CHEMICAL_CONVERSION,
                ((cdna, InputRole.TEMPLATE), (hydrolysis, InputRole.CONDITION)),
                (clean_cdna, degraded_rna))
    graph.discard(degraded_rna, "the full-length cDNA is the retained library substrate")
    size_rule = graph.input("cdna_size_selection", "full-length cDNA size-selection rule",
                            MaterialKind.SELECTION)
    sized = graph.product("sized_cdna", "size-selected full-length tRNA cDNA",
                          MaterialKind.MOLECULAR, state=clean_cdna.state)
    rejected = graph.product("short_cdna", "short or incomplete cDNA",
                             MaterialKind.WASTE)
    graph.react("select_cdna", "Size-select full-length cDNA", StepKind.SELECTION,
                ((clean_cdna, InputRole.TEMPLATE), (size_rule, InputRole.SELECTION)),
                (sized, rejected))
    graph.discard(rejected, "does not span the intended mature-tRNA product")
    circligase = graph.input("circligase", "CircLigase", MaterialKind.CATALYST)
    circle = graph.product(
        "circular_cdna", "circular full-length tRNA cDNA",
        MaterialKind.MOLECULAR,
        state=molecule("circular mim-tRNAseq cDNA", (
            role("circular full-length tRNA cDNA"),)))
    graph.react("circularize", "Circularize cDNA with CircLigase (60°C, 3 h)",
                StepKind.CIRCULARIZATION,
                ((sized, InputRole.TEMPLATE), (circligase, InputRole.CATALYST)),
                (circle,), geometry=CircularizationGeometry.known(
                    "single-strand intramolecular ligation",
                    statement="The protocol circularizes size-selected full-length cDNA with CircLigase."))
    reviewed = [
        "native_trna", "barcoded_adapter", "adapter_ligase", "adapter_ligation",
        "rt_oligo", "tgirt", "reverse_transcription", "rna_hydrolysis",
        "hydrolyse_rna", "cdna_size_selection", "select_cdna", "circligase",
        "circularize",
    ]
    _finish_pcr(graph, "library", circle, final, d.endpoint, reviewed)
    graph.steps["library_pcr"] = replace(graph.steps["library_pcr"], label="Amplify with KAPA HiFi (5–6 cycles)")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    return _mim_graph(_definition, FINAL_LIBRARY)
