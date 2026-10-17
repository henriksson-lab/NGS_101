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

TITLE = 'PANDORA-seq — pan-small-RNA sequencing'
DOI = '10.1038/s41556-021-00652-7'
CITATION = 'Shi et al., Nature Cell Biology (2021)'
SUMMARY = 'T4 PNK repairs blocked RNA termini and AlkB removes internal methylation barriers before conventional small-RNA adapter ligation.'
NOTES = '01_pandora_seq.html'
SEQUENCING_ENDING = 'NEBNext small-RNA libraries were read single-end for 50 or 75 cycles; kit sequencing-primer sequences are not reproduced in the paper.'
CAVEAT = 'The adapter and primer strings below come from a later NEB E7330 manual, not a verified historical 2021 kit revision. They describe that vendor kit family; the experiment-specific index and historical revision remain unassigned.'
SOURCE = 'Shi et al., Nature Cell Biology (2021); <a href="https://doi.org/10.1038/s41556-021-00652-7">doi:10.1038/s41556-021-00652-7</a>.'
FINAL_CAPTION = 'Role-level final library. Exact source-published oligos are listed above; undisclosed kit junctions remain role tokens.'
SEQUENCING_INTRO = ''
SEQUENCING_UNAVAILABLE = 'NEB confirms standard Illumina sequencing primers for E7330 libraries; this role-level final molecule has not yet reconstructed the exact binding site.'
SEQ_PRIMERS = ()
READ_LENGTHS = None
oligos = lambda: ('NEBNext E7330 3′ SR adapter — 5′-rAppAGATCGGAAGAGCACACGTCT-NH2-3′',
 'NEBNext E7330 5′ SR adapter — '
 '5′-rGrUrUrCrArGrArGrUrUrCrUrArCrArGrUrCrCrGrArCrGrArUrC-3′',
 'NEBNext E7330 SR RT primer — 5′-AGACGTGTGCTCTTCCGATCT-3′',
 'NEBNext E7330 SR PCR primer — '
 '5′-AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCG*A-3′',
 "NEBNext E7330 index 1 PCR primer (kit example; experiment's index choice unstated) — "
 '5′-CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC*T-3′')
_definition = SimpleNamespace(title=TITLE, doi=DOI, citation=CITATION, endpoint=SEQUENCING_ENDING)

PANDORA_RA3 = exact("NEBNext 3-prime SR adapter", "AGATCGGAAGAGCACACGTCT", "r2")
PANDORA_RA5 = exact("NEBNext 5-prime SR adapter", "GTTCAGAGTTCTACAGTCCGACGATC", "r1")
PANDORA_RT = exact("NEBNext SR RT primer", "AGACGTGTGCTCTTCCGATCT", "r2")

FINAL_LIBRARY = Construct((role("P5 flow-cell end"), *(role("small-RNA read entry"), role("repaired small-RNA cDNA")),
    role("sample index", tag="cbc", length_bp=6,
         molecular_feature=feature("sample_index", "sample_index", "unknown")),
    role("P7 flow-cell end")), name=TITLE)

def _pandora() -> Workflow:
    source = StructuredRNAInput("15–50-nt cellular small RNA", "heterogeneous blocked termini",
                                "AlkB demethylation", False)
    return chain(source.scene(), (
        ("Repair RNA termini with T4 polynucleotide kinase",
         molecule("end-repaired small RNA", (role("5-prime phosphate"),
                  role("small RNA"), role("3-prime hydroxyl"))),
         "PNK converts otherwise poorly ligatable terminal chemistries to 5′-P/3′-OH."),
        ("Demethylate internal RT-blocking bases with AlkB",
         molecule("PNK- and AlkB-treated small RNA", (role("ligatable demethylated small RNA"),)),
         "The two treatments address independent terminal and internal barriers."),
        ("Ligate the 3′ and 5′ small-RNA adapters sequentially",
         molecule("adapter-ligated small RNA", (PANDORA_RA5,
                  role("repaired small RNA"), PANDORA_RA3),
                  junctions=(("NEBNext 5-prime SR adapter", "repaired small RNA", "RNA ligation"),
                             ("repaired small RNA", "NEBNext 3-prime SR adapter", "RNA ligation"))),
         "The article names E7330S; the current NEB E7330 manual gives both adapter sequences."),
        ("Reverse-transcribe, PCR and PAGE-purify",
         molecule("PANDORA-seq library", (role("P5 small-RNA arm"),
                  role("small-RNA cDNA"), role("indexed P7 arm"))),
         "Libraries were sequenced single-end (50 or 75 cycles depending on sample)."),
    ))

workflow = _pandora

def _pandora_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    rna = graph.input(
        "small_rna", "15–50-nt cellular small RNA with heterogeneous blocked termini",
        MaterialKind.BIOLOGICAL,
        state=StructuredRNAInput(
            "15–50-nt cellular small RNA", "heterogeneous blocked termini",
            "AlkB demethylation", False).scene(), primary=True)
    pnk = graph.input("pnk", "T4 polynucleotide kinase", MaterialKind.CATALYST)
    repaired = graph.product(
        "repaired_rna", "small RNA with 5-prime phosphate and 3-prime hydroxyl",
        MaterialKind.MOLECULAR,
        state=molecule("end-repaired small RNA", (
            role("5-prime phosphate"), role("small RNA"), role("3-prime hydroxyl"))))
    graph.react("end_repair", "Repair RNA termini with T4PNK (37°C, 20 min)",
                StepKind.END_REPAIR,
                ((rna, InputRole.TEMPLATE), (pnk, InputRole.CATALYST)), (repaired,))
    alkb = graph.input("alkb", "AlkB demethylase", MaterialKind.CATALYST)
    demethylated = graph.product(
        "demethylated_rna", "PNK- and AlkB-treated small RNA",
        MaterialKind.MOLECULAR,
        state=molecule("PNK- and AlkB-treated small RNA", (
            role("ligatable demethylated small RNA"),)))
    graph.react("demethylate", "Demethylate with AlkB (37°C, 30 min)",
                StepKind.CHEMICAL_CONVERSION,
                ((repaired, InputRole.TEMPLATE), (alkb, InputRole.CATALYST)),
                (demethylated,))
    adapter3 = graph.input("three_adapter", "NEBNext 3-prime small-RNA adapter",
                           MaterialKind.ADAPTER, state=molecule("3-prime adapter", (
                               PANDORA_RA3,)))
    ligase3 = graph.input("three_ligase", "3-prime RNA adapter ligase",
                          MaterialKind.CATALYST)
    ligated3 = graph.product("three_ligated", "3-prime-adapted repaired small RNA",
                             MaterialKind.MOLECULAR, state=demethylated.state)
    graph.react("three_ligation", "Ligate the 3-prime small-RNA adapter",
                StepKind.LIGATION,
                ((demethylated, InputRole.TEMPLATE), (adapter3, InputRole.STRUCTURAL),
                 (ligase3, InputRole.CATALYST)), (ligated3,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The 2025 E7330 manual publishes the adapter sequence and terminal modifications; the model has not reconstructed the complete ligation substrate geometry."))
    adapter5 = graph.input("five_adapter", "NEBNext 5-prime small-RNA adapter",
                           MaterialKind.ADAPTER, state=molecule("5-prime adapter", (
                               PANDORA_RA5,)))
    ligase5 = graph.input("five_ligase", "5-prime RNA adapter ligase",
                          MaterialKind.CATALYST)
    ligated5 = graph.product(
        "adapted_rna", "two-adapter PANDORA small RNA",
        MaterialKind.MOLECULAR,
        state=molecule("adapter-ligated small RNA", (
            PANDORA_RA5, role("repaired small RNA"), PANDORA_RA3),
            junctions=(("NEBNext 5-prime SR adapter", "repaired small RNA", "RNA ligation"),
                       ("repaired small RNA", "NEBNext 3-prime SR adapter", "RNA ligation"))))
    graph.react("five_ligation", "Ligate the 5-prime small-RNA adapter",
                StepKind.LIGATION,
                ((ligated3, InputRole.TEMPLATE), (adapter5, InputRole.STRUCTURAL),
                 (ligase5, InputRole.CATALYST)), (ligated5,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The 2025 E7330 manual publishes the RNA adapter sequence; the model has not reconstructed the complete ligation geometry."))
    reviewed = [
        "small_rna", "pnk", "end_repair", "alkb", "demethylate",
        "three_adapter", "three_ligase", "three_ligation", "five_adapter",
        "five_ligase", "five_ligation",
    ]
    rt_primer = graph.input("library_rt_primer", "NEBNext small-RNA reverse-transcription primer",
                            MaterialKind.OLIGO, state=molecule("SR RT primer", (PANDORA_RT,)))
    reverse_transcriptase = graph.input("library_reverse_transcriptase", "reverse transcriptase",
                                       MaterialKind.CATALYST)
    cdna = graph.product("library_cdna", "PANDORA small-RNA cDNA",
                         MaterialKind.MOLECULAR,
                         state=molecule("PANDORA small-RNA cDNA", (
                             role("repaired small-RNA cDNA"),)))
    graph.react("library_rt", "Reverse-transcribe the adapted small RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((ligated5, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (cdna,))
    reviewed.extend(("library_rt_primer", "library_reverse_transcriptase", "library_rt"))
    library = _finish_pcr(
        graph, "library", cdna, final, d.endpoint, reviewed, sequence=False)
    page = graph.input("page_selection", "PAGE library-size selection rule",
                       MaterialKind.SELECTION)
    selected_library = graph.product("selected_library", "PAGE-purified PANDORA-seq library",
                                     MaterialKind.LIBRARY, state=final)
    rejected_library = graph.product("rejected_library", "off-size amplified products",
                                     MaterialKind.WASTE)
    graph.react("page_select", "PAGE-purify the amplified library", StepKind.SELECTION,
                ((library, InputRole.TEMPLATE), (page, InputRole.SELECTION)),
                (selected_library, rejected_library))
    graph.discard(rejected_library, "falls outside the selected small-RNA library size")
    graph.sequence(selected_library, SequencingPlan(limitation=d.endpoint))
    reviewed.extend(("page_selection", "page_select", "selected_library"))
    graph.steps["library_pcr"] = replace(graph.steps["library_pcr"], label="Amplify small-RNA cDNA (11–22 cycles)")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    return _pandora_graph(_definition, FINAL_LIBRARY)
