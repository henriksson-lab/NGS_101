"""Molecular model for MATQ-seq (Sheng et al., 2017)."""
from __future__ import annotations

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from batch_ngs import seg, truseq_library
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, MolecularState, Scene, Workflow, annotation_rows, oligo
from rna_amplification import poly_a_capture, terminal_tail
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role

GAT27 = "GTGAGTGATGGTTGAGGATGTGTGGAG"
GAT27_DT = GAT27 + "N" * 5 + "T" * 20
GAT27_5N3G = GAT27 + "N" * 5 + "GGG"
GAT27_5N3T = GAT27 + "N" * 5 + "TTT"
GAT21_6N3G = "GATGGTTGAGGATGTGTGGAG" + "N" * 6 + "GGG"
THREE_N_GAT24 = "N" * 3 + "AGTGATGGTTGAGGATGTGTGGAG"


def initial_scene():
    primer = [seg("GAT27 handle", GAT27, "tso"),
              seg("five random bases", "N" * 5, placeholder=True),
              seg("dT20", "T" * 20)]
    sc = poly_a_capture(primer, annealing_segment="dT20",
                        primer_label="GAT27dT primer")
    sc.note("primer", "the same RT also contains GAT27-5N3G and GAT27-5N3T random primers")
    return sc


FIRST_STRAND = Construct([
    seg("GAT27 handle", GAT27, "tso"), seg("primer N5", "N" * 5, placeholder=True),
    seg("first-strand cDNA", "X" * 38, placeholder=True)], name="MATQ first strand")


def first_strand_scene():
    sc = Scene(); sc.strand("cDNA", list(FIRST_STRAND), label="first-strand cDNA")
    sc.arrow("cDNA", "ten stepped-temperature annealing/extension cycles")
    sc.labels("cDNA")
    return sc


TAILED, _TAIL_SCENE = terminal_tail(FIRST_STRAND, "C", 12, enzyme="TdT")


def second_strand_scene():
    primer = [seg("GAT21 handle", "GATGGTTGAGGATGTGTGGAG", "tso"),
              seg("random N6", "N" * 6, placeholder=True), seg("3G", "GGG")]
    sc = Scene(); sc.strand("dC-tailed cDNA", list(TAILED), label="dC-tailed first strand")
    sc.anneal("GAT21-6N3G", primer, to="dC-tailed cDNA",
              pair=("3G", "poly(dC) tail"), label="second-strand primer",
              unpaired=("GAT21 handle", "random N6"))
    sc.arrow("GAT21-6N3G", "Deep Vent exo−; repeated annealing and extension")
    sc.labels("GAT21-6N3G")
    return sc


FINAL_LIBRARY, SEQ_PRIMERS = truseq_library(
    [seg("sheared MATQ cDNA", "X" * 42, placeholder=True)],
    "MATQ-seq paired-end library", inferred_adapters=True)


def final_rows():
    sc = Scene.duplex(list(FINAL_LIBRARY), label="MATQ-seq")
    sc.junction("top", "Read 1 arm", "sheared MATQ cDNA", "adapter ligation")
    sc.junction("top", "dA junction", "Index 1 / Read 2 arm", "adapter ligation")
    sc.labels("top")
    return [*sc.rows(), *annotation_rows(FINAL_LIBRARY)]


def workflow():
    wf = Workflow(MolecularState("Mixed oligo-dT and random priming",
                                 tuple(initial_scene().rows())))
    wf.react("SuperScript III reverse transcription", first_strand_scene().rows(),
             note="Ten low-to-high-temperature cycles repeatedly anneal the mixed primers across total RNA.")
    wf.react("T4 DNA polymerase primer digestion, optional Exonuclease I, then RNase H/RNase I",
             first_strand_scene().rows(),
             note="Unused primers and the RNA template are removed before homopolymer tailing.")
    wf.react("TdT adds a poly(dC) tail", _TAIL_SCENE.rows(),
             note="The common GAT27-bearing first strand receives a terminal dC tract.")
    wf.react("GAT21-6N3G priming and second-strand synthesis",
             second_strand_scene().rows(),
             note="The terminal GGG anchors at the dC tail; random bases promote repeated annealing/looping amplification.")
    wf.ligate_adapter("GAT27 PCR, shearing, end repair, dA tailing and TruSeq adapter ligation",
             AdapterLigationGeometry.unavailable("The defining article names TruSeq adapters but does not print their selected oligo or ligation-ready end geometry."), final_rows(),
             note="T4 Quick Ligase attaches TruSeq adapters; the exact adapter variant and molecular sequence are not given.")
    return wf


TITLE = "MATQ-seq — multiple annealing and dC-tailing single-cell RNA-seq"
NOTES = "01_matq-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/nmeth.4145">'
          'Sheng et al., <i>Nature Methods</i> (2017)</a>.')
SUMMARY = ("A mixture of oligo-dT and random primers repeatedly samples total RNA. After RT, "
           "TdT adds dC, a G-anchored primer makes the second strand, and one common handle amplifies the transcriptome.")
CAVEAT = ("INFERRED — the defining article names TruSeq adapters but does not identify the selected adapter/index oligo; terminal bases show canonical TruSeq geometry. "
          "Published read-length discrepancy: Online Methods say a NextSeq 500 2×75 midoutput kit (cost paragraph), while Supplementary Fig. 21 says paired-end reads of 85×2 bases. This model uses the explicitly named 2×75 kit; the origin of the 85×2 figure remains unresolved.")
FINAL_CAPTION = "INFERRED — representative sheared MATQ amplicon in a canonical paired-end Illumina library; the GAT27 amplification products themselves are source-defined."
SEQUENCING_INTRO = "The defining Online Methods name a NextSeq 500 2×75 run. Standard TruSeq primer sites are inferred because the selected adapter/index oligos are not printed. Supplementary Fig. 21 instead states 85×2 reads; see the discrepancy notice above."
READ_LENGTHS = {"Read 1": 75, "Read 2": 75}
SHOW_READ_CONTENT = True


def oligos():
    yield oligo("GAT27dT", [seg("GAT27", GAT27, "tso"),
                             seg("N5", "N" * 5, placeholder=True), seg("dT20", "T" * 20)])
    yield oligo("GAT27-5N3G", [seg("GAT27", GAT27, "tso"),
                                seg("N5", "N" * 5, placeholder=True), seg("3G", "GGG")])
    yield oligo("GAT27-5N3T", [seg("GAT27", GAT27, "tso"),
                                seg("N5", "N" * 5, placeholder=True), seg("3T", "TTT")])
    yield oligo("GAT21-6N3G", [seg("GAT21", "GATGGTTGAGGATGTGTGGAG", "tso"),
                                seg("N6", "N" * 6, placeholder=True), seg("3G", "GGG")])
    yield oligo("GAT27 PCR", [seg("GAT27", GAT27, "tso")])
    yield oligo("3NGAT24", [seg("N3", "N" * 3, placeholder=True),
                             seg("GAT24", "AGTGATGGTTGAGGATGTGTGGAG", "tso")])


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    rna = graph.input("rna", "single-cell total RNA", MaterialKind.BIOLOGICAL,
                      state=molecule("total RNA", (role("polyadenylated and non-polyadenylated transcripts"),)),
                      primary=True)
    primers = graph.input("mixed_rt_primers", "GAT27 oligo-dT and random-primer mixture",
                          MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "SuperScript III", MaterialKind.CATALYST)
    first = graph.product("first_strand", "GAT27-bearing first-strand cDNA",
                          MaterialKind.MOLECULAR, state=FIRST_STRAND)
    graph.react("reverse_transcribe", "Repeatedly anneal and extend the mixed RT primers",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna, InputRole.TEMPLATE), (primers, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (first,))
    cleanup = graph.input("cleanup_enzymes", "T4 DNA polymerase, optional Exonuclease I, and RNase H/If",
                          MaterialKind.CATALYST)
    cleaned = graph.product("cleaned", "primer- and RNA-depleted first-strand cDNA",
                            MaterialKind.MOLECULAR, state=FIRST_STRAND)
    graph.react("cleanup", "Remove unused primers and the RNA template", StepKind.DIGESTION,
                ((first, InputRole.TEMPLATE), (cleanup, InputRole.CATALYST)), (cleaned,))
    tdt = graph.input("tdt", "terminal deoxynucleotidyl transferase", MaterialKind.CATALYST)
    tailed = graph.product("tailed", "poly(dC)-tailed first-strand cDNA",
                           MaterialKind.MOLECULAR, state=TAILED)
    graph.react("tail", "Add a terminal poly(dC) tract", StepKind.TAILING,
                ((cleaned, InputRole.TEMPLATE), (tdt, InputRole.CATALYST)), (tailed,))
    second_primer = graph.input("second_strand_primer", "GAT21-6N3G primer",
                                MaterialKind.OLIGO)
    vent = graph.input("deep_vent", "Deep Vent exo-minus polymerase",
                       MaterialKind.CATALYST)
    duplex = graph.product("duplex", "MATQ second-strand amplification product",
                           MaterialKind.MOLECULAR, state=second_strand_scene())
    graph.react("second_strand", "Prime the dC tail and synthesize the second strand",
                StepKind.EXTENSION,
                ((tailed, InputRole.TEMPLATE), (second_primer, InputRole.STRUCTURAL),
                 (vent, InputRole.CATALYST)), (duplex,))
    handle_f = graph.input("gat27_primer", "GAT27 PCR primer", MaterialKind.OLIGO)
    handle_r = graph.input("gat24_primer", "3NGAT24 conversion primer", MaterialKind.OLIGO)
    wta_pol = graph.input("wta_polymerase", "whole-transcriptome PCR polymerase",
                          MaterialKind.CATALYST)
    amplified = graph.product("amplified", "GAT-handle-amplified MATQ cDNA",
                              MaterialKind.MOLECULAR, state=duplex.state)
    graph.react("wta", "24-cycle WTA then20-cycle3NGAT24 conversion (summarized)", StepKind.PCR,
                ((duplex, InputRole.TEMPLATE), (handle_f, InputRole.STRUCTURAL),
                 (handle_r, InputRole.STRUCTURAL), (wta_pol, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY, reason="The WTA product precedes shearing and sequencing-adapter installation."))
    shear = graph.input("shear_condition", "Covaris S220 shearing to about 200 bp",
                        MaterialKind.CONDITION)
    fragments = graph.product("fragments", "sheared MATQ cDNA",
                              MaterialKind.MOLECULAR, state=amplified.state)
    graph.react("shear", "Shear the whole-transcriptome amplicons", StepKind.FRAGMENTATION,
                ((amplified, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (fragments,))
    repair = graph.input("end_repair_mix", "NEBNext End Repair Module", MaterialKind.CATALYST)
    repaired = graph.product("repaired", "end-repaired MATQ fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    graph.react("end_repair", "Repair the sheared fragment ends", StepKind.END_REPAIR,
                ((fragments, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    klenow = graph.input("klenow_exo", "Klenow exo-minus DNA polymerase", MaterialKind.CATALYST)
    dA_tailed = graph.product("da_tailed", "dA-tailed MATQ fragments",
                              MaterialKind.MOLECULAR, state=repaired.state)
    graph.react("da_tail", "Add dA to the repaired fragment ends", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (klenow, InputRole.CATALYST)), (dA_tailed,))
    adapter = graph.input("illumina_adapter", "TruSeq adapter (selected variant unreported)",
                          MaterialKind.ADAPTER, source=SOURCE + " (Online Methods, library preparation)")
    ligase = graph.input("adapter_ligase", "T4 DNA Quick Ligase", MaterialKind.CATALYST,
                         source=SOURCE + " (Online Methods, library preparation)")
    ligated = graph.product("ligated", "adapter-ligated sheared MATQ cDNA",
                            MaterialKind.MOLECULAR, state=fragments.state, inferred=True,
                            source=SOURCE + " (adapter identity unavailable)")
    graph.react("ligate_adapter", "Ligate a TruSeq adapter",
                StepKind.LIGATION,
                ((dA_tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The defining article names TruSeq adapters but not the selected variant, its ordered bases, or its ligation-ready end."))
    pcr_f = graph.input("library_forward", "TruSeq PCR primer mix, forward site", MaterialKind.OLIGO,
                        inferred=True, source=SOURCE + " (TruSeq PCR primer named; constituent oligos not printed)")
    pcr_r = graph.input("library_reverse", "TruSeq PCR primer mix, reverse site", MaterialKind.OLIGO,
                        inferred=True, source=SOURCE + " (TruSeq PCR primer named; constituent oligos not printed)")
    library_pol = graph.input("library_polymerase", "KAPA HiFi HotStart Ready Mix",
                              MaterialKind.CATALYST)
    library = graph.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY, inferred=True,
                            source=SOURCE + " (canonical compatible shell inferred)")
    graph.react("library_pcr", "Amplify the paired-end sequencing library", StepKind.PCR,
                ((ligated, InputRole.TEMPLATE), (pcr_f, InputRole.STRUCTURAL),
                 (pcr_r, InputRole.STRUCTURAL), (library_pol, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY, reason="The article names the TruSeq PCR primer but does not print its ordered bases or selected index variant."))
    graph.sequence(library, SequencingPlan(library=FINAL_LIBRARY,
                                           primers=tuple(SEQ_PRIMERS),
                                           read_lengths=READ_LENGTHS))
    graph.declare_reviewed_scope([
        "rna", "mixed_rt_primers", "reverse_transcriptase", "reverse_transcribe",
        "cleanup_enzymes", "cleanup", "tdt", "tail", "second_strand_primer",
        "deep_vent", "second_strand", "gat27_primer", "gat24_primer",
        "wta_polymerase", "wta", "shear_condition", "shear", "end_repair_mix", "end_repair",
        "klenow_exo", "da_tail", "illumina_adapter",
        "adapter_ligase", "ligate_adapter", "library_forward", "library_reverse",
        "library_polymerase", "library_pcr", "library"], SOURCE)
    return graph.require_publishable()
