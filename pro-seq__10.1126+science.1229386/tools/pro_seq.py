"""Mahat et al. PRO-seq library, including its TruSeq Small RNA run geometry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Row, Scene, Segment, Workflow, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "PRO-seq — precision nuclear run-on sequencing"
NOTES = "01_pro-seq.html"
SOURCE = 'Defining paper: <a href="https://doi.org/10.1126/science.1229386">Kwak et al. (2013)</a>; detailed protocol: <a href="https://doi.org/10.1038/nprot.2016.086">Mahat et al. (2016)</a>.'
SUMMARY = "Engaged polymerases extend nascent RNA by one or a few biotin-NTPs. Biotin purification and sequential RNA-adapter ligations preserve the nucleotide-resolution 3′ end for sequencing."
CAVEAT = "The six-base sample index varies between RPI-n primers; the model leaves only those six bases variable."

# Mahat et al. 2016, Table 1. RNA U is represented as DNA T in the final cDNA.
VRA3 = "GATCGTCGGACTGTAGAACTCTGAAC"
VRA5 = "CCTTGGCACCCGAGAATTCCA"
RP1 = "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA"
RPI_BIND = "GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
RPI = il.P7 + "N" * 6 + RPI_BIND
SMALL_RNA_R1 = "GTTCAGAGTTCTACAGTCCGACGATC"
SMALL_RNA_I1 = "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC"

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"),
    seg("TruSeq Small RNA Read 1 site", SMALL_RNA_R1, "r1"),
    seg("reverse-complement nascent RNA", "X" * 34, placeholder=True),
    seg("VRA5-derived Read 2 side", "TGGAATTCTCGGGTGCCAAGG", "r2"),
    seg("Index 1 continuation", "AACTCCAGTCAC", "cbc"),
    seg("six-base i7 index reverse complement", "I" * 6, "cbc", placeholder=True,
        feature=feature("sample_i7", "sample_index", "whitelist", whitelist="published RPI-n primer set")),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="single-index PRO-seq library")
SEQ_PRIMERS = (
    sp.custom("Read 1", "TruSeq Small RNA Read 1", SMALL_RNA_R1,
              "Mahat et al. 2016 Table 1; TruSeq Small RNA geometry"),
    sp.custom("Index 1 (i7)", "TruSeq Small RNA Index 1", SMALL_RNA_I1,
              "Mahat et al. 2016 Table 1; RPI-n"),
)
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                         required_roles=("Read 1", "Index 1 (i7)")):
    raise ValueError("PRO-seq final library: " + "; ".join(problems))

FINAL_CAPTION = "Mahat et al. single-end, six-base single-index library. Read 1 sequences the reverse complement of nascent RNA from its informative 3′ end."
SEQUENCING_INTRO = "The TruSeq Small RNA Read 1 primer enters the reverse-complement nascent-RNA insert; the custom small-RNA index primer reads one six-base i7 index. There is no i5 read."

def _runon_scene():
    sc = Scene(); sc.strand("nascent RNA", [
        Segment("nascent RNA body", "X" * 42, placeholder=True),
        Segment("biotin-NTP run-on base", "N", "umi", placeholder=True)],
        label="biotin-labelled nascent RNA")
    sc.mark("nascent RNA", "biotin-NTP run-on base",
            "informative polymerase-active-site 3′ end")
    sc.labels("nascent RNA")
    return sc


def _vra3_geometry():
    adapter = Scene(); adapter.strand("VRA3", [Segment("VRA3", VRA3, "r1")],
                                      label="VRA3", mod5="/5Phos/",
                                      mod3="/3InvdT/")
    adapter.ligation_end("VRA3", "VRA3", kind="single-strand",
                         polarity="5-prime", owner="adapter",
                         text="5′-phosphorylated VRA3 ligating end")
    substrate = _runon_scene()
    substrate.ligation_end("nascent RNA", "biotin-NTP run-on base",
                           kind="single-strand", polarity="3-prime",
                           owner="substrate", text="nascent-RNA 3′-OH")
    return AdapterLigationGeometry.known(
        "single-strand", adapter, substrate,
        statement="Mahat et al. Table 1 prints VRA3 as 5′ phosphorylated and 3′ inverted-dT blocked; its 5′ end ligates to the informative nascent-RNA 3′-OH.")


def _three_prime_product_scene():
    sc = Scene(); sc.strand("product", [
        Segment("biotin-labelled nascent RNA", "X" * 42, placeholder=True),
        Segment("VRA3", VRA3, "r1")], label="3′-adapter-ligated nascent RNA",
        mod3="/3InvdT/")
    sc.junction("product", "biotin-labelled nascent RNA", "VRA3",
                "3′ RNA-adapter ligation")
    sc.labels("product")
    return sc


def _repaired_five_prime_scene():
    sc = Scene(); sc.strand("RNA", [
        Segment("selected nascent RNA", "X" * 42, placeholder=True),
        Segment("VRA3", VRA3, "r1")], label="decapped / 5′-repaired RNA",
        mod5="p", mod3="/3InvdT/")
    sc.mark("RNA", "selected nascent RNA", "5′ monophosphate after end repair")
    sc.labels("RNA")
    return sc


def _vra5_geometry():
    adapter = Scene(); adapter.strand("VRA5", [Segment("VRA5", VRA5, "r2")],
                                      label="VRA5")
    adapter.ligation_end("VRA5", "VRA5", kind="single-strand",
                         polarity="3-prime", owner="adapter",
                         text="VRA5 3′-OH ligating end")
    substrate = _repaired_five_prime_scene()
    substrate.ligation_end("RNA", "selected nascent RNA", kind="single-strand",
                           polarity="5-prime", owner="substrate",
                           text="repaired RNA 5′ monophosphate")
    return AdapterLigationGeometry.known(
        "single-strand", adapter, substrate,
        statement="The detailed protocol prints VRA5 and repairs the selected RNA 5′ end before T4 RNA-ligase-mediated 5′ adapter ligation.")


def _dual_adapter_scene():
    sc = Scene(); sc.strand("RNA", [Segment("VRA5", VRA5, "r2"),
        Segment("selected nascent RNA", "X" * 42, placeholder=True),
        Segment("VRA3", VRA3, "r1")], label="dual-adapter PRO-seq RNA",
        mod3="/3InvdT/")
    sc.junction("RNA", "VRA5", "selected nascent RNA", "5′ RNA-adapter ligation")
    sc.junction("RNA", "selected nascent RNA", "VRA3", "3′ RNA-adapter ligation")
    sc.labels("RNA")
    return sc


def _rt_scene():
    sc = Scene(); sc.strand("RNA", [Segment("VRA5", VRA5, "r2"),
        Segment("selected nascent RNA", "X" * 42, placeholder=True),
        Segment("VRA3", VRA3, "r1")], label="dual-adapter PRO-seq RNA")
    sc.anneal("RP1", [Segment("RP1 P5 tail", RP1[:-21], "p5"),
                       Segment("RP1 annealing suffix", RP1[-21:], "r1")],
              to="RNA", pair=("RP1 annealing suffix", "VRA3"), label="RP1",
              shift=5, unpaired=("RP1 P5 tail",))
    sc.arrow("RP1", "reverse transcriptase")
    sc.labels("RP1")
    return sc


def _pre_pcr_template():
    return Construct([Segment("RP1", RP1, "p5"),
                      Segment("reverse-complement nascent RNA", "X" * 42,
                              placeholder=True),
                      Segment("VRA5 reverse complement", revcomp(VRA5), "r2")],
                     name="PRO-seq first-strand cDNA")


def _pre_pcr_scene():
    return Scene.duplex(list(_pre_pcr_template()), label="PRO-seq first-strand cDNA")


def workflow():
    w = Workflow(MolecularState("Engaged polymerase and nascent RNA",
                                tuple(_runon_scene().rows())))
    w.react("Run on briefly with biotin-NTPs", _runon_scene().rows(),
            name="biotin-labelled nascent RNA",
            note="The last incorporated biotin nucleotide marks the active polymerase position and stalls further extension.")
    w.ligate_adapter("Ligate VRA3 to the nascent-RNA 3′ end", _vra3_geometry(),
                     _three_prime_product_scene().rows(),
                     name="3′-adapter-ligated biotin RNA",
                     note="The covalent junction is the information-bearing RNA 3′ boundary; inverted dT blocks VRA3 self-ligation.")
    w.react("Streptavidin-purify, decap and repair the RNA 5′ end",
            _repaired_five_prime_scene().rows(), name="5′-monophosphorylated selected RNA",
            note="Repeated biotin capture retains run-on products while end processing makes the opposite end ligatable.")
    w.ligate_adapter("Ligate VRA5 to the repaired RNA 5′ end", _vra5_geometry(),
                     _dual_adapter_scene().rows(), name="dual-adapter PRO-seq RNA",
                     note="VRA5's 3′ end joins the repaired RNA 5′ phosphate; both covalent RNA-adapter junctions remain marked.")
    w.react("Reverse-transcribe from RP1 annealed to VRA3", _rt_scene().rows(),
            name="RP1-primed first-strand cDNA",
            note="RP1 contributes P5 as a non-annealing 5′ tail and uses its 21-nt suffix to prime VRA3.")
    geometry = PCRGeometry.exact(
        _pre_pcr_template(), PCRPrimer("RP1", RP1), PCRPrimer("RPI-n", RPI),
        statement="Mahat et al. Table 1 prints RP1 and the RPI-n family; RP1 lands at the VRA3-derived end and the RPI-n 3′ suffix lands on VRA5-derived cDNA.")
    w.pcr("PCR-add the six-base i7 index and P7", geometry,
          Scene.duplex(list(FINAL_LIBRARY), label="indexed PRO-seq library").rows(),
          name="single-index PRO-seq library",
          note="The RPI-n member supplies one six-base i7 index; the protocol has no i5 read.")
    return w


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    nascent=graph.input("nascent_rna","engaged polymerase with nascent RNA",
                        MaterialKind.MOLECULAR,primary=True)
    biotin_ntp=graph.input("biotin_ntp","biotin-labelled run-on NTPs",
                           MaterialKind.MOLECULAR)
    polymerase=graph.input("engaged_polymerase","engaged RNA polymerase",
                           MaterialKind.CATALYST)
    labelled=graph.product("labelled","biotin-labelled nascent RNA",
                           MaterialKind.MOLECULAR,state=_runon_scene())
    graph.react("run_on","Run on briefly with biotin-NTPs",StepKind.RUN_ON,
                ((nascent,InputRole.TEMPLATE),(biotin_ntp,InputRole.STRUCTURAL),
                 (polymerase,InputRole.CATALYST)),(labelled,))
    vra3=graph.input("vra3","5-prime-phosphorylated 3-prime-blocked VRA3 adapter",
                     MaterialKind.ADAPTER,state=tuple(_vra3_geometry().rows()))
    ligase3=graph.input("ligase3","3-prime RNA-adapter ligase",MaterialKind.CATALYST)
    ligated3=graph.product("ligated3","VRA3-ligated biotin RNA",MaterialKind.MOLECULAR,
                           state=tuple(_three_prime_product_scene().rows()))
    graph.react("ligate_vra3","Ligate VRA3 to the informative RNA 3-prime end",
                StepKind.LIGATION,
                ((labelled,InputRole.TEMPLATE),(vra3,InputRole.STRUCTURAL),
                 (ligase3,InputRole.CATALYST)),(ligated3,),geometry=_vra3_geometry())
    streptavidin=graph.input("streptavidin","streptavidin capture reagent",
                            MaterialKind.SELECTION)
    captured=graph.product("captured","captured biotin run-on RNA",MaterialKind.MOLECULAR)
    uncaptured=graph.product("uncaptured","unlabelled RNA",MaterialKind.WASTE)
    graph.react("capture","Capture biotin-labelled run-on products",
                StepKind.AFFINITY_CAPTURE,
                ((ligated3,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION)),
                (captured,uncaptured))
    graph.discard(uncaptured,"unlabelled RNA is removed during streptavidin purification")
    repair=graph.input("end_repair_enzymes","RNA decapping and 5-prime end-repair enzymes",
                       MaterialKind.CATALYST)
    repaired=graph.product("repaired","5-prime-monophosphorylated selected RNA",
                           MaterialKind.MOLECULAR,
                           state=tuple(_repaired_five_prime_scene().rows()))
    graph.react("repair_five_prime","Decap and repair the selected RNA 5-prime end",
                StepKind.END_REPAIR,
                ((captured,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(repaired,))
    vra5=graph.input("vra5","VRA5 RNA adapter",MaterialKind.ADAPTER,
                     state=tuple(_vra5_geometry().rows()))
    ligase5=graph.input("ligase5","5-prime RNA-adapter ligase",MaterialKind.CATALYST)
    dual=graph.product("dual_adapter","dual-adapter PRO-seq RNA",MaterialKind.MOLECULAR,
                       state=tuple(_dual_adapter_scene().rows()))
    graph.react("ligate_vra5","Ligate VRA5 to the repaired RNA 5-prime end",
                StepKind.LIGATION,
                ((repaired,InputRole.TEMPLATE),(vra5,InputRole.STRUCTURAL),
                 (ligase5,InputRole.CATALYST)),(dual,),geometry=_vra5_geometry())
    rp1=graph.input("rp1","RP1 reverse-transcription/PCR primer",MaterialKind.OLIGO)
    rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
    cdna=graph.product("cdna","RP1-primed first-strand cDNA",MaterialKind.MOLECULAR,
                       state=_pre_pcr_template())
    graph.react("reverse_transcribe","Reverse-transcribe from RP1 annealed to VRA3",
                StepKind.REVERSE_TRANSCRIPTION,
                ((dual,InputRole.TEMPLATE),(rp1,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(cdna,))
    rpi=graph.input("rpi","RPI-n indexed PCR primer",MaterialKind.OLIGO)
    pcr_polymerase=graph.input("pcr_polymerase","library PCR polymerase",
                               MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY)
    geometry=PCRGeometry.exact(
        _pre_pcr_template(),PCRPrimer("RP1",RP1),PCRPrimer("RPI-n",RPI),
        statement="Mahat et al. Table 1 prints RP1 and the RPI-n family and their opposing sites.")
    graph.react("library_pcr","Add the six-base i7 index and P7",StepKind.PCR,
                ((cdna,InputRole.TEMPLATE),(rp1,InputRole.STRUCTURAL),
                 (rpi,InputRole.STRUCTURAL),(pcr_polymerase,InputRole.CATALYST)),
                (library,),geometry=geometry)
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "nascent_rna","biotin_ntp","engaged_polymerase","run_on","vra3",
        "ligase3","ligate_vra3","streptavidin","capture","end_repair_enzymes",
        "repair_five_prime","vra5","ligase5","ligate_vra5","rp1",
        "reverse_transcriptase","reverse_transcribe","rpi","pcr_polymerase",
        "library_pcr","library"),SOURCE)
    return graph.require_publishable()
