"""Original 2013 PRO-seq TruSeq option; later RPI sequence scope explicit."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_model import role
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "PRO-seq — precision nuclear run-on sequencing"
NOTES = "01_pro-seq.html"
SOURCE = 'Defining paper: <a href="https://doi.org/10.1126/science.1229386">Kwak et al. (2013)</a>; detailed protocol: <a href="https://doi.org/10.1038/nprot.2016.086">Mahat et al. (2016)</a>.'
SUMMARY = "Engaged polymerases extend nascent RNA by one or a few biotin-NTPs. Biotin purification and sequential RNA-adapter ligations preserve the nucleotide-resolution 3′ end for sequencing."
CAVEAT = ('The original 2013 SOM describes GX2/GX1 and TruSeq RP1/RPI alternatives. This page represents the explicitly named RP1/TruSeq option, not every original experiment. '
          'The original 2013 supplement covers both PRO-seq and PRO-cap; only the PRO-seq branch is represented. RNA U is written T in computational pairing views. '
          'INFERRED — the endpoint uses the later 2016 Table 1 RPI-n ordering family as a same-name comparator; its identity to a particular2013 ordered RPI oligo is unassigned. '
          'Exact loaded run primers are not assigned from the PCR sites. PRO-cap omits fragmentation, substitutes standard adapters, and uses AP/TAP without PNK; that branch is not imported here.')

# Original 2013 SOM prints VRA3, TruSeq VRA5 option and RP1. Later2016 Table 1 supplies the explicit RPI-n comparator; RNA U uses T for pairing arithmetic.
VRA3 = "GATCGTCGGACTGTAGAACTCTGAAC"
VRA5 = "CCTTGGCACCCGAGAATTCCA"
RP1 = "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA"
RPI_BIND = "GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
RPI = il.P7 + "N" * 6 + RPI_BIND
SMALL_RNA_R1 = revcomp(VRA3)

FINAL_LIBRARY = Construct([
    Segment("P5", il.P5, "p5", inferred=True),
    Segment("TruSeq Small RNA Read 1 site", SMALL_RNA_R1, "r1", inferred=True),
    role('reverse-complement nascent RNA'),
    Segment("VRA5-derived Read 2 side", revcomp(VRA5), "r2", inferred=True),
    Segment("Index 1 continuation", revcomp(RPI_BIND[:-len(VRA5)]), "r2", inferred=True),
    Segment("six-base i7 index reverse complement", "N" * 6, "cbc", placeholder=True,
        feature=feature("sample_i7", "sample_index", "whitelist", whitelist="published RPI-n primer set")),
    Segment("P7 reverse complement", il.P7_RC, "p7", inferred=True),
], name="INFERRED — RP1/RPI comparator PRO-seq endpoint")
# Derived binding sites are not evidence of physically loaded sequencing oligos.
SEQ_PRIMERS = ()
FINAL_CAPTION = 'INFERRED — original 2013 RP1/TruSeq adapter option with later 2016 RPI-n sequence comparator; selected historical index/loaded run primers unassigned. Read 1 starts at the reverse-complement nascent RNA 3-prime boundary.'
SEQUENCING_ENDING = 'Original 2013 SOM reports GenomeAnalyzer 2 or HiSeq 2000 sequencing; read sequences are reverse-complemented to recover nascent RNA sense orientation. Exact acquired cycles/loaded run primers are unassigned here.'
SEQUENCING_UNAVAILABLE = 'No complete source-specific loaded run-primer set is declared; later comparator PCR-site architecture does not establish the actual2013 loaded oligos.'
READ_SPECS = (
    sp.read_spec('Read 1', 'unassigned', 'reverse complement of nascent RNA, starting from its informative 3-prime end', 'include'),
    sp.read_spec('Index 1 (i7)', '6 in later 2016 comparator', 'RPI-n sample index; historical selected member unassigned', 'demultiplex'),
)

def oligos():
    rna = lambda sequence: ''.join('r' + base for base in sequence.replace('T', 'U'))
    return (
        'Original 2013 reverse 3-prime RNA adapter: 5′-p-' + rna(VRA3) + '-/3InvdT/-3′',
        'Original 2013 TruSeq reverse 5-prime RNA option: 5′-' + rna(VRA5) + '-3′',
        'Original 2013 RP1: 5′-' + RP1 + '-3′',
        'Later 2016 RPI-n comparator (not assigned to every 2013 library): 5′-' + RPI + '-3′',
    )

def _runon_scene():
    sc = Scene(); sc.strand("nascent RNA", [
        role('nascent RNA body'),
        Segment("biotin-NTP run-on base", "N", "r1", placeholder=True)],
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
        statement="Original 2013 SOM prints the reverse3-prime RNA adapter as 5′ phosphorylated and 3′ inverted-dT blocked; its 5′ end ligates to the informative nascent-RNA 3′-OH.")


def _three_prime_product_scene():
    sc = Scene(); sc.strand("product", [
        role('biotin-labelled nascent RNA'),
        Segment("VRA3", VRA3, "r1")], label="3′-adapter-ligated nascent RNA",
        mod3="/3InvdT/")
    sc.junction("product", "biotin-labelled nascent RNA", "VRA3",
                "3′ RNA-adapter ligation")
    sc.labels("product")
    return sc


def _repaired_five_prime_scene():
    sc = Scene(); sc.strand("RNA", [
        role('selected nascent RNA'),
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
        role('selected nascent RNA'),
        Segment("VRA3", VRA3, "r1")], label="dual-adapter PRO-seq RNA",
        mod3="/3InvdT/")
    sc.junction("RNA", "VRA5", "selected nascent RNA", "5′ RNA-adapter ligation")
    sc.junction("RNA", "selected nascent RNA", "VRA3", "3′ RNA-adapter ligation")
    sc.labels("RNA")
    return sc


def _rt_scene():
    sc = Scene(); sc.strand("RNA", [Segment("VRA5", VRA5, "r2"),
        role('selected nascent RNA'),
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
                      Segment("copied VRA3 prefix", revcomp(VRA3[:5]), "r1"),
                      role('reverse-complement nascent RNA'),
                      Segment("VRA5 reverse complement", revcomp(VRA5), "r2")],
                     name="PRO-seq first-strand cDNA")


def _pre_pcr_scene():
    sc = Scene(); sc.strand("cDNA", list(_pre_pcr_template()), label="PRO-seq first-strand cDNA"); return sc


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
    extraction = graph.input('rna_extraction', 'original Trizol extraction/precipitation then P30 desalting around hydrolysis', MaterialKind.CONDITION)
    fragments = graph.product('fragmented_rna', 'base-hydrolyzed RNA fragment pool', MaterialKind.MOLECULAR)
    graph.react('fragment_rna', 'Extract RNA, base-hydrolyze and desalt (PRO-seq only)', StepKind.FRAGMENTATION,
                ((labelled, InputRole.TEMPLATE), (extraction, InputRole.CONDITION)), (fragments,))
    streptavidin=graph.input('streptavidin', 'streptavidin capture reagent', MaterialKind.SELECTION)
    first = graph.product('first_capture', 'first-capture biotin-labelled fragments', MaterialKind.MOLECULAR, state=_runon_scene())
    first_wash=graph.product('first_wash', 'uncaptured RNA fragments', MaterialKind.WASTE)
    graph.react('capture_first', 'First streptavidin enrichment and extraction before adapter ligation', StepKind.AFFINITY_CAPTURE,
                ((fragments,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION)),(first,first_wash))
    graph.discard(first_wash,'Original first biotin purification removes uncaptured fragments.')
    vra3=graph.input("vra3","5-prime-phosphorylated 3-prime-blocked VRA3 adapter",
                     MaterialKind.ADAPTER,state=_vra3_geometry().scenes[0])
    ligase3=graph.input("ligase3","T4 RNA ligase I (source3-prime ligation)",MaterialKind.CATALYST)
    ligated3=graph.product("ligated3","VRA3-ligated biotin RNA",MaterialKind.MOLECULAR,
                           state=tuple(_three_prime_product_scene().rows()))
    graph.react("ligate_vra3","Ligate VRA3 to the informative RNA 3-prime end",
                StepKind.LIGATION,
                ((first,InputRole.TEMPLATE),(vra3,InputRole.STRUCTURAL),
                 (ligase3,InputRole.CATALYST)),(ligated3,),geometry=_vra3_geometry())
    captured=graph.product("captured","captured biotin run-on RNA",MaterialKind.MOLECULAR)
    uncaptured=graph.product("uncaptured","unlabelled RNA",MaterialKind.WASTE)
    graph.react("capture","Capture biotin-labelled run-on products",
                StepKind.AFFINITY_CAPTURE,
                ((ligated3,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION)),
                (captured,uncaptured))
    graph.discard(uncaptured,"unlabelled RNA is removed during streptavidin purification")
    repair=graph.input("end_repair_enzymes","TAP then PNK for PRO-seq; PRO-cap AP/TAP is a separate branch",
                       MaterialKind.CATALYST)
    repaired=graph.product("repaired","5-prime-monophosphorylated selected RNA",
                           MaterialKind.MOLECULAR,
                           state=tuple(_repaired_five_prime_scene().rows()))
    graph.react("repair_five_prime","Decap and repair the selected RNA 5-prime end",
                StepKind.END_REPAIR,
                ((captured,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(repaired,))
    vra5=graph.input("vra5","VRA5 RNA adapter",MaterialKind.ADAPTER,
                     state=_vra5_geometry().scenes[0])
    ligase5=graph.input("ligase5","5-prime RNA-adapter ligase",MaterialKind.CATALYST)
    dual=graph.product("dual_adapter","dual-adapter PRO-seq RNA",MaterialKind.MOLECULAR,
                       state=tuple(_dual_adapter_scene().rows()))
    graph.react("ligate_vra5","Ligate VRA5 to the repaired RNA 5-prime end",
                StepKind.LIGATION,
                ((repaired,InputRole.TEMPLATE),(vra5,InputRole.STRUCTURAL),
                 (ligase5,InputRole.CATALYST)),(dual,),geometry=_vra5_geometry())
    third = graph.product('third_capture', 'third-enriched dual-adapter biotin RNA', MaterialKind.MOLECULAR, state=_dual_adapter_scene())
    third_wash = graph.product('third_wash', 'third-capture unbound material', MaterialKind.WASTE)
    graph.react('capture_third', 'Third streptavidin enrichment after 5-prime ligation', StepKind.AFFINITY_CAPTURE,
                ((dual,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION)),(third,third_wash))
    graph.discard(third_wash,'Original SOM explicitly repeats streptavidin binding/extraction after 5-prime ligation.')
    rp1=graph.input("rp1","RP1 reverse-transcription/PCR primer",MaterialKind.OLIGO, state=Construct([Segment("RP1", RP1)]))
    rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
    cdna=graph.product("cdna","RP1-primed first-strand cDNA",MaterialKind.MOLECULAR,
                       state=_pre_pcr_template())
    graph.react("reverse_transcribe","Reverse-transcribe from RP1 annealed to VRA3",
                StepKind.REVERSE_TRANSCRIPTION,
                ((third,InputRole.TEMPLATE),(rp1,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(cdna,))
    rpi=graph.input("rpi","later 2016 RPI-n comparator; original selected oligo unassigned",MaterialKind.OLIGO, state=Construct([Segment("RPI 2016 comparator",RPI,placeholder=True,inferred=True)]))
    pcr_polymerase=graph.input("pcr_polymerase","Phusion polymerase (original 2013 SOM)",
                               MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY)
    geometry=PCRGeometry.exact(
        _pre_pcr_template(),PCRPrimer("RP1",RP1),PCRPrimer("RPI-n",RPI),
        statement="Original 2013 SOM prints RP1; later 2016 Table 1 prints the RPI-n comparator family and their opposing sites.")
    graph.react("library_pcr","INFERRED — PCR with source RP1 and later RPI-n comparator",StepKind.PCR,
                ((cdna,InputRole.TEMPLATE),(rp1,InputRole.STRUCTURAL),
                 (rpi,InputRole.STRUCTURAL),(pcr_polymerase,InputRole.CATALYST)),
                (library,),geometry=geometry)
    page = graph.input('page_selection', 'original PAGE library size purification', MaterialKind.SELECTION)
    selected = graph.product('page_library', FINAL_LIBRARY.name, MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True)
    excluded = graph.product('page_excluded', 'material outside excised library band', MaterialKind.WASTE)
    graph.react('purify_library', 'PAGE-purify the amplified library', StepKind.CLEANUP,
                ((library,InputRole.TEMPLATE),(page,InputRole.SELECTION)),(selected,excluded))
    graph.discard(excluded,'Original SOM PAGE-purifies amplified library; exact excluded species unassigned.')
    graph.sequence(selected,SequencingPlan(library=FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE))
    graph.declare_reviewed_scope((
        "nascent_rna","biotin_ntp","engaged_polymerase","run_on","vra3",
        "ligase3","ligate_vra3","streptavidin","capture","end_repair_enzymes",
        "repair_five_prime","vra5","ligase5","ligate_vra5","rp1",
        "reverse_transcriptase","reverse_transcribe","rpi","pcr_polymerase",
        "library_pcr","library", "rna_extraction","fragment_rna","capture_first","first_capture", "capture_third","third_capture","page_selection","page_library","purify_library"),SOURCE)
    return graph.require_publishable()
