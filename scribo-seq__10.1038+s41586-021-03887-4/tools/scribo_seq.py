"""VanInsberghe et al. single-cell ribosome profiling."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, MolecularState, Scene, Workflow, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "scRibo-seq — single-cell ribosome profiling"
NOTES = "01_scribo-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41586-021-03887-4">VanInsberghe et al. (2021)</a>.'
SUMMARY = "MNase generates ribosome-protected fragments in individual wells. One-pot end repair, 3′ and 5′ RNA-adapter ligation, reverse transcription and indexed PCR preserve a 10-nt UMI and a 10-nt cell index."
CAVEAT = ("Supplementary Table 1 (MOESM3 XLSX) prints the RNA/DNA adapters, RT primer and the indexed PCR primer families; "
          "the final construct below is assembled from those strings. Model correction (2026-10-17): earlier versions drew a "
          "TruSeq small-RNA RA5 Read 1 site, but the published 5′ RNA adapter and well-indexed forward primers carry the TruSeq "
          "Read 1 sequence. The sequencing primers themselves are not printed in Table 1: Read 1 and i5 are inferred to be the "
          "standard TruSeq primers, and the i7 read primer is inferred from the RPI 3′ region (the standard TruSeq index primer "
          "would pair over only 12 nt); the Methods name no custom sequencing primer. Reaction order and conditions were re-read in "
          "the published article Methods (open repository copy, 2026-10-17).")

# Supplementary Table 1 (MOESM3 XLSX, Sheet1) ordering strings, transcribed verbatim.
OMV630_3APP_ORDER = "/5rApp/TGGAATTCTCGGGTGCCAAGG/3ddC/"
OMV630_3APP = "TGGAATTCTCGGGTGCCAAGG"            # bases between the 5′ rApp and 3′ ddC
OMV572_RT = "GCCTTGGCACCCGAGAATTCCA"
OMV632_5P_ORDER = ("rCrUrArCrArCrGrArCrGrCrUrCrUrUrCrCrGrArUrCrUr(N1:25252525)"
                   "r(N1)r(N1)r(N1)r(N1)r(N1)r(N1)r(N1)r(N1)r(N1)")
OMV632_5P_FIXED = "CUACACGACGCUCUUCCGAUCU"          # fixed RNA before ten random rN
FWD_CBC_001 = "AATGATACGGCGACCACCGAGATCTACACGATTCACGACACACTCTTTCCCTACACGACGCTCTTCCG"  # well A1
RPI01 = "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"


def _assemble_sites():
    """Derive the library arms from the ordering strings; raise if they do not interlock."""
    if not FWD_CBC_001.startswith(il.P5):
        raise ValueError("forward primer does not begin with P5")
    fwd_tail = FWD_CBC_001[len(il.P5) + 10:]
    adapter_dna = OMV632_5P_FIXED.replace("U", "T")
    overlap = next(n for n in range(min(len(fwd_tail), len(adapter_dna)), 11, -1)
                   if fwd_tail.endswith(adapter_dna[:n]))
    read1_site = fwd_tail + adapter_dna[overlap:]
    if read1_site != il.TRUSEQ_READ1:
        raise ValueError("forward-primer tail plus 5′ adapter is not the TruSeq Read 1 site")
    if not RPI01.startswith(il.P7):
        raise ValueError("RPI primer does not begin with P7")
    i7_arm = revcomp(RPI01[len(il.P7) + 6:])
    if not i7_arm.startswith(OMV630_3APP):
        raise ValueError("RPI 3′ region does not copy the 3′ adapter")
    if revcomp(OMV572_RT) != OMV630_3APP + "C":
        raise ValueError("RT primer is not complementary to the 3′ adapter plus its 3′ ddC")
    return read1_site, i7_arm


READ1_SITE, SMALL_RNA_I1 = _assemble_sites()

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"),
    seg("ten-base cell index", "C" * 10, "cbc", placeholder=True,
        feature=feature("cell_index", "cell_barcode", "whitelist", whitelist="published well-specific forward primers")),
    seg("TruSeq Read 1 site (forward primer / 5′ RNA adapter)", READ1_SITE, "r1"),
    seg("ten-base UMI", "U" * 10, "umi", placeholder=True,
        feature=feature("umi", "umi", "random")),
    seg("ribosome-protected fragment", "R" * 35, placeholder=True, length_bp=35),
    seg("3-prime adapter / RPI Index 1 arm", SMALL_RNA_I1, "r2"),
    seg("six-base plate index reverse complement", "I" * 6, "cbc", placeholder=True,
        feature=feature("plate_i7", "sample_index", "whitelist", whitelist="published RPI primer set")),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="scRibo-seq library")
SEQ_PRIMERS = (
    sp.TRUSEQ["R1"],
    sp.custom("Index 1 (i7)", "small-RNA index primer (inferred from RPI 3′ region)", SMALL_RNA_I1,
              "derived from VanInsberghe et al. Supplementary Table 1 RPI primers",
              "Not printed in Table 1; the standard TruSeq Index 1 primer would pair over only 12 nt."),
    sp.TRUSEQ["I2"],
)
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS, required_roles=tuple(p.role for p in SEQ_PRIMERS)):
    raise ValueError("scRibo-seq final library: " + "; ".join(problems))
READ_LENGTHS = {"Read 1": 75, "Index 1 (i7)": 6, "Index 2 (i5)": 10}
SHOW_READ_CONTENT = True
FINAL_CAPTION = "Read 1 begins with the 10-nt UMI and then enters the 30–45-nt footprint. The six-cycle i7 read is the plate index; the ten-cycle i5 read is the well/cell index."
SEQUENCING_INTRO = "The defining paper specifies 75 cycles Read 1, six i7 cycles and ten i5 cycles on NextSeq 500. Identifier cycle spans are derived from the annotated final construct."

_initial = Scene()
_initial.strand("mRNA", [seg("translated mRNA", "R" * 48, placeholder=True)], label="cycloheximide-stalled translating mRNA")
_initial.mark("mRNA", "translated mRNA", "ribosome shields one footprint")
INITIAL_ROWS = tuple(_initial.rows())
INITIAL_NAME = "Single-cell translating mRNA"


def _three_prime_geometry():
    ready3=Scene()
    ready3.strand("substrate",[seg("RPF","R"*35,placeholder=True)],label="5′-phosphorylated, 3′-OH footprint")
    ready3.strand("adapter",[seg("3′ adapter",OMV630_3APP)],label="/5rApp/ adapter, /3ddC/ blocked")
    ready3.ligation_end("substrate","RPF",kind="single-strand",polarity="3-prime",owner="substrate")
    ready3.ligation_end("adapter","3′ adapter",kind="single-strand",polarity="5-prime",owner="adapter")
    return AdapterLigationGeometry.known("single-strand",ready3,statement="Supplementary Table 1 prints the 5′ pre-adenylated, 3′ dideoxy-C adapter; its 5′ end ligates to the RPF 3′-OH.")


def _five_prime_geometry():
    ready5=Scene()
    ready5.strand("substrate",[seg("RPF","R"*35,placeholder=True),seg("3′ adapter",OMV630_3APP)],label="3′-adapter-ligated footprint")
    ready5.strand("adapter",[seg("5′ adapter fixed RNA",OMV632_5P_FIXED),seg("ten-base UMI","U"*10,placeholder=True)],label="5′ RNA adapter with 3′ UMI")
    ready5.ligation_end("substrate","RPF",kind="single-strand",polarity="5-prime",owner="substrate")
    ready5.ligation_end("adapter","ten-base UMI",kind="single-strand",polarity="3-prime",owner="adapter")
    return AdapterLigationGeometry.known("single-strand",ready5,statement="Supplementary Table 1 prints the 5′ RNA adapter with ten random 3′ nucleotides; its 3′ end ligates to the phosphorylated RPF 5′ end.")


def workflow():
    footprint = Scene(); footprint.strand("RPF", [seg("ribosome-protected RNA", "R" * 35, placeholder=True, length_bp=35)], label="MNase-resistant RPF")
    footprint.mark("RPF", "ribosome-protected RNA", "30–45 nt retained after MNase")
    ligated = Scene(); ligated.strand("RNA", [
        seg("5-prime adapter fixed RNA", OMV632_5P_FIXED, "r1"),
        seg("ten-base UMI", "U" * 10, "umi", placeholder=True,
            feature=feature("umi", "umi", "random")),
        seg("RPF", "R" * 35, placeholder=True), seg("3-prime pre-adenylated adapter", OMV630_3APP, "r2")], label="doubly ligated footprint RNA")
    ligated.junction("RNA", "ten-base UMI", "RPF", "RNA ligation")
    ligated.junction("RNA", "RPF", "3-prime pre-adenylated adapter", "RNA ligation")
    ligated.labels("RNA")
    repaired = Scene(); repaired.strand("RPF", [seg("end-repaired RPF", "R" * 35, placeholder=True)], label="5′-phosphate / 3′-hydroxyl RPF")
    repaired.mark("RPF", "end-repaired RPF", "T4 PNK prepares both ends")
    indexed = Scene.duplex(list(FINAL_LIBRARY), label="indexed single-cell RPF library")
    # Size selection removes molecules, not sequence features from a retained molecule.
    # Keep the complete indexed construct so the final Workflow state still carries the
    # cell index, UMI and plate index installed in the preceding reaction.
    selected = Scene.duplex(list(FINAL_LIBRARY), label="175–185-bp PAGE-selected library")
    selected.mark("top", "ribosome-protected fragment", "30–40-nt footprint-size selection")
    w = Workflow(MolecularState(INITIAL_NAME, INITIAL_ROWS))
    w.react("Digest exposed RNA with MNase", footprint.rows(), name="Ribosome-protected fragment",
            note="MNase (37 °C, 30 min) leaves the ribosome-shielded footprint in each single-cell well; Proteinase K/EGTA/EDTA/GuSCN stop mix ends digestion.")
    w.react("Repair footprint ends", repaired.rows(), name="Ligatable RPF",
            note="T4 PNK prepares both ends for the two RNA ligations.")
    g3=_three_prime_geometry()
    three=Scene(); three.strand("RNA",[seg("RPF","R"*35,placeholder=True),seg("3′ adapter",OMV630_3APP)],label="3′-adapter-ligated RPF"); three.junction("RNA","RPF","3′ adapter","RNA ligation")
    w.ligate_adapter("Ligate the pre-adenylated 3′ RNA adapter",g3,three.rows(),name="3′-adapter-ligated RPF",note="T4 Rnl2(tr) KQ performs the ATP-independent single-strand ligation (4 °C, 18 h). OMV572_miRNA4_RT is then pre-annealed (65 °C → 37 °C → 25 °C) before the 5′ ligation.")
    g5=_five_prime_geometry()
    w.ligate_adapter("Ligate the UMI-bearing 5′ RNA adapter",g5,ligated.rows(),name="Adapter-ligated RPF",note="T4 Rnl1 joins the second single-strand RNA adapter (37 °C, 2 h) to the RT-primer-annealed molecule.")
    w.react("Reverse-transcribe and index each well", indexed.rows(), name="Indexed cDNA library",
            note="Maxima H-minus RT (50 °C, 1 h) extends the pre-annealed RT primer; well forward and RPI reverse primers install cell and plate indices in 10 Q5 PCR cycles.")
    w.react("Pool and footprint-size select", selected.rows(), name="Footprint-selected pool",
            note="PAGE retains the 175–185-bp library band, corresponding to a 30–40-nt insert.")
    return w


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    mrna = graph.input("translated_mrna", "cycloheximide-stalled translating mRNA",
                       MaterialKind.BIOLOGICAL, state=INITIAL_ROWS, primary=True)
    mnase = graph.input("mnase", "micrococcal nuclease", MaterialKind.CATALYST)
    rpf = graph.product("rpf", "ribosome-protected RNA fragment", MaterialKind.MOLECULAR)
    exposed = graph.product("digested_rna", "MNase-digested exposed RNA", MaterialKind.WASTE)
    graph.react("mnase_digest", "Digest exposed RNA around each ribosome",
                StepKind.DIGESTION, ((mrna, InputRole.TEMPLATE), (mnase, InputRole.CATALYST)),
                (rpf, exposed))
    graph.discard(exposed, "unprotected RNA is digested and does not enter the library")
    pnk = graph.input("pnk", "T4 polynucleotide kinase", MaterialKind.CATALYST)
    repaired = graph.product("repaired_rpf", "5′-phosphate/3′-hydroxyl RPF",
                             MaterialKind.MOLECULAR)
    graph.react("end_repair", "Repair both footprint ends", StepKind.END_REPAIR,
                ((rpf, InputRole.TEMPLATE), (pnk, InputRole.CATALYST)), (repaired,))
    adapter3 = graph.input("adapter_3p", "OMV630_miRNA4_3App pre-adenylated, 3′-ddC-blocked DNA adapter",
                           MaterialKind.ADAPTER,
                           state=Construct([seg("OMV630_miRNA4_3App", OMV630_3APP, "r2")],
                                           name="OMV630_miRNA4_3App " + OMV630_3APP_ORDER),
                           source=SOURCE + " Supplementary Table 1")
    rnl2 = graph.input("rnl2", "T4 RNA ligase 2 truncated KQ", MaterialKind.CATALYST)
    lig3 = graph.product("ligated_3p", "3′-adapter-ligated footprint", MaterialKind.MOLECULAR)
    graph.react("ligate_3p", "Ligate the pre-adenylated 3′ adapter", StepKind.LIGATION,
                ((repaired, InputRole.TEMPLATE), (adapter3, InputRole.STRUCTURAL),
                 (rnl2, InputRole.CATALYST)), (lig3,),
                geometry=_three_prime_geometry())
    rt_primer = graph.input("rt_primer", "OMV572_miRNA4_RT primer (complements the 3′ adapter including its ddC)",
                            MaterialKind.OLIGO,
                            state=Construct([seg("OMV572_miRNA4_RT", OMV572_RT)], name="OMV572_miRNA4_RT"),
                            source=SOURCE + " Supplementary Table 1")
    primed = graph.product("primed_3p", "RT-primer-annealed 3′-adapter-ligated footprint",
                           MaterialKind.MOLECULAR)
    graph.react("anneal_rt_primer", "Pre-anneal the RT primer to the 3′ adapter (65 °C 1 min, 37 °C 2 min, 25 °C 2 min)",
                StepKind.HYBRIDIZATION,
                ((lig3, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL)), (primed,))
    adapter5 = graph.input("adapter_5p", "OMV632_miRNA5_5A_10U RNA adapter with ten random 3′ rN (UMI)",
                           MaterialKind.ADAPTER,
                           state=Construct([seg("5′ adapter fixed RNA", OMV632_5P_FIXED, "r1"),
                                            seg("ten-base UMI", "U" * 10, "umi", placeholder=True,
                                                feature=feature("umi", "umi", "random"))],
                                           name="OMV632_miRNA5_5A_10U"),
                           source=SOURCE + " Supplementary Table 1")
    rnl1 = graph.input("rnl1", "T4 RNA ligase 1", MaterialKind.CATALYST)
    doubly = graph.product("doubly_ligated", "doubly adapter-ligated footprint",
                           MaterialKind.MOLECULAR)
    graph.react("ligate_5p", "Ligate the UMI-bearing 5′ adapter", StepKind.LIGATION,
                ((primed, InputRole.TEMPLATE), (adapter5, InputRole.STRUCTURAL),
                 (rnl1, InputRole.CATALYST)), (doubly,),
                geometry=_five_prime_geometry())
    rt = graph.input("reverse_transcriptase", "Maxima H-minus reverse transcriptase",
                     MaterialKind.CATALYST)
    cdna = graph.product("cdna", "footprint cDNA", MaterialKind.MOLECULAR)
    graph.react("reverse_transcribe", "Reverse-transcribe from the pre-annealed RT primer",
                StepKind.REVERSE_TRANSCRIPTION,
                ((doubly, InputRole.TEMPLATE), (doubly, InputRole.SELF_PRIMER),
                 (rt, InputRole.CATALYST)), (cdna,))
    fwd = graph.input("well_forward_primer", "miRv6-PCR_F-cbc well-indexed forward primer family (P5–10-nt cell index–TruSeq Read 1 tail)",
                      MaterialKind.OLIGO, source=SOURCE + " Supplementary Table 1")
    rev = graph.input("rpi_reverse_primer", "RPI plate-indexed reverse primer family (P7–6-nt index–3′-adapter complement)",
                      MaterialKind.OLIGO, source=SOURCE + " Supplementary Table 1")
    pol = graph.input("pcr_polymerase", "Q5 Hot Start High-Fidelity 2× Master Mix (10 cycles)", MaterialKind.CATALYST)
    indexed = graph.product("indexed_library", "indexed scRibo-seq library",
                            MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react("index_pcr", "PCR with well and plate indexing primers", StepKind.PCR,
                ((cdna, InputRole.TEMPLATE), (fwd, InputRole.STRUCTURAL),
                 (rev, InputRole.STRUCTURAL), (pol, InputRole.CATALYST)), (indexed,),
                geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY, reason="Each well and plate uses a different member of the published forward and RPI primer families; no single member is drawn."))
    page = graph.input("page_selection", "175–185-bp PAGE selection", MaterialKind.SELECTION)
    selected = graph.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY,
                             state=FINAL_LIBRARY)
    offsize = graph.product("offsize_library", "off-size library molecules", MaterialKind.WASTE)
    graph.react("size_select", "Select the footprint-sized library band", StepKind.SELECTION,
                ((indexed, InputRole.TEMPLATE), (page, InputRole.SELECTION)),
                (selected, offsize))
    graph.discard(offsize, "molecules outside the 175–185-bp PAGE band are excluded")
    graph.sequence(selected, SequencingPlan(library=FINAL_LIBRARY,
                                            primers=tuple(SEQ_PRIMERS),
                                            read_lengths=READ_LENGTHS))
    graph.declare_reviewed_scope((
        "translated_mrna", "mnase", "mnase_digest", "pnk", "end_repair", "adapter_3p",
        "rnl2", "ligate_3p", "rt_primer", "anneal_rt_primer", "adapter_5p", "rnl1", "ligate_5p",
        "reverse_transcriptase", "reverse_transcribe", "well_forward_primer",
        "rpi_reverse_primer", "pcr_polymerase", "index_pcr", "page_selection",
        "size_select", "library"), SOURCE)
    return graph.require_publishable()
