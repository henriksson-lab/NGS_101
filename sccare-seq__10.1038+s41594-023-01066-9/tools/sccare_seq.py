"""scCARE-seq molecular model."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Row, Scene, Segment, feature, oligo, revcomp
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry, PCRPrimer
from protocol_extensions import molecule
from protocol_graph import Disposition, Terminal
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            StepKind, Terminal)

RT_HANDLE = "ACGACGCTCTTCCGATCT"
RT_BARCODES = (
    "CTTAGGAC", "GTTGGAGC", "TCTTACGC", "GATGGTCA",
    "GAGAATCC", "CGCAACGA", "CGGACCTA", "GAATTAAC",
)
UMI = "N" * 8
POLYT = "T" * 30 + "VN"

# Representative members of the three complete index-primer sets printed in
# Supplementary Table 1. Their invariant arms come from the shared Illumina model.
N501_INDEX = "TAGATCGC"
N701_INDEX = "TCGCCTTA"
T501_INDEX = "TATAGCCT"
N501 = nx.n5xx_primer(N501_INDEX)
N701 = nx.n7xx_primer(N701_INDEX)
T501 = il.nebnext_i5_primer(T501_INDEX)
COMMON_PRIMER = nx.S7 + nx.ME[:13]
DNA_PRIMER = nx.S5 + nx.ME[:13]
RNA_PRIMER = "CCCTAC" + RT_HANDLE  # A5-PreA, Supplementary Table 1 (corrected 2026-10-17: was missing the C)
TN5_N5 = nx.ADAPTOR_S5
TN5_N7 = nx.ADAPTOR_S7
TN5_REV = nx.ME_RC
CARRIER = "TCAGGTTTTCCTGAA"

# Article Methods ("scCARE-seq experimental procedure") and GEO GSE211395: RT, MboI
# digestion, ligation and second-strand synthesis happen in intact nuclei *before*
# single nuclei are sorted into wells, so the RT barcode is shared by all
# nuclei of one RT reaction; per-cell identity comes from the well's i5/i7 index pair.
CELL = feature(
    "sccare_rt_barcode", "sample_index", "whitelist",
    whitelist="eight RT barcodes printed in Supplementary Table 1 (one per RT reaction; RT precedes single-nucleus sorting per article Methods and GEO GSE211395)",
)
UMI_FEATURE = feature("sccare_umi", "umi", "random")
I5 = feature("sccare_i5", "sample_index", "whitelist",
             whitelist="N501–N516 or T501–T516 in Supplementary Table 1")
I7 = feature("sccare_i7", "sample_index", "whitelist",
             whitelist="N701–N724 in Supplementary Table 1")


def _dna_library() -> Construct:
    """N501/N701 representative, assembled only from source-printed parts."""
    return Construct([
        Segment("P5", il.P5, "p5"),
        Segment("i5 (N501)", N501_INDEX, "cbc", feature=I5),
        Segment("S5", nx.S5, "s5"),
        Segment("left mosaic end", nx.ME, "me"),
        Segment("contact-junction DNA", "X" * 46, placeholder=True),
        Segment("right mosaic end reverse complement", nx.ME_RC, "me"),
        Segment("S7 reverse complement", nx.S7_RC, "s7"),
        Segment("i7 (N701) reverse complement", revcomp(N701_INDEX), "cbc", feature=I7),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="scCARE-seq DNA-contact library (N501/N701 representative)")


def _rna_library() -> Construct:
    """Hybrid TruSeq-i5/Nextera-i7 RNA branch shown in Extended Data Fig. 1."""
    return Construct([
        Segment("P5", il.P5, "p5"),
        Segment("i5 (T501)", T501_INDEX, "cbc", feature=I5),
        Segment("TruSeq Read 1 arm", il.TRUSEQ_READ1, "r1"),
        Segment("UMI", UMI, "umi", placeholder=True, feature=UMI_FEATURE),
        Segment("RT barcode", RT_BARCODES[0], "cbc", feature=CELL),
        Segment("oligo(dT)VN / cDNA boundary", POLYT, placeholder=True),
        Segment("cDNA insert", "X" * 42, placeholder=True),
        Segment("mosaic end reverse complement", nx.ME_RC, "me"),
        Segment("S7 reverse complement", nx.S7_RC, "s7"),
        Segment("i7 (N701) reverse complement", revcomp(N701_INDEX), "cbc", feature=I7),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="scCARE-seq RNA library (T501/N701 representative)")


DNA = _dna_library()
DNA_P = tuple(sp.NEXTERA[k] for k in ("R1", "I1", "I2", "R2"))
RNA = _rna_library()
RNA_P = (sp.TRUSEQ["R1"], sp.NEXTERA["I1"], sp.TRUSEQ["I2"], sp.NEXTERA["R2"])
for _library, _primers in ((DNA, DNA_P), (RNA, RNA_P)):
    _problems = sp.verify(_library, _primers, required_roles=tuple(p.role for p in _primers))
    if _problems:
        raise ValueError(f"{_library.name}: " + "; ".join(_problems))

DNA_GEOMETRY = PCRGeometry.exact(
    DNA, PCRPrimer("N501", N501), PCRPrimer("N701", N701),
    statement="Supplementary Table 1 prints the complete N501 and N701 primers.",
)
RNA_GEOMETRY = PCRGeometry.exact(
    RNA, PCRPrimer("T501", T501), PCRPrimer("N701", N701),
    statement=("Supplementary Table 1 prints the complete T501 and N701 primers; "
               "article Methods 'mRNA library preparation' names TruSeq i5 + Nextera i7."),
)

TITLE = "scCARE-seq — chromatin architecture and RNA expression"
NOTES = "01_sccare-seq.html"
SOURCE = '<a href="https://doi.org/10.1038/s41594-023-01066-9">Qu, Sun et al., <i>Nature Structural &amp; Molecular Biology</i> (2023)</a>.'
SUMMARY = "In fixed nuclei, a UMI- and RT-barcoded oligo(dT) primes cDNA, MboI digestion and T4 ligation join chromatin contacts and the cDNA is made double-stranded; single nuclei are then sorted into wells, lysed, tagmented, pre-amplified with a three-primer mix and split into Nextera-indexed DNA and TruSeq-i5/Nextera-i7 RNA PCRs."
CAVEAT = ("Step order, enzymes and conditions come from the published article's Methods "
          "('scCARE-seq experimental procedure', Qu et al. 2023) and oligos from Supplementary Table 1. "
          "Unlike bulk CARE-seq, the single-cell protocol uses the plain oligo(dT) primer, has no "
          "biotin fill-in of the MboI ends and no streptavidin enrichment: the Methods ligate the "
          "digested ends directly with T4 DNA ligase, so the drawn contact junction is a reconstituted "
          "GATC (derived from MboI specificity, not printed). The role of the LIANTI carrier ssDNA in "
          "the lysis buffer is not explained. Pre-amplification uses the printed common (N7-PreA), DNA "
          "(N5-PreA) and RNA (A5-PreA) primers together; which primer pair copies which molecule is "
          "inferred from their sequences, and the pre-amplified products are not drawn base by base. "
          "The Tn5 preparation text says the transposase was loaded with an 'annealed barcoded "
          "adapter', whereas Supplementary Table 1 prints unbarcoded Tn5-N5/Tn5-N7; the model uses the "
          "printed adapters. Read lengths and index cycles are not stated (Methods: 'sequenced on the "
          "Illumina platform'; GEO lists NovaSeq 6000). ") + ("Supplementary Table 1 prints all eight RT barcodes and the complete N501–N516, N701–N724 and T501–T516 primer sets. The page uses Oligo-dT-1, N501, N701 and T501 as exact representatives; X and N remain variable biological or randomized sequence.")
FINAL_LIBRARIES = (
    ("DNA-contact library", DNA, DNA_P, "N501/N701-amplified, proximity-ligated genomic fragments.", "Paired genomic reads report chromatin contacts; the two sample indices identify the indexed library."),
    ("RNA library", RNA, RNA_P, "T501/N701-amplified well-barcoded UMI cDNA.", "Read 1 encounters the UMI and RT barcode; per-cell identity is inferred to come from the well i5/i7 index pair; Read 2 enters cDNA from the Nextera end."),
)


def oligos():
    yield oligo("Tn5-N5", [Segment("S5", nx.S5, "s5"), Segment("ME", nx.ME, "me")])
    yield oligo("Tn5-N7", [Segment("S7", nx.S7, "s7"), Segment("ME", nx.ME, "me")])
    yield oligo("Tn5-Rev", [Segment("ME reverse complement", TN5_REV, "me")], mods="/5Phos/")
    yield oligo("Carrier ssDNA", [Segment("carrier", CARRIER)])
    yield oligo("Common primer (N7-PreA)", [Segment("common primer", COMMON_PRIMER)])
    yield oligo("DNA primer (N5-PreA)", [Segment("DNA primer", DNA_PRIMER)])
    yield oligo("RNA primer (A5-PreA)", [Segment("RNA primer", RNA_PRIMER)])
    yield oligo("Oligo-dT-1", [
        Segment("RNA handle", RT_HANDLE, "r1"),
        Segment("UMI", UMI, "umi", placeholder=True, feature=UMI_FEATURE),
        Segment("RT barcode 1", RT_BARCODES[0], "cbc", feature=CELL),
        Segment("oligo(dT)VN", POLYT, placeholder=True),
    ], mods="none in scCARE-seq (article Methods RT primer); bulk CARE-seq uses Oligo-biotin-dT with one biotin-dT at position 8, 5'-ACGACGC[BiodT]CTTCCGATCT…")
    yield oligo("N501", [Segment("complete N501 primer", N501)])
    yield oligo("N701", [Segment("complete N701 primer", N701)])
    yield oligo("T501", [Segment("complete T501 primer", T501)])


def _contact_scene() -> Scene:
    """MboI ends re-ligated without fill-in reconstitute one GATC (derived, not printed)."""
    sc = Scene()
    sc.strand("contact", [molecule("locus A", "X" * 22, placeholder=True),
                          molecule("reconstituted MboI site", "GATC"),
                          molecule("locus B", "Y" * 22, placeholder=True)],
              label="contact chimera")
    sc.junction("contact", "locus A", "reconstituted MboI site", "proximity ligation")
    sc.labels("contact")
    return sc


def sections():
    # Order: article Methods "scCARE-seq experimental procedure" (fixation, RT, chromatin
    # digestion, proximity ligation, second strand synthesis, FACS and cell lysis,
    # tagmentation, preamplification, Hi-C / mRNA library preparation), consistent with
    # the GEO GSE211395 sample extract protocol.
    return [
        ("Reverse-transcribe RNA in fixed nuclei", [Row(chunks=[("RNA handle — UMI — RT barcode — oligo(dT) → poly(A) RNA", "umi", False)])], "1% formaldehyde-fixed nuclei; SuperScript IV RT with the printed UMI/RT-barcode oligo(dT) primers."),
        ("Digest with MboI and proximity-ligate in situ", _contact_scene().rows(), "MboI 37 °C 3 h after 0.6% SDS/Triton; T4 DNA ligase 16 °C 4 h without fill-in."),
        ("Sort single nuclei into wells, lyse and tagment", [Row(chunks=[("Tn5-N5 ** ligated DNA / ds cDNA ** Tn5-N7", "me", False)])], "Second-strand synthesis in nuclei, FACS into wells, reverse-crosslinking lysis and Tn5 tagmentation per well."),
        ("Pre-amplify, split and amplify", [Row(chunks=[("pre-amplification → half: N5xx/N7xx DNA-contact PCR | half: T5xx/N7xx RNA PCR", None, False)])], "Three-primer Q5 pre-amplification, ExoI, split 12 µl + 12 µl into indexed PCRs."),
    ]


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    nuclei = g.input("nuclei", "cells fixed in 1% formaldehyde 10 min (glycine-quenched), lysed in IGEPAL buffer to nuclei; 1.6 × 10^4 nuclei per RT tube", MaterialKind.BIOLOGICAL, state=sections()[0][1], primary=True)
    rt_primer = g.input("rt_primer", "UMI/RT-barcode oligo(dT) primer Oligo-dT-1..8, 8 µl of 25 µM (Supplementary Table 1; Methods)", MaterialKind.OLIGO)
    rt = g.input("rt", "SuperScript IV reverse transcriptase; 4→55 °C step programme, 55 °C 10 min", MaterialKind.CATALYST)
    copied = g.product("copied", "nuclei containing barcoded first-strand cDNA and fixed chromatin", MaterialKind.MOLECULAR, state=sections()[0][1])
    g.react("reverse_transcribe", "Reverse-transcribe RNA in intact fixed nuclei", StepKind.REVERSE_TRANSCRIPTION, ((nuclei, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (copied,))
    mboi = g.input("mboi", "MboI 25 U/µl in NEBuffer 2, 37 °C 3 h, after 0.6% SDS 62 °C 10 min and Triton X-100 quench", MaterialKind.CATALYST)
    digested = g.product("digested", "nuclei with MboI-digested chromatin (5′-GATC overhangs) and cDNA", MaterialKind.MOLECULAR, state=sections()[0][1])
    g.react("digest", "Digest chromatin with MboI in intact nuclei", StepKind.DIGESTION, ((copied, InputRole.TEMPLATE), (mboi, InputRole.CATALYST)), (digested,))
    ligase = g.input("contact_ligase", "T4 DNA ligase in 1× T4 ligase buffer, Triton X-100, BSA; 16 °C 4 h (no end fill-in in the single-cell Methods)", MaterialKind.CATALYST)
    contacts = g.product("contacts", "nuclei containing GATC-junction proximity ligation products and cDNA", MaterialKind.MOLECULAR, state=sections()[1][1])
    g.react("proximity_ligate", "Ligate digested chromatin ends in intact nuclei", StepKind.PROXIMITY_LIGATION, ((digested, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (contacts,))
    ss_enzyme = g.input("second_strand_enzyme", "NEBNext mRNA Second Strand Synthesis module (E6111), 16 °C 180 min", MaterialKind.CATALYST)
    ds = g.product("ds_cdna", "nuclei containing ligated chromatin and double-stranded cDNA", MaterialKind.MOLECULAR, state=sections()[1][1])
    g.react("second_strand", "Synthesize second-strand cDNA in intact nuclei", StepKind.EXTENSION, ((contacts, InputRole.TEMPLATE), (ss_enzyme, InputRole.CATALYST)), (ds,))
    sorter = g.input("facs", "40 µm strainer, DAPI stain, FACS of single nuclei into 96-well plates", MaterialKind.SELECTION)
    well = g.product("well", "one nucleus per well", MaterialKind.MOLECULAR, state=sections()[1][1])
    unsorted = g.product("unsorted", "nuclei not sorted into wells", MaterialKind.WASTE)
    g.react("sort", "Sort single nuclei into wells", StepKind.PARTITION, ((ds, InputRole.TEMPLATE), (sorter, InputRole.SELECTION)), (well, unsorted))
    g.discard(unsorted, "not deposited in a well")
    prok = g.input("reverse_crosslink", "3 µl scReverse Crosslink Buffer with 1 mg/ml proteinase K; 55 °C 2 h, 65 °C 12 h, 70 °C 30 min", MaterialKind.CATALYST)
    carrier = g.input("carrier_ssdna", "carrier ssDNA, 0.5 µM in the lysis buffer ('equally to LIANTI method used'; role not explained)", MaterialKind.OLIGO,
                      state=Construct([Segment("carrier ssDNA", CARRIER)], name="Carrier ssDNA"), source=SOURCE + " Supplementary Table 1")
    lysed = g.product("lysed", "reverse-crosslinked, deproteinized single-nucleus DNA and cDNA", MaterialKind.MOLECULAR, state=sections()[1][1])
    g.react("lyse", "Lyse and reverse-crosslink each sorted nucleus", StepKind.DEPROTEINIZATION, ((well, InputRole.TEMPLATE), (prok, InputRole.CATALYST), (carrier, InputRole.STRUCTURAL)), (lysed,))
    transposome = g.input("transposome", "Tn5-N5/Tn5-N7 + 5′-phosphorylated Tn5-Rev loaded in-house Tn5 (Supplementary Table 1)", MaterialKind.ADAPTER)
    tn5 = g.input("tn5", "0.025 µM in-house Tn5 in TAPS–DMF, 55 °C 10 min, 72 °C 5 min; then proteinase K 55/65/70 °C", MaterialKind.CATALYST)
    tagged = g.product("tagged", "Tn5-tagged ligation products and cDNA in one well", MaterialKind.MOLECULAR, state=sections()[2][1])
    twaste = g.product("tn5_waste", "candidate non-library tagmentation products", MaterialKind.MOLECULAR)
    g.react("tagment", "Tagment each well with Tn5", StepKind.TAGMENTATION, ((lysed, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)), (tagged, twaste))
    g.mark_terminal(Terminal(twaste.id, Disposition.UNRESOLVED, 'Candidate non-library tagmentation products lack the intended selected end combination. No class-specific physical removal is represented; their amplification and recovery fates are unresolved.'))
    src_t1 = SOURCE + " Supplementary Table 1"
    common = g.input("preamp_common", "Common primer N7-PreA (Supplementary Table 1)", MaterialKind.OLIGO, state=Construct([Segment("N7-PreA", COMMON_PRIMER)], name="N7-PreA"), source=src_t1)
    dna_pre = g.input("preamp_dna", "DNA primer N5-PreA (Supplementary Table 1)", MaterialKind.OLIGO, state=Construct([Segment("N5-PreA", DNA_PRIMER)], name="N5-PreA"), source=src_t1)
    rna_pre = g.input("preamp_rna", "RNA primer A5-PreA (Supplementary Table 1)", MaterialKind.OLIGO, state=Construct([Segment("A5-PreA", RNA_PRIMER)], name="A5-PreA"), source=src_t1)
    pre_pol = g.input("preamp_polymerase", "Q5 High-Fidelity 2× Master Mix; 72 °C 5 min, 12–14 cycles (65 °C anneal)", MaterialKind.CATALYST)
    pre_state = Construct([Segment("pre-amplified DNA-contact and cDNA fragments", "X" * 40, placeholder=True)], name="scCARE-seq pre-amplified well library (role level)")
    pre = g.product("preamplified", "pre-amplified single-cell library", MaterialKind.MOLECULAR, state=pre_state)
    g.react("preamplify", "Pre-amplify each well with the three-primer mix", StepKind.PCR, ((tagged, InputRole.TEMPLATE), (common, InputRole.STRUCTURAL), (dna_pre, InputRole.STRUCTURAL), (rna_pre, InputRole.STRUCTURAL), (pre_pol, InputRole.CATALYST)), (pre,),
            geometry=PCRGeometry.unavailable(pre_state, reason="Methods add common, DNA and RNA primers (40 µM each) as one mix; the per-molecule primer pairing (N5-PreA/N7-PreA on contact fragments, A5-PreA/N7-PreA on cDNA) is inferred from sequence and the products are not drawn base by base."))
    exo = g.input("exo_i", "Exonuclease I 2 µl of 5 U/µl, 37 °C 60 min, 80 °C 20 min", MaterialKind.CATALYST)
    cleaned = g.product("exo_treated", "pre-amplified library with excess primers degraded", MaterialKind.MOLECULAR, state=pre_state)
    spent = g.product("degraded_primers", "ExoI-degraded pre-amplification primers", MaterialKind.MOLECULAR)
    g.react("exo_digest", "Remove pre-amplification primers with Exonuclease I", StepKind.DIGESTION, ((pre, InputRole.TEMPLATE), (exo, InputRole.CATALYST)), (cleaned, spent))
    g.mark_terminal(Terminal(spent.id, Disposition.UNRESOLVED, "nucleotides remain in the well; no separation until the pooled column/bead cleanups after indexing PCR"))
    dna = g.product("dna_branch", "12 µl half for Hi-C (DNA-contact) library amplification", MaterialKind.MOLECULAR, state=pre_state)
    rna = g.product("rna_branch", "12 µl half for mRNA library amplification", MaterialKind.MOLECULAR, state=pre_state)
    g.react("split", "Split each well into two new 96-well plates", StepKind.SPLIT, ((cleaned, InputRole.TEMPLATE),), (dna, rna))
    reviewed = ["nuclei", "rt_primer", "rt", "reverse_transcribe", "mboi", "digest", "contact_ligase", "proximity_ligate", "second_strand_enzyme", "second_strand", "facs", "sort", "reverse_crosslink", "carrier_ssdna", "lyse", "transposome", "tn5", "tagment", "preamp_common", "preamp_dna", "preamp_rna", "preamp_polymerase", "preamplify", "exo_i", "exo_digest", "split"]
    finish_pcr_branch(g, "dna", dna, DNA, DNA_P, reviewed, label="scCARE-seq DNA-contact", geometry=DNA_GEOMETRY,
                      step_label="Nextera i5/i7 PCR, 4 cycles; pool, Zymo DCC-5 and 0.8× AMPure XP",
                      forward_label="Nextera i5 primer (N501 shown; N501–N516)", reverse_label="Nextera i7 primer (N701 shown; N701–N724)",
                      polymerase_label="Q5 High-Fidelity 2× Master Mix (60 °C anneal)",
                      forward_state=Construct([Segment("N501", N501)], name="N501"), reverse_state=Construct([Segment("N701", N701)], name="N701"), primer_source=src_t1)
    finish_pcr_branch(g, "rna", rna, RNA, RNA_P, reviewed, label="scCARE-seq RNA", geometry=RNA_GEOMETRY,
                      step_label="TruSeq i5 / Nextera i7 PCR, 10 cycles; pool, Zymo DCC-5 and 0.8× AMPure XP",
                      forward_label="TruSeq i5 primer (T501 shown; T501–T516)", reverse_label="Nextera i7 primer (N701 shown; N701–N724)",
                      polymerase_label="Q5 High-Fidelity 2× Master Mix (65 °C anneal)",
                      forward_state=Construct([Segment("T501", T501)], name="T501"), reverse_state=Construct([Segment("N701", N701)], name="N701"), primer_source=src_t1)
    g.declare_reviewed_scope(reviewed, SOURCE)
    return g.require_publishable()
