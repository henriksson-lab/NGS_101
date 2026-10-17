"""CRISPR-sciATAC — combinatorially indexed ATAC with a matched guide-RNA readout.

Every structural oligo below is transcribed from Liscovitch-Brauer et al. (2021)
Supplementary Table 2, which publishes the common oligos and all six barcode plates in
full.  The read layout follows Supplementary Fig. 1d,e and the reaction order follows the
Methods ("Transposome assembly", "CRISPR-sciATAC species mixing experiment").
"""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, Scene, duplex_rows, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from tagmentation import loaded_tn5_adapter

TITLE = "CRISPR-sciATAC — pooled perturbations with accessibility readout"
NOTES = "01_crispr-sciatac.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/s41587-021-00902-x">'
          'Liscovitch-Brauer et al., <i>Nature Biotechnology</i> (2021)</a>.')
SUMMARY = ("Accessible chromatin is tagmented by a barcoded TnY transposome and the "
           "polyadenylated CROP-seq guide transcript is reverse-transcribed with the "
           "matching barcode in the same well. Two further barcodes are added by PCR, so "
           "the accessibility and guide libraries carry the same three-part cell identity.")
CAVEAT = ("Every oligo drawn here is verbatim from Supplementary Table 2. Three regions "
          "are length-preserving placeholders: the accessible genomic insert, the 20-nt "
          "guide spacer (every gRNA in Supplementary Table 2 is 20 nt) and the 42-nt "
          "sgRNA-scaffold stretch between the spacer and the reverse-transcription site, "
          "whose length comes from the CROPseq-Guide-Puro map (Addgene 86708) named in "
          "the Methods rather than from the defining paper. The guide library is drawn "
          "with the unstaggered P5 PCR2 primer; the published primer set also carries "
          "+1, +2, +3 and +8 diversity bases, which shift every later Read 1 cycle. "
          "<b>Published sequence conflict:</b> Supplementary Table 2 prints the U6 outer "
          "PCR1 primer as TGCATATACGATCACAAGGCTGTTAG, while the CROPseq-Guide-Puro "
          "vector map (Addgene 86708), named in Methods, has "
          "TGCATATACGATACAAGGCTGTTAG at the corresponding hU6 promoter site. "
          "The model preserves the printed primer and uses regional PCR geometry "
          "because its exact site is absent on that vector map. The lab's actual U6 "
          "variant and whether the table has a typo remain unresolved.")

SUPP_T2 = "Liscovitch-Brauer et al. 2021, Supplementary Table 2"
SUPP_F1 = "Liscovitch-Brauer et al. 2021, Supplementary Fig. 1d,e"
WHITELIST = "CRISPR-sciATAC Supplementary Table 2 barcode plates"

# --------------------------------------------------------------- Supplementary Table 2
# "Common oligos", verbatim.
PMENT = "CTGTCTCTTATACACATCT"                                       # supplied /5Phos/
T5_TAGMENT = "TCGTCGGCAGCGTCTCCACGCAACCAACGCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG"
READ1_ATACSEQ = "CACCGTCTCCGCCTCAGATGTGTATAAGAGACAG"
INDEX1_ATACSEQ = "CTGTCTCTTATACACATCTGCCGTCCTCGATCGC"
READ2_ATACSEQ = "GCGATCGAGGACGGCAGATGTGTATAAGAGACAG"
U6_OUTER = "TGCATATACGATCACAAGGCTGTTAG"
# Plate-oligo flanks, read off the six barcode plates (identical in all 96/144 rows).
T7_SPACER = "CTGTCCCTGTCC"                 # "Spacer on T7"
RT_HANDLE = "CCGGGAGCTGCATGT"              # 5' arm of every RTi sgRNA primer
RT_ANNEAL = "TTTTTCAAGTTGATAA"             # 3' arm; anneals inside the sgRNA scaffold
GUIDE_HANDLE = "CGCCGTCCAGCTGA"            # 5' arm of every P5i sgRNA PCR1 primer
U6_INNER_SITE = "TTGTGGAAAGGACGAAACAC"     # 3' arm of every P7i sgRNA PCR2 primer

BARCODE_NT = 8                 # BC1, BC2 and BC3 are all 8 nt (every plate row)
BC1_COUNT = 96                 # barcoded T7 tagmentation / RTi sgRNA plates
BC2_COUNT = 96                 # P5i ATAC / P5i sgRNA PCR1 plates
BC3_COUNT = 144                # P7i ATAC PCR1 / P7i sgRNA PCR2 plates
GUIDE_SPACER_NT = 20           # every gRNA in Supplementary Table 2
SCAFFOLD_TO_RT_NT = 42         # spacer end -> RT site, CROPseq-Guide-Puro (Addgene 86708)
STAGGER_NT = 8                 # longest published P5 sgRNA PCR2 diversity stagger

# ------------------------------------------------------------------ derived decomposition
# The mosaic end, the Read 1/Read 2 entry points and the reverse complements below are
# derived from the verbatim oligos above rather than retyped.
READ1_SITE = READ1_ATACSEQ[:-len(nx.ME)]        # CACCGTCTCCGCCTC
READ2_SITE = READ2_ATACSEQ[:-len(nx.ME)]        # GCGATCGAGGACGGC
T5_SPACER = T5_TAGMENT[len(nx.S5):-len(READ1_ATACSEQ)]
SCAFFOLD_RT_SITE = revcomp(RT_ANNEAL)           # TTATCAACTTGAAAAA, inside the scaffold

if PMENT != nx.ME_RC:
    raise ValueError("pMENT is the canonical Tn5 mosaic-end complement")
if not (READ1_ATACSEQ.endswith(nx.ME) and READ2_ATACSEQ.endswith(nx.ME)):
    raise ValueError("both ATACseq read primers must end in the canonical mosaic end")
if nx.S5 + T5_SPACER + READ1_SITE + nx.ME != T5_TAGMENT:
    raise ValueError("the T5 tagmentation oligo does not decompose into s5/spacer/ME")
if INDEX1_ATACSEQ != nx.ME_RC + revcomp(READ2_SITE):
    raise ValueError("Index 1 ATACseq must read out of the mosaic end across the T7 arm")
if not il.TRUSEQ_READ1.startswith(il.P5[-4:]):
    raise ValueError("the TruSeq Read 1 arm must overlap the 3' ACAC of P5")

# ---------------------------------------------------------------------- ordered oligos
# Barcode windows are written N, exactly as Supplementary Fig. 1d,e writes them.
BC = "N" * BARCODE_NT
T7_TAGMENT = nx.S7 + T7_SPACER + BC + READ2_SITE + nx.ME
RT_PRIMER = RT_HANDLE + BC + RT_ANNEAL
P5I_ATAC = il.P5 + BC + nx.S5
P7I_ATAC = il.P7 + BC + nx.S7
P5I_SGRNA = GUIDE_HANDLE + BC + RT_HANDLE
P7I_SGRNA = il.P7 + BC + U6_INNER_SITE
TRUSEQ_R1_ARM = il.TRUSEQ_READ1[4:]             # P5 already supplies the leading ACAC
P5_SGRNA_PCR2 = il.P5 + TRUSEQ_R1_ARM + GUIDE_HANDLE
P5_SGRNA_PCR2_STAGGERED = il.P5 + TRUSEQ_R1_ARM + "N" * STAGGER_NT + GUIDE_HANDLE


def _barcode_feature(part: str):
    return feature(f"cell_barcode_{part.lower()}", "cell_barcode", "combinatorial",
                   group="cell_barcode", part=part, whitelist=WHITELIST)


def _guide_feature():
    return feature("guide_spacer", "guide_barcode", "whitelist",
                   whitelist="CRISPR-sciATAC Supplementary Table 2 gRNA library")


_BC_LETTER = {"BC1": "B", "BC2": "C", "BC3": "D"}


def barcode(part: str, *, rc: bool = False):
    """One drawn cell-barcode window; the three parts stay visually distinct."""
    name = part + (" reverse complement" if rc else "")
    return seg(name, _BC_LETTER[part] * BARCODE_NT, "cbc", placeholder=True,
               feature=_barcode_feature(part))


def ordered_barcode(part: str):
    """The same window as the vendor writes it in an ordered oligo (N, not a stand-in)."""
    return seg(part, BC, "cbc", feature=_barcode_feature(part))


def guide_spacer(*, rc: bool = False):
    name = "guide spacer" + (" reverse complement" if rc else "")
    return seg(name, "G" * GUIDE_SPACER_NT, "cbc", placeholder=True,
               length_bp=GUIDE_SPACER_NT, feature=_guide_feature())


def scaffold_stretch(*, rc: bool = False):
    name = "sgRNA scaffold" + (" reverse complement" if rc else "")
    return seg(name, "X" * SCAFFOLD_TO_RT_NT, placeholder=True,
               length_bp=SCAFFOLD_TO_RT_NT,
               note="length from CROPseq-Guide-Puro, Addgene 86708")


# -------------------------------------------------------------------- ordered-oligo states
def t5_oligo() -> tuple:
    return (seg("s5", nx.S5, "s5"), seg("T5 spacer", T5_SPACER),
            seg("Read 1 ATACseq site", READ1_SITE, "r1"),
            seg("mosaic end", nx.ME, "me"))


def t7_oligo() -> tuple:
    return (seg("s7", nx.S7, "s7"), seg("T7 spacer", T7_SPACER),
            ordered_barcode("BC1"),
            seg("Read 2 ATACseq site", READ2_SITE, "r2"),
            seg("mosaic end", nx.ME, "me"))


def pment_oligo() -> tuple:
    return (seg("pMENT", PMENT, "me"),)


def rt_oligo() -> tuple:
    return (seg("RT handle", RT_HANDLE), ordered_barcode("BC1"),
            seg("sgRNA scaffold anneal", RT_ANNEAL))


def p5i_atac_oligo() -> tuple:
    return (seg("P5", il.P5, "p5"), ordered_barcode("BC2"), seg("s5", nx.S5, "s5"))


def p7i_atac_oligo() -> tuple:
    return (seg("P7", il.P7, "p7"), ordered_barcode("BC3"), seg("s7", nx.S7, "s7"))


def p5i_sgrna_oligo() -> tuple:
    return (seg("guide PCR handle", GUIDE_HANDLE), ordered_barcode("BC2"),
            seg("RT handle", RT_HANDLE))


def p7i_sgrna_oligo() -> tuple:
    return (seg("P7", il.P7, "p7"), ordered_barcode("BC3"),
            seg("U6 promoter 3' end", U6_INNER_SITE))


def u6_outer_oligo() -> tuple:
    return (seg("U6 outer primer", U6_OUTER),)


def p5_sgrna_pcr2_oligo() -> tuple:
    return (seg("P5", il.P5, "p5"), seg("TruSeq Read 1 arm", TRUSEQ_R1_ARM, "r1"),
            seg("guide PCR handle", GUIDE_HANDLE))


# ------------------------------------------------------------------------- MEDS adapters
def meds_a() -> Scene:
    """MEDS A: the unbarcoded T5 tagmentation oligo annealed to pMENT."""
    return loaded_tn5_adapter(list(t5_oligo()), mosaic_end="mosaic end",
                              label="T5 tagment sciATAC", bottom_name="pMENT",
                              bottom_label="5'-phosphorylated pMENT", mod5="p")


def meds_b() -> Scene:
    """MEDS B: one of 96 barcoded T7 tagmentation oligos annealed to pMENT."""
    top = [seg("s7", nx.S7, "s7"), seg("T7 spacer", T7_SPACER), barcode("BC1"),
           seg("Read 2 ATACseq site", READ2_SITE, "r2"),
           seg("mosaic end", nx.ME, "me")]
    return loaded_tn5_adapter(top, mosaic_end="mosaic end",
                              label="Barcoded T7 tagment sciATAC", bottom_name="pMENT",
                              bottom_label="5'-phosphorylated pMENT", mod5="p")


# --------------------------------------------------------------------- ATAC constructs
def tagged_fragment(insert_nt: int = 34) -> Construct:
    """Gap-filled T5/T7 transposition product: only this end pair can bridge P5 to P7."""
    return Construct([
        seg("s5", nx.S5, "s5"), seg("T5 spacer", T5_SPACER),
        seg("Read 1 ATACseq site", READ1_SITE, "r1"),
        seg("left mosaic end", nx.ME, "me"),
        seg("accessible genomic DNA", "X" * insert_nt, placeholder=True),
        seg("right mosaic end reverse complement", nx.ME_RC, "me"),
        seg("Read 2 ATACseq site reverse complement", revcomp(READ2_SITE), "r2"),
        barcode("BC1", rc=True),
        seg("T7 spacer reverse complement", revcomp(T7_SPACER)),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
    ], name="BC1-tagmented accessible fragment")


def atac_library(insert_nt: int = 34) -> Construct:
    body = list(tagged_fragment(insert_nt))
    return Construct([seg("P5", il.P5, "p5"), barcode("BC2"), *body,
                      barcode("BC3", rc=True),
                      seg("P7 reverse complement", il.P7_RC, "p7")],
                     name="CRISPR-sciATAC accessibility library")


# -------------------------------------------------------------------- guide constructs
def guide_transcript() -> list:
    """The CROP-seq Pol II transcript that also carries the U6 cassette and the gRNA."""
    return [seg("puromycin-resistance transcript and WPRE", "X" * 24, placeholder=True),
            seg("U6 promoter 3' end", U6_INNER_SITE),
            guide_spacer(), scaffold_stretch(),
            seg("RT-primer site in the scaffold", SCAFFOLD_RT_SITE),
            seg("scaffold 3' end and poly(A) tail", "X" * 28, placeholder=True)]


def rt_scene() -> Scene:
    """The barcoded RT primer annealed inside the sgRNA scaffold of the transcript."""
    sc = Scene()
    sc.strand("transcript", guide_transcript(),
              label="gRNA-bearing polyadenylated transcript")
    sc.anneal("rt", [seg("RT handle", RT_HANDLE), barcode("BC1"),
                     seg("sgRNA scaffold anneal", RT_ANNEAL)],
              to="transcript", pair=("sgRNA scaffold anneal", "RT-primer site in the scaffold"),
              label="RTi sgRNA primer", unpaired=("RT handle", "BC1"))
    sc.arrow("rt", "reverse transcription")
    sc.labels("transcript")
    return sc


def guide_cdna() -> Construct:
    """First-strand cDNA (drawn in its own 5'->3' sense) with BC1 already attached."""
    return Construct([
        seg("RT handle", RT_HANDLE), barcode("BC1"),
        seg("sgRNA scaffold anneal", RT_ANNEAL),
        scaffold_stretch(rc=True), guide_spacer(rc=True),
        seg("U6 promoter 3' end reverse complement", revcomp(U6_INNER_SITE)),
        seg("U6 promoter region reverse complement", "[U6 promoter and vector]",
            placeholder=True),
    ], name="BC1-barcoded guide cDNA")


def guide_amplicon() -> Construct:
    """PCR1 product: the P5i sgRNA primer has added the guide handle and BC2."""
    body = list(guide_cdna())
    return Construct([seg("guide PCR handle", GUIDE_HANDLE), barcode("BC2"), *body],
                     name="first-round guide amplicon")


def guide_library() -> Construct:
    body = list(guide_amplicon())[:-1]            # PCR2 ends at the U6 inner primer site
    return Construct([seg("P5", il.P5, "p5"),
                      seg("TruSeq Read 1 arm", TRUSEQ_R1_ARM, "r1"), *body,
                      barcode("BC3", rc=True),
                      seg("P7 reverse complement", il.P7_RC, "p7")],
                     name="CRISPR-sciATAC guide library")


ATAC = atac_library()
GUIDE = guide_library()
TAGGED = tagged_fragment()
GUIDE_CDNA = guide_cdna()
GUIDE_AMPLICON = guide_amplicon()

# ------------------------------------------------------------------ sequencing, Supp. Fig. 1d,e
ATAC_PRIMERS = (
    sp.custom("Read 1", "Read 1 ATACseq", READ1_ATACSEQ, SUPP_T2),
    sp.custom("Index 1 (i7)", "Index 1 ATACseq", INDEX1_ATACSEQ, SUPP_T2,
              "One 43-cycle index read spanning BC1, the T7 spacer, s7 and BC3; "
              "it carries no sample index."),
    sp.custom("Index 2 (i5)",
              "flow-cell P5 oligo (MiSeq / HiSeq 2500 forward-strand i5 read)", il.P5,
              'Illumina "Indexed Sequencing Overview Guide" #15057455',
              "Supplementary Table 2 lists no i5 primer: on a forward-strand instrument "
              "the flow-cell P5 oligo primes the 8-cycle BC2 read."),
    sp.custom("Read 2", "Read 2 ATACseq", READ2_ATACSEQ, SUPP_T2),
)
GUIDE_PRIMERS = (sp.TRUSEQ["R1"],)

ATAC_READ_LENGTHS = {"Read 1": 45, "Index 1 (i7)": 43, "Index 2 (i5)": 8, "Read 2": 45}
GUIDE_READ_LENGTHS = {"Read 1": 162}
ATAC_READ_USES = {"Read 1": "include", "Index 1 (i7)": "extract",
                  "Index 2 (i5)": "extract", "Read 2": "include"}
GUIDE_READ_USES = {"Read 1": "extract"}

FINAL_LIBRARIES = (
    ("Accessibility library", ATAC, ATAC_PRIMERS,
     "BC1 arrives on the T7 transposome arm, BC2 and BC3 on the two PCR1 primers; "
     "the three together are the cell identity. There is no sample index.",
     "All three ATACseq primers are custom (Supplementary Table 2); BC2 is read off the "
     "flow-cell P5 oligo on the forward-strand instruments used in the paper."),
    ("Guide library", GUIDE, GUIDE_PRIMERS,
     "A TruSeq-compatible amplicon: the same BC1/BC2/BC3 combination precedes the guide "
     "spacer, which is read directly rather than through a separate guide barcode.",
     "One 162-cycle TruSeq Read 1 covers all three barcodes and the guide spacer, so the "
     "run takes no index reads at all."),
)
READ_LENGTHS = {"Accessibility library": ATAC_READ_LENGTHS,
                "Guide library": GUIDE_READ_LENGTHS}
READ_USES = {"Accessibility library": ATAC_READ_USES, "Guide library": GUIDE_READ_USES}

for _name, _lib, _primers in (("accessibility", ATAC, ATAC_PRIMERS),
                              ("guide", GUIDE, GUIDE_PRIMERS)):
    _problems = sp.verify(_lib, _primers,
                          required_roles=tuple(p.role for p in _primers))
    if _problems:
        raise ValueError(f"CRISPR-sciATAC {_name} library: " + "; ".join(_problems))


def oligos():
    return (
        f"<b>pMENT</b> — 5′-/5Phos/{PMENT}-3′",
        f"<b>T5 tagment sciATAC (MEDS A top)</b> — 5′-{T5_TAGMENT}-3′",
        f"<b>Barcoded T7 tagment sciATAC (MEDS B top, {BC1_COUNT} wells)</b> — "
        f"5′-{nx.S7}{T7_SPACER}[8-nt BC1]{READ2_SITE}{nx.ME}-3′",
        f"<b>RTi sgRNA ({BC1_COUNT} wells, BC1 matched to MEDS B)</b> — "
        f"5′-{RT_HANDLE}[8-nt BC1]{RT_ANNEAL}-3′",
        f"<b>P5i ATAC PCR1 ({BC2_COUNT} wells)</b> — 5′-{il.P5}[8-nt BC2]{nx.S5}-3′",
        f"<b>P7i ATAC PCR1 ({BC3_COUNT} wells)</b> — 5′-{il.P7}[8-nt BC3]{nx.S7}-3′",
        f"<b>P5i sgRNA PCR1 ({BC2_COUNT} wells, BC2 matched to P5i ATAC)</b> — "
        f"5′-{GUIDE_HANDLE}[8-nt BC2]{RT_HANDLE}-3′",
        f"<b>U6 outer sgRNA PCR1</b> — 5′-{U6_OUTER}-3′",
        f"<b>P5 sgRNA PCR2</b> — 5′-{P5_SGRNA_PCR2}-3′",
        f"<b>P5 sgRNA PCR2 +1/+2/+3/+8 stagger</b> — "
        f"5′-{il.P5}{TRUSEQ_R1_ARM}[0–{STAGGER_NT}-nt diversity stagger]{GUIDE_HANDLE}-3′",
        f"<b>P7i sgRNA PCR2 ({BC3_COUNT} wells, BC3 matched to P7i ATAC)</b> — "
        f"5′-{il.P7}[8-nt BC3]{U6_INNER_SITE}-3′",
        f"<b>Read 1 ATACseq</b> — 5′-{READ1_ATACSEQ}-3′",
        f"<b>Index 1 ATACseq</b> — 5′-{INDEX1_ATACSEQ}-3′",
        f"<b>Read 2 ATACseq</b> — 5′-{READ2_ATACSEQ}-3′",
    )


def _atac_pcr_geometry() -> PCRGeometry:
    return PCRGeometry.exact(
        TAGGED,
        PCRPrimer("P5i ATAC PCR1", P5I_ATAC, anneal_nt=len(nx.S5)),
        PCRPrimer("P7i ATAC PCR1", P7I_ATAC, anneal_nt=len(nx.S7)),
        statement=("Supplementary Table 2 publishes both indexed ATAC PCR1 plates in "
                   "full; their 3' ends are the s5 and s7 entry points of the two "
                   "transposome arms, so they face inward on every tagmented fragment."))


def _guide_pcr1_geometry() -> PCRGeometry:
    return PCRGeometry.regional(
        GUIDE_CDNA,
        PCRPrimer("P5i sgRNA PCR1", P5I_SGRNA, anneal_nt=len(RT_HANDLE)),
        reverse_pool="U6 outer sgRNA PCR1",
        reverse_region="U6 promoter region reverse complement",
        statement=("Supplementary Table 2 gives the U6 outer primer verbatim, but the "
                   "defining paper does not reproduce the CROPseq-Guide-Puro U6 promoter "
                   "it anneals to, so its exact coordinate on this cDNA is not reported."))


def _guide_pcr2_geometry() -> PCRGeometry:
    return PCRGeometry.exact(
        GUIDE_AMPLICON,
        PCRPrimer("P5 sgRNA PCR2", P5_SGRNA_PCR2, anneal_nt=len(GUIDE_HANDLE)),
        PCRPrimer("P7i sgRNA PCR2", P7I_SGRNA, anneal_nt=len(U6_INNER_SITE)),
        statement=("Supplementary Table 2 publishes both second-round guide primers; the "
                   "P7i primer is nested inside the U6 outer primer, which is what trims "
                   "the amplicon to the U6 3' end."))


def protocol_graph() -> ProtocolGraph:
    g = ProtocolGraph(TITLE, authority=SOURCE)

    cells = g.input("cells", "Cas9-expressing cells before perturbation",
                    MaterialKind.BIOLOGICAL, primary=True)
    guide_vector = g.input(
        "guide_vector", "CROPseq-Guide-Puro lentiviral transfer vector (Addgene 86708)",
        MaterialKind.MOLECULAR)
    lentivirus = g.input("lentiviral_machinery", "lentiviral packaging and integration",
                         MaterialKind.CATALYST)
    perturbed = g.product(
        "perturbed_cells",
        "cells with an integrated U6-gRNA cassette and a polyadenylated gRNA transcript",
        MaterialKind.BIOLOGICAL)
    g.react("transduce", "Transduce the pooled gRNA library and select with puromycin",
            StepKind.INTEGRATION,
            ((cells, InputRole.TEMPLATE), (guide_vector, InputRole.STRUCTURAL),
             (lentivirus, InputRole.CATALYST)), (perturbed,))

    fixation = g.input("fixation", "glyoxal fixation and gentle lysis buffers",
                       MaterialKind.CONDITION)
    nuclei = g.product("nuclei", "fixed nuclei retaining chromatin and gRNA transcripts",
                       MaterialKind.BIOLOGICAL)
    g.react("isolate_nuclei", "Fix with glyoxal and lyse to recover intact nuclei",
            StepKind.PERMEABILIZATION,
            ((perturbed, InputRole.TEMPLATE), (fixation, InputRole.CONDITION)),
            (nuclei,))

    t5 = g.input("t5_oligo", "T5 tagment sciATAC oligo", MaterialKind.OLIGO,
                 state=t5_oligo())
    t7 = g.input("t7_oligo", f"barcoded T7 tagment sciATAC oligo ({BC1_COUNT} wells)",
                 MaterialKind.OLIGO, state=t7_oligo())
    pment = g.input("pment", "5'-phosphorylated pMENT mosaic-end oligo",
                    MaterialKind.OLIGO, state=pment_oligo())
    meds = g.product("meds", "MEDS A and barcoded MEDS B mosaic-end duplexes",
                     MaterialKind.ADAPTER,
                     state=[*meds_a().rows(), *meds_b().rows()])
    g.react("anneal_meds", "Anneal pMENT to the T5 and barcoded T7 oligos (MEDS A/B)",
            StepKind.HYBRIDIZATION,
            ((t5, InputRole.TEMPLATE), (pment, InputRole.STRUCTURAL),
             (t7, InputRole.STRUCTURAL)), (meds,))

    tny = g.input("tny", "TnY transposase (Vibrio parahaemolyticus ViPar, P50K + M53Q)",
                  MaterialKind.CATALYST)
    tagged = g.product(
        "tagged_nuclei",
        "nuclei holding BC1-tagmented accessible DNA and intact gRNA transcripts",
        MaterialKind.MOLECULAR,
        state=[*duplex_rows(TAGGED, label="accessibility branch"),
               Row(chunks=[("the polyadenylated gRNA transcript is still intact",
                            None, False)])])
    tagment_waste = g.product("tagment_waste",
                              "T5/T5 and T7/T7 transposition products", MaterialKind.WASTE)
    g.react("tagment", "Tagment accessible chromatin well by well with barcoded TnY",
            StepKind.TAGMENTATION,
            ((nuclei, InputRole.TEMPLATE), (meds, InputRole.STRUCTURAL),
             (tny, InputRole.CATALYST)), (tagged, tagment_waste))
    g.discard(tagment_waste,
              "carries the same transposome arm at both ends and so cannot bridge P5 to P7")

    rt_primer = g.input("rt_primer",
                        f"RTi sgRNA primer ({BC1_COUNT} wells, BC1 matched to MEDS B)",
                        MaterialKind.OLIGO, state=rt_oligo())
    rtase = g.input("rtase", "RevertAid reverse transcriptase", MaterialKind.CATALYST)
    copied = g.product(
        "rt_nuclei",
        "nuclei holding the tagmented fragment and BC1-barcoded guide cDNA",
        MaterialKind.MOLECULAR,
        state=[*duplex_rows(TAGGED, label="accessibility branch"), *rt_scene().rows()])
    g.react("reverse_transcribe",
            "Reverse-transcribe the gRNA transcript in situ with the matching BC1 primer",
            StepKind.REVERSE_TRANSCRIPTION,
            ((tagged, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
             (rtase, InputRole.CATALYST)), (copied,))

    redistribution = g.input("redistribution",
                             "pooling and ~20-nuclei-per-well redistribution (2 ul at 10 nuclei/ul)",
                             MaterialKind.CONDITION)
    redistributed = g.product("redistributed", "nuclei assigned a second plate coordinate",
                              MaterialKind.MOLECULAR, state=copied.state)
    collisions = g.product("collisions", "empty or colliding wells", MaterialKind.WASTE)
    g.react("redistribute", "Pool all wells and redistribute into the PCR plate",
            StepKind.PARTITION,
            ((copied, InputRole.TEMPLATE), (redistribution, InputRole.CONDITION)),
            (redistributed, collisions), changes_molecule=False)
    g.discard(collisions, "does not yield an interpretable three-part cell identity")

    protease = g.input("protease", "proteinase K", MaterialKind.CATALYST)
    extracted = g.product("extracted", "released tagmented DNA and guide cDNA",
                          MaterialKind.MOLECULAR, state=copied.state)
    g.react("extract", "Digest nuclei with proteinase K and quench with PMSF",
            StepKind.DEPROTEINIZATION,
            ((redistributed, InputRole.TEMPLATE), (protease, InputRole.CATALYST)),
            (extracted,))

    atac_template = g.product("atac_template", "accessibility PCR1 template",
                              MaterialKind.MOLECULAR, state=TAGGED)
    guide_template = g.product("guide_template", "guide PCR1 template",
                               MaterialKind.MOLECULAR, state=GUIDE_CDNA)
    g.react("split", "Direct each released molecule to its own primer pair",
            StepKind.SPLIT, ((extracted, InputRole.TEMPLATE),),
            (atac_template, guide_template), changes_molecule=False)

    polymerase = g.input("pfux7", "PfuX7 DNA polymerase", MaterialKind.CATALYST)

    p5i_atac = g.input("p5i_atac", f"P5i ATAC PCR1 primer ({BC2_COUNT} wells)",
                       MaterialKind.OLIGO, state=p5i_atac_oligo())
    p7i_atac = g.input("p7i_atac", f"P7i ATAC PCR1 primer ({BC3_COUNT} wells)",
                       MaterialKind.OLIGO, state=p7i_atac_oligo())
    atac_amplicon = g.product("atac_amplicon", "indexed accessibility amplicon",
                              MaterialKind.MOLECULAR, state=ATAC)
    g.react("atac_pcr",
            "Fill the transposition gaps and amplify with the BC2/BC3 ATAC primers",
            StepKind.PCR,
            ((atac_template, InputRole.TEMPLATE), (p5i_atac, InputRole.STRUCTURAL),
             (p7i_atac, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
            (atac_amplicon,), geometry=_atac_pcr_geometry())
    atac_lib = g.product("atac_library", "CRISPR-sciATAC accessibility library",
                         MaterialKind.LIBRARY, state=ATAC)
    atac_excess = g.product("atac_excess", "free oligos and off-size products",
                            MaterialKind.WASTE)
    g.react("atac_cleanup", "Column-purify and size-select with 0.9X Ampure XP beads",
            StepKind.CLEANUP, ((atac_amplicon, InputRole.TEMPLATE),),
            (atac_lib, atac_excess))
    g.discard(atac_excess, "removed by the 0.9X bead size selection")
    g.sequence(atac_lib, SequencingPlan(
        library=ATAC, primers=ATAC_PRIMERS,
        acquired_roles=tuple(p.role for p in ATAC_PRIMERS),
        read_lengths=ATAC_READ_LENGTHS, read_uses=ATAC_READ_USES))

    p5i_sgrna = g.input("p5i_sgrna", f"P5i sgRNA PCR1 primer ({BC2_COUNT} wells)",
                        MaterialKind.OLIGO, state=p5i_sgrna_oligo())
    u6_outer = g.input("u6_outer", "U6 outer sgRNA PCR1 primer", MaterialKind.OLIGO,
                       state=u6_outer_oligo())
    guide_round1 = g.product("guide_round1", "first-round guide amplicon",
                             MaterialKind.MOLECULAR, state=GUIDE_AMPLICON)
    g.react("guide_pcr1",
            "Amplify the guide cDNA with the U6 outer and BC2 primers",
            StepKind.PCR,
            ((guide_template, InputRole.TEMPLATE), (p5i_sgrna, InputRole.STRUCTURAL),
             (u6_outer, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
            (guide_round1,), geometry=_guide_pcr1_geometry())

    p5_pcr2 = g.input("p5_sgrna_pcr2", "P5 sgRNA PCR2 primer (unstaggered and +1/+2/+3/+8)",
                      MaterialKind.OLIGO, state=p5_sgrna_pcr2_oligo())
    p7i_sgrna = g.input("p7i_sgrna", f"P7i sgRNA PCR2 primer ({BC3_COUNT} wells)",
                        MaterialKind.OLIGO, state=p7i_sgrna_oligo())
    guide_round2 = g.product("guide_round2", "second-round guide amplicon",
                             MaterialKind.MOLECULAR, state=GUIDE)
    g.react("guide_pcr2",
            "Add P5, the TruSeq Read 1 arm and BC3 in a nested second PCR",
            StepKind.PCR,
            ((guide_round1, InputRole.TEMPLATE), (p5_pcr2, InputRole.STRUCTURAL),
             (p7i_sgrna, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
            (guide_round2,), geometry=_guide_pcr2_geometry())
    guide_lib = g.product("guide_library", "CRISPR-sciATAC guide library",
                          MaterialKind.LIBRARY, state=GUIDE)
    guide_excess = g.product("guide_excess", "off-size guide PCR products",
                             MaterialKind.WASTE)
    g.react("guide_cleanup", "Gel-extract the expected ~250-bp guide band",
            StepKind.CLEANUP, ((guide_round2, InputRole.TEMPLATE),),
            (guide_lib, guide_excess))
    g.discard(guide_excess, "runs away from the expected ~250-bp band on the 2% E-gel")
    g.sequence(guide_lib, SequencingPlan(
        library=GUIDE, primers=GUIDE_PRIMERS,
        acquired_roles=tuple(p.role for p in GUIDE_PRIMERS),
        read_lengths=GUIDE_READ_LENGTHS, read_uses=GUIDE_READ_USES))

    g.declare_reviewed_scope((
        "cells", "guide_vector", "lentiviral_machinery", "transduce",
        "fixation", "isolate_nuclei", "t5_oligo", "t7_oligo", "pment", "anneal_meds",
        "tny", "tagment", "rt_primer", "rtase", "reverse_transcribe",
        "redistribution", "redistribute", "protease", "extract", "split",
        "pfux7", "p5i_atac", "p7i_atac", "atac_pcr", "atac_cleanup", "atac_library",
        "p5i_sgrna", "u6_outer", "guide_pcr1", "p5_sgrna_pcr2", "p7i_sgrna",
        "guide_pcr2", "guide_cleanup", "guide_library",
    ), SOURCE)
    return g.require_publishable()
