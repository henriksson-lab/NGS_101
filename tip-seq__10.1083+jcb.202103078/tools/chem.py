"""TIP-seq / sciTIP-seq molecular model from Bartlett et al. Table S1."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature, oligo
from pcr import PCRGeometry, PCRPrimer
from tagmentation import loaded_tn5_adapter

TITLE = "TIP-seq — targeted insertion of promoters"
NOTES = "01_tip_seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1083/jcb.202103078">'
          'Bartlett et al., <i>Journal of Cell Biology</i> (2021)</a>.')
SUMMARY = ("Antibody-tethered pA–Tn5 inserts a T7 promoter beside targeted chromatin. "
           "T7 transcription linearly amplifies that one-ended insertion before random-hexamer "
           "reverse transcription, second-strand synthesis, MeB tagmentation and indexed PCR.")
CAVEAT = ("The workflow draws the published representative sciTIP r5001/S502/N701 branch. "
          "Table S1 supplies 384 r5 transposons and the full i5/i7 PCR-index families.")

ME_T7_BULK = "GAATTTAATACGACTCACTATAGGGTAGATCGCTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
ME_SC_R5001 = "GAATTTAATACGACTCACTATAGGGAGATCCACGCGCTCTAATTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
ME_REV = nx.ME_RC
ME_B = nx.S7 + nx.ME
SSS_SCI = "GGGAGATCCACGCGC"
R5_001 = "TCTAAT"
S502 = il.P5 + "CTCTCTAT" + SSS_SCI
N701 = il.P7 + "TCGCCTTA" + nx.S7
NXT_R1 = nx.S5 + nx.ME
NXT_I7 = nx.ME_RC + nx.S7_RC
NXT_I5 = il.P5


def r5_feature():
    return feature("scitip_r5", "cell_barcode", "whitelist",
                   group="scitip_cell", part="targeted Tn5 well",
                   whitelist="Bartlett et al. Table S1 r5001-r5384")


def i5_feature():
    return feature("scitip_i5", "cell_barcode", "combinatorial",
                   group="scitip_cell", part="PCR well i5",
                   whitelist="Bartlett et al. Table S1 modified S5 indexes")


def i7_feature():
    return feature("scitip_i7", "cell_barcode", "combinatorial",
                   group="scitip_cell", part="PCR well i7",
                   whitelist="Bartlett et al. Table S1 N7 indexes")


def me_t7_parts():
    return [
        seg("T7 upstream stabilizer/promoter", ME_SC_R5001[:22]),
        seg("T7 start + second-strand site", SSS_SCI, "r1"),
        seg("r5-001 cell barcode", R5_001, "cbc", feature=r5_feature()),
        seg("S5", nx.S5, "s5"),
        seg("mosaic end", nx.ME, "me"),
    ]


def loaded_me_t7():
    sc = loaded_tn5_adapter(me_t7_parts(), mosaic_end="mosaic end",
                            label="ME-scT7-r5001",
                            bottom_name="ME-rev",
                            bottom_label="5′-phosphorylated ME-rev",
                            mod5="p")
    sc.mark("Tn5 top", "r5-001 cell barcode", "first combinatorial cell index")
    return sc


def loaded_me_b():
    return loaded_tn5_adapter(
        [seg("S7", nx.S7, "s7"), seg("mosaic end", nx.ME, "me")],
        mosaic_end="mosaic end", label="ME-B-only Tn5",
        bottom_name="ME-rev", bottom_label="5′-phosphorylated ME-rev", mod5="p")


TAGGED = Construct([*me_t7_parts(), seg("target-proximal genomic DNA", "X" * 38,
                                        placeholder=True)],
                   name="ME-scT7-tagged chromatin")

RNA_PRODUCT = [
    seg("transcribed T7 start + second-strand site", SSS_SCI, "r1"),
    seg("transcribed r5-001 barcode", R5_001, "cbc", feature=r5_feature()),
    seg("transcribed S5", nx.S5, "s5"),
    seg("transcribed mosaic end", nx.ME, "me"),
    seg("target-proximal RNA copy", "X" * 38, placeholder=True),
]

PRE_PCR = Construct([
    seg("second-strand site", SSS_SCI, "r1"),
    seg("r5-001 cell barcode", R5_001, "cbc", feature=r5_feature()),
    seg("S5", nx.S5, "s5"), seg("left mosaic end", nx.ME, "me"),
    seg("target-proximal cDNA", "X" * 38, placeholder=True),
    seg("right mosaic end reverse complement", nx.ME_RC, "me"),
    seg("S7 reverse complement", nx.S7_RC, "s7"),
], name="twice-tagmented sciTIP cDNA")

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"),
    seg("PCR-well i5", "I" * 8, "cbc", placeholder=True, feature=i5_feature()),
    seg("second-strand site", SSS_SCI, "r1"),
    seg("r5-001 cell barcode", R5_001, "cbc", feature=r5_feature()),
    seg("S5", nx.S5, "s5"), seg("left mosaic end", nx.ME, "me"),
    seg("target-proximal cDNA", "X" * 38, placeholder=True),
    seg("right mosaic end reverse complement", nx.ME_RC, "me"),
    seg("S7 reverse complement", nx.S7_RC, "s7"),
    seg("PCR-well i7 reverse complement", "I" * 8, "cbc", placeholder=True,
        feature=i7_feature()),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="sciTIP-seq r5001/S502/N701 library")

SEQ_PRIMERS = (
    sp.custom("Read 1", "nXTv2 Read 1", NXT_R1, "Bartlett et al. Table S1"),
    sp.custom("Index 1 (i7)", "nXTv2 i7 index read", NXT_I7,
              "Bartlett et al. Table S1"),
    sp.custom("Index 2 (i5)", "nXTv2 i5 index read", NXT_I5,
              "Bartlett et al. Table S1"),
    sp.NEXTERA["R2"],
)
_problems = sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                      required_roles=tuple(p.role for p in SEQ_PRIMERS))
if _problems:
    raise ValueError("TIP-seq final library: " + "; ".join(_problems))

FINAL_CAPTION = ("Representative sciTIP molecule: r5-001 records the targeted-Tn5 well; "
                 "the S502 i5 and N701 i7 indexes record the second combinatorial well.")
SEQUENCING_INTRO = ("The published modified run acquires 51-cycle paired genomic reads, "
                    "8-cycle i7 and 29-cycle i5+r5 structure.")
READ_LENGTHS = {"Read 1": 51, "Index 1 (i7)": 8,
                "Index 2 (i5)": 29, "Read 2": 51}
SHOW_READ_CONTENT = True


def oligos():
    yield oligo("ME-scT7-r5001", me_t7_parts())
    yield oligo("ME-rev", [Segment("ME reverse", ME_REV, "me")], mods="/5Phos/")
    yield oligo("ME-B", [seg("S7", nx.S7, "s7"), seg("ME", nx.ME, "me")])
    yield oligo("sss_sci-nXTv2", [seg("second-strand primer", SSS_SCI, "r1")])
    yield oligo("sciTIP S502", [seg("P5", il.P5, "p5"),
                                 seg("i5", "CTCTCTAT", "cbc"),
                                 seg("annealing region", SSS_SCI, "r1")])
    yield oligo("N701", [seg("P7", il.P7, "p7"), seg("i7", "TCGCCTTA", "cbc"),
                          seg("annealing region", nx.S7, "s7")])


def _rna_rows():
    sc = Scene(); sc.strand("RNA", RNA_PRODUCT, label="T7 RNA copy")
    sc.mark("RNA", "transcribed r5-001 barcode", "cell barcode retained in every transcript")
    return sc.rows()


def _random_hexamer_rt_rows():
    sc = Scene(); sc.strand("RNA", RNA_PRODUCT, label="one T7 RNA copy")
    sc.anneal("random hexamer", [seg("random N6", "N" * 6, placeholder=True)],
              to="RNA", pair=("random N6", "target-proximal RNA copy"), shift=24,
              unpaired=("random N6",), label="random-hexamer RT primer")
    sc.arrow("random hexamer", "MMLV reverse transcription")
    return sc.rows()


def _second_strand_rows():
    sc = Scene.duplex(RNA_PRODUCT, label="RNA/cDNA after random-hexamer RT")
    sc.anneal("sss_sci-nXTv2", [seg("sss_sci-nXTv2", SSS_SCI, "r1")],
              to="bottom", pair=("sss_sci-nXTv2", "transcribed T7 start + second-strand site'"),
              label="published second-strand primer", above=True)
    sc.arrow("sss_sci-nXTv2", "Taq second-strand synthesis")
    return sc.rows()


def workflow():
    start = Scene.duplex([seg("antibody-targeted chromatin", "X" * 38,
                              placeholder=True)], label="antibody-bound chromatin")
    w = Workflow(MolecularState("Antibody-bound chromatin", tuple(start.rows())))
    w.react("Anneal ME-scT7-r5001 to ME-rev and load pA–Tn5",
            loaded_me_t7().rows(), name="barcoded promoter transposome")
    w.react("Tagment target-proximal DNA and fill the staggered gap",
            Scene.duplex(list(TAGGED), label="promoter-tagged chromatin").rows(),
            name="ME-scT7-tagged genomic DNA")
    w.react("Transcribe outward from the inserted T7 promoter",
            _rna_rows(), name="linearly amplified RNA")
    w.react("Prime first-strand cDNA with a random hexamer",
            _random_hexamer_rt_rows(), name="RNA with RT primer")
    w.react("Reverse-transcribe, remove RNA and anneal sss_sci-nXTv2",
            _second_strand_rows(), name="second-strand priming state")
    w.react("Synthesize the second strand",
            Scene.duplex(RNA_PRODUCT, label="double-stranded cDNA").rows(),
            name="double-stranded cDNA")
    w.react("Anneal ME-B to ME-rev and load a MeB-only Tn5",
            loaded_me_b().rows(), name="MeB-only transposome")
    w.react("Fragment cDNA and install the opposite MeB end",
            Scene.duplex(list(PRE_PCR), label="PCR-ready sciTIP cDNA").rows(),
            name="twice-tagmented cDNA")
    geometry = PCRGeometry.exact(
        PRE_PCR, PCRPrimer("sciTIP S502", S502, anneal_nt=len(SSS_SCI)),
        PCRPrimer("N701", N701, anneal_nt=len(nx.S7)),
        statement="Table S1 prints both representative index primers and their annealing regions.")
    w.pcr("Amplify with the S502 and N701 indexed primers", geometry,
          Scene.duplex(list(FINAL_LIBRARY), label="indexed sciTIP library").rows(),
          name="sciTIP sequencing library")
    return w
