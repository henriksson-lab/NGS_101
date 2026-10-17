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
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from tagmentation import loaded_tn5_adapter

TITLE = "TIP-seq — targeted insertion of promoters"
NOTES = "01_tip_seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1083/jcb.202103078">'
          'Bartlett et al., <i>Journal of Cell Biology</i> (2021)</a>.')
SUMMARY = ("Antibody-tethered pA–Tn5 inserts a T7 promoter beside targeted chromatin. "
           "T7 transcription linearly amplifies that one-ended insertion before random-hexamer "
           "reverse transcription, second-strand synthesis, MeB tagmentation and indexed PCR.")
CAVEAT = ("The workflow draws the published representative sciTIP r5001/S502/N701 branch. "
          "Table S1 supplies 384 r5 transposons and the full i5/i7 PCR-index families. "
          "Bulk TIP-seq (ME-T7-bulk, sss_bulk-nXT, standard Nextera index PCR) is not drawn. "
          "Published discrepancy in the Index 2 read: Table S1 lists nXTv2_i5_index_read as the "
          "P5 sequence AATGATACGGCGACCACCGAGATCTACAC, a forward-strand primer that would read i5 "
          "first (cycles 1-8) and r5 last (24-29). The Methods sciTIP-seq library note instead "
          "states r5 at index-2 positions 1-6 and the 8-bp i5 at 21-29, separated by a 15-bp "
          "spacer, for the NovaSeq 6000 SP v1.5 run. The page models the acquired run with a "
          "reverse-complement-workflow Index 2 primer, which places r5 at 1-6 and i5 at 22-29; "
          "the printed 21-29 span (9 positions for an 8-bp index) and the role of the Table S1 "
          "P5-sequence primer remain unexplained. "
          "SPRI cleanups, GuHCl Tn5 stripping and the final 0.85x size selection are not drawn.")

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
    sp.custom("Index 2 (i5)", "reverse-complement-workflow Index 2 primer (NovaSeq v1.5)",
              nx.INDEX2_PRIMER, "Bartlett et al. Methods sciTIP-seq: r5 at index-2 cycles 1-6",
              "Table S1 instead prints the P5-sequence nXTv2_i5_index_read; see the "
              "published-discrepancy notice."),
    sp.NEXTERA["R2"],
)
_problems = sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                      required_roles=tuple(p.role for p in SEQ_PRIMERS))
if _problems:
    raise ValueError("TIP-seq final library: " + "; ".join(_problems))

FINAL_CAPTION = ("Representative sciTIP molecule: r5-001 records the targeted-Tn5 well; "
                 "the S502 i5 and N701 i7 indexes record the second combinatorial well.")
SEQUENCING_INTRO = ("The published modified NovaSeq 6000 SP v1.5 run (Table S2 RunInfo) "
                    "acquires 51-cycle paired genomic reads, an 8-cycle i7 read and a "
                    "29-cycle index 2 read covering r5, the 15-bp spacer and i5.")
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
    yield oligo("nXTv2_i5_index_read (Table S1; P5 sequence, forward-strand)",
                [seg("P5", NXT_I5, "p5")])


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


def protocol_graph():
    source = ("Bartlett et al., Journal of Cell Biology (2021), Table S1; "
              "doi:10.1083/jcb.202103078")
    graph = ProtocolGraph("TIP-seq / sciTIP-seq", authority=source)
    chromatin = graph.input(
        "target_chromatin", "antibody-bound target chromatin",
        MaterialKind.BIOLOGICAL,
        state=MolecularState("Antibody-bound chromatin", tuple(
            Scene.duplex([seg("antibody-targeted chromatin", "X" * 38,
                              placeholder=True)], label="antibody-bound chromatin").rows())),
        primary=True)
    promoter_transposome = graph.input(
        "promoter_transposome", "loaded ME-scT7-r5001/ME-rev pA–Tn5 transposome",
        MaterialKind.ADAPTER, state=loaded_me_t7())
    pa_tn5 = graph.input("pa_tn5", "antibody-tethered pA–Tn5 transposase",
                         MaterialKind.CATALYST)
    tagged = graph.product(
        "promoter_tagged", "ME-scT7-tagged target-proximal genomic DNA",
        MaterialKind.MOLECULAR,
        state=Scene.duplex(list(TAGGED), label="promoter-tagged chromatin"))
    first_losses = graph.product("first_tn5_losses", "non-library first-Tn5 products",
                                 MaterialKind.WASTE)
    graph.react("promoter_tagmentation", "Tagment target-proximal DNA with promoter-loaded pA–Tn5",
                StepKind.TAGMENTATION,
                ((chromatin, InputRole.TEMPLATE), (promoter_transposome, InputRole.STRUCTURAL),
                 (pa_tn5, InputRole.CATALYST)), (tagged, first_losses))
    graph.discard(first_losses, "does not retain the promoter-bearing target-proximal insertion")
    fill = graph.input("gap_fill_polymerase", "Taq 5X Master Mix (72 C, 3 min gap fill)",
                       MaterialKind.CATALYST)
    filled = graph.product("filled_tagged_dna", "gap-filled promoter-tagged genomic DNA",
                           MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("gap_fill", "Fill the staggered Tn5 gap", StepKind.EXTENSION,
                ((tagged, InputRole.TEMPLATE), (fill, InputRole.CATALYST)), (filled,))
    t7 = graph.input("t7_polymerase", "HiScribe T7 high-yield kit T7 RNA polymerase mix",
                     MaterialKind.CATALYST)
    rna = graph.product("t7_rna", "linearly amplified target-proximal RNA",
                        MaterialKind.MOLECULAR,
                        state=MolecularState("T7 RNA copy", tuple(_rna_rows())))
    graph.react("t7_transcription", "Transcribe outward from the inserted T7 promoter",
                StepKind.TRANSCRIPTION,
                ((filled, InputRole.TEMPLATE), (t7, InputRole.CATALYST)), (rna,))
    random_hexamer = graph.input("random_hexamer", "random-hexamer RT primer",
                                 MaterialKind.OLIGO)
    mmlv = graph.input("mmlv", "SMART MMLV reverse transcriptase", MaterialKind.CATALYST)
    first_cdna = graph.product("first_cdna", "random-hexamer first-strand cDNA",
                               MaterialKind.MOLECULAR,
                               state=MolecularState("RNA/cDNA", tuple(
                                   Scene.duplex(RNA_PRODUCT, label="RNA/cDNA after RT").rows())))
    graph.react("reverse_transcription", "Prime and reverse-transcribe T7 RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna, InputRole.TEMPLATE), (random_hexamer, InputRole.STRUCTURAL),
                 (mmlv, InputRole.CATALYST)), (first_cdna,))
    hydrolysis = graph.input("rna_removal", "RNase H (diluted), 37 C 20 min",
                             MaterialKind.CATALYST)
    clean_cdna = graph.product("clean_cdna", "RNA-free first-strand cDNA",
                               MaterialKind.MOLECULAR, state=first_cdna.state)
    degraded_rna = graph.product("degraded_rna", "removed RNA fragments", MaterialKind.WASTE)
    graph.react("remove_rna", "Degrade RNA in RNA–cDNA hybrids with RNase H", StepKind.DIGESTION,
                ((first_cdna, InputRole.TEMPLATE), (hydrolysis, InputRole.CATALYST)),
                (clean_cdna, degraded_rna))
    graph.discard(degraded_rna, "the cDNA strand is retained for second-strand synthesis")
    second_primer = graph.input("second_strand_primer", "published sss_sci-nXTv2 primer",
                                MaterialKind.OLIGO)
    taq = graph.input("second_strand_polymerase", "Taq DNA polymerase",
                      MaterialKind.CATALYST)
    double = graph.product("double_stranded_cdna", "double-stranded target-proximal cDNA",
                           MaterialKind.MOLECULAR,
                           state=Scene.duplex(RNA_PRODUCT, label="double-stranded cDNA"))
    graph.react("second_strand", "Synthesize the second cDNA strand",
                StepKind.EXTENSION,
                ((clean_cdna, InputRole.TEMPLATE), (second_primer, InputRole.STRUCTURAL),
                 (taq, InputRole.CATALYST)), (double,))
    meb = graph.input("meb_transposome", "loaded MeB-only Tn5 transposome",
                      MaterialKind.ADAPTER, state=loaded_me_b())
    tn5 = graph.input("second_tn5", "MeB-only Tn5 transposase", MaterialKind.CATALYST)
    pre_pcr = graph.product("pre_pcr", "twice-tagmented sciTIP cDNA",
                            MaterialKind.MOLECULAR, state=PRE_PCR)
    second_losses = graph.product("second_tn5_losses", "non-amplifiable second-Tn5 products",
                                  MaterialKind.WASTE)
    graph.react("second_tagmentation", "Fragment cDNA and install the opposite MeB end",
                StepKind.TAGMENTATION,
                ((double, InputRole.TEMPLATE), (meb, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (pre_pcr, second_losses))
    graph.discard(second_losses, "does not retain the amplifiable two-ended sciTIP configuration")
    s502 = graph.input("s502", "published sciTIP S502 indexed primer", MaterialKind.OLIGO)
    n701 = graph.input("n701", "published N701 indexed primer", MaterialKind.OLIGO)
    pcr_polymerase = graph.input("pcr_polymerase", "NEBNext High-Fidelity 2X PCR Master Mix",
                                 MaterialKind.CATALYST)
    library = graph.product("sequencing_library", "sciTIP sequencing library",
                            MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    geometry = PCRGeometry.exact(
        PRE_PCR, PCRPrimer("sciTIP S502", S502, anneal_nt=len(SSS_SCI)),
        PCRPrimer("N701", N701, anneal_nt=len(nx.S7)),
        statement="Table S1 prints both representative index primers and annealing regions.")
    graph.react("indexed_pcr", "Amplify with S502 and N701 indexed primers",
                StepKind.PCR,
                ((pre_pcr, InputRole.TEMPLATE), (s502, InputRole.STRUCTURAL),
                 (n701, InputRole.STRUCTURAL), (pcr_polymerase, InputRole.CATALYST)),
                (library,), geometry=geometry)
    graph.sequence(library, SequencingPlan(
        FINAL_LIBRARY, tuple(SEQ_PRIMERS), tuple(p.role for p in SEQ_PRIMERS),
        READ_LENGTHS,
        {"Read 1": "include", "Index 1 (i7)": "demultiplex",
         "Index 2 (i5)": "extract_and_demultiplex", "Read 2": "include"}))
    graph.declare_reviewed_scope((
        "target_chromatin", "promoter_transposome", "pa_tn5", "promoter_tagmentation",
        "gap_fill_polymerase", "gap_fill", "t7_polymerase", "t7_transcription",
        "random_hexamer", "mmlv", "reverse_transcription", "rna_removal", "remove_rna",
        "second_strand_primer", "second_strand_polymerase", "second_strand",
        "meb_transposome", "second_tn5", "second_tagmentation", "s502", "n701",
        "pcr_polymerase", "indexed_pcr", "sequencing_library",
    ), source)
    return graph.require_publishable()
