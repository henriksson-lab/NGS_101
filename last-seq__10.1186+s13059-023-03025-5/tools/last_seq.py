"""Molecular model for LAST-seq (Lyu and Chen, 2023)."""
from __future__ import annotations

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import (Construct, MolecularState, Row, Scene, Segment, Workflow,
                      annotation_rows, feature, oligo)
from rna_amplification import direct_ssrna_ivt, poly_a_capture
from pcr import PCRGeometry, PCRPrimer
from chemdraw import revcomp
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role
from topology import HairpinTopology


def seg(name, bases, tag=None, **kw):
    return Segment(name, bases, tag, **kw)


T7_PROMOTER = "TAATACGACTCACTATAGGG"
RANDOM_RT_HANDLE = "TACACGACGCTCTTCCGATCT"
RANDOM_RT = RANDOM_RT_HANDLE + "NNNNNN"
DTRU_LINKER = "CCACCTTTCATTCACCCT"
DTRU_DNA = "TTTTTT"
DTRU_RNA = "U" * 19
DTRU_TAIL = DTRU_LINKER + DTRU_DNA + DTRU_RNA
HAIRPIN_LINKER_ARM = "GGGTGAATGAAAGGTGG"
HAIRPIN_CORE = "CTCTAATACGACTCACTATAGGGAGACGTGTGCTCTTCCGATCT"
HAIRPIN_CLOSING_ARM = "TCGTATTAGAG"
CELL_BARCODE_LEN = 6
UMI_LEN = 8
BC2 = "GTATCC"
BC2_HAIRPIN = HAIRPIN_LINKER_ARM + HAIRPIN_CORE + "N" * UMI_LEN + BC2 + HAIRPIN_CLOSING_ARM
BC2_HAIRPIN_TOPOLOGY = HairpinTopology(
    "LAST-seq BC2 hairpin module", BC2_HAIRPIN,
    first_stem_start=len(HAIRPIN_LINKER_ARM),
    second_stem_start=len(BC2_HAIRPIN) - len(HAIRPIN_CLOSING_ARM),
    stem_length=len(HAIRPIN_CLOSING_ARM), colour="var(--c-cbc)")

# Additional file 2, Table S3 library primers (NEBNext1 shown; index window CGTGAT).
NEBNEXT_UNIVERSAL = il.NEBNEXT_UNIVERSAL_PRIMER          # printed with a 3' phosphorothioate
NEBNEXT_INDEX_5 = "CAAGCAGAAGACGGCATACGAGAT"
NEBNEXT_INDEX_3 = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"  # printed ...GATC*T
NEBNEXT_INDEX_LEN = 6
LIB_P5 = "AATGATACGGCGACCACCGAGATCTACAC"
LIB_P7 = "CAAGCAGAAGACGGCATACGAGAT"


def _overlap_3prime_site(primer: str, site: str) -> int:
    """Longest primer 3' end found verbatim in ``site``."""
    for k in range(len(primer), 0, -1):
        if primer[-k:] in site:
            return k
    return 0


# Computed annealing of the Table S3 primers on the RT-product ends: the NEBNext
# primers have 21-nt 3' sites on the RandomRT handle / hairpin PCR handle, whereas the
# short P5/P7 primers have only 5 and 3 nt there.  Methods name "P5/P7 primers" for
# the first 8-cycle PCR; see CAVEAT.
NEB_SITE_NT = (_overlap_3prime_site(NEBNEXT_UNIVERSAL, RANDOM_RT_HANDLE),
               _overlap_3prime_site(NEBNEXT_INDEX_3, HAIRPIN_CORE))
SHORT_SITE_NT = (_overlap_3prime_site(LIB_P5, RANDOM_RT_HANDLE),
                 _overlap_3prime_site(LIB_P7, HAIRPIN_CORE))
if min(NEB_SITE_NT) < 18 or max(SHORT_SITE_NT) >= 10:
    raise ValueError("LAST-seq PCR primer site computation changed: %r %r"
                     % (NEB_SITE_NT, SHORT_SITE_NT))

CELL_FEATURE = feature("last_cell_barcode", "cell_barcode", "whitelist",
                       whitelist="Supplementary Table S3")
UMI_FEATURE = feature("last_umi", "umi", "random")


def capture_primer():
    return [
        seg("hairpin module", "[hairpin: PCR handle + T7 promoter + UMI + cell barcode]",
            "r2", placeholder=True),
        seg("complementary linker", DTRU_LINKER, "tso"),
        seg("dT7", DTRU_DNA),
        seg("rU19", DTRU_RNA),
    ]


def capture_scene():
    sc = poly_a_capture(capture_primer(), annealing_segment="rU19",
                        primer_label="3′-ddC-blocked LAST-seq primer")
    sc.junction("primer", "hairpin module", "complementary linker", "primer assembly ligation")
    sc.note("primer", "the dT/rU tract is nicked by RNase H; 3′ ddC blocks extension from the primer end")
    return sc


def primer_ligation_geometry():
    sc = Scene(); sc.strand("hairpin", [seg("hairpin linker arm", HAIRPIN_LINKER_ARM),
                                         seg("folded barcode/promoter hairpin", "[folded hairpin]", placeholder=True)],
                                    label="hairpin module")
    sc.anneal("rU-dT module", [seg("phosphorylated paired linker", DTRU_LINKER[:-1]),
                                seg("linker 3-prime dT overhang", DTRU_LINKER[-1])],
              to="hairpin", pair=("phosphorylated paired linker", "hairpin linker arm"),
              label="rU-dT module", mod5="p")
    sc.ligation_end("rU-dT module", "phosphorylated paired linker", kind="splinted",
                    polarity="5-prime", owner="adapter", text="5′-phosphorylated nick end")
    sc.overhang("rU-dT module", "linker 3-prime dT overhang", polarity="3-prime",
                owner="adapter", text="one unpaired dT")
    return AdapterLigationGeometry.known(
        "splinted", sc,
        statement="Supplementary Table S2 gives a 17-base hairpin arm complementary to the first 17 bases of the 18-base rU-dT linker; E. coli DNA ligase seals their annealed nick and one linker dT remains unpaired.")


def promoter_scene():
    promoter = seg("T7 promoter", T7_PROMOTER, "t7")
    template = [seg("cell barcode", "B" * CELL_BARCODE_LEN, "cbc", placeholder=True,
                    feature=CELL_FEATURE),
                seg("UMI", "U" * UMI_LEN, "umi", placeholder=True, feature=UMI_FEATURE),
                seg("mRNA template", "X" * 34, placeholder=True)]
    return direct_ssrna_ivt(template, promoter)


def arna_scene():
    sc = Scene(); sc.strand("aRNA", [
        seg("aRNA transcript", "X" * 34, placeholder=True),
        seg("copied UMI", "U" * UMI_LEN, "umi", placeholder=True, feature=UMI_FEATURE),
        seg("copied cell barcode", "B" * CELL_BARCODE_LEN, "cbc", placeholder=True,
            feature=CELL_FEATURE)], label="antisense RNA copies")
    sc.note("aRNA", "hundreds of linear-amplification products from one original ssRNA")
    return sc


def random_rt_scene():
    arna = [seg("aRNA body", "X" * 34, placeholder=True)]
    primer = [seg("Read 1 handle", RANDOM_RT_HANDLE, "r1"),
              seg("random N6", "N" * 6, placeholder=True)]
    sc = Scene(); sc.strand("aRNA", arna, label="amplified antisense RNA")
    sc.anneal("random RT primer", primer, to="aRNA", pair=("random N6", "aRNA body"),
              shift=len(arna[0]) - 6, label="random RT primer",
              unpaired=("Read 1 handle", "random N6", "aRNA body"))
    sc.arrow("random RT primer", "SuperScript IV")
    sc.labels("random RT primer")
    return sc


FINAL_LIBRARY = Construct([
    seg("P5", LIB_P5, "p5"),
    seg("Read 1 arm", NEBNEXT_UNIVERSAL[len(LIB_P5):], "r1"),
    seg("random-primer N6", "N" * 6, placeholder=True, length_bp=6),
    seg("aRNA-derived cDNA", "X" * 38, placeholder=True),
    seg("mRNA poly(A) tract (variable)", "A" * 8, placeholder=True),
    seg("rU-dT linker complement", revcomp(DTRU_LINKER), "tso"),
    seg("hairpin closing-arm complement", revcomp(HAIRPIN_CLOSING_ARM)),
    seg("cell barcode", "B" * CELL_BARCODE_LEN, "cbc", placeholder=True,
        feature=CELL_FEATURE),
    seg("UMI", "U" * UMI_LEN, "umi", placeholder=True, feature=UMI_FEATURE),
    seg("Read 2 arm (NEBNext index primer)", revcomp(NEBNEXT_INDEX_3), "r2"),
    seg("i7 reverse complement", "I" * NEBNEXT_INDEX_LEN, "cbc", placeholder=True,
        feature=feature("last_sample_index", "sample_index", "whitelist",
                        whitelist="NEBNext1-9 and NEBNext12 index primers, Table S3")),
    seg("P7 reverse complement", revcomp(LIB_P7), "p7"),
], name="LAST-seq library")

SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
_errors = sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                    required_roles=tuple(p.role for p in SEQ_PRIMERS))
if _errors:
    raise ValueError("LAST-seq final library: " + "; ".join(_errors))


def final_rows():
    sc = Scene.duplex(list(FINAL_LIBRARY), label="LAST-seq")
    sc.labels("top")
    return [*sc.rows(), *annotation_rows(FINAL_LIBRARY)]


def workflow():
    wf = Workflow(MolecularState("Separate hairpin and rU-dT primer modules",
                                 (Row(chunks=[("hairpin module + 5′-phosphorylated rU-dT module", None, False)]),)))
    wf.ligate_adapter("Anneal and nick-ligate the LAST-seq capture primer",
                      primer_ligation_geometry(), capture_scene().rows(),
                      name="assembled hairpin/rU-dT primer",
                      note="The gel-purified product is then annealed to mRNA poly(A).")
    wf.react("RNase H nicking and Klenow exo− extension", promoter_scene().rows(),
             note="A short rA/dT hybrid is nicked and extended so the original RNA gains a duplex T7 promoter without reverse transcription.")
    wf.react("T7 in-vitro transcription directly from ssRNA", arna_scene().rows(),
             note="The original single-stranded RNA is the transcription template; amplification is linear.")
    wf.react("Random-prime the pooled aRNA and reverse-transcribe", random_rt_scene().rows(),
             note="Reverse transcription occurs only after direct RNA amplification.")
    wf.react("INFERRED — 8-cycle Q5 PCR with NEBNext primers, then further P5/P7 PCR", final_rows(),
             note=("Methods name 'P5/P7 primers' for the 8-cycle PCR; only the Table S3 NEBNext "
                   "primers have full sites on the RT-product ends (see caveat). PCR completes the "
                   "Illumina library; the inline cell barcode and UMI originate in the capture hairpin."))
    return wf


TITLE = "LAST-seq — direct linear amplification of single-stranded RNA"
NOTES = "01_last-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1186/s13059-023-03025-5">'
          'Lyu and Chen, <i>Genome Biology</i> (2023)</a>.')
SUMMARY = ("A ligated hairpin/rU-dT primer installs a short duplex T7 promoter at an mRNA’s "
           "3′ end. T7 polymerase then amplifies the original RNA directly, before any reverse transcription.")
CAVEAT = ("Source ambiguity in library PCR: Methods 'Library preparation of the amplified "
          "products' say the RT products were amplified for 8 cycles with Q5 and 'P5/P7 "
          "primers', then 'further PCR amplified'. Additional file 2, Table S3 prints both "
          "short P5 (29 nt) and P7 (24 nt) primers and NEBNext Universal/Index primers. "
          "Computed from those sequences, the short P5/P7 3′ ends match only %d and %d nt of "
          "the RandomRT handle and hairpin PCR handle, whereas the NEBNext primers have "
          "%d- and %d-nt sites. The model therefore assigns the first PCR to the NEBNext "
          "primers and the unnamed second PCR to short P5/P7 re-amplification; both "
          "assignments are inferred, and the second PCR's cycle number is not stated. "
          % (SHORT_SITE_NT + NEB_SITE_NT) +
          "Published barcode collision: Additional file 2, Table S3, BC5 and BC7 rows "
          "(extracted lines 48 and 50) both print CAGGAA, although the article's LAST-seq "
          "assay section describes 80 selected cellular barcodes. The model displays BC2 "
          "(GTATCC), a separate printed member, as an example and does not silently "
          "assign a replacement to BC5 or BC7. Which entry was intended to differ "
          "remains unresolved.")
FINAL_CAPTION = ("Representative PCR-completed LAST-seq molecule. Read 2 reads the hairpin-derived "
                 "UMI, then the cell barcode, then the 11-nt hairpin closing arm; Read 1 enters "
                 "the randomly primed aRNA-derived insert. The P7 side is drawn from the "
                 "NEBNext index primer (6-nt index) printed in Table S3; the poly(A) length is variable.")
SEQUENCING_INTRO = ("NextSeq 550, High Output v2.5 kit: Read 1, 50 cycles; Read 2, 25 cycles. "
                    "The paper reports 6-nt cellular barcodes and 8-nt UMIs; index-read cycles "
                    "are not stated.")
READ_LENGTHS = {"Read 1": 50, "Read 2": 25}
SHOW_READ_CONTENT = True


def oligos():
    yield oligo("rU-dT module", [seg("linker", DTRU_LINKER, "tso"),
                                  seg("dT7", DTRU_DNA), seg("rU19", DTRU_RNA)],
                mods="5′ phosphate; 3′ ddC")
    yield oligo("NEBNext Universal (Table S3)",
                [seg("P5", LIB_P5, "p5"),
                 seg("Read 1 arm", NEBNEXT_UNIVERSAL[len(LIB_P5):], "r1")],
                mods="3′-terminal phosphorothioate (printed C*T)")
    yield oligo("NEBNext index primer (Table S3)",
                [seg("P7", NEBNEXT_INDEX_5, "p7"),
                 seg("i7 index", "I" * NEBNEXT_INDEX_LEN, "cbc", placeholder=True),
                 seg("Read 2 arm", NEBNEXT_INDEX_3, "r2")],
                mods="3′-terminal phosphorothioate (printed C*T)")
    yield oligo("P5 (Table S3)", [seg("P5", LIB_P5, "p5")])
    yield oligo("P7 (Table S3)", [seg("P7", LIB_P7, "p7")])
    yield oligo("RandomRT primer6", [seg("Read 1 handle", RANDOM_RT_HANDLE, "r1"),
                                      seg("N6", "N" * 6, placeholder=True)])


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    rna = graph.input("rna", "polyadenylated cellular RNA", MaterialKind.BIOLOGICAL,
                      state=molecule("mRNA", (role("transcript body"), role("poly(A) tail"))),
                      primary=True)
    hairpin = graph.input("hairpin", "barcode/UMI/T7-promoter hairpin module",
                          MaterialKind.MOLECULAR)
    rudt = graph.input("rudt", "5-prime-phosphorylated rU-dT module",
                       MaterialKind.ADAPTER)
    ligase = graph.input("primer_ligase", "E. coli DNA ligase", MaterialKind.CATALYST)
    assembled = graph.product("capture_primer", "assembled hairpin/rU-dT capture primer",
                              MaterialKind.MOLECULAR, state=capture_scene())
    lig_condition = graph.input("ligation_condition",
                                "self-fold hairpin by heating/cooling; 16 °C, 30 min, 20 µl",
                                MaterialKind.CONDITION)
    ligated = graph.product("ligated_primer", "nick-ligated hairpin/rU-dT primer",
                            MaterialKind.MOLECULAR, state=capture_scene())
    graph.react("assemble_primer", "Nick-ligate the LAST-seq capture primer",
                StepKind.LIGATION,
                ((hairpin, InputRole.TEMPLATE), (rudt, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST), (lig_condition, InputRole.CONDITION)),
                (ligated,), geometry=primer_ligation_geometry())
    page_gel = graph.input("urea_page", "10% urea PAGE, gel soak in NaOAc, ethanol precipitation",
                           MaterialKind.CONDITION)
    unligated = graph.product("unligated_modules", "unligated hairpin and rU-dT modules",
                              MaterialKind.WASTE)
    graph.react("gel_purify_primer", "Gel-purify the ligated primer", StepKind.CLEANUP,
                ((ligated, InputRole.TEMPLATE), (page_gel, InputRole.CONDITION)),
                (assembled, unligated))
    graph.discard(unligated, "outside the excised ligation-product gel band")
    lysis = graph.input("lysis", "sort into 2 µl lysis buffer (0.1% Triton X-100, 2 nM "
                        "primer, dNTP, RNA carrier, spike-in); 65 °C 5 min, then ice",
                        MaterialKind.CONDITION)
    hybrid = graph.product("captured_rna", "capture primer annealed to mRNA poly(A)",
                           MaterialKind.MOLECULAR, state=capture_scene())
    graph.react("capture", "Lyse single cells and anneal the primer to mRNA poly(A)",
                StepKind.HYBRIDIZATION,
                ((rna, InputRole.TEMPLATE), (assembled, InputRole.STRUCTURAL),
                 (lysis, InputRole.CONDITION)), (hybrid,))
    rnase = graph.input("rnase_h", "RNase H", MaterialKind.CATALYST)
    klenow = graph.input("klenow", "Klenow exo-minus polymerase", MaterialKind.CATALYST)
    nick_condition = graph.input("nick_condition", "evaporate lysate; 1X NEBuffer2, "
                                 "0.02 U/µl RNase H, 0.3 U/µl Klenow exo−, 37 °C 30 min",
                                 MaterialKind.CONDITION)
    promoter_rna = graph.product("promoter_rna", "ssRNA carrying a duplex T7 promoter",
                                 MaterialKind.MOLECULAR, state=promoter_scene())
    graph.react("install_promoter", "Nick the rA/dT tract and extend the duplex promoter",
                StepKind.EXTENSION,
                ((hybrid, InputRole.TEMPLATE), (rnase, InputRole.CATALYST),
                 (klenow, InputRole.CATALYST), (nick_condition, InputRole.CONDITION)),
                (promoter_rna,))
    beads_05 = graph.input("beads_05x", "pool 10 or 80 lysates; 0.5X beads + 1.3X buffer",
                           MaterialKind.CONDITION)
    bound = graph.product("bead_bound_rna", "pooled promoter-bearing RNA on beads",
                          MaterialKind.MOLECULAR, state=promoter_scene())
    waste_05 = graph.product("pool_cleanup_waste", "material not retained by 0.5X beads",
                             MaterialKind.WASTE)
    graph.react("pool_cleanup", "Pool cell lysates and bead-purify", StepKind.CLEANUP,
                ((promoter_rna, InputRole.TEMPLATE), (beads_05, InputRole.CONDITION)),
                (bound, waste_05))
    graph.discard(waste_05, "not retained by the 0.5X bead cleanup")
    ivt_condition = graph.input("ivt_condition", "beads rehydrated in 25 µl T7 IVT mix; "
                                "37 °C overnight, up to 16 h", MaterialKind.CONDITION)
    t7 = graph.input("t7_polymerase", "T7 RNA polymerase", MaterialKind.CATALYST)
    arna = graph.product("arna", "linearly amplified antisense RNA copies",
                         MaterialKind.MOLECULAR, state=arna_scene())
    graph.react("transcribe", "Transcribe antisense RNA directly from the original ssRNA",
                StepKind.TRANSCRIPTION,
                ((bound, InputRole.TEMPLATE), (t7, InputRole.CATALYST),
                 (ivt_condition, InputRole.CONDITION)), (arna,))
    beads_08 = graph.input("beads_08x", "0.8X beads + 1X buffer; elute in 7 µl",
                           MaterialKind.CONDITION)
    arna_clean = graph.product("arna_clean", "purified aRNA", MaterialKind.MOLECULAR,
                               state=arna_scene())
    waste_08 = graph.product("arna_cleanup_waste", "material not retained by 0.8X beads",
                             MaterialKind.WASTE)
    graph.react("arna_cleanup", "Bead-purify aRNA", StepKind.CLEANUP,
                ((arna, InputRole.TEMPLATE), (beads_08, InputRole.CONDITION)),
                (arna_clean, waste_08))
    graph.discard(waste_08, "not retained by the 0.8X bead cleanup")
    rt_condition = graph.input("rt_condition", "70 °C 2 min with dNTP and primer; "
                               "23 °C 10 min, 50 °C 15 min, 80 °C 10 min",
                               MaterialKind.CONDITION)
    random_primer = graph.input("random_rt_primer", "RandomRT primer6", MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "SuperScript IV", MaterialKind.CATALYST)
    cdna = graph.product("cdna", "random-primed cDNA from pooled aRNA",
                         MaterialKind.MOLECULAR, state=random_rt_scene())
    graph.react("reverse_transcribe", "Random-prime and reverse-transcribe pooled aRNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((arna_clean, InputRole.TEMPLATE), (random_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST), (rt_condition, InputRole.CONDITION)), (cdna,))
    forward = graph.input("pcr_forward", "NEBNext Universal primer (Table S3)",
                          MaterialKind.OLIGO, inferred=True,
                          source=SOURCE + " (assignment to the first PCR inferred; see caveat)")
    reverse = graph.input("pcr_reverse", "NEBNext index primer (Table S3)",
                          MaterialKind.OLIGO, inferred=True,
                          source=SOURCE + " (assignment to the first PCR inferred; see caveat)")
    polymerase = graph.input("pcr_polymerase", "Q5 High-Fidelity Master Mix",
                             MaterialKind.CATALYST)
    pcr1_condition = graph.input("pcr1_condition", "8 cycles", MaterialKind.CONDITION)
    amplified = graph.product("pcr1_product", "first-PCR library", MaterialKind.MOLECULAR,
                              state=FINAL_LIBRARY)
    graph.react("library_pcr", "Amplify RT products with indexing primers (8 cycles)",
                StepKind.PCR,
                ((cdna, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST),
                 (pcr1_condition, InputRole.CONDITION)),
                (amplified,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY, reason=("The pooled random-RT precursor is heterogeneous; "
                                           "computed NEBNext 3′ sites on the RandomRT handle "
                                           "and hairpin PCR handle are %d and %d nt."
                                           % NEB_SITE_NT)))
    beads_09 = graph.input("beads_09x", "0.9X beads", MaterialKind.CONDITION)
    pcr1_clean = graph.product("pcr1_clean", "purified first-PCR library",
                               MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    waste_09 = graph.product("pcr1_cleanup_waste", "primers and short products not retained",
                             MaterialKind.WASTE)
    graph.react("pcr1_cleanup", "Bead-purify the first PCR", StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (beads_09, InputRole.CONDITION)),
                (pcr1_clean, waste_09))
    graph.discard(waste_09, "not retained by the 0.9X bead cleanup")
    p5 = PCRPrimer("P5 (Table S3)", LIB_P5)
    p7 = PCRPrimer("P7 (Table S3)", LIB_P7)
    p5_in = graph.input("pcr2_p5", "P5 primer (Table S3)", MaterialKind.OLIGO,
                        inferred=True, source=SOURCE + " (second-PCR primers not named)")
    p7_in = graph.input("pcr2_p7", "P7 primer (Table S3)", MaterialKind.OLIGO,
                        inferred=True, source=SOURCE + " (second-PCR primers not named)")
    polymerase2 = graph.input("pcr2_polymerase", "second-PCR polymerase (not stated)",
                              MaterialKind.CATALYST)
    reamplified = graph.product("pcr2_product", "re-amplified library",
                                MaterialKind.MOLECULAR, state=FINAL_LIBRARY, inferred=True)
    graph.react("library_pcr2", "Further PCR amplification (primers and cycles not stated)",
                StepKind.PCR,
                ((pcr1_clean, InputRole.TEMPLATE), (p5_in, InputRole.STRUCTURAL),
                 (p7_in, InputRole.STRUCTURAL), (polymerase2, InputRole.CATALYST)),
                (reamplified,), geometry=PCRGeometry.exact(
                    FINAL_LIBRARY, p5, p7,
                    statement=("Inferred assignment: Table S3 'library enrichment' P5/P7 "
                               "primers on the completed library ends.")))
    beads_final = graph.input("beads_final", "bead purification; elute in 11 µl 0.1× TE",
                              MaterialKind.CONDITION)
    library = graph.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY)
    waste_final = graph.product("final_cleanup_waste", "material not retained by beads",
                                MaterialKind.WASTE)
    graph.react("final_cleanup", "Bead-purify the final library", StepKind.CLEANUP,
                ((reamplified, InputRole.TEMPLATE), (beads_final, InputRole.CONDITION)),
                (library, waste_final))
    graph.discard(waste_final, "not retained by the final bead cleanup")
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=tuple(SEQ_PRIMERS),
        read_lengths=READ_LENGTHS))
    graph.declare_reviewed_scope([
        "rna", "hairpin", "rudt", "primer_ligase", "ligation_condition",
        "assemble_primer", "urea_page", "gel_purify_primer", "lysis", "capture",
        "rnase_h", "klenow", "nick_condition", "install_promoter", "beads_05x",
        "pool_cleanup", "t7_polymerase", "ivt_condition", "transcribe", "beads_08x",
        "arna_cleanup", "random_rt_primer", "reverse_transcriptase", "rt_condition",
        "reverse_transcribe", "pcr_forward", "pcr_reverse", "pcr_polymerase",
        "pcr1_condition", "library_pcr", "beads_09x", "pcr1_cleanup", "pcr2_p5", "pcr2_p7",
        "pcr2_polymerase", "library_pcr2", "beads_final", "final_cleanup", "library"], SOURCE)
    return graph.require_publishable()
