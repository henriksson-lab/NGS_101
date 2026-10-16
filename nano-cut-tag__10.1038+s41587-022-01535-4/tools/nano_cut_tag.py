"""nano-CUT&Tag molecular model from Bartosovic et al. Supplementary Table 1."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature, oligo, revcomp
from pcr import PCRGeometry
from tagmentation import loaded_tn5_adapter
from protocol_graph import (AnnotationSlot, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind)
from protocol_model import molecule, role

TITLE = "nano-CUT&Tag — multimodal nanobody–Tn5 chromatin profiling"
NOTES = "01_nano-cut-tag.html"
SOURCE = '<a href="https://doi.org/10.1038/s41587-022-01535-4">Bartosovic et al., <i>Nature Biotechnology</i> (2023)</a>.'
SUMMARY = ("Nanobody–Tn5 fusions first transfer a target-specific MeA adapter. Chromium "
           "linear amplification adds cell identity; a second MeB-only Tn5 installs the "
           "opposite end before PCR.")
CAVEAT = ("The target-specific adapters and Tn5 loading geometry are source-defined. "
          "The linked 10x ATAC v1.1 guide prints gel-bead and SI-PCR oligo "
          "scaffolds; the 16-nt cell barcode and 8-nt sample index vary.")

OUTER = "TCCACGC"
BARCODE_A = "TATAGCCT"
BARCODE_B = "ATAGAGGC"
BARCODE_C = "CCTATCCT"
CONNECTOR = "GCGATCGAGGACGGC"
ME_REV = nx.ME_RC
CUSTOM_R1 = CONNECTOR + nx.ME
CUSTOM_I2 = nx.ME_RC + revcomp(CONNECTOR)


ATAC_GUIDE = "10x Chromium Next GEM Single Cell ATAC v1.1 guide CG000209 Rev G, Appendix A1"
GEL_BEAD_FORWARD = il.P5 + "N" * 16 + nx.S5
SI_PCR_FORWARD = "AATGATACGGCGACCACCGAGA"
SI_PCR_REVERSE = "CAAGCAGAAGACGGCATACGAGAT" + "N" * 8 + "GTCTCGTGGGCTCGG"


def modality_feature():
    return feature("nanoct_modality", "feature_barcode", "whitelist",
                   whitelist="Bartosovic et al. Supplementary Table 1")


def mea_parts(barcode=BARCODE_A, name="A", stagger=0):
    parts = [seg("S5", nx.S5, "s5"), seg("outer connector", OUTER),
             seg(f"Bcd{name} modality barcode", barcode, "cbc",
                 feature=modality_feature())]
    if stagger:
        parts.append(seg("sequencing stagger", "N" * stagger, placeholder=True))
    return [*parts, seg("inner connector", CONNECTOR),
            seg("MeA mosaic end", nx.ME, "me")]


def meb_parts():
    return [seg("S7", nx.S7, "s7"), seg("MeB mosaic end", nx.ME, "me")]


def loaded(parts, end, label):
    return loaded_tn5_adapter(parts, mosaic_end=end, label=label,
                              bottom_name="Tn5_Rev", bottom_label="5′-phosphorylated Tn5_Rev",
                              bottom_sequence=ME_REV, mod5="p")


FIRST_TAGGED = Construct([*mea_parts(), seg("targeted genomic DNA", "X" * 38,
                                            placeholder=True)],
                         name="BcdA-tagged target fragment")

PRE_PCR = Construct([
    seg("10x gel-bead P5", il.P5, "p5"),
    seg("10x cell barcode", "B" * 16, "cbc", placeholder=True,
        feature=feature("cell_barcode", "cell_barcode", "whitelist",
                        whitelist="10x Chromium ATAC whitelist")),
    *mea_parts(), seg("targeted genomic DNA", "X" * 38, placeholder=True),
    seg("MeB mosaic end reverse complement", nx.ME_RC, "me"),
    seg("S7 reverse complement", nx.S7_RC, "s7"),
], name="twice-tagmented nano-CT molecule before SI-PCR")

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"),
    seg("10x cell barcode", "B" * 16, "cbc", placeholder=True,
        feature=feature("cell_barcode", "cell_barcode", "whitelist",
                        whitelist="10x Chromium ATAC whitelist")),
    *mea_parts(), seg("targeted genomic DNA", "X" * 38, placeholder=True),
    seg("MeB mosaic end reverse complement", nx.ME_RC, "me"),
    seg("S7 reverse complement", nx.S7_RC, "s7"),
    seg("i7 sample index reverse complement", "I" * 8, "cbc", placeholder=True,
        feature=feature("sample_index_i7", "sample_index", "unknown")),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="nano-CUT&Tag BcdA library")

SEQ_PRIMERS = (
    sp.custom("Read 1", "nano-CT custom Read 1", CUSTOM_R1,
              "Bartosovic et al. Supplementary Table 1"),
    sp.NEXTERA["I1"],
    sp.custom("Index 2 (i5)", "nano-CT custom Index 2", CUSTOM_I2,
              "Bartosovic et al. Supplementary Table 1"),
    sp.NEXTERA["R2"],
)
_problems = sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                      required_roles=tuple(p.role for p in SEQ_PRIMERS))
if _problems:
    raise ValueError("nano-CT final library: " + "; ".join(_problems))

FINAL_CAPTION = ("Representative BcdA product. Custom R1 traverses MeA into genomic DNA; "
                 "custom Index 2 traverses the modality adapter toward the 10x cell barcode.")
SEQUENCING_INTRO = ("The source specifies a 36-8-48-36 R1-I1-I2-R2 run with custom R1 and "
                    "Index 2 primers; R2 and Index 1 are standard.")
READ_LENGTHS = {"Read 1": 36, "Index 1 (i7)": 8,
                "Index 2 (i5)": 48, "Read 2": 36}
SHOW_READ_CONTENT = True


def oligos():
    for name, barcode in (("A", BARCODE_A), ("B", BARCODE_B), ("C", BARCODE_C)):
        for stagger in range(4):
            yield oligo(f"Tn5_P5_MeA_Bcd{name}_{stagger}N",
                        mea_parts(barcode, name, stagger))
    yield oligo("Tn5_P7_MeB", meb_parts())
    yield oligo("Tn5_Rev", [Segment("ME reverse", ME_REV, "me")], mods="/5Phos/")
    yield oligo("R1_seq", [seg("binding region", CUSTOM_R1, "r1")])
    yield oligo("I2_seq", [seg("binding region", CUSTOM_I2, "r2")])

    yield oligo("10x gel-bead oligo (CG000209 A1)", [
        seg("P5", il.P5, "p5"),
        seg("10x cell barcode", "N" * 16, "cbc", placeholder=True,
            feature=feature("cell_barcode", "cell_barcode", "whitelist",
                            whitelist="10x Chromium ATAC whitelist")),
        seg("partial Read1N", nx.S5, "s5")])
    yield oligo("10x SI-PCR forward (CG000209 A1)",
                [seg("partial P5", SI_PCR_FORWARD, "p5")])
    yield oligo("10x indexed PCR reverse (CG000209 A1)", [
        seg("P7", il.P7, "p7"),
        seg("sample index", "N" * 8, "cbc", placeholder=True,
            feature=feature("sample_index_i7", "sample_index", "unknown")),
        seg("partial Read2N", "GTCTCGTGGGCTCGG", "r2")])


def workflow():
    start = Scene.duplex([seg("antibody-targeted chromatin", "X" * 38,
                              placeholder=True)], label="primary-antibody-bound chromatin")
    w = Workflow(MolecularState("Antibody-bound chromatin", tuple(start.rows())))
    w.react("Anneal BcdA MeA to Tn5_Rev and load nanobody–Tn5",
            loaded(mea_parts(), "MeA mosaic end", "BcdA nano-Tn5").rows(),
            name="target-coded nano-Tn5")
    w.react("Tagment antibody-proximal DNA and fill the staggered gap",
            Scene.duplex(list(FIRST_TAGGED), label="BcdA-tagged genomic template").rows(),
            name="BcdA-tagged genomic DNA")
    linear = Construct([seg("10x gel-bead P5", il.P5, "p5"),
                        seg("10x cell barcode", "B" * 16, "cbc", placeholder=True,
                            feature=feature("cell_barcode", "cell_barcode", "whitelist",
                                            whitelist="10x Chromium ATAC whitelist")),
                        *list(FIRST_TAGGED)], name="cell-barcoded linear-amplification product")
    w.react("Linearly amplify in Chromium droplets",
            Scene.duplex(list(linear), label="cell-barcoded copies").rows(),
            name="cell-barcoded linear products",
            note="The 10x ATAC v1.1 guide Appendix A1 prints P5–N16–partial Read1N for the gel-bead oligo.")
    w.react("Anneal MeB to Tn5_Rev and load standard Tn5",
            loaded(meb_parts(), "MeB mosaic end", "MeB-only Tn5").rows(),
            name="MeB-only transposome")
    w.react("Tagment the linear products with MeB-only Tn5",
            Scene.duplex(list(PRE_PCR), label="twice-tagmented nano-CT template").rows(),
            name="PCR-ready nano-CT molecule")
    geometry = PCRGeometry.unavailable(
        PRE_PCR,
        reason="The 10x guide prints primer scaffolds; this representative workflow has not yet reconstructed their strand-resolved PCR binding geometry.")
    w.pcr("Amplify with the 10x scATAC SI-PCR reagents", geometry,
          Scene.duplex(list(FINAL_LIBRARY), label="nano-CT sequencing library").rows(),
          name="nano-CT sequencing library")
    return w


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    nuclei = graph.input(
        "nuclei", "permeabilized nuclei", MaterialKind.BIOLOGICAL,
        state=molecule("permeabilized nuclei", (role("chromatin epitopes"),)),
        primary=True)
    antibody = graph.input("primary_antibodies", "species-resolved primary antibodies",
                           MaterialKind.SELECTION)
    mea_adapter = graph.input(
        "loaded_mea", "target-coded MeA-loaded nanobody–Tn5",
        MaterialKind.ADAPTER,
        state=loaded(mea_parts(), "MeA mosaic end", "BcdA nano-Tn5"))
    nano_tn5 = graph.input("nanobody_tn5", "species-selective nanobody–Tn5 fusion",
                           MaterialKind.CATALYST)
    targeted = graph.product("targeted_nuclei", "antibody-bound nano-Tn5 chromatin",
                             MaterialKind.MOLECULAR, state=nuclei.state)
    wash = graph.product("unbound_targeting_reagents", "unbound antibodies and nano-Tn5",
                         MaterialKind.WASTE)
    graph.react("target_chromatin", "Bind target-coded nano-Tn5 through primary antibodies",
                StepKind.AFFINITY_CAPTURE,
                ((nuclei, InputRole.TEMPLATE), (antibody, InputRole.SELECTION),
                 (mea_adapter, InputRole.STRUCTURAL), (nano_tn5, InputRole.CATALYST)),
                (targeted, wash))
    graph.discard(wash, "is removed before target-proximal tagmentation")
    first_tagged = graph.product("mea_tagged", FIRST_TAGGED.name,
                                 MaterialKind.MOLECULAR, state=FIRST_TAGGED)
    first_waste = graph.product("mea_tn5_waste", "non-amplifiable first-tagmentation products",
                                MaterialKind.WASTE)
    graph.react("mea_tagment", "Transfer the target-coded MeA adapter",
                StepKind.TAGMENTATION,
                ((targeted, InputRole.TEMPLATE), (mea_adapter, InputRole.STRUCTURAL),
                 (nano_tn5, InputRole.CATALYST)), (first_tagged, first_waste))
    graph.discard(first_waste, "does not retain the target-coded MeA end")
    gap_polymerase = graph.input("gap_fill_polymerase", "Tn5-gap-fill polymerase",
                                 MaterialKind.CATALYST)
    filled = graph.product("gap_filled_mea", "gap-filled MeA-tagged genomic DNA",
                           MaterialKind.MOLECULAR, state=FIRST_TAGGED)
    graph.react("gap_fill", "Fill the first staggered Tn5 gap", StepKind.EXTENSION,
                ((first_tagged, InputRole.TEMPLATE),
                 (gap_polymerase, InputRole.CATALYST)), (filled,))
    gel_bead = graph.input(
        "gel_bead_oligo", "10x scATAC gel-bead amplification oligo",
        MaterialKind.OLIGO, state=Construct([seg("P5", il.P5, "p5"),
            seg("10x cell barcode", "N" * 16, "cbc", placeholder=True,
                feature=feature("cell_barcode", "cell_barcode", "whitelist",
                                whitelist="10x Chromium ATAC whitelist")),
            seg("partial Read1N", nx.S5, "s5")]), source=ATAC_GUIDE)
    linear_polymerase = graph.input("linear_polymerase", "Chromium linear-amplification polymerase",
                                    MaterialKind.CATALYST)
    linear = graph.product(
        "cell_barcoded", "cell-barcoded MeA-tagged linear products",
        MaterialKind.MOLECULAR,
        state=molecule("cell-barcoded linear products", (
            role("10x cell barcode and amplification handle"),
            role("BcdA-coded MeA genomic fragment"))),
        source=SOURCE + "; " + ATAC_GUIDE)
    graph.react("linear_amplification", "Linearly amplify and cell-barcode in Chromium droplets",
                StepKind.EXTENSION,
                ((filled, InputRole.TEMPLATE), (gel_bead, InputRole.STRUCTURAL),
                 (linear_polymerase, InputRole.CATALYST)), (linear,))
    meb_adapter = graph.input("loaded_meb", "MeB-loaded standard Tn5 adapter",
                              MaterialKind.ADAPTER,
                              state=loaded(meb_parts(), "MeB mosaic end", "MeB-only Tn5"))
    standard_tn5 = graph.input("standard_tn5", "standard Tn5 transposase",
                               MaterialKind.CATALYST)
    twice_tagged = graph.product("twice_tagged", PRE_PCR.name,
                                 MaterialKind.MOLECULAR, state=PRE_PCR)
    second_waste = graph.product("meb_tn5_waste", "non-amplifiable second-tagmentation products",
                                 MaterialKind.WASTE)
    graph.react("meb_tagment", "Install the opposite MeB end on linear products",
                StepKind.TAGMENTATION,
                ((linear, InputRole.TEMPLATE), (meb_adapter, InputRole.STRUCTURAL),
                 (standard_tn5, InputRole.CATALYST)), (twice_tagged, second_waste))
    graph.discard(second_waste, "does not retain the MeA/MeB PCR-ready configuration")
    forward = graph.input("si_pcr_forward", "10x scATAC SI-PCR forward primer",
                          MaterialKind.OLIGO, source=ATAC_GUIDE)
    reverse = graph.input("si_pcr_reverse", "10x scATAC SI-PCR reverse primer",
                          MaterialKind.OLIGO, source=ATAC_GUIDE)
    polymerase = graph.input("pcr_polymerase", "SI-PCR polymerase", MaterialKind.CATALYST)
    library = graph.product("sequencing_library", FINAL_LIBRARY.name,
                            MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react("library_pcr", "Amplify with the 10x scATAC SI-PCR reagents",
                StepKind.PCR,
                ((twice_tagged, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    PRE_PCR,
                    reason="The 10x guide prints primer scaffolds; this representative workflow has not yet reconstructed their strand-resolved PCR binding geometry."))
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=tuple(SEQ_PRIMERS),
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS), read_lengths=READ_LENGTHS))
    graph.declare_reviewed_scope((
        "nuclei", "primary_antibodies", "loaded_mea", "nanobody_tn5",
        "target_chromatin", "mea_tagment", "gap_fill_polymerase", "gap_fill",
        "gel_bead_oligo", "linear_polymerase", "linear_amplification", "loaded_meb",
        "standard_tn5", "meb_tagment", "si_pcr_forward", "si_pcr_reverse",
        "pcr_polymerase", "library_pcr", "sequencing_library"), SOURCE)
    graph.annotate(SourceDiscrepancy(
        target="loaded_mea", title="Corrected MeA barcode-adapter sequence",
        claims=(
            SourceClaim("Author protocols.io v2, materials table lines 123–138",
                        "BcdA/B/C oligos omit CCACGC between S5 and each 8-base modality barcode."),
            SourceClaim("Author protocols.io v2, correction lines 271–293; published Supplementary Table 1 rows 3–14",
                        "BcdA/B/C oligos include CCACGC before the modality barcode.")),
        model_choice="The drawn BcdA adapter includes CCACGC before TATAGCCT.",
        evidence="The protocol labels the earlier strings erroneous, and the published oligo workbook prints the corrected family.",
        uncertainty="The source does not establish whether any earlier experiment used the erroneous strings.",
        slot=AnnotationSlot.BEFORE_WORKFLOW))
    return graph.require_publishable()
