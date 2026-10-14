"""SNARE-seq RNA and accessible-chromatin molecular states."""
import nextera as nx

from chemdraw import Construct, MolecularState, Row, Scene, Segment, Workflow, feature
from multibranch import illumina_branch
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (Disposition, Input, InputRole, Material, MaterialKind,
                            ProtocolGraph, SequencingPlan, SourceCoverage, SourceItem,
                            Step, StepKind, Terminal)
from protocol_extensions import analyte_split

TITLE = "SNARE-seq — transcriptome and chromatin accessibility"
NOTES = "01_snare-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/s41587-019-0290-0">'
          'Chen et al., <i>Nature Biotechnology</i> (2019)</a>.')
SUMMARY = ("A phosphorylated poly(A) splint anneals across the complete Nextera S5+ME "
           "end of pre-tagmented chromatin. A Drop-seq bead oligo captures that tail and "
           "nuclear mRNA, transferring the same cell barcode to both branches.")

BEAD_SMART = "AAGCAGTGGTATCAACGCAGAGTAC"
BEAD_SPACER = "T" * 8
CELL_NT, UMI_NT = 12, 8
BEAD_DT, SPLINT_POLYA = "T" * 30, "A" * 31
SPLINT = nx.ME_RC + nx.S5_RC + SPLINT_POLYA
NEXTERA_R1_BK = nx.S5 + nx.ME
NEXTERA_AD2_BK = nx.S7_RC
TSO = "AAGCAGTGGTATCAACGCAGAGTGAAT"  # followed by rGrGrG as ordered
TSO_PCR = "AAGCAGTGGTATCAACGCAGAGT"
NEXTERA_R2 = nx.S7 + nx.ME


def _cell_segment():
    return Segment("12-nt cell barcode", "N" * CELL_NT, "cbc",
                   feature=feature("cell_barcode", "cell_barcode", "whitelist",
                                   whitelist="Drop-seq Macosko-2011-10 bead whitelist"))


def _umi_segment():
    return Segment("8-nt UMI", "N" * UMI_NT, "umi",
                   feature=feature("umi", "umi", "random"))


def _bead_segments():
    return [Segment("bead spacer", BEAD_SPACER),
            Segment("bead SMART handle", BEAD_SMART, "r1"),
            _cell_segment(), _umi_segment(), Segment("bead oligo-dT", BEAD_DT)]


def _tagged_segments():
    return [Segment("S5 handle", nx.S5, "s5"),
            Segment("left mosaic end", nx.ME, "me"),
            Segment("accessible genomic fragment", "X" * 42, placeholder=True),
            Segment("right mosaic-end complement", nx.ME_RC, "me"),
            Segment("S7-handle complement", nx.S7_RC, "s7")]


def tagmented_scene():
    sc = Scene.duplex(_tagged_segments(), label="productive TDE1-tagged chromatin")
    sc.junction("top", "left mosaic end", "accessible genomic fragment", "Tn5 transfer")
    sc.junction("top", "accessible genomic fragment", "right mosaic-end complement",
                "Tn5 transfer")
    sc.note("top", "Only an S5 end binds the capture splint; later PCR selects the asymmetric S5/S7 product.")
    sc.labels("top")
    return sc


def capture_scene():
    """Exact three-strand S5-end/splint/bead capture intermediate."""
    sc = Scene()
    sc.strand("tagged DNA", [Segment("S5 handle", nx.S5, "s5"),
                             Segment("mosaic end", nx.ME, "me"),
                             Segment("accessible genomic DNA", "X" * 36,
                                     placeholder=True)],
              label="tagmented accessible DNA")
    sc.anneal("phosphorylated capture splint", [
        Segment("mosaic-end reverse complement", nx.ME_RC, "me"),
        Segment("S5-handle reverse complement", nx.S5_RC, "s5"),
        Segment("31-nt poly(A) capture tail", SPLINT_POLYA),
    ], to="tagged DNA", pair=("mosaic-end reverse complement", "mosaic end"),
              label="Nextera-R1-rc-polyA", mod5="/5Phos/")
    sc.anneal("bead oligo", _bead_segments(), to="phosphorylated capture splint",
              pair=("bead oligo-dT", "31-nt poly(A) capture tail"),
              label="Drop-seq bead oligo")
    sc.mark("phosphorylated capture splint", "31-nt poly(A) capture tail",
            "hybridization transfers this chromatin end to the barcode bead")
    sc.labels("phosphorylated capture splint")
    sc.labels("bead oligo")
    return sc


def atac_template():
    return Construct([
        *_bead_segments(),
        Segment("copied S5 handle", nx.S5, "s5"),
        Segment("copied left mosaic end", nx.ME, "me"),
        Segment("accessible genomic fragment", "X" * 42, placeholder=True),
        Segment("right mosaic-end complement", nx.ME_RC, "me"),
        Segment("S7-handle complement", nx.S7_RC, "s7"),
    ], name="gap-filled and ligated SNARE-seq chromatin template")


def gap_filled_scene():
    sc = Scene.duplex(list(atac_template()), label="barcoded chromatin PCR template")
    sc.junction("bottom", "copied left mosaic end'", "accessible genomic fragment'",
                "T7 ligation after gap fill")
    sc.mark("top", "12-nt cell barcode", "same bead identity as the RNA branch")
    sc.labels("top")
    return sc


def preamp_product():
    # The PCR starts inside the SMART handle, so the bead's eight-T spacer is absent.
    return Construct(list(atac_template())[1:], name="SNARE-seq common ATAC preamp product")


def preamp_scene():
    sc = Scene.duplex(list(preamp_product()), label="Nextera-R2 / TSO-PCR amplicon")
    sc.mark("top", "12-nt cell barcode", "cell identity")
    sc.mark("top", "8-nt UMI", "used for RNA molecules")
    sc.labels("top")
    return sc


RNA, RP = illumina_branch("SNARE-seq RNA", "3-prime cDNA", cell=12, umi=8)
ATAC, AP = illumina_branch("SNARE-seq accessibility", "tagmented genomic DNA",
                           architecture="nextera", cell=12)
FINAL_LIBRARIES = (
    ("RNA library", RNA, RP, "Cell barcode and UMI precede transcript-derived sequence.",
     "Declared Illumina primers on the Nextera-XT-finished RNA branch."),
    ("Accessibility library", ATAC, AP,
     "The matching cell barcode accompanies accessible genomic DNA.",
     "Declared Nextera primers on the DNA branch."),
)


def oligos():
    return [f"<b>Nextera-R1-rc-polyA</b> — 5′-/5Phos/{SPLINT}-3′",
            f"<b>Nextera-R1-bk</b> — 5′-/5Phos/{NEXTERA_R1_BK}/3InvdT/-3′",
            f"<b>Nextera-Ad2-bk</b> — 5′-/5Phos/{NEXTERA_AD2_BK}/3InvdT/-3′",
            f"<b>TSO</b> — 5′-{TSO}rGrGrG-3′",
            f"<b>TSO-PCR</b> — 5′-{TSO_PCR}-3′",
            f"<b>Nextera-R2</b> — 5′-{NEXTERA_R2}-3′",
            "<b>P501-Tso, example i5 primer</b> — 5′-AATGATACGGCGACCACCGAGATCTACACTAGATCGCGCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGT*A*C-3′",
            "<b>Ad2.1_TAAGGCGA, example i7 primer</b> — 5′-CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGAGATGT-3′",
            "<b>Read1CustomSeqB</b> — 5′-GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC-3′"]


def workflow():
    start = (Row(chunks=[("permeabilized nucleus: nuclear poly(A) RNA + accessible chromatin",
                          None, False)]),)
    w = Workflow(MolecularState("Permeabilized nucleus", start))
    w.react("Tagment accessible chromatin with commercial TDE1",
            tagmented_scene().rows(), name="productive asymmetric tagmented fragment",
            note="The productive S5/ME–insert–MErc/S7rc structure is one of the three possible Nextera end combinations.")
    w.react("Hybridize the phosphorylated poly(A) splint and barcode-bead oligo",
            capture_scene().rows(), name="barcode-bead-captured chromatin",
            note="The exact 33-nt Nextera-complementary arm spans ME and S5; its 31-A tail binds the bead's 30-T tract.")
    w.react("Fill the staggered gaps and seal the phosphorylated splint-side nick",
            gap_filled_scene().rows(), name="gap-filled, barcoded chromatin template",
            note="Polymerase extends the captured molecule and T7 ligase seals the phosphorylated splint/genomic junction; the opposite S5-side nick is not claimed as sealed.")
    geometry = PCRGeometry.exact(
        atac_template(), PCRPrimer("TSO-PCR", TSO_PCR),
        PCRPrimer("Nextera-R2", NEXTERA_R2),
        statement="Supplementary Table 4 prints both common-amplification primers; their inward-facing sites select the asymmetric S5/S7 product.")
    w.pcr("Co-amplify captured chromatin with Nextera-R2 and TSO-PCR",
          geometry, preamp_scene().rows(), name="common chromatin preamp",
          note="RNA molecules amplify through TSO-PCR sites at both ends; the exact ATAC pair shown selects the asymmetric chromatin product.")
    w.react("Split into RNA and accessibility finishing branches",
            analyte_split("TSO enrichment, Nextera XT and indexed RNA PCR",
                          "P5xx-TSO / Ad2.x indexed chromatin PCR"),
            name="separate sequencing libraries",
            note="Supplementary Table 4 prints the finishing-primer families; index choice varies by sample.")
    return w


def protocol_graph() -> ProtocolGraph:
    """RNA and chromatin branches joined by one physical bead barcode."""
    source = "Chen et al. 2019, doi:10.1038/s41587-019-0290-0"
    graph = ProtocolGraph(TITLE, authority=source)

    def add_input(key, label, kind, *, primary=False):
        graph.add_input(Material(key, label, kind, primary=primary, source=source))

    def add_step(key, label, kind, ports, products, *, geometry=None):
        products = tuple(products)
        graph.add_step(Step(key, label, kind, tuple(ports),
                            tuple(x.id for x in products), geometry=geometry,
                            source=source), products)

    add_input("nucleus", "permeabilized nucleus containing poly(A) RNA and chromatin",
              MaterialKind.BIOLOGICAL, primary=True)
    add_input("loaded_tn5_adapters", "TDE1-loaded Nextera adapters", MaterialKind.ADAPTER)
    add_input("tde1", "TDE1 transposase", MaterialKind.CATALYST)
    tagged = Material("tagmented_nucleus", "nucleus with productively tagmented chromatin",
                      MaterialKind.MOLECULAR, state=tagmented_scene(), source=source)
    tag_waste = Material("non_library_tn5_products",
                         "symmetric, untagged and otherwise non-library Tn5 products",
                         MaterialKind.WASTE, source=source)
    add_step("tagment", "Tagment accessible chromatin", StepKind.TAGMENTATION,
             (Input("nucleus", InputRole.TEMPLATE),
              Input("loaded_tn5_adapters", InputRole.STRUCTURAL),
              Input("tde1", InputRole.CATALYST)), (tagged, tag_waste))
    graph.mark_terminal(Terminal(tag_waste.id, Disposition.DISCARD,
                                 "not the asymmetric S5/S7 product selected by PCR"))

    add_input("capture_splint", "phosphorylated Nextera-R1-rc-polyA splint",
              MaterialKind.OLIGO)
    add_input("barcode_bead_oligo", "Drop-seq barcode-bead oligo", MaterialKind.OLIGO)
    captured = Material("captured_nucleus", "RNA and chromatin captured on one barcode bead",
                        MaterialKind.MOLECULAR, state=capture_scene(), source=source)
    add_step("capture", "Hybridize the poly(A) splint and barcode-bead oligo",
             StepKind.HYBRIDIZATION,
             (Input(tagged.id, InputRole.TEMPLATE),
              Input("capture_splint", InputRole.STRUCTURAL),
              Input("barcode_bead_oligo", InputRole.STRUCTURAL)), (captured,))

    rna_template = Material(
        "rna_bead_template", "bead-captured nuclear poly(A) RNA",
        MaterialKind.MOLECULAR,
        state=(Row(chunks=[("bead barcode + UMI + oligo(dT) … poly(A) RNA", None, False)]),),
        source=source)
    atac_template_material = Material(
        "atac_bead_template", "bead-captured tagmented chromatin",
        MaterialKind.MOLECULAR, state=capture_scene(), source=source)
    add_step("separate_analytes", "Separate the RNA and chromatin molecular branches",
             StepKind.SPLIT, (Input(captured.id, InputRole.TEMPLATE),),
             (rna_template, atac_template_material))

    add_input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    add_input("tso", "template-switching oligo", MaterialKind.OLIGO)
    rna_cdna = Material(
        "barcoded_cdna", "cell-barcoded, UMI-bearing cDNA",
        MaterialKind.MOLECULAR,
        state=(Row(chunks=[("TSO — cDNA — UMI — cell barcode — SMART handle", None, False)]),),
        source=source)
    add_step("rna_rt", "Reverse-transcribe and template-switch nuclear RNA",
             StepKind.REVERSE_TRANSCRIPTION,
             (Input(rna_template.id, InputRole.TEMPLATE),
              Input("barcode_bead_oligo", InputRole.STRUCTURAL),
              Input("tso", InputRole.STRUCTURAL),
              Input("reverse_transcriptase", InputRole.CATALYST)), (rna_cdna,))

    add_input("gap_fill_polymerase", "gap-fill polymerase", MaterialKind.CATALYST)
    gap_filled = Material("gap_filled_atac", "gap-filled barcoded chromatin template",
                          MaterialKind.MOLECULAR, state=gap_filled_scene(), source=source)
    add_step("gap_fill", "Fill the staggered chromatin gaps", StepKind.EXTENSION,
             (Input(atac_template_material.id, InputRole.TEMPLATE),
              Input("gap_fill_polymerase", InputRole.CATALYST)), (gap_filled,))
    add_input("t7_ligase", "T7 DNA ligase", MaterialKind.CATALYST)
    sealed = Material("sealed_atac", "splint-side-nick-sealed chromatin template",
                      MaterialKind.MOLECULAR, state=gap_filled_scene(), source=source)
    add_step("seal_splint_nick", "Seal the phosphorylated splint-side nick",
             StepKind.NICK_SEALING,
             (Input(gap_filled.id, InputRole.TEMPLATE),
              Input("t7_ligase", InputRole.CATALYST)), (sealed,))

    add_input("tso_pcr", "TSO-PCR", MaterialKind.OLIGO)
    add_input("nextera_r2", "Nextera-R2", MaterialKind.OLIGO)
    preamp_geometry = PCRGeometry.exact(
        atac_template(), PCRPrimer("TSO-PCR", TSO_PCR),
        PCRPrimer("Nextera-R2", NEXTERA_R2),
        statement="Supplementary Table 4 prints both common-amplification primers; their inward-facing sites select the asymmetric S5/S7 chromatin product.")
    atac_preamp = Material("atac_preamp", "common chromatin preamplification product",
                           MaterialKind.MOLECULAR, state=preamp_product(), source=source)
    add_step("atac_preamp_pcr", "Co-amplify captured chromatin", StepKind.PCR,
             (Input(sealed.id, InputRole.TEMPLATE),
              Input("tso_pcr", InputRole.STRUCTURAL),
              Input("nextera_r2", InputRole.STRUCTURAL)),
             (atac_preamp,), geometry=preamp_geometry)

    for prefix, label in (("rna", "RNA finishing primer pool"),
                          ("atac", "accessibility finishing primer pool")):
        add_input(prefix + "_finish_forward", label + " forward", MaterialKind.OLIGO)
        add_input(prefix + "_finish_reverse", label + " reverse", MaterialKind.OLIGO)
    rna_geometry = PCRGeometry.unavailable(
        RNA, reason="Supplementary Table 4 publishes an index-varying finishing-primer family rather than one immutable selected pair.")
    atac_geometry = PCRGeometry.unavailable(
        ATAC, reason="Supplementary Table 4 publishes an index-varying finishing-primer family rather than one immutable selected pair.")
    rna_library = Material("rna_library", "SNARE-seq RNA library", MaterialKind.LIBRARY,
                           state=RNA, source=source)
    atac_library = Material("atac_library", "SNARE-seq accessibility library",
                            MaterialKind.LIBRARY, state=ATAC, source=source)
    add_step("finish_rna", "Finish the indexed RNA library", StepKind.PCR,
             (Input(rna_cdna.id, InputRole.TEMPLATE),
              Input("rna_finish_forward", InputRole.STRUCTURAL),
              Input("rna_finish_reverse", InputRole.STRUCTURAL)),
             (rna_library,), geometry=rna_geometry)
    add_step("finish_atac", "Finish the indexed accessibility library", StepKind.PCR,
             (Input(atac_preamp.id, InputRole.TEMPLATE),
              Input("atac_finish_forward", InputRole.STRUCTURAL),
              Input("atac_finish_reverse", InputRole.STRUCTURAL)),
             (atac_library,), geometry=atac_geometry)
    graph.mark_terminal(Terminal(rna_library.id, Disposition.SEQUENCING),
                        SequencingPlan(RNA, tuple(RP), tuple(p.role for p in RP)))
    graph.mark_terminal(Terminal(atac_library.id, Disposition.SEQUENCING),
                        SequencingPlan(ATAC, tuple(AP), tuple(p.role for p in AP)))

    # This explicit list is the reviewed source inventory.  Validation separately checks
    # that every external input, step and sequencing output occurs in it.
    reviewed_targets = (
        "nucleus", "loaded_tn5_adapters", "tde1", "tagment", "capture_splint",
        "barcode_bead_oligo", "capture", "separate_analytes", "reverse_transcriptase",
        "tso", "rna_rt", "gap_fill_polymerase", "gap_fill", "t7_ligase",
        "seal_splint_nick", "tso_pcr", "nextera_r2", "atac_preamp_pcr",
        "rna_finish_forward", "rna_finish_reverse", "atac_finish_forward",
        "atac_finish_reverse", "finish_rna", "finish_atac", "rna_library",
        "atac_library",
    )
    graph.declare_reviewed_scope(reviewed_targets, source)
    return graph.require_publishable()
