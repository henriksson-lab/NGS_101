"""SNARE-seq RNA and accessible-chromatin molecular states."""
import nextera as nx

from chemdraw import Construct, MolecularState, Row, Scene, Segment, Workflow, feature
from multibranch import illumina_branch
from pcr import PCRGeometry, PCRPrimer
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
            f"<b>Nextera-R2</b> — 5′-{NEXTERA_R2}-3′"]


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
