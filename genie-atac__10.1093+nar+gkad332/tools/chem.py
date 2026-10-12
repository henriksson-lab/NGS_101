"""GenIE-ATAC single-adapter targeted accessibility readout."""
import nextera as nx
from batch_ngs import seg, truseq_amplicon_library
from chemdraw import MolecularState, Scene, Segment, Workflow
from tagmentation import loaded_tn5_adapter


TITLE = "GenIE-ATAC — targeted allelic chromatin accessibility"
NOTES = "01_genie-atac.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1093/nar/gkad332">'
          'Cooper et al., <i>Nucleic Acids Research</i> (2023)</a>.')
SUMMARY = ("Tn5 loaded only with MEDS-A creates one universal transposon-primer end. "
           "A locus-specific biotinylated linear extension and streptavidin capture enrich "
           "one accessible locus before nested amplification and MiSeq sequencing.")
CAVEAT = ("MEDS-A/MEDS-REV and the rs12269414 example primers are source-exact. Locus "
          "primers are experiment-specific. The paper refers to a published MiSeq barcoding "
          "PCR without printing the selected index-primer pair, so only that terminal shell is inferred.")

MEDS_HANDLE = nx.S5
MEDS_ME = nx.ME
MEDS_A = MEDS_HANDLE + MEDS_ME
MEDS_REV = nx.ME_RC
MEDS_REV_WRITTEN = "p" + MEDS_REV + "/intT"
TN5A_PRIMER = "TCGTCGGCAGCGTCAG"

# A representative target from Supplementary Table 2. Every genomic part is a parameter
# in a real experiment; retaining one complete example makes the molecular path concrete.
TARGET = "rs12269414"
BIOTIN_LINEAR = "CCAGGAGTTGGAACCTGC"
NESTED_TAIL = "TCGGCATTCCTGCTGAACCGCTCTTCCGATCT"
NESTED_LOCUS = "CCTCTTGCCTTGGCTGTG"
NESTED_FORWARD = NESTED_TAIL + NESTED_LOCUS

FINAL_LIBRARY, SEQ_PRIMERS = truseq_amplicon_library([
    seg("nested rs12269414 locus primer-derived end", NESTED_LOCUS),
    seg("targeted accessible-DNA amplicon", "X" * 80, placeholder=True),
    seg("MEDS mosaic end reverse complement", nx.ME_RC, "me"),
    seg("MEDS-A primer-entry complement", nx.S5_RC, "s5"),
], "GenIE-ATAC targeted accessibility amplicon", inferred_adapters=True)
FINAL_CAPTION = ("Representative rs12269414 amplicon. One biological boundary is fixed by "
                 "the nested locus primer and the other by a nearby MEDS-A insertion. The "
                 "terminal MiSeq index shell is inferred because the selected barcoding "
                 "primer pair is not printed in the defining paper.")
READ_LENGTHS = {"Read 1": 150, "Index 1 (i7)": 8,
                "Index 2 (i5)": 8, "Read 2": 150}
SHOW_READ_CONTENT = True
SEQUENCING_INTRO = ("The reported 2 × 150-base MiSeq reads approach the locus-defined "
                    "amplicon from both ends and were merged when their inserts overlapped.")


def oligos():
    return (
        f"<b>MEDS-A top</b> — 5′-{MEDS_A}-3′",
        f"<b>MEDS-REV bottom</b> — 5′-{MEDS_REV_WRITTEN}-3′",
        f"<b>Tn5A PCR primer</b> — 5′-{TN5A_PRIMER}-3′",
        f"<b>{TARGET} biotinylated linear-PCR primer</b> — 5′-/5Biosg/{BIOTIN_LINEAR}-3′",
        f"<b>{TARGET} nested forward oligo</b> — 5′-{NESTED_FORWARD}-3′",
    )


def transposome_scene() -> Scene:
    return loaded_tn5_adapter(
        [Segment("MEDS-A primer-entry handle", MEDS_HANDLE, "s5"),
         Segment("MEDS-A mosaic end", MEDS_ME, "me")],
        mosaic_end="MEDS-A mosaic end", label="MEDS-A",
        bottom_name="MEDS-REV mosaic-end bottom", bottom_label="MEDS-REV",
        bottom_sequence=MEDS_REV, mod5="p", mod3="/intT/")


def _transferred_scene(*, gap_filled: bool) -> Scene:
    parts = [Segment("left MEDS-A handle", MEDS_HANDLE, "s5"),
             Segment("left mosaic end", MEDS_ME, "me"),
             Segment("accessible genomic fragment", "X" * 46, placeholder=True),
             Segment("right mosaic-end complement", nx.ME_RC, "me"),
             Segment("right MEDS-A handle complement", nx.S5_RC, "s5")]
    sc = Scene.duplex(parts, label=("gap-filled symmetric MEDS-A fragment"
                                    if gap_filled else "newly transferred MEDS-A fragment"))
    sc.junction("top", "left mosaic end", "accessible genomic fragment",
                "Tn5 transfer")
    sc.junction("top", "accessible genomic fragment", "right mosaic-end complement",
                "Tn5 transfer")
    if gap_filled:
        sc.mark("top", "accessible genomic fragment", "9-base stagger copied on both strands")
    else:
        sc.note("top", "Tn5 leaves the characteristic 9-base staggered single-strand gaps; the next state fills them before amplification.")
    sc.labels("top")
    return sc


def _linear_product_scene() -> Scene:
    sc = Scene()
    sc.strand("linear product", [
        Segment("5-prime biotin", "[5′ biotin]", "w1", placeholder=True),
        Segment("linear locus primer", BIOTIN_LINEAR),
        Segment("copied target and SNP", "X" * 42, placeholder=True),
        Segment("copied mosaic-end complement", nx.ME_RC, "me"),
        Segment("copied MEDS-A handle complement", nx.S5_RC, "s5"),
    ], label=f"{TARGET} linear extension product")
    sc.arrow("linear product", "extension from locus to a MEDS-A insertion")
    sc.mark("linear product", "5-prime biotin", "streptavidin capture")
    sc.labels("linear product")
    return sc


def _nested_scene() -> Scene:
    sc = Scene.duplex([
        Segment("nested-primer sequencing tail", NESTED_TAIL, "r1"),
        Segment("nested locus primer", NESTED_LOCUS),
        Segment("targeted accessible-DNA interval", "X" * 54, placeholder=True),
        Segment("MEDS mosaic end reverse complement", nx.ME_RC, "me"),
        Segment("Tn5A primer-entry complement", nx.S5_RC, "s5"),
    ], label=f"semi-nested {TARGET} amplicon")
    sc.mark("top", "nested locus primer", "locus-specific boundary")
    sc.mark("top", "Tn5A primer-entry complement", "tagmentation boundary")
    sc.labels("top")
    return sc


def workflow() -> Workflow:
    w = Workflow(MolecularState("Exact single-adapter Tn5 complex",
                                tuple(transposome_scene().rows())))
    w.react("Tagment accessible chromatin with MEDS-A-only Tn5",
            _transferred_scene(gap_filled=False).rows(),
            name="Gapped symmetric MEDS-A-tagged fragments",
            note="Both transferred ends use the same MEDS-A sequence; there is no MEDS-B end in GenIE-ATAC.")
    w.react("Fill the 9-base Tn5 gaps",
            _transferred_scene(gap_filled=True).rows(),
            name="Amplifiable MEDS-A-tagged fragments",
            note="The filled duplex makes the common Tn5A primer site available at either insertion boundary.")
    w.react(f"Linearly extend from the biotinylated {TARGET} primer",
            _linear_product_scene().rows(), name="Biotinylated locus-to-MEDS-A products",
            note="Repeated single-primer cycles copy only molecules that span the chosen locus and a nearby open-chromatin insertion.")
    w.react("Capture biotinylated strands on streptavidin beads",
            _linear_product_scene().rows(), name="Enriched targeted strands",
            note="Unrelated genomic MEDS-A fragments lack the locus-primer-derived biotin and wash away.")
    w.react("Amplify with the nested locus primer and common Tn5A primer",
            _nested_scene().rows(), name="Semi-nested targeted amplicon",
            note="The genomic primer fixes one boundary; MEDS-A fixes the accessibility-dependent boundary.")
    final = Scene.duplex(list(FINAL_LIBRARY), label="indexed MiSeq amplicon",
                         unpaired=tuple(s.name + "'" for s in FINAL_LIBRARY if s.inferred))
    w.react("Add MiSeq barcodes by a second PCR",
            final.rows(), name="Indexed GenIE-ATAC library",
            note="The selected index-primer pair is not printed in the paper, so its standard sequencing shell is visibly inferred.")
    return w
