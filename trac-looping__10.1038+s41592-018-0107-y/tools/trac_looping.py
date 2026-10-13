"""TrAC-looping molecular model (Lai et al. 2018)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import nextera as nx
from chemdraw import (MolecularState, Scene, Segment, Workflow,
                      annotation_rows, oligo)
from chromatin_epigenetics import tagmented_library
from contact_capture import BivalentTransposaseLinker
from circular import CircularizationGeometry, ClosedCircle
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role

BIO67F = "CTGTCTCTTATACACATCTCCGAGCCCACGAGACTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
BP67R = "CTGTCTCTTATACACATCTGACGCTGCCGACGAGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
SPACER = "CCGAGCCCACGAGACTCGTCGGCAGCGTC"
LINKER = BivalentTransposaseLinker(nx.ME_RC, SPACER, "67-bp bivalent ME linker")
CONTACT_LINEAR = LINKER.bridged_product()
CONTACT_CIRCLE = ClosedCircle(CONTACT_LINEAR,
                              CONTACT_LINEAR.segments[-1].name,
                              CONTACT_LINEAR.segments[0].name)


def rca_rows():
    """Three model-derived traversals of the random-primed contact circle."""
    repeats = []
    for number in range(1, 4):
        repeats.extend(Segment(f"copy {number}: {part.name}", part.top, part.tag,
                               placeholder=part.placeholder, inferred=part.inferred,
                               length_bp=part.length_bp)
                       for part in CONTACT_LINEAR)
    sc = Scene(); sc.strand("concatemer", repeats,
                            label="TempliPhi rolling-circle concatemer")
    return sc.rows()


def bridge_scene():
    sc = Scene.duplex(list(LINKER.bridged_product()), label="two accessible loci")
    sc.junction("top", "contact locus A", "ME A", "Tn5 insertion")
    sc.junction("top", "ME B", "contact locus B", "Tn5 insertion")
    sc.note("top", "tetrameric Tn5 inserts the same short linker into two chromatin regions in trans")
    sc.labels("top")
    return sc


FINAL_LIBRARY, SEQ_PRIMERS = tagmented_library("RCA-derived TrAC contact insert")


def final_rows():
    sc = Scene.duplex(list(FINAL_LIBRARY), label="paired-end library")
    sc.labels("top")
    return [*sc.rows(), *annotation_rows(FINAL_LIBRARY)]


def workflow():
    wf = Workflow(MolecularState("Bivalent transposition in intact chromatin", tuple(bridge_scene().rows())))
    wf.react("Digest genomic arms with a four-base cutter and capture biotinylated products",
             bridge_scene().rows(), note="The bridge, not proximity ligation, records the contact.")
    wf.react("Intramolecular circularization with T7 DNA ligase", bridge_scene().rows(),
             note="Dilute ligation closes each enriched contact molecule into a circle.")
    wf.react("Rolling-circle amplification", bridge_scene().rows(),
             note="TempliPhi amplification copies the circle before indexed library PCR.")
    wf.react("N501/N701 indexed PCR", final_rows(),
             note="The printed primers complete canonical Nextera-compatible P5 and P7 ends.")
    return wf


TITLE = "TrAC-looping — contacts captured by bivalent Tn5 insertion"
NOTES = "01_trac-looping.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/s41592-018-0107-y">'
          'Lai et al., <i>Nature Methods</i> (2018)</a>.')
SUMMARY = ("A tetrameric Tn5 complex inserts a biotinylated two-ME linker across two nearby "
           "accessible chromatin regions. Restriction digestion, affinity capture, circle formation "
           "and rolling-circle amplification recover the linked pair without Hi-C proximity ligation.")
FINAL_CAPTION = "PCR-completed Nextera-compatible molecule copied from a circularized TrAC contact product."
SEQUENCING_INTRO = "Paired reads enter the two genomic ends represented in the amplified circle; dual indexes come from N5xx/N7xx PCR primers."


def oligos():
    yield oligo("Bio67F", [Segment("bivalent top strand", BIO67F, "me")], mods="5′ phosphate; internal biotin-dT")
    yield oligo("67bpR", [Segment("bivalent bottom strand", BP67R, "me")], mods="5′ phosphate")


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    chromatin = graph.input(
        "chromatin", "intact chromatin", MaterialKind.BIOLOGICAL,
        state=molecule("intact chromatin", (role("spatially proximal accessible loci"),)),
        primary=True)
    linker = graph.input("bivalent_linker", "published biotinylated bivalent-ME linker",
                         MaterialKind.ADAPTER, state=LINKER.bridged_product())
    tn5 = graph.input("tetrameric_tn5", "tetrameric Tn5 transposase complex",
                      MaterialKind.CATALYST)
    bridged = graph.product("bridged_contacts", "bivalently tagmented chromatin contacts",
                            MaterialKind.MOLECULAR, state=bridge_scene())
    tn5_waste = graph.product("tn5_waste", "non-contact or non-amplifiable Tn5 products",
                              MaterialKind.WASTE)
    graph.react("bivalent_tagmentation", "Insert one bivalent linker across nearby loci",
                StepKind.TAGMENTATION,
                ((chromatin, InputRole.TEMPLATE), (linker, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (bridged, tn5_waste))
    graph.discard(tn5_waste, "does not retain the bivalent trans contact configuration")
    cutter = graph.input("four_base_cutter", "protocol four-base restriction endonuclease",
                         MaterialKind.CATALYST)
    digested = graph.product("digested_contacts", "restriction-trimmed bivalent contact products",
                             MaterialKind.MOLECULAR, state=bridged.state)
    graph.react("digest", "Digest genomic arms with a four-base cutter",
                StepKind.DIGESTION,
                ((bridged, InputRole.TEMPLATE), (cutter, InputRole.CATALYST)),
                (digested,))
    streptavidin = graph.input("streptavidin", "streptavidin capture beads",
                               MaterialKind.SELECTION)
    captured = graph.product("captured_contacts", "biotinylated bivalent contact products",
                             MaterialKind.MOLECULAR, state=digested.state)
    unbound = graph.product("unbound_fragments", "non-biotinylated restriction fragments",
                            MaterialKind.WASTE)
    graph.react("capture", "Capture the biotinylated bivalent products",
                StepKind.AFFINITY_CAPTURE,
                ((digested, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)),
                (captured, unbound))
    graph.discard(unbound, "does not retain the biotinylated bivalent linker")
    t7_ligase = graph.input("t7_ligase", "T7 DNA ligase", MaterialKind.CATALYST)
    circles = graph.product(
        "contact_circles", "intramolecularly circularized TrAC contact molecules",
        MaterialKind.MOLECULAR, state=CONTACT_CIRCLE)
    graph.react("circularize", "Intramolecularly circularize captured products",
                StepKind.CIRCULARIZATION,
                ((captured, InputRole.TEMPLATE), (t7_ligase, InputRole.CATALYST)),
                (circles,), geometry=CircularizationGeometry.known(
                    "double-strand intramolecular ligation",
                    statement="The defining protocol dilute-ligates captured contact molecules with T7 DNA ligase."))
    phi29 = graph.input("phi29", "TempliPhi phi29 DNA polymerase",
                        MaterialKind.CATALYST)
    rca = graph.product("rca_product", "rolling-circle-amplified TrAC contact DNA",
                        MaterialKind.MOLECULAR, state=tuple(rca_rows()))
    graph.react("rca", "Rolling-circle amplify the contact circles", StepKind.EXTENSION,
                ((circles, InputRole.TEMPLATE), (phi29, InputRole.CATALYST)), (rca,))
    n501 = graph.input("n501", "published N501 indexed PCR primer",
                       MaterialKind.OLIGO)
    n701 = graph.input("n701", "published N701 indexed PCR primer",
                       MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "indexed-library PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("sequencing_library", FINAL_LIBRARY.name,
                            MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react("indexed_pcr", "Amplify the RCA product with N501 and N701",
                StepKind.PCR,
                ((rca, InputRole.TEMPLATE), (n501, InputRole.STRUCTURAL),
                 (n701, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,
                    reason="The paper establishes the indexed Nextera components, but the variable RCA traversal is not represented as one exact pre-PCR construct."))
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=tuple(SEQ_PRIMERS),
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "chromatin", "bivalent_linker", "tetrameric_tn5", "bivalent_tagmentation",
        "four_base_cutter", "digest", "streptavidin", "capture", "t7_ligase",
        "circularize", "phi29", "rca", "n501", "n701", "pcr_polymerase",
        "indexed_pcr", "sequencing_library"), SOURCE)
    return graph.require_publishable()
