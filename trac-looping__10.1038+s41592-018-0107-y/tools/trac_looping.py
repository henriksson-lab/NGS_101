"""TrAC-looping molecular model (Lai et al. 2018)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow, feature, revcomp,
                      annotation_rows, oligo)
from circular import CircularizationGeometry, ClosedCircle
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind, Terminal)
from protocol_model import molecule, role

BIO67F = "CTGTCTCTTATACACATCTCCGAGCCCACGAGACTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
BP67R = "CTGTCTCTTATACACATCTGACGCTGCCGACGAGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
SPACER = BIO67F[len(nx.ME_RC):-len(nx.ME)]
if BP67R != revcomp(BIO67F):
    raise ValueError("source Bio67F/67bpR duplex does not pair")
LINKER_SEGMENTS = (Segment("ME A", nx.ME_RC, "me"),
                   Segment("linker spacer", SPACER, "linker"),
                   Segment("ME B", nx.ME, "me"))
if "".join(s.top for s in LINKER_SEGMENTS) != BIO67F:
    raise ValueError("local linker differs from source Bio67F")
CONTACT_LINEAR = Construct([Segment("contact locus A", "X"*24, placeholder=True),
                            *LINKER_SEGMENTS,
                            Segment("contact locus B", "X"*24, placeholder=True)],
                           name="representative bivalent contact")


def linker_scene():
    sc = Scene()
    sc.strand("Bio67F", [Segment("source Bio67F", BIO67F)], mod5="phosphate")
    sc.anneal("67bpR", [Segment("source67bpR", BP67R)], to="Bio67F",
              pair=("source67bpR", "source Bio67F"), mod5="phosphate")
    sc.note("Bio67F", "Source internal biotin-dT at nucleotide 35 (one-based); not a terminal biotin.")
    return sc

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
    sc = Scene.duplex(list(CONTACT_LINEAR), label="two accessible loci")
    sc.junction("top", "contact locus A", "ME A", "Tn5 insertion")
    sc.junction("top", "ME B", "contact locus B", "Tn5 insertion")
    sc.note("top", "tetrameric Tn5 inserts the same short linker into two chromatin regions in trans")
    sc.labels("top")
    return sc


N501 = il.P5 + "TAGATCGC" + nx.READ1_PRIMER
N701 = il.P7 + "TCGCCTTA" + nx.READ2_PRIMER
FINAL_LIBRARY = Construct([
    Segment("P5", il.P5, "p5"),
    Segment("N501 index", "TAGATCGC", "cbc", feature=feature("sample_i5", "sample_index", "whitelist", whitelist="TrAC Supplementary Table2")),
    Segment("Read1 site", nx.READ1_PRIMER, "r1"),
    Segment("RCA-derived contact insert", "X"*36, placeholder=True),
    Segment("Read2 site reverse complement", revcomp(nx.READ2_PRIMER), "r2"),
    Segment("N701 index reverse complement", revcomp("TCGCCTTA"), "cbc", feature=feature("sample_i7", "sample_index", "whitelist", whitelist="TrAC Supplementary Table2")),
    Segment("P7 reverse complement", il.P7_RC, "p7")], name="TrAC-looping N501/N701 library")
SEQ_PRIMERS = tuple(sp.NEXTERA[k] for k in ("R1", "I1", "I2", "R2"))
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS, required_roles=tuple(p.role for p in SEQ_PRIMERS)):
    raise ValueError("invalid TrAC endpoint: " + "; ".join(problems))


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
CAVEAT = ("The local linker duplex follows printed Bio67F/67bpR; the former repeated-ME orientation was a repository error. This graph traces the NlaIII aliquot; the source's parallel MluCI aliquot is not reconstructed. Contact-circle flanks and RCA traversal are schematic, not exact heterogeneous experimental products.")
FINAL_CAPTION = "PCR-completed Nextera-compatible molecule copied from a circularized TrAC contact product."
SEQUENCING_INTRO = "Paired reads enter the two genomic ends represented in the amplified circle; dual indexes come from N5xx/N7xx PCR primers."


def oligos():
    yield oligo("ME3ddCTPtop", [Segment("modified ME top", "CTGTCTCTTATACACATCTCTGATGGCGCGAGGGA", "me")], mods="5′ phosphate; 3′ /3ddCTP/")
    yield oligo("5NH2MEbottom", [Segment("modified ME bottom", "TCCCTCGCGCCATCAGAGATGTGTATAAGAGACAG", "me")], mods="5′ /5AmC6/")
    yield oligo("Bio67F", [Segment("bivalent top strand", BIO67F, "me")], mods="5′ phosphate; internal biotin-dT")
    yield oligo("67bpR", [Segment("bivalent bottom strand", BP67R, "me")], mods="5′ phosphate")
    yield oligo("N501, example i5 PCR primer", [Segment("indexed P5 primer", "AATGATACGGCGACCACCGAGATCTACACTAGATCGCTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG", "p5")])
    yield oligo("N701, example i7 PCR primer", [Segment("indexed P7 primer", "CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG", "p7")])


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    chromatin = graph.input(
        "chromatin", "intact chromatin", MaterialKind.BIOLOGICAL,
        state=molecule("intact chromatin", (role("spatially proximal accessible loci"),)),
        primary=True)
    linker = graph.input("bivalent_linker", "published biotinylated bivalent-ME linker",
                         MaterialKind.ADAPTER, state=linker_scene())
    tn5 = graph.input("tetrameric_tn5", "tetrameric Tn5 transposase complex",
                      MaterialKind.CATALYST)
    bridged = graph.product("bridged_contacts", "bivalently tagmented chromatin contacts",
                            MaterialKind.MOLECULAR, state=bridge_scene())
    tn5_waste = graph.product("tn5_waste", "candidate non-contact or library-ineligible Tn5 products",
                              MaterialKind.MOLECULAR, inferred=True)
    graph.react("bivalent_tagmentation", "Insert one bivalent linker across nearby loci",
                StepKind.TAGMENTATION,
                ((chromatin, InputRole.TEMPLATE), (linker, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (bridged, tn5_waste))
    graph.mark_terminal(Terminal(tn5_waste.id, Disposition.UNRESOLVED, "Candidate products, not an established measured population; loss of contact configuration or PCR eligibility does not prove physical removal."))
    proteinase = graph.input("proteinase", "EDTA/SDS/Proteinase K quench and DNA extraction", MaterialKind.CATALYST)
    purified = graph.product("purified", "purified transposition products", MaterialKind.MOLECULAR)
    graph.react("release_dna", "Quench, reverse crosslinks and recover DNA by extraction/precipitation", StepKind.MAGIC,
                ((bridged, InputRole.TEMPLATE), (proteinase, InputRole.CATALYST)), (purified,),
                reason="Source-stated release and recovery summarized; individual products/fates are not reconstructed.")
    repair = graph.input("gap_polymerase", "T4 DNA polymerase with dNTPs", MaterialKind.CATALYST)
    repaired = graph.product("repaired", "gap-repaired transposition products", MaterialKind.MOLECULAR)
    graph.react("gap_repair", "Repair transposition gaps with T4 DNA polymerase", StepKind.END_REPAIR,
                ((purified, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    ampure = graph.input("ampure", "AMPure XP recovery of fragments greater than 300 bp", MaterialKind.SELECTION)
    selected = graph.product("selected", "size-retained gap-repaired DNA", MaterialKind.MOLECULAR)
    short = graph.product("short", "unrecovered short DNA", MaterialKind.WASTE)
    graph.react("initial_size_selection", "Recover greater-than 300 bp DNA before restriction digestion", StepKind.SELECTION,
                ((repaired, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)), (selected, short))
    graph.discard(short, "not recovered in the selected DNA fraction")
    nla = graph.product("nla_aliquot", "NlaIII digestion aliquot", MaterialKind.MOLECULAR)
    mlu = graph.product("mlu_aliquot", "parallel MluCI aliquot", MaterialKind.MOLECULAR)
    graph.react("split_digests", "Split eluate into NlaIII and MluCI aliquots", StepKind.SPLIT,
                ((selected, InputRole.TEMPLATE),), (nla, mlu), changes_molecule=False)
    graph.mark_terminal(Terminal(mlu.id, Disposition.UNRESOLVED, "The source processes this parallel aliquot with MluCI; its full downstream molecular path is outside the selected NlaIII trace."))
    cutter = graph.input("four_base_cutter", "NlaIII restriction enzyme (selected aliquot)",
                         MaterialKind.CATALYST)
    digested = graph.product("digested_contacts", "restriction-trimmed bivalent contact products",
                             MaterialKind.MOLECULAR, state=bridged.state)
    graph.react("digest", "Digest the selected aliquot with NlaIII",
                StepKind.DIGESTION,
                ((nla, InputRole.TEMPLATE), (cutter, InputRole.CATALYST)),
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
    release = graph.input("capture_elution", "SDS/Proteinase K elution, phenol-chloroform and ethanol recovery", MaterialKind.CONDITION)
    eluted = graph.product("eluted", "DNA eluted from streptavidin beads for circularization", MaterialKind.MOLECULAR)
    graph.react("elute_capture", "Elute captured DNA and purify before dilute ligation", StepKind.MAGIC,
                ((captured, InputRole.TEMPLATE), (release, InputRole.CONDITION)), (eluted,),
                reason="Published elution/recovery steps are summarized; exact released molecular end states are not reconstructed.")
    t7_ligase = graph.input("t7_ligase", "T7 DNA ligase", MaterialKind.CATALYST)
    circles = graph.product(
        "contact_circles", "intramolecularly circularized TrAC contact molecules",
        MaterialKind.MOLECULAR, state=CONTACT_CIRCLE)
    graph.react("circularize", "Intramolecularly circularize captured products",
                StepKind.CIRCULARIZATION,
                ((eluted, InputRole.TEMPLATE), (t7_ligase, InputRole.CATALYST)),
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
                       MaterialKind.OLIGO, state=Construct([Segment("N501", N501)]))
    n701 = graph.input("n701", "published N701 indexed PCR primer",
                       MaterialKind.OLIGO, state=Construct([Segment("N701", N701)]))
    polymerase = graph.input("pcr_polymerase", "NEB Phusion HF master mix",
                             MaterialKind.CATALYST)
    library = graph.product("amplified_library", FINAL_LIBRARY.name,
                            MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    graph.react("indexed_pcr", "Amplify RCA product with N501/N701 for 11 cycles",
                StepKind.PCR,
                ((rca, InputRole.TEMPLATE), (n501, InputRole.STRUCTURAL),
                 (n701, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,
                    reason="The paper establishes the indexed Nextera components, but the variable RCA traversal is not represented as one exact pre-PCR construct."))
    gel = graph.input("size_gel", "E-gel 250–650 bp recovery interval", MaterialKind.SELECTION)
    final = graph.product("sequencing_library", FINAL_LIBRARY.name, MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    offsize = graph.product("offsize", "products outside recovered E-gel interval", MaterialKind.WASTE)
    graph.react("final_size_selection", "Recover 250–650 bp PCR products from E-gel", StepKind.SELECTION,
                ((library, InputRole.TEMPLATE), (gel, InputRole.SELECTION)), (final, offsize))
    graph.discard(offsize, "outside recovered E-gel size interval")
    graph.sequence(final, SequencingPlan(
        library=FINAL_LIBRARY, primers=tuple(SEQ_PRIMERS),
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "chromatin", "bivalent_linker", "tetrameric_tn5", "bivalent_tagmentation",
        "four_base_cutter", "digest", "streptavidin", "capture", "t7_ligase",
        "circularize", "phi29", "rca", "n501", "n701", "pcr_polymerase",
        "indexed_pcr", "sequencing_library", "proteinase", "release_dna", "gap_polymerase", "gap_repair",
        "ampure", "initial_size_selection", "split_digests", "capture_elution", "elute_capture", "size_gel", "final_size_selection"), SOURCE)
    return graph.require_publishable()
