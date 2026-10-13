"""Li et al. 2017 long-read ChIA-PET molecular model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from batch_ngs import nextera_library, seg
from chemdraw import Construct, Row, Scene, Segment
from proximity_linker import ThreePrimeBridgeLinker
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role

TITLE = "Long-read ChIA-PET — bridge-linked chromatin contacts"
NOTES = "01_long-read-chia-pet.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/nprot.2017.012">'
          'Li et al., <i>Nature Protocols</i> (2017)</a>.')
SUMMARY = ("Protein-associated chromatin fragments are repaired and dA-tailed, joined in "
           "proximity through a double-ended biotin bridge, then tagmented; streptavidin "
           "retains only fragments that still contain the contact junction.")
CAVEAT = ("The bridge-linker bases and modification are published. The final sequencing "
          "ends use the canonical Nextera components supplied by the named Illumina kit.")

BRIDGE_F = "CGCGATATCTTATCTGACT"
BRIDGE_R = "GTCAGATAAGATATCGCGT"
BRIDGE = ThreePrimeBridgeLinker("biotin bridge linker", BRIDGE_F, BRIDGE_R,
                                overhang="T", modified_top_position=9)

def contact_segments() -> list[Segment]:
    core = BRIDGE.top_core
    return [
        seg("left contact tag", "X" * 28, placeholder=True),
        seg("left A:T ligation junction", "A", bottom="T"),
        seg("bridge core left", core[:9], "cbc"),
        seg("internal biotin-dT", core[9], "umi"),
        seg("bridge core right", core[10:], "cbc"),
        seg("right T:A ligation junction", "T", bottom="A"),
        seg("right contact tag", "X" * 28, placeholder=True),
    ]

CONTACT = Construct(contact_segments(), name="DNA–bridge–DNA proximity product")
FINAL_LIBRARY, SEQ_PRIMERS = nextera_library(contact_segments(),
                                              "long-read ChIA-PET library")
FINAL_CAPTION = ("A bridge-containing proximity product after Nextera tagmentation, "
                 "streptavidin selection and indexed PCR.")
SEQUENCING_INTRO = ("Paired reads enter from the two Tn5 ends. A read reaching either "
                    "published bridge-linker orientation proves that the selected fragment "
                    "contains a proximity-ligation junction.")
TAGMENTATION = (55, 5)
MAX_PCR_CYCLES = 13


def contact_scene() -> Scene:
    sc = Scene.duplex(list(CONTACT), label="proximity product")
    sc.junction("top", "left contact tag", "left A:T ligation junction", "ligation")
    sc.junction("top", "right T:A ligation junction", "right contact tag", "ligation")
    sc.mark("top", "internal biotin-dT", "streptavidin handle")
    return sc


def sections():
    return [
        ("Enrich protein-associated chromatin", [
            Row(chunks=[("fragment A — target protein / antibody — fragment B", None, False)])
        ], "Dual-cross-linked chromatin is sonicated and immunoprecipitated with an "
           "antibody to the protein or histone mark of interest."),
        ("Anneal the double-ended bridge linker", BRIDGE.scene().rows(),
         "Both oligos are 5′ phosphorylated. Their reverse-complementary 18-nt cores leave "
         "one 3′ T overhang at each end; the highlighted internal dT carries biotin."),
        ("Join two dA-tailed chromatin fragments", contact_scene().rows(),
         "Each terminal T pairs with a genomic 3′ dA. ** marks the two covalent genomic-"
         "fragment/bridge ligations made in one proximity-ligation reaction."),
        ("Tagment, select the junction and amplify", [
            Row(chunks=[("S5—ME ** genomic tag — biotin bridge — genomic tag ** ME—S7",
                         "me", False)]),
            Row(chunks=[("                         ↓ streptavidin capture", "umi", False)]),
        ], "Nextera Tn5 cuts the ligation product and transfers adapters. Streptavidin "
           "retains only tagmented pieces containing the biotin bridge before indexed PCR."),
    ]


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    chromatin = graph.input(
        "crosslinked_chromatin", "dual-cross-linked chromatin", MaterialKind.BIOLOGICAL,
        state=molecule("dual-cross-linked chromatin", (
            role("protein-associated genomic locus A"),
            role("protein-associated genomic locus B"))), primary=True)
    sonication = graph.input("sonication", "chromatin-sonication condition",
                             MaterialKind.CONDITION)
    fragments = graph.product("chromatin_fragments", "sonicated cross-linked chromatin",
                              MaterialKind.MOLECULAR, state=chromatin.state)
    graph.react("fragment", "Sonicate dual-cross-linked chromatin",
                StepKind.FRAGMENTATION,
                ((chromatin, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)),
                (fragments,))
    antibody = graph.input("target_antibody", "target-specific antibody and capture beads",
                           MaterialKind.SELECTION)
    enriched = graph.product("enriched_chromatin", "protein-associated chromatin fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    flowthrough = graph.product("ip_flowthrough", "unbound chromatin fragments",
                                MaterialKind.WASTE)
    graph.react("immunoprecipitate", "Immunoprecipitate protein-associated chromatin",
                StepKind.AFFINITY_CAPTURE,
                ((fragments, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                (enriched, flowthrough))
    graph.discard(flowthrough, "is not retained by the protein-directed immunoprecipitation")
    repair_mix = graph.input("end_repair_mix", "DNA end-repair enzyme mix",
                             MaterialKind.CATALYST)
    repaired = graph.product("end_repaired", "end-repaired chromatin fragments",
                             MaterialKind.MOLECULAR, state=enriched.state)
    graph.react("end_repair", "Repair immunoprecipitated DNA ends", StepKind.END_REPAIR,
                ((enriched, InputRole.TEMPLATE), (repair_mix, InputRole.CATALYST)),
                (repaired,))
    tailing_polymerase = graph.input("tailing_polymerase", "dA-tailing polymerase",
                                     MaterialKind.CATALYST)
    dntp = graph.input("datp", "dATP tailing nucleotide", MaterialKind.MOLECULAR)
    tailed = graph.product("da_tailed", "3′-dA-tailed chromatin fragments",
                           MaterialKind.MOLECULAR, state=repaired.state)
    graph.react("da_tail", "Add one 3′ dA to repaired chromatin ends", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (dntp, InputRole.STRUCTURAL),
                 (tailing_polymerase, InputRole.CATALYST)), (tailed,))
    bridge = graph.input("bridge_linker", "published double-ended biotin bridge linker",
                         MaterialKind.ADAPTER, state=BRIDGE.scene())
    ligase = graph.input("bridge_ligase", "bridge-linker DNA ligase",
                         MaterialKind.CATALYST)
    contact = graph.product("bridge_contact", CONTACT.name, MaterialKind.MOLECULAR,
                            state=contact_scene())
    graph.react("bridge_ligation", "Join two dA-tailed fragments through the bridge linker",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (bridge, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (contact,),
                geometry=AdapterLigationGeometry.known(
                    "cohesive", BRIDGE.scene(),
                    statement="The published phosphorylated bridge duplex exposes one 3′ dT at each end for ligation to genomic 3′ dA ends."))
    loaded_adapter = graph.input("loaded_nextera_adapter", "Nextera adapter-loaded transposome",
                                 MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagmented_contact", "tagmented bridge-containing contact fragments",
                           MaterialKind.MOLECULAR, state=contact.state)
    tn5_waste = graph.product("tn5_waste", "non-amplifiable tagmentation products",
                              MaterialKind.WASTE)
    graph.react("tagment", "Tagment the bridge-linked contact products",
                StepKind.TAGMENTATION,
                ((contact, InputRole.TEMPLATE), (loaded_adapter, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tn5_waste))
    graph.discard(tn5_waste, "does not retain a library-eligible pair of Tn5 ends")
    streptavidin = graph.input("streptavidin", "streptavidin capture beads",
                               MaterialKind.SELECTION)
    captured = graph.product("captured_contact", "biotin-bridge-bearing tagmented fragments",
                             MaterialKind.MOLECULAR, state=tagged.state)
    nonbridge = graph.product("nonbridge_fragments", "tagmented fragments lacking the biotin bridge",
                              MaterialKind.WASTE)
    graph.react("capture_bridge", "Capture fragments retaining the internal biotin bridge",
                StepKind.AFFINITY_CAPTURE,
                ((tagged, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)),
                (captured, nonbridge))
    graph.discard(nonbridge, "does not contain the biotinylated contact bridge")
    forward = graph.input("pcr_forward", "Nextera indexed PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input("pcr_reverse", "Nextera indexed PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "on-bead PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("sequencing_library", FINAL_LIBRARY.name,
                            MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react("library_pcr", "Index and amplify the captured contact library",
                StepKind.PCR,
                ((captured, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,
                    reason="The selected population contains variable Tn5 insertion sites rather than one exact pre-PCR construct."))
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=tuple(SEQ_PRIMERS),
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "crosslinked_chromatin", "sonication", "fragment", "target_antibody",
        "immunoprecipitate", "end_repair_mix", "end_repair", "tailing_polymerase",
        "datp", "da_tail", "bridge_linker", "bridge_ligase", "bridge_ligation",
        "loaded_nextera_adapter", "tn5", "tagment", "streptavidin",
        "capture_bridge", "pcr_forward", "pcr_reverse", "pcr_polymerase",
        "library_pcr", "sequencing_library"), SOURCE)
    return graph.require_publishable()
