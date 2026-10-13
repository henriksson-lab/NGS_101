"""Henikoff et al. low-salt CUTAC molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Row
from chromatin_epigenetics import antibody_enzyme_rows, tagmented_library, tagmentation_scene
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role

TITLE = "CUTAC — tethered tagmentation of nearby accessible DNA"
NOTES = "01_cutac.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1186/s13059-020-02266-4">Henikoff et al., <i>Genome Biology</i> (2020)</a>.'
SUMMARY = "An antibody tethers pA/G–Tn5 to a transcription-associated chromatin mark; low-salt activation shifts tagmentation into the adjacent nucleosome-depleted accessible region."
CAVEAT = "The link from the antibody-bound nucleosome to the neighboring accessible site is spatial, not covalent. The final adapter structure is canonical Nextera."
FINAL_LIBRARY, SEQ_PRIMERS = tagmented_library("adjacent accessible DNA")
FINAL_CAPTION = "Paired-end library from low-salt, antibody-tethered tagmentation of an adjacent nucleosome-depleted region."
SEQUENCING_INTRO = "Paired reads enter the short accessible-DNA fragment from the two targeted Tn5 insertion sites."

def sections():
    return [
        ("Target a transcription-associated nucleosome", antibody_enzyme_rows("pA/G–Tn5", "loaded Tn5 tethered beside an NDR"),
         "An H3K4me2/3 or RNA-polymerase antibody positions adapter-loaded transposase on chromatin flanking accessible DNA."),
        ("Activate in low salt", [Row(chunks=[("marked nucleosome | low-salt reach → [nucleosome-depleted region] ← low-salt reach | marked nucleosome", None, False)])],
         "Reduced ionic strength promotes tagmentation of the nearby accessible gap rather than ordinary nucleosomal CUT&Tag fragments."),
        ("Tag accessible DNA and amplify", tagmentation_scene("adjacent accessible DNA").rows(),
         "Tn5 cleavage and adapter transfer are coupled; ** marks both adapter–insert boundaries."),
    ]


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    nuclei = graph.input(
        "nuclei", "permeabilized cells or nuclei", MaterialKind.BIOLOGICAL,
        state=molecule("permeabilized chromatin", (
            role("transcription-associated marked nucleosome"),
            role("adjacent nucleosome-depleted accessible DNA"))), primary=True)
    antibody = graph.input("target_antibody", "H3K4me2/3 or transcription-target antibody",
                           MaterialKind.SELECTION)
    transposome = graph.input("loaded_pag_tn5", "adapter-loaded pA/G–Tn5",
                              MaterialKind.ADAPTER)
    pag_tn5 = graph.input("pag_tn5", "pA/G–Tn5 transposase", MaterialKind.CATALYST)
    targeted = graph.product("targeted_chromatin", "antibody-tethered pA/G–Tn5 chromatin",
                             MaterialKind.MOLECULAR, state=nuclei.state)
    wash = graph.product("unbound_targeting_reagents", "unbound targeting reagents",
                         MaterialKind.WASTE)
    graph.react("target_chromatin", "Target pA/G–Tn5 to a transcription-associated nucleosome",
                StepKind.AFFINITY_CAPTURE,
                ((nuclei, InputRole.TEMPLATE), (antibody, InputRole.SELECTION),
                 (transposome, InputRole.STRUCTURAL), (pag_tn5, InputRole.CATALYST)),
                (targeted, wash))
    graph.discard(wash, "is removed before low-salt tagmentation activation")
    low_salt = graph.input("low_salt_magnesium", "low-salt magnesium activation condition",
                           MaterialKind.CONDITION)
    tagged = graph.product("tagmented_accessible_dna", "tagmented adjacent accessible DNA",
                           MaterialKind.MOLECULAR, state=targeted.state)
    tn5_waste = graph.product("tn5_waste", "non-accessible or non-amplifiable Tn5 products",
                              MaterialKind.WASTE)
    graph.react("tagment", "Tagment the adjacent nucleosome-depleted region",
                StepKind.TAGMENTATION,
                ((targeted, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (pag_tn5, InputRole.CATALYST), (low_salt, InputRole.CONDITION)),
                (tagged, tn5_waste))
    graph.discard(tn5_waste, "does not retain the productive accessible-DNA configuration")
    release = graph.input("dna_release", "DNA release and deproteinization condition",
                          MaterialKind.CONDITION)
    recovered = graph.product("recovered_dna", "released tagmented accessible DNA",
                              MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("release_dna", "Release and deproteinize tagmented DNA",
                StepKind.DEPROTEINIZATION,
                ((tagged, InputRole.TEMPLATE), (release, InputRole.CONDITION)),
                (recovered,))
    forward = graph.input("pcr_forward", "Nextera indexed PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input("pcr_reverse", "Nextera indexed PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "indexed-library PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("sequencing_library", FINAL_LIBRARY.name,
                            MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react("library_pcr", "Complete and index the CUTAC library", StepKind.PCR,
                ((recovered, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,
                    reason="Variable Tn5 insertion positions do not define one exact pre-PCR construct."))
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=tuple(SEQ_PRIMERS),
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "nuclei", "target_antibody", "loaded_pag_tn5", "pag_tn5",
        "target_chromatin", "low_salt_magnesium", "tagment", "dna_release",
        "release_dna", "pcr_forward", "pcr_reverse", "pcr_polymerase",
        "library_pcr", "sequencing_library"), SOURCE)
    return graph.require_publishable()
