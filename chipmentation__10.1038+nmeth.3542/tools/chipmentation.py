"""Schmidl et al. ChIPmentation molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Row
from chromatin_epigenetics import antibody_enzyme_rows, tagmented_library, tagmentation_scene
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role

TITLE = "ChIPmentation — tagmentation of bead-bound ChIP DNA"
NOTES = "01_chipmentation.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nmeth.3542">Schmidl et al., <i>Nature Methods</i> (2015)</a>.'
SUMMARY = "Immunoprecipitate crosslinked chromatin first, then let Tn5 fragment and adapter-tag the bead-bound ChIP DNA in a single reaction before indexed PCR."
CAVEAT = "The defining chemistry uses Nextera reagents; the final molecule is assembled from the repository’s canonical Nextera components."
FINAL_LIBRARY, SEQ_PRIMERS = tagmented_library("ChIP-selected DNA")
FINAL_CAPTION = "Indexed Nextera library whose insert was selected by chromatin immunoprecipitation before tagmentation."
SEQUENCING_INTRO = ("Read primers can enter the immunoprecipitated fragment from either Tn5 insertion site; "
                    "the published sequencing summary includes both single-read and paired-end runs. "
                    "Sample indexes are installed during PCR.")

def sections():
    return [
        ("Fragment and immunoprecipitate chromatin", antibody_enzyme_rows("bead capture", "retain target-associated chromatin"),
         "Crosslinked, sonicated chromatin is enriched with a target-specific antibody."),
        ("Tagment on beads", tagmentation_scene("ChIP-selected DNA").rows(),
         "Tn5 replaces conventional end repair, dA addition and adapter ligation; ** marks both adapter-transfer boundaries."),
        ("Reverse crosslinks and amplify", [Row(chunks=[("bead-bound tagged DNA → proteinase/reverse crosslink → indexed PCR", None, False)])],
         "PCR adds complete flow-cell adapters and sample indexes to productive two-ended fragments."),
    ]


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    chromatin = graph.input(
        "crosslinked_chromatin", "cross-linked chromatin", MaterialKind.BIOLOGICAL,
        state=molecule("cross-linked chromatin", (role("target-associated chromatin"),)),
        primary=True)
    sonication = graph.input("sonication", "chromatin-sonication condition",
                             MaterialKind.CONDITION)
    fragments = graph.product("chromatin_fragments", "sonicated cross-linked chromatin",
                              MaterialKind.MOLECULAR, state=chromatin.state)
    graph.react("fragment", "Sonicate cross-linked chromatin", StepKind.FRAGMENTATION,
                ((chromatin, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)),
                (fragments,))
    antibody = graph.input("target_antibody", "target-specific antibody and capture beads",
                           MaterialKind.SELECTION)
    selected = graph.product("selected_chromatin", "bead-bound ChIP-selected chromatin",
                             MaterialKind.MOLECULAR, state=fragments.state)
    flowthrough = graph.product("ip_flowthrough", "unbound chromatin fragments",
                                MaterialKind.WASTE)
    graph.react("immunoprecipitate", "Immunoprecipitate target-associated chromatin",
                StepKind.AFFINITY_CAPTURE,
                ((fragments, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                (selected, flowthrough))
    graph.discard(flowthrough, "is not retained by the target-directed immunoprecipitation")
    transposome = graph.input("loaded_nextera_adapter", "Nextera adapter-loaded transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Nextera Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagmented_chromatin", "bead-bound tagmented ChIP DNA",
                           MaterialKind.MOLECULAR, state=selected.state)
    tn5_waste = graph.product("tn5_waste", "non-amplifiable tagmentation products",
                              MaterialKind.WASTE)
    graph.react("tagment", "Tagment bead-bound ChIP DNA", StepKind.TAGMENTATION,
                ((selected, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tn5_waste))
    graph.discard(tn5_waste, "does not retain a library-eligible pair of Nextera ends")
    proteinase = graph.input("proteinase", "proteinase and cross-link reversal condition",
                             MaterialKind.CONDITION)
    released = graph.product("released_tagged_dna", "de-cross-linked tagmented ChIP DNA",
                             MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("release_dna", "Digest protein and reverse cross-links",
                StepKind.DEPROTEINIZATION,
                ((tagged, InputRole.TEMPLATE), (proteinase, InputRole.CONDITION)),
                (released,))
    forward = graph.input("pcr_forward", "Nextera indexed PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input("pcr_reverse", "Nextera indexed PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "indexed-library PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("sequencing_library", FINAL_LIBRARY.name,
                            MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react("library_pcr", "Complete and index the ChIPmentation library",
                StepKind.PCR,
                ((released, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,
                    reason="Variable Tn5 insertion positions do not define one exact pre-PCR construct."))
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=tuple(SEQ_PRIMERS),
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "crosslinked_chromatin", "sonication", "fragment", "target_antibody",
        "immunoprecipitate", "loaded_nextera_adapter", "tn5", "tagment",
        "proteinase", "release_dna", "pcr_forward", "pcr_reverse",
        "pcr_polymerase", "library_pcr", "sequencing_library"), SOURCE)
    return graph.require_publishable()
