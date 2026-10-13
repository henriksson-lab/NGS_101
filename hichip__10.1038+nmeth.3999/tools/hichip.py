"""Mumbach et al. 2016 HiChIP molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Row
from chromatin_epigenetics import proximity_rows, tagmented_library, tagmentation_scene
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role

TITLE = "HiChIP — protein-directed chromatin contacts"
NOTES = "01_hichip.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nmeth.3999">Mumbach et al., <i>Nature Methods</i> (2016)</a>.'
SUMMARY = "Build biotin-marked in situ Hi-C contacts, immunoprecipitate contacts associated with a chosen protein, and tagment the bead-bound DNA into an Illumina library."
CAVEAT = "The contact panel derives the original MboI fill and GATCGATC ligation scar. The opposite-strand biotin is present in the duplex but only the top-strand biotin is separately coloured; the final adapter structure is canonical Nextera."
FINAL_LIBRARY, SEQ_PRIMERS = tagmented_library("protein-enriched contact fragment")
FINAL_CAPTION = "Tn5-built paired-end library from immunoprecipitated, biotin-bearing proximity-ligation products."
SEQUENCING_INTRO = "Paired reads recover the two genomic ends of each protein-enriched contact molecule; an internal ligation scar identifies informative contacts."
SHOW_READ_CONTENT = True

def sections():
    return [
        ("Create in situ proximity ligations", proximity_rows(),
         "MboI ends are filled with biotin-dATP and ligated in intact nuclei; ** marks the exact boundary in the derived GATCGATC contact scar."),
        ("Immunoprecipitate protein-associated contacts", [Row(chunks=[("contact chromatin + target antibody → protein-directed bead capture", "cbc", False)])],
         "Chromatin is sheared after proximity ligation and enriched with an antibody to the factor or histone mark of interest."),
        ("Tagment bead-bound DNA", tagmentation_scene("protein-enriched contact fragment").rows(),
         "Tn5 installs sequencing ends directly on immunoprecipitated chromatin; ** marks the two adapter-transfer boundaries."),
    ]


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    cells = graph.input("cells", "input cells", MaterialKind.BIOLOGICAL,
                        state=molecule("cells", (role("native chromatin contacts"),)),
                        primary=True)
    formaldehyde = graph.input("formaldehyde", "formaldehyde cross-linking condition",
                               MaterialKind.CONDITION)
    fixed = graph.product("fixed_chromatin", "cross-linked chromatin",
                          MaterialKind.MOLECULAR, state=cells.state)
    graph.react("crosslink", "Cross-link chromatin contacts", StepKind.CROSS_LINKING,
                ((cells, InputRole.TEMPLATE), (formaldehyde, InputRole.CONDITION)),
                (fixed,))
    mboi = graph.input("mboi", "MboI restriction endonuclease", MaterialKind.CATALYST)
    digested = graph.product("digested_chromatin", "MboI-digested chromatin",
                             MaterialKind.MOLECULAR, state=fixed.state)
    graph.react("digest", "Digest fixed chromatin with MboI", StepKind.DIGESTION,
                ((fixed, InputRole.TEMPLATE), (mboi, InputRole.CATALYST)), (digested,))
    biotin_datp = graph.input("biotin_datp", "biotin-dATP fill-in nucleotide",
                              MaterialKind.MOLECULAR)
    fill_polymerase = graph.input("fill_polymerase", "restriction-end fill-in polymerase",
                                  MaterialKind.CATALYST)
    filled = graph.product("filled_ends", "biotin-filled MboI restriction ends",
                           MaterialKind.MOLECULAR, state=digested.state)
    graph.react("fill_ends", "Fill restriction ends with biotin-dATP",
                StepKind.EXTENSION,
                ((digested, InputRole.TEMPLATE), (biotin_datp, InputRole.STRUCTURAL),
                 (fill_polymerase, InputRole.CATALYST)), (filled,))
    ligase = graph.input("proximity_ligase", "in-situ proximity DNA ligase",
                         MaterialKind.CATALYST)
    contacts = graph.product(
        "contacts", "biotin-bearing in-situ proximity-ligation products",
        MaterialKind.MOLECULAR,
        state=molecule("HiChIP contact chromatin", (
            role("protein-associated chromatin"), role("biotin-marked GATCGATC junction"))))
    graph.react("proximity_ligation", "Ligate spatially adjacent filled restriction ends",
                StepKind.PROXIMITY_LIGATION,
                ((filled, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (contacts,))
    sonication = graph.input("sonication", "chromatin-sonication condition",
                             MaterialKind.CONDITION)
    sheared = graph.product("sheared_contacts", "sheared proximity-ligated chromatin",
                            MaterialKind.MOLECULAR, state=contacts.state)
    graph.react("shear", "Shear proximity-ligated chromatin", StepKind.FRAGMENTATION,
                ((contacts, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)),
                (sheared,))
    antibody = graph.input("target_antibody", "target-specific antibody and capture beads",
                           MaterialKind.SELECTION)
    enriched = graph.product("protein_contacts", "protein-associated contact chromatin",
                             MaterialKind.MOLECULAR, state=sheared.state)
    flowthrough = graph.product("ip_flowthrough", "unbound chromatin fragments",
                                MaterialKind.WASTE)
    graph.react("immunoprecipitate", "Immunoprecipitate protein-associated contacts",
                StepKind.AFFINITY_CAPTURE,
                ((sheared, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                (enriched, flowthrough))
    graph.discard(flowthrough, "is not retained by the target-directed immunoprecipitation")
    loaded_adapter = graph.input("loaded_nextera_adapter", "Nextera adapter-loaded transposome",
                                 MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagmented_contacts", "on-bead tagmented contact fragments",
                           MaterialKind.MOLECULAR, state=enriched.state)
    tn5_waste = graph.product("tn5_waste", "non-amplifiable tagmentation products",
                              MaterialKind.WASTE)
    graph.react("tagment", "Tagment bead-bound protein-associated contacts",
                StepKind.TAGMENTATION,
                ((enriched, InputRole.TEMPLATE), (loaded_adapter, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tn5_waste))
    graph.discard(tn5_waste, "does not carry the two library-eligible Tn5 ends")
    streptavidin = graph.input("streptavidin", "streptavidin capture beads",
                               MaterialKind.SELECTION)
    captured = graph.product("captured_contacts", "biotin-junction-bearing tagmented contacts",
                             MaterialKind.MOLECULAR, state=tagged.state)
    nonjunction = graph.product("nonjunction_fragments", "tagmented fragments without retained biotin junctions",
                                MaterialKind.WASTE)
    graph.react("capture_biotin", "Capture biotin-bearing contact products",
                StepKind.AFFINITY_CAPTURE,
                ((tagged, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)),
                (captured, nonjunction))
    graph.discard(nonjunction, "does not retain a biotin-marked proximity junction")
    forward = graph.input("pcr_forward", "Nextera index PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input("pcr_reverse", "Nextera index PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "index PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("sequencing_library", FINAL_LIBRARY.name,
                            MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react("library_pcr", "Complete and index the tagmented contact library",
                StepKind.PCR,
                ((captured, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,
                    reason="The pre-PCR on-bead contact population is not represented as one exact construct."))
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=tuple(SEQ_PRIMERS),
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "cells", "formaldehyde", "crosslink", "mboi", "digest", "biotin_datp",
        "fill_polymerase", "fill_ends", "proximity_ligase", "proximity_ligation",
        "sonication", "shear", "target_antibody", "immunoprecipitate",
        "loaded_nextera_adapter", "tn5", "tagment", "streptavidin",
        "capture_biotin", "pcr_forward", "pcr_reverse", "pcr_polymerase",
        "library_pcr", "sequencing_library"), SOURCE)
    return graph.require_publishable()
