"""Flyamer et al. 2017 single-nucleus Hi-C chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from restriction import DPNII
from single_cell_hic import ContactWorkflow, unresolved_illumina_library
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role
WORKFLOW = ContactWorkflow("snHi-C", DPNII, "before fixation", "none", "MDA")
TITLE = "snHi-C — DpnII contact ligation followed by MDA"
NOTES = "01_snhi-c.html"
CITATION = 'Defining source: <a href="https://doi.org/10.1038/nature21711">Flyamer et al., <i>Nature</i> (2017)</a>.'
SUMMARY = "This view follows the isolated-zygote-pronucleus branch: pronuclei are extracted before fixation; DpnII contacts are ligated without biotin enrichment and whole molecules are amplified by phi29 MDA."
STEPS = (("Isolate and fix nuclei", "Prepare intact single nuclei before the restriction reaction."),
         ("Digest and ligate", "DpnII cuts GATC; compatible chromatin ends are proximity-ligated without biotin fill-in."),
         ("Reverse crosslinks", "Purify the complete single-nucleus ligation mixture rather than selecting streptavidin-bound junctions."),
         ("MDA then library prep", "Amplify with phi29 multiple-displacement amplification and convert amplified DNA to an indexed sequencing library."))
READOUT = "Paired reads identify restriction fragments; computational filtering distinguishes contact ligations from other MDA products."
JUNCTION_CAPTION = "A DpnII-compatible contact product. Unlike biotin Hi-C, this protocol does not select a filled, marked junction. ** marks ligation."
LIBRARY_CAPTION = "The contact product passes through MDA, AMPure cleanup, sonication and a second AMPure or QIAquick cleanup before VBCF library preparation with 300–1300 bp size selection. Exact adapter assignment remains unresolved."
LIBRARY_CAVEAT = "INFERRED — the outer library is unresolved in this model; the detailed VBCF preparation and cited GenomiPhi method remain deferred. Whole-zygote preparation is a separate source branch, with nuclei extracted after SDS treatment."
SEQUENCING_UNAVAILABLE = "The source specifies 2 × 125-bp HiSeq 2500 sequencing ; the exact preparation and primer assignment have not been established in this bounded pass."
READ_LENGTHS = {"Read 1": 125, "Read 2": 125}
REQUIRED_ROLES = ("Read 1", "Read 2")
def contact_product(): return WORKFLOW.junction()
def final_library(): return unresolved_illumina_library(contact_product(), "snHi-C MDA-derived library", indexed=None)


def protocol_graph() -> ProtocolGraph:
    """Single nuclei are isolated before DpnII digestion and phi29 amplification."""
    g = ProtocolGraph(TITLE, authority=CITATION)
    nuclei = g.input("nuclei", "prepared individual nuclei", MaterialKind.BIOLOGICAL,
                     state=molecule("chromatin nuclei", (role("chromosomal DNA"),)),
                     primary=True)
    isolation = g.input("isolation", "single-nucleus isolation", MaterialKind.CONDITION)
    single = g.product("single", "one isolated nucleus", MaterialKind.MOLECULAR,
                       state=nuclei.state)
    excluded = g.product("excluded", "other nuclei", MaterialKind.WASTE)
    g.react("isolate", "Isolate one nucleus before fixation and digestion",
            StepKind.SAMPLING,
            ((nuclei, InputRole.TEMPLATE), (isolation, InputRole.CONDITION)),
            (single, excluded), changes_molecule=False)
    g.discard(excluded, "outside the retained single-nucleus reaction")
    dpnii = g.input("dpnii", "DpnII", MaterialKind.CATALYST)
    digested = g.product("digested", "DpnII-cut chromatin", MaterialKind.MOLECULAR,
                         state=single.state, chemical_state=WORKFLOW.digest())
    g.react("digest", "Digest single-nucleus chromatin with DpnII",
            StepKind.DIGESTION,
            ((single, InputRole.TEMPLATE), (dpnii, InputRole.CATALYST)), (digested,))
    ligase = g.input("contact_ligase", "proximity-ligation DNA ligase",
                     MaterialKind.CATALYST)
    contacts = g.product("contacts", "unselected DpnII contact products",
                         MaterialKind.MOLECULAR, state=contact_product())
    g.react("contact_ligation", "Proximity-ligate compatible DpnII ends",
            StepKind.PROXIMITY_LIGATION,
            ((digested, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)),
            (contacts,))
    phi29 = g.input("phi29", "phi29 strand-displacing polymerase",
                    MaterialKind.CATALYST)
    decrosslinking = g.input("decrosslinking", "MDA workflow modified to include overnight 65°C decrosslinking (14–16 h); internal order unassigned", MaterialKind.CONDITION)
    amplified = g.product("amplified", "whole-genome MDA contact mixture",
                          MaterialKind.MOLECULAR, state=contacts.state)
    g.react("mda", "Amplify the entire single-nucleus ligation mixture",
            StepKind.EXTENSION,
            ((contacts, InputRole.TEMPLATE), (phi29, InputRole.CATALYST),
             (decrosslinking, InputRole.CONDITION)),
            (amplified,))
    ampure = g.input("ampure", "AMPure XP post-MDA cleanup", MaterialKind.SELECTION)
    recovered = g.product("recovered", "cleaned MDA products", MaterialKind.MOLECULAR, state=amplified.state)
    mda_waste = g.product("mda_waste", "post-MDA cleanup flow-through and washes; composition unassigned", MaterialKind.WASTE)
    g.react("mda_cleanup", "Clean MDA products with AMPure XP", StepKind.CLEANUP,
            ((amplified, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)), (recovered, mda_waste), changes_molecule=False)
    g.discard(mda_waste, "removed during AMPure cleanup")
    sonication = g.input("sonication", "sonication", MaterialKind.CONDITION)
    fragmented = g.product("fragmented", "sonicated MDA DNA", MaterialKind.MOLECULAR)
    g.react("sonicate", "Sonicate cleaned MDA DNA", StepKind.FRAGMENTATION,
            ((recovered, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)), (fragmented,))
    cleanup = g.input("cleanup", "AMPure XP or QIAquick cleanup (source alternatives)", MaterialKind.SELECTION)
    cleaned = g.product("cleaned", "cleaned sonicated DNA", MaterialKind.MOLECULAR)
    fragment_waste = g.product("fragment_waste", "post-sonication cleanup flow-through and washes; composition unassigned", MaterialKind.WASTE)
    g.react("fragment_cleanup", "Clean sonicated DNA", StepKind.CLEANUP,
            ((fragmented, InputRole.TEMPLATE), (cleanup, InputRole.SELECTION)), (cleaned, fragment_waste), changes_molecule=False)
    g.discard(fragment_waste, "removed during cleanup")
    final = final_library()
    library = g.product("library", final.name, MaterialKind.LIBRARY,
                        state=final, inferred=True, source=LIBRARY_CAVEAT)
    g.react("library_prep", "Convert MDA material to paired-end library",
            StepKind.MAGIC, ((cleaned, InputRole.TEMPLATE),), (library,),
            reason="VBCF library preparation with 300–1300 bp selection is source-described; exact transition order, adapter and primer assignment remain unverified in this model.")
    g.sequence(library, SequencingPlan(acquired_roles=("Read 1", "Read 2"),
                                       read_lengths=READ_LENGTHS, limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(tuple(g.steps) + tuple(g.materials), CITATION)
    return g.require_publishable()
