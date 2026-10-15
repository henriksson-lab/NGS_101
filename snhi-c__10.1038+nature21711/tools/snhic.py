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
SUMMARY = "Nuclei are isolated before fixation; DpnII contacts are ligated without biotin enrichment and whole molecules are amplified by phi29 MDA."
STEPS = (("Isolate and fix nuclei", "Prepare intact single nuclei before the restriction reaction."),
         ("Digest and ligate", "DpnII cuts GATC; compatible chromatin ends are proximity-ligated without biotin fill-in."),
         ("Reverse crosslinks", "Purify the complete single-nucleus ligation mixture rather than selecting streptavidin-bound junctions."),
         ("MDA then library prep", "Amplify with phi29 multiple-displacement amplification and convert amplified DNA to an indexed sequencing library."))
READOUT = "Paired reads identify restriction fragments; computational filtering distinguishes contact ligations from other MDA products."
JUNCTION_CAPTION = "A DpnII-compatible contact product. Unlike biotin Hi-C, this protocol does not select a filled, marked junction. ** marks ligation."
LIBRARY_CAPTION = "The contact product passes through phi29 MDA before paired-end library conversion. The supplementary methods do not identify the library kit or adapter sequences."
LIBRARY_CAVEAT = "INFERRED — only the presence of an unspecified Illumina library is established; adapter sequence and index architecture are not reported."
SEQUENCING_UNAVAILABLE = "The source specifies 2 × 125-bp HiSeq 2500 sequencing but not the library kit or sequencing-primer sequences, so primer binding cannot be placed base by base."
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
    amplified = g.product("amplified", "whole-genome MDA contact mixture",
                          MaterialKind.MOLECULAR, state=contacts.state)
    g.react("mda", "Amplify the entire single-nucleus ligation mixture",
            StepKind.EXTENSION,
            ((contacts, InputRole.TEMPLATE), (phi29, InputRole.CATALYST)),
            (amplified,))
    final = final_library()
    library = g.product("library", final.name, MaterialKind.LIBRARY,
                        state=final, inferred=True, source=LIBRARY_CAVEAT)
    g.react("library_prep", "Convert MDA material to paired-end library",
            StepKind.MAGIC, ((amplified, InputRole.TEMPLATE),), (library,),
            reason="The supplementary methods do not identify the library kit or adapter and primer sequences, so the exact library-construction transitions cannot be reconstructed.")
    g.sequence(library, SequencingPlan(acquired_roles=("Read 1", "Read 2"),
                                       limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(tuple(g.steps) + tuple(g.materials), CITATION)
    return g.require_publishable()
