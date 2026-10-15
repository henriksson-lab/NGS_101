"""Nagano et al. 2013 single-cell Hi-C chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from restriction import BGLII
import illumina as il
import seqprimers as sp
from single_cell_hic import (ContactWorkflow, HISTORICAL_PE_READ2,
                             historical_inline_pe_library)
from chemdraw import Construct, Segment
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, MolecularRelation, ProtocolGraph,
                            SequencingPlan, StepKind)
from protocol_model import molecule, role

WORKFLOW = ContactWorkflow("Nagano single-cell Hi-C", BGLII,
    "after bulk proximity ligation", "biotin-fill", "PCR", capture=True)
TITLE = "Nagano single-cell Hi-C — ligate first, isolate nuclei second"
NOTES = "01_nagano-single-cell-hi-c.html"
CITATION = 'Defining source: <a href="https://doi.org/10.1038/nature12593">Nagano et al., <i>Nature</i> (2013)</a>; detailed protocol: <a href="https://doi.org/10.1038/nprot.2015.127">Nagano et al. (2015)</a>.'
SUMMARY = "BglII contact ends are biotin-filled and ligated in fixed nuclei; individual nuclei are isolated only after ligation."
STEPS = (("Crosslink and digest", "Digest fixed chromatin in intact nuclei with BglII."),
         ("Fill and ligate", "Fill the 5′ overhang with biotin-14-dATP, then proximity-ligate."),
         ("Isolate one nucleus", "Dilute and sort individual ligated nuclei before reversal of crosslinks."),
         ("Capture and amplify", "Fragment DNA, enrich biotin-bearing junctions on streptavidin and PCR-amplify an Illumina library."))
READOUT = "Paired reads start in the two genomic fragments on either side of a captured contact junction."
JUNCTION_CAPTION = "BglII-compatible ends form the contact boundary; the experimental junction is selected through the biotin-bearing fill-in. ** marks ligation."
LIBRARY_CAPTION = "PCR-completed paired-end library with the paper's 3-bp identification tag at each read end. The CAA member of the twelve-adapter set is shown."
JUNCTIONS = (("left filled end", "right filled end", "proximity ligation"),)
def contact_product():
    j = WORKFLOW.digest().fill_in(biotin_base="A").junction
    return Construct([
        Segment("locus A", "X"*22, placeholder=True),
        Segment("junction before top biotin", j[:2], "me"),
        Segment("top-strand biotin-dA", j[2], "w1"),
        Segment("left filled end", j[3:5], "me"),
        Segment("right filled end", j[5:7], "me"),
        Segment("opposite-strand biotin-dA", j[7], "w1"),
        Segment("junction after opposite biotin", j[8:], "me"),
        Segment("locus B", "X"*22, placeholder=True),
    ], name="BglII biotin-filled contact")
SEQ_PRIMERS = (
    sp.custom("Read 1", "historical Illumina Read 1", il.TRUSEQ_READ1,
              "Nagano et al. 2013 Supplementary Methods"),
    sp.custom("Read 2", "historical paired-end Read 2", HISTORICAL_PE_READ2,
              "Nagano et al. 2013 Supplementary Methods"),
)
REQUIRED_ROLES = ("Read 1", "Read 2")
def final_library(): return historical_inline_pe_library(contact_product(), "Nagano scHi-C library")


def protocol_graph() -> ProtocolGraph:
    """Bulk contact ligation precedes the single-nucleus sampling boundary."""
    g = ProtocolGraph(TITLE, authority=CITATION)
    sample = g.input("sample", "fixed intact nuclei", MaterialKind.BIOLOGICAL,
                     state=molecule("fixed chromatin", (role("chromosomal DNA"),)),
                     primary=True)
    enzyme = g.input("bgl2", "BglII", MaterialKind.CATALYST)
    cut = g.product("cut", "BglII-cut chromatin", MaterialKind.MOLECULAR,
                    state=sample.state, chemical_state=WORKFLOW.digest())
    g.react("digest", "Digest fixed chromatin with BglII", StepKind.DIGESTION,
            ((sample, InputRole.TEMPLATE), (enzyme, InputRole.CATALYST)), (cut,))
    biotin = g.input("biotin_datp", "biotin-14-dATP", MaterialKind.MOLECULAR)
    polymerase = g.input("fill_polymerase", "restriction-end fill polymerase",
                         MaterialKind.CATALYST)
    marked = g.product("marked", "biotin-filled BglII ends", MaterialKind.MOLECULAR,
                       state=cut.state, chemical_state=WORKFLOW.digest().fill_in(biotin_base="A"))
    g.react("fill", "Fill BglII overhangs with biotin-dATP", StepKind.EXTENSION,
            ((cut, InputRole.TEMPLATE), (biotin, InputRole.STRUCTURAL),
             (polymerase, InputRole.CATALYST)), (marked,),
            relations=(MolecularRelation.RESTRICTION_FILL,))
    ligase = g.input("contact_ligase", "proximity-ligation DNA ligase",
                     MaterialKind.CATALYST)
    contacts = g.product("contacts", "bulk proximity-ligated contacts",
                         MaterialKind.MOLECULAR, state=contact_product())
    g.react("contact_ligation", "Ligate contacts before isolating individual nuclei",
            StepKind.PROXIMITY_LIGATION,
            ((marked, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)),
            (contacts,), relations=(MolecularRelation.CONTACT_JUNCTION,))
    isolation = g.input("isolation", "single-nucleus dilution and sorting",
                        MaterialKind.CONDITION)
    nucleus = g.product("nucleus", "one ligated nucleus", MaterialKind.MOLECULAR,
                        state=contacts.state)
    excluded = g.product("excluded", "other ligated nuclei", MaterialKind.WASTE)
    g.react("isolate", "Isolate one nucleus after contact ligation", StepKind.SAMPLING,
            ((contacts, InputRole.TEMPLATE), (isolation, InputRole.CONDITION)),
            (nucleus, excluded), changes_molecule=False)
    g.discard(excluded, "outside the retained single nucleus")
    shear = g.input("shear", "contact-DNA fragmentation", MaterialKind.CONDITION)
    fragments = g.product("fragments", "fragmented contact DNA",
                          MaterialKind.MOLECULAR, state=nucleus.state)
    g.react("fragment", "Reverse crosslinks and fragment contact DNA",
            StepKind.FRAGMENTATION,
            ((nucleus, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (fragments,))
    beads = g.input("streptavidin", "streptavidin beads", MaterialKind.SELECTION)
    captured = g.product("captured", "biotin-selected contacts", MaterialKind.MOLECULAR,
                         state=fragments.state)
    unbound = g.product("unbound", "fragments without a captured biotin contact",
                        MaterialKind.WASTE)
    g.react("capture", "Enrich biotin-bearing contact junctions",
            StepKind.AFFINITY_CAPTURE,
            ((fragments, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
            (captured, unbound))
    g.discard(unbound, "lacks a retained biotin-marked junction")
    adapter = g.input("adapter", "3-bp-tagged historical Illumina adapter",
                      MaterialKind.ADAPTER)
    adapter_ligase = g.input("adapter_ligase", "library-adapter ligase",
                             MaterialKind.CATALYST)
    ligated = g.product("ligated", "adapter-ligated contact fragments",
                        MaterialKind.MOLECULAR, state=captured.state)
    g.react("adapter_ligation", "Install historical paired-end adapters",
            StepKind.LIGATION,
            ((captured, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
             (adapter_ligase, InputRole.CATALYST)), (ligated,),
            geometry=AdapterLigationGeometry.unavailable(
                "The printed historical primer sites do not determine the ligation-ready adapter end."))
    forward = g.input("pcr_forward", "library PCR forward primer", MaterialKind.OLIGO)
    reverse = g.input("pcr_reverse", "library PCR reverse primer", MaterialKind.OLIGO)
    pcr = g.input("pcr_polymerase", "library PCR polymerase", MaterialKind.CATALYST)
    final = final_library()
    library = g.product("library", final.name, MaterialKind.LIBRARY, state=final)
    g.react("library_pcr", "Amplify the historical paired-end contact library",
            StepKind.PCR,
            ((ligated, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
             (reverse, InputRole.STRUCTURAL), (pcr, InputRole.CATALYST)),
            (library,), geometry=PCRGeometry.unavailable(
                final, reason="The captured pre-PCR molecule is a heterogeneous contact population."))
    g.sequence(library, SequencingPlan(library=final, primers=SEQ_PRIMERS,
                                       acquired_roles=REQUIRED_ROLES))
    g.declare_reviewed_scope(tuple(g.steps) + tuple(g.materials), CITATION)
    return g.require_publishable()
