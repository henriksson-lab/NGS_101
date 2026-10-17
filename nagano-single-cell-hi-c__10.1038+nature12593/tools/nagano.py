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
from chemdraw import Construct, Scene, Segment, bases_pair, revcomp
from adapter_ligation import AdapterLigationGeometry
from endprep import dA_tailed_scene, repair_and_dA_tail
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, MolecularRelation, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal)
from protocol_model import molecule, role

WORKFLOW = ContactWorkflow("Nagano single-cell Hi-C", BGLII,
    "after bulk proximity ligation", "biotin-fill", "PCR", capture=True)
TITLE = "Nagano single-cell Hi-C — ligate first, isolate nuclei second"
NOTES = "01_nagano-single-cell-hi-c.html"
CITATION = 'Defining source: <a href="https://doi.org/10.1038/nature12593">Nagano et al., <i>Nature</i> (2013)</a>; <a href="https://www.nature.com/nature/journal/v502/n7469/extref/nature12593-s1.pdf">original supplementary methods</a>; detailed protocol: <a href="https://doi.org/10.1038/nprot.2015.127">Nagano et al. (2015)</a>.'
SUMMARY = "BglII contact ends are biotin-filled and ligated in fixed nuclei; individual nuclei are isolated only after ligation."
STEPS = (("Crosslink and digest", "Digest fixed chromatin in intact nuclei with BglII."),
         ("Fill and ligate", "Fill the 5′ overhang with biotin-14-dATP, then proximity-ligate."),
         ("Isolate one nucleus", "Pick individual ligated nuclei under a microscope with mouth-controlled glass capillaries before reversal of crosslinks."),
         ("Capture and amplify", "Capture biotin-bearing junctions, cut with AluI, add a 3′ dA, ligate the printed Y-adapter, and PCR-amplify."))
READOUT = "Each paired read starts with a 3-bp inline identification tag, then reads one genomic fragment flanking the captured contact junction."
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
READ_LENGTHS = {"Read 1": 40, "Read 2": 40}
PCR_FORWARD = il.P5 + il.TRUSEQ_READ1[4:]
PCR_REVERSE = il.P7 + HISTORICAL_PE_READ2


def pcr_primer_scene(name, sequence):
    scene = Scene()
    scene.strand(name, [Segment(name, sequence)], mod3="terminal phosphorothioate")
    return scene

# Reuse only the endpoint: original 2013 SI pp3–5 independently establishes these
# CAA-tagged arms and both historical PCR primers; no later reaction settings are inherited.
def final_library(): return historical_inline_pe_library(contact_product(), "Nagano scHi-C library")


ALUI_CONTACT = contact_product()
ALUI_TAILED = repair_and_dA_tail(ALUI_CONTACT)


def adapter_geometry() -> AdapterLigationGeometry:
    """The CAA member of the source-written adapter, before ligation."""
    if not bases_pair("T", ALUI_TAILED.three_prime_overhang):
        raise ValueError("published adapter 3′ dT cannot pair with the prepared contact end")
    top_arm = "ACACTCTTTCCCTACACGA"
    paired = "CGCTCTTCCGATCTCAA"
    bottom_arm = "GTTCAGCAGGAATGCCGAG"
    sc = Scene()
    sc.strand("F", [Segment("F Y arm", top_arm),
                    Segment("paired core", paired),
                    Segment("3′ dT overhang", "T")], label="CAA adapter F",
              mod3="terminal phosphorothioate")
    sc.anneal("R", [Segment("R paired core", revcomp(paired)),
                    Segment("R Y arm", bottom_arm)], to="F",
              pair=("R paired core", "paired core"), label="CAA adapter R",
              unpaired=("R Y arm",), mod5="phosphate")
    sc.ligation_end("F", "3′ dT overhang", kind="cohesive",
                    polarity="3-prime", owner="adapter",
                    text="3′ dT ligates to the dA-tailed contact fragment")
    return AdapterLigationGeometry.known(
        "cohesive", dA_tailed_scene(ALUI_TAILED, label="AluI-cut contact fragment"), sc,
        statement="2013 Supplementary Information pp. 2–5 prints dA tailing and both CAA adapter strands, including the 3′ T overhang, 5′ phosphate and terminal phosphorothioate.")


LIBRARY_ADAPTER_GEOMETRY = adapter_geometry()
LIBRARY_ADAPTER_CAPTION = ("Source-printed CAA Y-adapter before ligation: the 17-bp core pairs, "
                           "the sequencing arms do not pair, and the terminal 3′ dT meets the fragment's 3′ dA.")


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
    polymerase = g.input("fill_polymerase", "DNA Polymerase I Large (Klenow) Fragment",
                         MaterialKind.CATALYST)
    marked = g.product("marked", "biotin-filled BglII ends", MaterialKind.MOLECULAR,
                       state=cut.state, chemical_state=WORKFLOW.digest().fill_in(biotin_base="A"))
    g.react("fill", "Fill BglII overhangs with biotin-dATP", StepKind.EXTENSION,
            ((cut, InputRole.TEMPLATE), (biotin, InputRole.STRUCTURAL),
             (polymerase, InputRole.CATALYST)), (marked,),
            relations=(MolecularRelation.RESTRICTION_FILL,))
    ligase = g.input("contact_ligase", "T4 DNA ligase",
                     MaterialKind.CATALYST)
    contacts = g.product("contacts", "bulk proximity-ligated contacts",
                         MaterialKind.MOLECULAR, state=contact_product())
    g.react("contact_ligation", "Ligate contacts before isolating individual nuclei",
            StepKind.PROXIMITY_LIGATION,
            ((marked, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)),
            (contacts,), relations=(MolecularRelation.CONTACT_JUNCTION,))
    isolation = g.input("isolation", "microscope picking with mouth-controlled glass capillaries",
                        MaterialKind.CONDITION)
    nucleus = g.product("nucleus", "one ligated nucleus", MaterialKind.MOLECULAR,
                        state=contacts.state)
    excluded = g.product("excluded", "other ligated nuclei", MaterialKind.MOLECULAR)
    g.react("isolate", "Isolate one nucleus after contact ligation", StepKind.SAMPLING,
            ((contacts, InputRole.TEMPLATE), (isolation, InputRole.CONDITION)),
            (nucleus, excluded), changes_molecule=False)
    g.mark_terminal(Terminal(excluded.id, Disposition.UNRESOLVED,
        "Other nuclei lie outside this picked-nucleus lineage; their experimental fate is not represented."))
    reversal = g.input("reversal", "crosslink reversal at 65 °C", MaterialKind.CONDITION)
    fragments = g.product("fragments", "crosslink-reversed contact DNA",
                          MaterialKind.MOLECULAR, state=nucleus.state)
    g.react("reverse_crosslinks", "Reverse crosslinks after single-nucleus isolation",
            StepKind.DEPROTEINIZATION,
            ((nucleus, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)), (fragments,))
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
    alui = g.input("alui", "AluI", MaterialKind.CATALYST)
    cut_contacts = g.product("alui_cut", "AluI-cut captured contacts",
                             MaterialKind.MOLECULAR, state=ALUI_CONTACT)
    g.react("alui_digest", "Cut captured contacts with AluI", StepKind.DIGESTION,
            ((captured, InputRole.TEMPLATE), (alui, InputRole.CATALYST)),
            (cut_contacts,))
    datp = g.input("datp", "dATP", MaterialKind.MOLECULAR)
    klenow = g.input("klenow_exo", "Klenow fragment (3′→5′ exo−)", MaterialKind.CATALYST)
    tailed = g.product("tailed", "3′ dA-tailed contact fragments",
                       MaterialKind.MOLECULAR, state=ALUI_TAILED)
    g.react("da_tail", "Add a 3′ dA to AluI-cut contact fragments", StepKind.TAILING,
            ((cut_contacts, InputRole.TEMPLATE), (datp, InputRole.STRUCTURAL),
             (klenow, InputRole.CATALYST)), (tailed,))
    adapter = g.input("adapter", "3-bp-tagged historical Illumina Y-adapter",
                      MaterialKind.ADAPTER)
    adapter_ligase = g.input("adapter_ligase", "T4 DNA ligase for adapter ligation",
                             MaterialKind.CATALYST)
    ligated = g.product("ligated", "adapter-ligated contact fragments",
                        MaterialKind.MOLECULAR, state=tailed.state)
    g.react("adapter_ligation", "Install historical paired-end adapters",
            StepKind.LIGATION,
            ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
             (adapter_ligase, InputRole.CATALYST)), (ligated,),
            geometry=LIBRARY_ADAPTER_GEOMETRY)
    forward = g.input("pcr_forward", "source-printed library PCR forward primer", MaterialKind.OLIGO,
                      state=pcr_primer_scene("PCR forward", PCR_FORWARD))
    reverse = g.input("pcr_reverse", "source-printed library PCR reverse primer", MaterialKind.OLIGO,
                      state=pcr_primer_scene("PCR reverse", PCR_REVERSE))
    pcr = g.input("pcr_polymerase", "Platinum Pfx DNA polymerase", MaterialKind.CATALYST)
    final = final_library()
    library = g.product("amplified", final.name, MaterialKind.MOLECULAR, state=final)
    g.react("library_pcr", "Amplify for 25 cycles with Platinum Pfx",
            StepKind.PCR,
            ((ligated, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
             (reverse, InputRole.STRUCTURAL), (pcr, InputRole.CATALYST)),
            (library,), geometry=PCRGeometry.unavailable(
                final, reason="The captured pre-PCR molecule is a heterogeneous contact population."))
    ampure = g.input("ampure", "AMPure XP magnetic beads after removal of streptavidin beads", MaterialKind.SELECTION)
    cleaned = g.product("cleaned", "AMPure-purified amplified contact library", MaterialKind.MOLECULAR, state=final)
    cleanup_waste = g.product("cleanup_waste", "material not recovered from post-PCR AMPure purification", MaterialKind.WASTE)
    g.react("post_pcr_cleanup", "Remove streptavidin beads and purify PCR library with AMPure XP",
            StepKind.CLEANUP, ((library, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)),
            (cleaned, cleanup_waste))
    g.discard(cleanup_waste, "not recovered in the library eluate; exact molecular composition not assigned")
    gel = g.input("size_gel", "2% agarose gel, 300–700 bp recovery window", MaterialKind.SELECTION)
    selected = g.product("library", final.name, MaterialKind.LIBRARY, state=final)
    offsize = g.product("offsize", "products outside the recovered gel interval", MaterialKind.WASTE)
    g.react("size_selection", "Recover 300–700 bp libraries from 2% agarose gel",
            StepKind.SELECTION, ((cleaned, InputRole.TEMPLATE), (gel, InputRole.SELECTION)),
            (selected, offsize))
    g.discard(offsize, "outside the recovered 300–700 bp gel interval")
    g.sequence(selected, SequencingPlan(library=final, primers=SEQ_PRIMERS,
                                       acquired_roles=REQUIRED_ROLES, read_lengths=READ_LENGTHS))
    g.declare_reviewed_scope(tuple(g.steps) + tuple(g.materials), CITATION)
    return g.require_publishable()
