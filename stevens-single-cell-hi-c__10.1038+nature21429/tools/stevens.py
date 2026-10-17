"""Stevens et al. 2017 single-cell Hi-C chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from restriction import MBOI
from single_cell_hic import (ContactWorkflow, HISTORICAL_PE_READ2,
                             historical_inline_pe_library)
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Row, Scene, Workflow, revcomp
from endprep import dA_tailed_scene, repair_and_dA_tail
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, MolecularRelation, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal)
from adapter_ligation import AdapterLigationGeometry
from protocol_model import molecule, role
WORKFLOW = ContactWorkflow("Stevens single-cell Hi-C", MBOI,
    "before digestion and ligation", "biotin-fill", "PCR", capture=True)
TITLE = "Stevens single-cell Hi-C — single nuclei before contact ligation"
NOTES = "01_stevens-single-cell-hi-c.html"
CITATION = 'Defining source: <a href="https://doi.org/10.1038/nature21429">Stevens et al., <i>Nature</i> (2017)</a>.'
SUMMARY = "Individual nuclei are isolated first; MboI contact ends are biotin-filled and proximity-ligated. AluI is used later to fragment captured DNA for library preparation."
STEPS = (("Isolate nuclei", "Lyse fixed cells and distribute individual nuclei before library construction."),
         ("Contact digest", "MboI cuts GATC in each isolated nucleus."),
         ("Mark and ligate", "Fill MboI ends with biotin-14-dATP, then ligate spatial neighbours."),
         ("Capture and library prep", "Capture contacts, digest captured DNA with AluI, A-tail, ligate a 3-bp-tagged Illumina adapter and PCR. Nextera XT was a published alternative."))
READOUT = "Paired reads report the two MboI fragments joined within one nucleus; the first three bases identify the library adapter."
JUNCTION_CAPTION = "Biotin-filled MboI ends form the contact boundary. AluI acts only after contact capture. ** marks ligation."
LIBRARY_CAPTION = "The source-printed adapter-ligation route, with the CAA 3-bp identification tag shown. The paper also reports a dual-index Nextera XT alternative."

ALUI_INSERT = Construct([seg("AluI contact fragment", "X"*38, placeholder=True)])
ALUI_TAILED = repair_and_dA_tail(ALUI_INSERT)
ADAPTER_F = il.TRUSEQ_READ1 + "CAA" + "T"
ADAPTER_R = revcomp(ADAPTER_F[-18:-1]) + "GTTCAGCAGGAATGCCGAG"


def adapter_scene():
    # Independently printed in Stevens2017 SI pp9–10; no reaction choices inherited.
    paired = ADAPTER_F[-18:-1]
    sc = Scene()
    sc.strand("F", [seg("unpaired F arm", ADAPTER_F[:-18]), seg("paired CAA core", paired),
                    seg("3-prime dT", "T")], mod3="terminal phosphorothioate")
    sc.anneal("R", [seg("paired R core", revcomp(paired)), seg("unpaired R arm", ADAPTER_R[17:])],
              to="F", pair=("paired R core", "paired CAA core"),
              unpaired=("unpaired R arm",), mod5="phosphate")
    sc.ligation_end("F", "3-prime dT", kind="cohesive", polarity="3-prime",
                    owner="adapter", text="source-printed 3-prime dT meets dA-tailed fragment")
    return sc


def adapter_geometry():
    return AdapterLigationGeometry.known("cohesive", dA_tailed_scene(ALUI_TAILED), adapter_scene(),
        statement="Stevens2017 SI prints the CAA-tagged adapter pair, including its 17 bp core, 5-prime phosphate and terminal phosphorothioate.")

def adapter_workflow():
    w=Workflow(MolecularState("AluI-digested captured contacts",tuple(Scene.duplex(list(ALUI_INSERT)).rows())))
    w.react("Add 3′ dA with Klenow exo-minus",dA_tailed_scene(ALUI_TAILED).rows(),name="dA-tailed contact fragments",note="The source-printed ligation route dA-tails the AluI fragments.")
    geometry=adapter_geometry()
    linked = [Row(chunks=[("CAA-adapter-linked contact DNA; forked pre-PCR ends not reconstructed", None, False)])]
    w.ligate_adapter("Ligate the 3-bp-tagged Illumina adapter",geometry,linked,name="adapter-ligated contact DNA",note="PCR-completed arms appear only in the final-library panel. Nextera XT is an alternative route.")
    return w
def contact_product(): return WORKFLOW.junction()
SEQ_PRIMERS = (
    sp.custom("Read 1", "historical Illumina Read 1", il.TRUSEQ_READ1,
              "Stevens et al. 2017 Supplementary Methods"),
    sp.custom("Read 2", "historical paired-end Read 2", HISTORICAL_PE_READ2,
              "Stevens et al. 2017 Supplementary Methods"),
)
REQUIRED_ROLES = ("Read 1", "Read 2")
READ_LENGTHS = {"Read 1": 75, "Read 2": 75}
PCR_FORWARD = il.P5 + il.TRUSEQ_READ1[4:]
PCR_REVERSE = il.P7 + HISTORICAL_PE_READ2


def primer_scene(name, sequence):
    sc = Scene()
    sc.strand(name, [seg(name, sequence)], mod3="terminal phosphorothioate")
    return sc

def final_library(): return historical_inline_pe_library(contact_product(), "Stevens scHi-C library")


def protocol_graph() -> ProtocolGraph:
    """Stevens' single-nucleus MboI route with AluI only after contact capture."""
    g = ProtocolGraph(TITLE, authority=CITATION)
    sample = g.input("sample", "fixed chromatin nuclei", MaterialKind.BIOLOGICAL,
                     state=molecule("fixed chromatin", (role("chromosomal DNA"),)),
                     primary=True)
    isolation = g.input("isolation", "single-nucleus distribution", MaterialKind.CONDITION)
    nucleus = g.product("nucleus", "one isolated nucleus", MaterialKind.MOLECULAR,
                        state=sample.state)
    excluded = g.product("excluded", "other nuclei", MaterialKind.MOLECULAR)
    g.react("isolate", "Distribute individual nuclei before library construction",
            StepKind.SAMPLING,
            ((sample, InputRole.TEMPLATE), (isolation, InputRole.CONDITION)),
            (nucleus, excluded), changes_molecule=False)
    g.mark_terminal(Terminal(excluded.id, Disposition.UNRESOLVED, "Other nuclei are outside this retained single-nucleus lineage; their disposition is unrepresented."))
    agarose = g.input("agarose", "1% low-melting-point agarose overlay after imaging", MaterialKind.CONDITION)
    mboi = g.input("mboi", "MboI", MaterialKind.CATALYST)
    cut = g.product("cut", "MboI-cut chromatin", MaterialKind.MOLECULAR,
                    state=nucleus.state, chemical_state=WORKFLOW.digest())
    g.react("digest", "Digest single-nucleus chromatin with MboI",
            StepKind.DIGESTION,
            ((nucleus, InputRole.TEMPLATE), (mboi, InputRole.CATALYST), (agarose, InputRole.CONDITION)), (cut,))
    biotin = g.input("biotin_datp", "biotin-14-dATP", MaterialKind.MOLECULAR)
    fill_polymerase = g.input("fill_polymerase", "DNA Polymerase I Large (Klenow) Fragment",
                             MaterialKind.CATALYST)
    marked = g.product("marked", "biotin-filled MboI ends", MaterialKind.MOLECULAR,
                       state=cut.state, chemical_state=WORKFLOW.digest().fill_in(biotin_base="A"))
    g.react("fill", "Fill MboI overhangs with biotin-dATP", StepKind.EXTENSION,
            ((cut, InputRole.TEMPLATE), (biotin, InputRole.STRUCTURAL),
             (fill_polymerase, InputRole.CATALYST)), (marked,),
            relations=(MolecularRelation.RESTRICTION_FILL,))
    ligase = g.input("contact_ligase", "T4 DNA ligase",
                     MaterialKind.CATALYST)
    contacts = g.product("contacts", "single-nucleus contact molecules",
                         MaterialKind.MOLECULAR, state=contact_product())
    g.react("contact_ligation", "Ligate biotin-filled MboI ends",
            StepKind.PROXIMITY_LIGATION,
            ((marked, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)),
            (contacts,), relations=(MolecularRelation.CONTACT_JUNCTION,))
    heat = g.input("crosslink_reversal", "65 C overnight crosslink reversal", MaterialKind.CONDITION)
    reversed_dna = g.product("reversed_contacts", "crosslink-reversed single-cell contacts", MaterialKind.MOLECULAR, state=contacts.state)
    g.react("reverse_crosslinks", "Reverse crosslinks before streptavidin capture", StepKind.DEPROTEINIZATION,
            ((contacts, InputRole.TEMPLATE), (heat, InputRole.CONDITION)), (reversed_dna,))
    beads = g.input("streptavidin", "streptavidin beads", MaterialKind.SELECTION)
    captured = g.product("captured", "captured contact DNA", MaterialKind.MOLECULAR,
                         state=contacts.state)
    unbound = g.product("unbound", "contact-poor unbound DNA", MaterialKind.WASTE)
    g.react("capture", "Capture biotin-marked contact junctions",
            StepKind.AFFINITY_CAPTURE,
            ((reversed_dna, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
            (captured, unbound))
    g.discard(unbound, "lacks a captured biotin junction")
    alui = g.input("alui", "AluI fragmentation enzyme", MaterialKind.CATALYST)
    fragments = g.product("fragments", "AluI-fragmented captured contacts",
                          MaterialKind.MOLECULAR, state=ALUI_INSERT)
    g.react("fragment", "Digest captured contacts with AluI", StepKind.DIGESTION,
            ((captured, InputRole.TEMPLATE), (alui, InputRole.CATALYST)),
            (fragments,))
    prep = g.input("end_prep", "Klenow Fragment (3-prime to 5-prime exo-minus)",
                   MaterialKind.CATALYST)
    tailed = g.product("tailed", "dA-tailed AluI fragments", MaterialKind.MOLECULAR,
                       state=ALUI_TAILED)
    g.react("tail", "dA-tail AluI fragments with Klenow exo-minus", StepKind.TAILING,
            ((fragments, InputRole.TEMPLATE), (prep, InputRole.CATALYST)), (tailed,))
    adapter = g.input("adapter", "source-printed CAA-tagged Illumina adapter",
                      MaterialKind.ADAPTER, state=adapter_scene())
    adapter_ligase = g.input("adapter_ligase", "T4 DNA ligase for adapter ligation",
                             MaterialKind.CATALYST)
    ligated = g.product("ligated", "adapter-ligated contact fragments",
                        MaterialKind.MOLECULAR)
    g.react("adapter_ligation", "Ligate the 3-bp-tagged adapter",
            StepKind.LIGATION,
            ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
             (adapter_ligase, InputRole.CATALYST)), (ligated,),
            geometry=adapter_geometry())
    forward = g.input("pcr_forward", "source-printed PCR forward primer", MaterialKind.OLIGO, state=primer_scene("PCR forward", PCR_FORWARD))
    reverse = g.input("pcr_reverse", "source-printed PCR reverse primer", MaterialKind.OLIGO, state=primer_scene("PCR reverse", PCR_REVERSE))
    pcr = g.input("pcr_polymerase", "Platinum Pfx DNA polymerase", MaterialKind.CATALYST)
    final = final_library()
    library = g.product("amplified", final.name, MaterialKind.MOLECULAR, state=final)
    g.react("library_pcr", "Amplify adapter-ligated contacts for 25 cycles with Platinum Pfx", StepKind.PCR,
            ((ligated, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
             (reverse, InputRole.STRUCTURAL), (pcr, InputRole.CATALYST)),
            (library,), geometry=PCRGeometry.unavailable(
                final, reason="The captured contact population is not one exact pre-PCR construct."))
    ampure = g.input("ampure", "AMPure XP beads, 0.7x post-PCR recovery", MaterialKind.SELECTION)
    cleaned = g.product("cleaned", "AMPure-purified single-cell library", MaterialKind.MOLECULAR, state=final)
    waste = g.product("cleanup_waste", "material not recovered in post-PCR cleanup", MaterialKind.WASTE)
    g.react("post_pcr_cleanup", "Purify adapter-route library with 0.7x AMPure XP", StepKind.CLEANUP,
            ((library, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)), (cleaned, waste))
    g.discard(waste, "not retained in library eluate; exact composition unrepresented")
    selection = g.input("size_selection", "concentrated AMPure double selection, 300–700 bp", MaterialKind.SELECTION)
    selected = g.product("library", final.name, MaterialKind.LIBRARY, state=final)
    offsize = g.product("offsize", "products outside pooled library-size interval", MaterialKind.WASTE)
    g.react("recover_size_window", "Recover pooled 300–700 bp libraries by double AMPure selection", StepKind.SELECTION,
            ((cleaned, InputRole.TEMPLATE), (selection, InputRole.SELECTION)), (selected, offsize))
    g.discard(offsize, "outside recovered size interval")
    g.sequence(selected, SequencingPlan(library=final, primers=SEQ_PRIMERS,
                                       acquired_roles=REQUIRED_ROLES, read_lengths=READ_LENGTHS))
    g.declare_reviewed_scope(tuple(g.steps) + tuple(g.materials), CITATION)
    return g.require_publishable()
