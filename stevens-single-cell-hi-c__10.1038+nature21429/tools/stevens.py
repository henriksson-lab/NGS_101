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
from chemdraw import Construct, MolecularState, Scene, Workflow
from endprep import (ThreePrimeOverhangAdapter, cohesive_ligation_geometry,
                     dA_tailed_scene, ligation_product_scene, repair_and_dA_tail)
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, MolecularRelation, ProtocolGraph,
                            SequencingPlan, StepKind)
from protocol_model import molecule, role
WORKFLOW = ContactWorkflow("Stevens single-cell Hi-C", MBOI,
    "before digestion and ligation", "biotin-fill", "PCR", capture=True)
TITLE = "Stevens single-cell Hi-C — single nuclei before AluI ligation"
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
ALUI_ADAPTER = ThreePrimeOverhangAdapter(
    "3-bp-tagged Illumina adapter ligating stem", "GCTCTTCCGATCT", "GATCGGAAGAGC")


def adapter_geometry():
    return cohesive_ligation_geometry(
        ALUI_TAILED, ALUI_ADAPTER,
        statement="The source-printed customized Illumina adapter retains the standard 3′-T ligating end, complementary to the fragment 3′-dA; CAA is the inline identification tag.")

def adapter_workflow():
    w=Workflow(MolecularState("AluI-digested captured contacts",tuple(Scene.duplex(list(ALUI_INSERT)).rows())))
    w.react("End-repair and add 3′ dA",dA_tailed_scene(ALUI_TAILED).rows(),name="dA-tailed contact fragments",note="The source-printed ligation route dA-tails the AluI fragments.")
    geometry=adapter_geometry()
    product=ligation_product_scene(final_library(),"3-bp identification tag","opposite 3-bp identification tag",label="Stevens scHi-C library")
    w.ligate_adapter("Ligate the 3-bp-tagged Illumina adapter",geometry,product.rows(),name="adapter-ligated contact library",note="Nextera XT was a published alternative.")
    return w
def contact_product(): return WORKFLOW.junction()
SEQ_PRIMERS = (
    sp.custom("Read 1", "historical Illumina Read 1", il.TRUSEQ_READ1,
              "Stevens et al. 2017 Supplementary Methods"),
    sp.custom("Read 2", "historical paired-end Read 2", HISTORICAL_PE_READ2,
              "Stevens et al. 2017 Supplementary Methods"),
)
REQUIRED_ROLES = ("Read 1", "Read 2")
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
    excluded = g.product("excluded", "other nuclei", MaterialKind.WASTE)
    g.react("isolate", "Distribute individual nuclei before library construction",
            StepKind.SAMPLING,
            ((sample, InputRole.TEMPLATE), (isolation, InputRole.CONDITION)),
            (nucleus, excluded), changes_molecule=False)
    g.discard(excluded, "outside the retained single-nucleus reaction")
    mboi = g.input("mboi", "MboI", MaterialKind.CATALYST)
    cut = g.product("cut", "MboI-cut chromatin", MaterialKind.MOLECULAR,
                    state=nucleus.state, chemical_state=WORKFLOW.digest())
    g.react("digest", "Digest single-nucleus chromatin with MboI",
            StepKind.DIGESTION,
            ((nucleus, InputRole.TEMPLATE), (mboi, InputRole.CATALYST)), (cut,))
    biotin = g.input("biotin_datp", "biotin-14-dATP", MaterialKind.MOLECULAR)
    fill_polymerase = g.input("fill_polymerase", "restriction-end fill polymerase",
                             MaterialKind.CATALYST)
    marked = g.product("marked", "biotin-filled MboI ends", MaterialKind.MOLECULAR,
                       state=cut.state, chemical_state=WORKFLOW.digest().fill_in(biotin_base="A"))
    g.react("fill", "Fill MboI overhangs with biotin-dATP", StepKind.EXTENSION,
            ((cut, InputRole.TEMPLATE), (biotin, InputRole.STRUCTURAL),
             (fill_polymerase, InputRole.CATALYST)), (marked,),
            relations=(MolecularRelation.RESTRICTION_FILL,))
    ligase = g.input("contact_ligase", "proximity-ligation DNA ligase",
                     MaterialKind.CATALYST)
    contacts = g.product("contacts", "single-nucleus contact molecules",
                         MaterialKind.MOLECULAR, state=contact_product())
    g.react("contact_ligation", "Ligate biotin-filled MboI ends",
            StepKind.PROXIMITY_LIGATION,
            ((marked, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)),
            (contacts,), relations=(MolecularRelation.CONTACT_JUNCTION,))
    beads = g.input("streptavidin", "streptavidin beads", MaterialKind.SELECTION)
    captured = g.product("captured", "captured contact DNA", MaterialKind.MOLECULAR,
                         state=contacts.state)
    unbound = g.product("unbound", "contact-poor unbound DNA", MaterialKind.WASTE)
    g.react("capture", "Capture biotin-marked contact junctions",
            StepKind.AFFINITY_CAPTURE,
            ((contacts, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
            (captured, unbound))
    g.discard(unbound, "lacks a captured biotin junction")
    alui = g.input("alui", "AluI fragmentation enzyme", MaterialKind.CATALYST)
    fragments = g.product("fragments", "AluI-fragmented captured contacts",
                          MaterialKind.MOLECULAR, state=ALUI_INSERT)
    g.react("fragment", "Digest captured contacts with AluI", StepKind.DIGESTION,
            ((captured, InputRole.TEMPLATE), (alui, InputRole.CATALYST)),
            (fragments,))
    prep = g.input("end_prep", "end-repair and dA-tailing enzymes",
                   MaterialKind.CATALYST)
    tailed = g.product("tailed", "dA-tailed AluI fragments", MaterialKind.MOLECULAR,
                       state=ALUI_TAILED)
    g.react("tail", "End-repair and dA-tail the AluI fragments", StepKind.TAILING,
            ((fragments, InputRole.TEMPLATE), (prep, InputRole.CATALYST)), (tailed,))
    adapter = g.input("adapter", "3-bp-tagged historical Illumina adapter",
                      MaterialKind.ADAPTER)
    adapter_ligase = g.input("adapter_ligase", "library-adapter ligase",
                             MaterialKind.CATALYST)
    ligated = g.product("ligated", "adapter-ligated contact fragments",
                        MaterialKind.MOLECULAR, state=final_library())
    g.react("adapter_ligation", "Ligate the 3-bp-tagged adapter",
            StepKind.LIGATION,
            ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
             (adapter_ligase, InputRole.CATALYST)), (ligated,),
            geometry=adapter_geometry())
    forward = g.input("pcr_forward", "library PCR forward primer", MaterialKind.OLIGO)
    reverse = g.input("pcr_reverse", "library PCR reverse primer", MaterialKind.OLIGO)
    pcr = g.input("pcr_polymerase", "library PCR polymerase", MaterialKind.CATALYST)
    final = final_library()
    library = g.product("library", final.name, MaterialKind.LIBRARY, state=final)
    g.react("library_pcr", "Amplify the adapter-ligated contact library", StepKind.PCR,
            ((ligated, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
             (reverse, InputRole.STRUCTURAL), (pcr, InputRole.CATALYST)),
            (library,), geometry=PCRGeometry.unavailable(
                final, reason="The captured contact population is not one exact pre-PCR construct."))
    g.sequence(library, SequencingPlan(library=final, primers=SEQ_PRIMERS,
                                       acquired_roles=REQUIRED_ROLES))
    g.declare_reviewed_scope(tuple(g.steps) + tuple(g.materials), CITATION)
    return g.require_publishable()
