"""Source-owned TruePrime Single Cell WGA workflow from Picher 2016."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1]), str(HERE.parents[1] / "lib")]

from batch_page import render
from chemdraw import Construct
from checks import Check
from protocol_model import molecule, role
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind


PAGE_FILE = "trueprime-single-cell.html"
TITLE = "TruePrime Single Cell WGA"
NOTES = "01_trueprime_single_cell.html"
SOURCE = '<a href="https://doi.org/10.1038/ncomms13296">Picher et al., <i>Nature Communications</i> (2016)</a>.'
SUMMARY = ("An intact HEK293 cell is lysed in Buffer L2, neutralized, then TthPrimPol makes short DNA primers "
           "that phi29 polymerase extends by strand displacement.")
CAVEAT = ("The model follows Picher et al., Methods ‘Single-cell WGA’: 10 min on ice with L2, then "
          "3 or 6 h amplification at 30 °C. The contemporaneous SYGNIS "
          "<a href=\"https://www.pacificscience.co.th/userfiles2/fckeditor/121912/file/Sygnis-TruePrime-Single-Cell-WGA-Handbook(1).pdf\"><i>TruePrime Single Cell WGA Kit v2.0 Handbook</i></a>, "
          "p. 11, instead gives 5 min at room temperature plus 1 min at 95 °C "
          "after L2, then 3 h with 6 h optional. The paper does not state which handbook revision "
          "was used, so the difference remains unresolved. The downstream NEBNext kit and adapter "
          "sequence are not resolved in this model. The <a href=\"../trueprime__10.1038+ncomms13296/trueprime.html\">purified-DNA WGA branch</a> "
          "starts instead with Buffer D alkaline denaturation.")
FINAL_LIBRARY = Construct([role("TruePrime-amplified single-cell genomic DNA")],
                          name="single-cell WGA product")
FINAL_CAPTION = "Amplified single-cell genomic DNA, before sequencing adapters."
SEQ_PRIMERS = ()
READ_LENGTHS = None
SEQUENCING_INTRO = ""
SEQUENCING_ENDING = "The paper reports downstream Covaris fragmentation, NEBNext library preparation and paired-end HiSeq 2500 sequencing."
SEQUENCING_UNAVAILABLE = "The NEBNext kit generation, adapter bases and sequencing-primer landing sites are unresolved."


def oligos():
    return ()


# Source-owned Picher 2016 amplification and downstream experiment.
def add_amplification(graph, template, *, authority: str, duration: str):
    """Prime exposed DNA de novo, then extend it by strand displacement."""
    primpol = graph.input("tth_primpol", "TthPrimPol", MaterialKind.CATALYST)
    primed = graph.product(
        "primed", "templates bearing de novo 7–9-nt DNA primers",
        MaterialKind.MOLECULAR,
        state=molecule("de novo primed template", (
            role("exposed genomic template"), role("7–9-nt DNA primers"))),
    )
    graph.react("prime", "Synthesize short DNA primers directly from dNTPs",
                StepKind.EXTENSION,
                ((template, InputRole.TEMPLATE), (primpol, InputRole.CATALYST)),
                (primed,), source=authority)
    phi29 = graph.input("phi29", "phi29 DNA polymerase", MaterialKind.CATALYST)
    amplified = graph.product(
        "amplified", "branched multiple-displacement WGA product",
        MaterialKind.MOLECULAR,
        state=molecule("TruePrime-amplified genomic DNA", (
            role("long genomic copies"), role("displaced single strands"))),
    )
    graph.react("amplify", f"Extend primers by strand displacement at 30 °C for {duration}",
                StepKind.EXTENSION,
                ((primed, InputRole.TEMPLATE), (phi29, InputRole.CATALYST)),
                (amplified,), source=authority)
    heat = graph.input("inactivation_condition", "65 °C for 10 min", MaterialKind.CONDITION)
    product = graph.product(
        "inactivated_wga", "TruePrime-amplified genomic DNA after phi29 inactivation",
        MaterialKind.MOLECULAR, state=amplified.state)
    graph.react("inactivate", "Inactivate phi29 at 65 °C for 10 min",
                StepKind.THERMAL_TREATMENT,
                ((amplified, InputRole.TEMPLATE), (heat, InputRole.CONDITION)),
                (product,), source=authority)
    return product


def add_reported_sequencing(graph, amplified, *, authority: str):
    """Show the reported downstream library without assigning unknown NEBNext oligos."""
    condition = graph.input("fragmentation_condition", "downstream fragmentation condition",
                            MaterialKind.CONDITION, inferred=True, source=authority)
    fragments = graph.product("fragments", "fragmented TruePrime-amplified DNA",
                              MaterialKind.MOLECULAR, state=amplified.state,
                              inferred=True, source=authority)
    graph.react("fragment", "Fragment the amplified DNA for library preparation",
                StepKind.FRAGMENTATION,
                ((amplified, InputRole.TEMPLATE), (condition, InputRole.CONDITION)),
                (fragments,), source=authority + " (reported Covaris experiment; library model unresolved)")
    library = graph.product("sequencing_library", "NEBNext library with unresolved kit identity",
                            MaterialKind.LIBRARY, inferred=True, source=authority)
    graph.react("library_prep", "Prepare the reported NEBNext sequencing library",
                StepKind.MAGIC, ((fragments, InputRole.TEMPLATE),), (library,),
                reason="The reported NEBNext preparation is not yet resolved into kit-specific reactions or adapter-end geometry in this model.",
                source=authority)
    graph.sequence(library, SequencingPlan(
        acquired_roles=("Read 1", "Read 2"),
        limitation="Paired-end HiSeq 2500 was reported; this model has not resolved the NEBNext kit version or primer sites."))
    return library


def protocol_graph() -> ProtocolGraph:
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    cell = graph.input("single_cell", "one visually inspected HEK293 cell in 2.5 µl PBS",
                       MaterialKind.BIOLOGICAL, primary=True,
                       state=molecule("intact HEK293 cell", (role("genomic DNA in intact cell"),)))
    l2 = graph.input("buffer_l2", "2.5 µl lysis Buffer L2", MaterialKind.CONDITION)
    lysate = graph.product("cell_lysate", "lysed single cell", MaterialKind.MOLECULAR,
                           state=molecule("single-cell lysate", (role("released genomic DNA"),)))
    graph.react("lyse", "Lyse the cell for 10 min on ice with Buffer L2",
                StepKind.EXTRACTION,
                ((cell, InputRole.TEMPLATE), (l2, InputRole.CONDITION)), (lysate,), source=SOURCE)
    neutralizer = graph.input("buffer_n", "2.5 µl neutralization Buffer N", MaterialKind.CONDITION)
    template = graph.product("neutralized_template", "neutralized single-cell DNA",
                             MaterialKind.MOLECULAR,
                             state=molecule("neutralized cell lysate", (role("exposed genomic template"),)))
    graph.react("neutralize", "Neutralize the lysate with Buffer N",
                StepKind.QUENCHING,
                ((lysate, InputRole.TEMPLATE), (neutralizer, InputRole.CONDITION)),
                (template,), source=SOURCE)
    amplified = add_amplification(graph, template, authority=SOURCE, duration="3 or 6 h")
    add_reported_sequencing(graph, amplified, authority=SOURCE)
    graph.declare_reviewed_scope([
        "single_cell", "buffer_l2", "lyse", "buffer_n", "neutralize",
        "tth_primpol", "prime", "phi29", "amplify", "inactivation_condition", "inactivate", "fragmentation_condition",
        "fragment", "library_prep",
        "sequencing_library"], SOURCE)
    return graph.require_publishable()


def render_page():
    return render(sys.modules[__name__])


def selftest():
    c = Check()
    graph = protocol_graph()
    c("single-cell input and L2 lysis, without purified-DNA Buffer D",
      ("single_cell" in graph.materials, "buffer_l2" in graph.materials,
       "denaturation_condition" in graph.materials), (True, True, False))
    c.report()
