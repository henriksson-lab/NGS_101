"""Source-owned TruePrime Single Cell WGA workflow from Picher 2016."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1]), str(HERE.parents[1] / "lib")]

from batch_page import render
from chemdraw import Construct
from checks import Check
from protocol_model import molecule, role
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_families.primpol_mda import add_amplification, add_reported_sequencing


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
          "sequence are unidentified. The <a href=\"../trueprime__10.1038+ncomms13296/trueprime.html\">purified-DNA WGA branch</a> "
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
