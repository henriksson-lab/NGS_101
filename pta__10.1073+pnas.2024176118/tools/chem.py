"""Primary template-directed amplification and its defining library route."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import pta
from batch_ngs import seg, truseq_library
from chemdraw import Construct, MolecularState, Scene, Workflow
from endprep import dA_tailed_scene, repair_and_dA_tail


TITLE = "PTA — primary template-directed whole-genome amplification"
NOTES = "01_pta.html"
SOURCE = ('<a href="https://doi.org/10.1073/pnas.2024176118">Gonzalez-Pena et al., '
          '<i>PNAS</i> (2021)</a>, Methods and Supplementary Figure S1.')
SUMMARY = ("Random-primed phi29 synthesis is interrupted by exonuclease-resistant "
           "alpha-thio-ddNTPs. Shorter products undergo less daughter-on-daughter "
           "amplification, biasing yield toward repeated copying of the primary genome.")
CAVEAT = ("The defining paper does not publish the random-primer or IDT UDI adapter "
          "sequences. Their roles are preserved; canonical TruSeq adapter bases are "
          "shown as inferred only to locate standard sequencing primers.")

FINAL_LIBRARY, SEQ_PRIMERS = truseq_library(
    [seg("PTA-amplified genomic insert", "X" * 48, placeholder=True)],
    "PTA KAPA HyperPlus library", dual_index=True, inferred_adapters=True)
READ_LENGTHS = {"Index 1 (i7)": 8, "Index 2 (i5)": 8}
FINAL_CAPTION = ("INFERRED adapter bases: the published KAPA HyperPlus plus IDT UDI route "
                 "fixes a standard dual-index Illumina topology but does not print its "
                 "oligos. The genomic insert was already shortened by PTA, so HyperPlus "
                 "fragmentation is omitted.")
SEQUENCING_INTRO = ("The study used NovaSeq 6000 but did not report genomic read lengths "
                    "or exact UDI oligos. Canonical primer sites follow the stated KAPA "
                    "HyperPlus/IDT dual-index library geometry.")


def oligos():
    return ()


def _initial_rows():
    sc = Scene.duplex([seg("single-cell genomic DNA", "X" * 48, placeholder=True)],
                      label="single-cell genome")
    return sc.rows()


def workflow():
    wf = Workflow(MolecularState("Denatured single-cell genome", tuple(_initial_rows())))
    wf.react("Anneal the undisclosed exonuclease-resistant random primer",
             pta.priming_rows(), name="Random-primed primary template",
             note="The paper gives 500 µM primer and 15 minutes at 65 °C, but does not disclose its bases, length or modification pattern.")
    wf.react("Extend with phi29 and strand displacement at 30 °C for 8 hours",
             pta.product_rows(), name="Terminated primary-template copies",
             note="Equal alpha-thio-ddATP, -ddCTP, -ddGTP and -ddTTP stop synthesis irreversibly; shorter displaced products are poor secondary templates, producing quasi-linear amplification.")
    body = seg("extendable PTA product", "X" * 48, placeholder=True)
    prepared = repair_and_dA_tail(Construct([body], name="PTA product"))
    wf.react("End-repair and dA-tail the purified PTA products without fragmentation",
             dA_tailed_scene(prepared, label="ligation-ready PTA product").rows(),
             name="End-repaired PTA product",
             note="The defining KAPA HyperPlus workflow explicitly omits fragmentation because PTA already creates shorter products.")
    sc = Scene.duplex(list(FINAL_LIBRARY), label="INFERRED adapter-completed library")
    sc.labels("top")
    wf.react("INFERRED — ligate IDT UDI adapters and amplify 15 PCR cycles",
             sc.rows(), name="Dual-index Illumina library",
             note="Exact adapter and primer bases are unavailable. Gel extraction retains 400–700-bp products for NovaSeq 6000 sequencing.")
    return wf


def final_rows():
    sc = Scene.duplex(list(FINAL_LIBRARY), label="PTA library")
    sc.labels("top")
    return sc.rows()
