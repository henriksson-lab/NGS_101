"""TrAC-looping molecular model (Lai et al. 2018)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import nextera as nx
from chemdraw import MolecularState, Scene, Segment, Workflow, annotation_rows, oligo
from chromatin_epigenetics import tagmented_library
from contact_capture import BivalentTransposaseLinker

BIO67F = "CTGTCTCTTATACACATCTCCGAGCCCACGAGACTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
BP67R = "CTGTCTCTTATACACATCTGACGCTGCCGACGAGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
SPACER = "CCGAGCCCACGAGACTCGTCGGCAGCGTC"
LINKER = BivalentTransposaseLinker(nx.ME_RC, SPACER, "67-bp bivalent ME linker")


def bridge_scene():
    sc = Scene.duplex(list(LINKER.bridged_product()), label="two accessible loci")
    sc.junction("top", "contact locus A", "ME A", "Tn5 insertion")
    sc.junction("top", "ME B", "contact locus B", "Tn5 insertion")
    sc.note("top", "tetrameric Tn5 inserts the same short linker into two chromatin regions in trans")
    sc.labels("top")
    return sc


FINAL_LIBRARY, SEQ_PRIMERS = tagmented_library("RCA-derived TrAC contact insert")


def final_rows():
    sc = Scene.duplex(list(FINAL_LIBRARY), label="paired-end library")
    sc.labels("top")
    return [*sc.rows(), *annotation_rows(FINAL_LIBRARY)]


def workflow():
    wf = Workflow(MolecularState("Bivalent transposition in intact chromatin", tuple(bridge_scene().rows())))
    wf.react("Digest genomic arms with a four-base cutter and capture biotinylated products",
             bridge_scene().rows(), note="The bridge, not proximity ligation, records the contact.")
    wf.react("Intramolecular circularization with T7 DNA ligase", bridge_scene().rows(),
             note="Dilute ligation closes each enriched contact molecule into a circle.")
    wf.react("Rolling-circle amplification", bridge_scene().rows(),
             note="TempliPhi amplification copies the circle before indexed library PCR.")
    wf.react("N501/N701 indexed PCR", final_rows(),
             note="The printed primers complete canonical Nextera-compatible P5 and P7 ends.")
    return wf


TITLE = "TrAC-looping — contacts captured by bivalent Tn5 insertion"
NOTES = "01_trac-looping.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/s41592-018-0107-y">'
          'Lai et al., <i>Nature Methods</i> (2018)</a>.')
SUMMARY = ("A tetrameric Tn5 complex inserts a biotinylated two-ME linker across two nearby "
           "accessible chromatin regions. Restriction digestion, affinity capture, circle formation "
           "and rolling-circle amplification recover the linked pair without Hi-C proximity ligation.")
FINAL_CAPTION = "PCR-completed Nextera-compatible molecule copied from a circularized TrAC contact product."
SEQUENCING_INTRO = "Paired reads enter the two genomic ends represented in the amplified circle; dual indexes come from N5xx/N7xx PCR primers."


def oligos():
    yield oligo("Bio67F", [Segment("bivalent top strand", BIO67F, "me")], mods="5′ phosphate; internal biotin-dT")
    yield oligo("67bpR", [Segment("bivalent bottom strand", BP67R, "me")], mods="5′ phosphate")
