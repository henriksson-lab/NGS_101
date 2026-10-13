"""PrimeFlow RNA adjacent-probe branched-DNA signal-amplification architecture."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from branched_dna import BranchModel, SequenceExample, adjacent_target_example_scene
from chemdraw import Segment

TITLE = "PrimeFlow RNA Assay"
NOTES = "01_primeflow-rna.html"
SOURCE = ('Commercial protocol: <a href="https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0019788_PrimeFlowRNAAssay_UG.pdf">'
          'Thermo Fisher MAN0019788 Rev B.0</a>; defining flow-RNA application: '
          '<a href="https://doi.org/10.1038/ncomms6641">Porichis et al. (2014)</a>.')
SUMMARY = "Adjacent target-probe pairs detect RNA inside fixed, permeabilized suspension cells; a branched-DNA tree amplifies fluorescence for single-cell flow cytometry, optionally alongside antibody staining."
MODEL = BranchModel(20, 40, 20, 20)
PREPARATION = "suspension cells ± surface antibodies → fix → permeabilize → intracellular antibodies (optional)"
PREPARATION_CAPTION = "The manual preserves RNA and cell identity through fixation and permeabilization so RNA fluorescence can be combined with immunophenotyping."
TARGET_CAPTION = "A gene-specific set contains 20–40 adjacent probe pairs. Stable preamplifier binding requires both halves of a pair. Porichis et al. disclose the target regions used in their sets, but not the universal tails or amplification-tree oligos."
SEQUENCE_BOUNDARY = ("Porichis et al. Supplementary Table 1 publishes the target-RNA "
                     "intervals and sequences used to design its IL-2, IL-21, IFN-γ, IDO, "
                     "HMBS and miRNA probe sets. These determine the complementary "
                     "target-binding arms. The manufactured probes' universal tails, "
                     "including their lengths and docking bases, and the preamplifiers, "
                     "amplifiers and label probes remain undisclosed.")
DETECTION = "fluorophore-conjugated label probes → fluorescence intensity per intact cell"
DETECTION_CAPTION = "Types 1, 4, 6 and 10 use distinct amplification structures and fluorophores, enabling up to four RNA targets before standard flow-cytometric compensation and analysis."
READOUT = "PrimeFlow produces per-cell fluorescence measured by a flow cytometer, not a sequencing library. There are therefore no sequencing adapters, indices, sequencing-primer sites or base-call reads."

# Porichis et al. Supplementary Table 1, human IL-2 NM_021803 bases 2–24 and
# 25–52. The table reports target regions, not complete manufactured oligos.
IL2_TARGET_2_24 = "CCTTGGTCTCGTTTTCACTTCAG"
IL2_TARGET_25_52 = "TCTCATAAGTACCAACAGTAGAGCTAGA"


def sequence_example() -> SequenceExample:
    return SequenceExample(
        "Paper-disclosed adjacent target regions",
        "The first adjacent IL-2 pair is shown as a concrete example. The two target "
        "sequences are transcribed from Supplementary Table 1; the paired probe arms "
        "are their computed reverse complements. They are not complete orderable "
        "PrimeFlow probes because neither the bases nor lengths of the universal tails "
        "were published.",
        (("IL-2 target region, NM_021803 bases 2–24",
          (Segment("published target region", IL2_TARGET_2_24, "me"),)),
         ("IL-2 target region, NM_021803 bases 25–52",
          (Segment("published target region", IL2_TARGET_25_52, "me"),))),
        adjacent_target_example_scene(IL2_TARGET_2_24, IL2_TARGET_25_52),
        "The target regions are contiguous. Complementarity and probe orientation are "
        "derived by the shared model; the non-nucleotide tail placeholders mark the "
        "unpublished commercial part without inventing a length. The downstream tree "
        "is therefore a topological view, not a claimed oligo sequence.")
