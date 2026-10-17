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
MODEL = BranchModel(20, 40, labels_per_tree=400)
TREE_MULTIPLICITY_SOURCE = ("The 2022 user guide (p.5) states only that a complete tree has "
                            "400 label-probe binding sites; it gives no amplifier-per-"
                            "preamplifier or label-per-amplifier split, so none is modelled.")
PREPARATION = "suspension cells ± surface antibodies → Fixation Buffer 1 → permeabilize → intracellular antibodies (optional) → Fixation Buffer 2"
PREPARATION_CAPTION = ("The 2022 guide (Day 1 steps 4–19) fixes surface-stained cells in Fixation Buffer 1 "
                       "(30 min, 2–8 °C), permeabilizes them with RNase inhibitors, allows "
                       "intracellular antibody staining, then fixes again in Fixation Buffer 2 "
                       "(60 min, room temperature) before washes. Porichis et al. (2014, Methods) "
                       "likewise fixed, permeabilized with cold methanol and fixed a second time, "
                       "so RNA fluorescence can be combined with immunophenotyping.")
TARGET_CAPTION = ("A gene-specific set contains 20–40 adjacent probe pairs (2022 guide). Preamplifier "
                  "molecules hybridize only when both halves of a pair have bound. Target probes "
                  "are diluted 1:20 and hybridized for 2 h at 40 °C, then washed (guide and "
                  "Porichis et al. Methods). Porichis et al. disclose the target regions used in "
                  "their sets, but not the universal tails or amplification-tree oligos.")
SEQUENCE_BOUNDARY = ("Porichis et al. Supplementary Table 1 publishes the target-RNA "
                     "intervals and sequences used to design its IL-2, IL-21, IFN-γ, IDO, "
                     "HMBS and miRNA probe sets. The article's Methods says probe "
                     "sequences are detailed there, but the table columns contain target "
                     "regions and coordinates, not complete manufactured probes. These "
                     "determine the complementary "
                     "target-binding arms. The manufactured probes' universal tails, "
                     "including their lengths and docking bases, and the preamplifier, "
                     "amplifier and label-probe sequences are not given in the reviewed "
                     "table or 2022 vendor guide. Linked ViewRNA sources remain under review. "
                     "Porichis et al. (2014, Methods and Results) report three available fluorochrome "
                     "probe sets (Alexa 488, 647 and 750); the model follows the later "
                     "PrimeFlow user guide (2022, probe-type table), which lists four "
                     "types (1, 4, 6 and 10). This is a product-version difference, not "
                     "evidence that the 2014 experiment used Type 10. "
                     "The same 2022 guide is internally inconsistent about structures: its "
                     "product description (p.5) says “three different amplification structures” "
                     "allow up to four RNA targets, whereas its assay guidelines (p.9) and "
                     "workflow figure list four probe types (1, 4, 6 and 10), each with its own "
                     "dye. The model shows four types and does not claim how many distinct "
                     "amplification structures underlie them. The guide states 400 label-probe "
                     "sites per complete tree but not how they divide between amplifiers and "
                     "label probes; step 3 therefore models only the total and draws "
                     "representative arms without per-level counts.")
DETECTION = "fluorophore-conjugated label probes → fluorescence intensity per intact cell"
DETECTION_CAPTION = ("PreAmp Mix and Amp Mix are each hybridized for 1.5 h at 40 °C with washes "
                     "between, then label probes (1:100) for 1 h at 40 °C (guide Day 2; Porichis "
                     "et al. Methods). Type 1 (Alexa Fluor 647), Type 4 (488), Type 6 (750) and "
                     "Type 10 (568) label probes allow up to four RNA targets before standard "
                     "flow-cytometric compensation and analysis.")
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
