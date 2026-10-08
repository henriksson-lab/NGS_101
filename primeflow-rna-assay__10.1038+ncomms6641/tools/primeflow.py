"""PrimeFlow RNA adjacent-probe branched-DNA signal-amplification architecture."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from branched_dna import BranchModel

TITLE = "PrimeFlow RNA Assay"
NOTES = "01_primeflow-rna.html"
SOURCE = ('Commercial protocol: <a href="https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0019788_PrimeFlowRNAAssay_UG.pdf">'
          'Thermo Fisher MAN0019788 Rev B.0</a>; defining flow-RNA application: '
          '<a href="https://doi.org/10.1038/ncomms6641">Porichis et al. (2014)</a>.')
SUMMARY = "Adjacent target-probe pairs detect RNA inside fixed, permeabilized suspension cells; a branched-DNA tree amplifies fluorescence for single-cell flow cytometry, optionally alongside antibody staining."
MODEL = BranchModel(20, 40, 20, 20)
PREPARATION = "suspension cells ± surface antibodies → fix → permeabilize → intracellular antibodies (optional)"
PREPARATION_CAPTION = "The manual preserves RNA and cell identity through fixation and permeabilization so RNA fluorescence can be combined with immunophenotyping."
TARGET_CAPTION = "A gene-specific set contains 20–40 adjacent probe pairs. Stable preamplifier binding requires both halves of a pair; exact target-probe and tail sequences are proprietary."
DETECTION = "fluorophore-conjugated label probes → fluorescence intensity per intact cell"
DETECTION_CAPTION = "Types 1, 4, 6 and 10 use distinct amplification structures and fluorophores, enabling up to four RNA targets before standard flow-cytometric compensation and analysis."
READOUT = "PrimeFlow produces per-cell fluorescence measured by a flow cytometer, not a sequencing library. There are therefore no sequencing adapters, indices, sequencing-primer sites or base-call reads."
