"""RNAscope double-Z branched-DNA signal-amplification architecture."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from branched_dna import BranchModel

TITLE = "RNAscope"
NOTES = "01_rnascope.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1016/j.jmoldx.2011.08.002">'
          'Wang et al., <i>Journal of Molecular Diagnostics</i> (2012)</a>; architecture '
          'cross-checked against current Bio-Techne RNAscope documentation.')
SUMMARY = "Adjacent double-Z probes recognize RNA in fixed cells or tissue and nucleate a branched-DNA signal tree for fluorescent or chromogenic microscopy. The RNA itself is not amplified."
MODEL = BranchModel(20, 20, 20, 20)
PREPARATION = "fixed tissue or cells → target retrieval / permeabilization → accessible RNA"
PREPARATION_CAPTION = "Fixation preserves morphology; pretreatment exposes RNA to the probe set. The exact pretreatment depends on specimen and RNAscope kit."
TARGET_CAPTION = "Each probe has an 18–25-nt target-complementary region, spacer and 14-nt Z tail. Only two probes bound adjacently create the complete 28-nt preamplifier site. A standard set contains 20 pairs across roughly 1 kb."
DETECTION = "label probes → fluorescent punctum or enzyme-generated chromogenic dot"
DETECTION_CAPTION = "Microscopy localizes individual RNA molecules in their tissue context; label chemistry depends on the fluorescent or bright-field assay."
READOUT = "RNAscope produces spatial microscopy signal, not a sequencing library. There are therefore no flow-cell adapters, sample indices, sequencing primers or nucleotide reads to draw."
