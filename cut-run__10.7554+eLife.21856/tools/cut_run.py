"""Original Skene & Henikoff CUT&RUN molecular model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import duplex_fragment, seg, targeting_rows
from chemdraw import Construct, Row

TITLE = "CUT&RUN — antibody-tethered MNase chromatin profiling"
NOTES = "01_cut-run.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.7554/eLife.21856">Skene & Henikoff, <i>eLife</i> (2017)</a>.'
SUMMARY = "An antibody recruits Protein A–MNase to native chromatin. Calcium activates local cleavage; fragments cut on both sides of the bound particle diffuse into the supernatant and receive conventional Illumina adapters afterward."
CAVEAT = "CUT&RUN does not install adapters during targeted cleavage. The defining paper used the KAPA DNA polymerase library-preparation chemistry, but KAPA supplied no adapters in that kit and the paper does not identify the separately chosen adapter/index oligos. Their sequences, index count and sequencing-primer set are therefore left unavailable instead of being replaced by generic TruSeq."

FINAL_LIBRARY = Construct([
    seg("unreported left adapter", "[adapter]", placeholder=True, inferred=True),
    seg("released target-proximal DNA", "X" * 38, placeholder=True),
    seg("unreported right adapter / index", "[adapter/index]", placeholder=True,
        inferred=True),
], name="CUT&RUN library (adapter identity unavailable)")
SEQ_PRIMERS = ()
FINAL_CAPTION = "INFERRED — the released MNase fragment is end-repaired, dA-tailed, adapter-ligated and PCR-amplified. Dotted role tokens mark adapter information that the defining paper does not report."
SEQUENCING_ENDING = "The paper reports 25-cycle paired-end HiSeq 2500 sequencing, but does not name the adapter oligos or sequencing primers."
SEQUENCING_UNAVAILABLE = "adapter and sequencing-primer identities are absent from the defining paper"

def sections():
    return [
        ("Immobilize nuclei and tether pA–MNase",
         targeting_rows("MNase", "inactive nuclease positioned at the epitope"),
         "Concanavalin-A beads retain permeabilized nuclei during antibody and pA–MNase binding and washing."),
        ("Activate cleavage and release the selected fragment",
         [Row(chunks=[("chromatin ---- ^ [antibody-bound particle] ^ ---- chromatin", None, False)]),
          Row(chunks=[("              Ca²⁺-activated cuts; two-cut fragment enters supernatant", "me", False)])],
         "Calcium activates MNase at 0 °C. Chelators stop digestion; fragments cleaved on both sides of a particle are recovered from the soluble fraction."),
        ("Recover DNA for conventional library preparation",
         duplex_fragment("released MNase fragment").rows(),
         "Proteinase treatment and extraction yield an unadapted duplex. Adapter ligation is a separate downstream operation, unlike CUT&Tag."),
    ]
