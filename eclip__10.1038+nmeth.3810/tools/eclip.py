"""Enhanced CLIP molecular workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg, truseq_library
from chemdraw import Row
from rna_special import adapter_ligation_scene

TITLE = "eCLIP — enhanced crosslinking and immunoprecipitation"
NOTES = "01_eclip.html"
SOURCE = 'Defining paper and protocol: <a href="https://doi.org/10.1038/nmeth.3810">Van Nostrand et al. (2016)</a>; <a href="https://www.encodeproject.org/documents/fa2a3246-6039-46ba-b960-17ef06e7876a/">ENCODE eCLIP SOP</a>.'
SUMMARY = "UV-crosslinked RBP–RNA complexes are immunoprecipitated and RNase-trimmed. A 3′ RNA adapter precedes reverse transcription; a second adapter is ligated to cDNA so both read-through and crosslink-truncated cDNAs can amplify."
CAVEAT = "Barcode identities vary by experiment. Adapter roles are source-supported; the final flow-cell arms are inferred canonical TruSeq geometry."
FINAL_LIBRARY, SEQ_PRIMERS = truseq_library([seg("RBP-bound RNA cDNA", "X"*34, placeholder=True)], "eCLIP library", inferred_adapters=True)
FINAL_CAPTION = "Final eCLIP amplicon. Its insert can end at the reverse-transcription stop immediately downstream of the protein–RNA crosslink."
SEQUENCING_INTRO = "Read start and cDNA termination together locate RBP-associated RNA; the sample-matched input controls background."

def sections():
    return [
        ("Crosslink, trim and immunoprecipitate", [Row(chunks=[("RNA ————— ×RBP× ————— RNA", None, False)]), Row(chunks=[("UV crosslink   + partial RNase   + antibody beads", "umi", False)]), Row(chunks=[("                 retained RBP–RNA fragment", None, False)])], "UV fixes direct contacts; limited RNase and RBP immunoprecipitation define the captured RNA fragment."),
        ("Ligate the 3′ RNA adapter on beads", adapter_ligation_scene(fragment_name="RBP-bound RNA fragment").rows(), "The indexed RNA adapter is ligated before gel purification and reverse transcription."),
        ("Reverse transcription can stop at the crosslink", [Row(chunks=[("RNA  5′ ————— × peptide ————— adapter — 3′", None, False)]), Row(chunks=[("cDNA 3′ ————— ← RT stop", "r1", False)]), Row(chunks=[("or   3′ ——————————————————— ← read-through cDNA", "r2", False)])], "Both truncated and read-through cDNAs are retained."),
        ("Ligate a second adapter to cDNA", [Row(chunks=[("ssDNA adapter ** cDNA ending at crosslink", "me", False)]), Row(chunks=[("ssDNA adapter ** full-length cDNA — RNA-adapter complement", "me", False)])], "The cDNA 3′ ligation bypasses the circularization step used by iCLIP and makes both cDNA classes amplifiable."),
    ]
