"""Precision nuclear run-on sequencing workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg, truseq_library
from chemdraw import Row
from rna_special import adapter_ligation_scene

TITLE = "PRO-seq — precision nuclear run-on sequencing"
NOTES = "01_pro-seq.html"
SOURCE = 'Defining paper: <a href="https://doi.org/10.1101/gr.150027.112">Kwak et al. (2013)</a>; detailed protocol: <a href="https://doi.org/10.1038/nprot.2016.086">Mahat et al. (2016)</a>.'
SUMMARY = "Engaged polymerases extend nascent RNA by one or a few biotin-NTPs. Biotin purification and sequential RNA-adapter ligations preserve the nucleotide-resolution 3′ end for sequencing."
CAVEAT = "The publication defines adapters and barcodes in its oligo tables; the final panel represents their run roles with inferred canonical TruSeq arms rather than silently substituting bases."
FINAL_LIBRARY, SEQ_PRIMERS = truseq_library([seg("nascent RNA-derived cDNA", "X"*34, placeholder=True)], "PRO-seq library", inferred_adapters=True)
FINAL_CAPTION = "Directional library whose insert boundary reports the biotin-run-on RNA 3′ end."
SEQUENCING_INTRO = "The insert-facing read reports the position and orientation of engaged RNA polymerase."

def sections():
    return [
        ("Biotin nuclear run-on", [Row(chunks=[("DNA template ————— RNAP — nascent RNA—3′", None, False)]), Row(chunks=[("                            + biotin-NTP → RNA—●—3′ (stalled)", "umi", False)])], "Permeabilized cells undergo a short run-on with biotinylated NTPs; the bulky incorporated nucleotide both marks and limits extension."),
        ("Ligate the 3′ RNA adapter", adapter_ligation_scene(fragment_name="biotinylated nascent RNA").rows(), "The adapter is joined at the informative nascent-RNA 3′ end."),
        ("Capture, decap and ligate the 5′ adapter", [Row(chunks=[("biotin-RNA—3′ adapter  → streptavidin capture → 5′ end repair", None, False)]), Row(chunks=[("5′ adapter ** nascent RNA ** 3′ adapter", "r1", False)])], "Repeated streptavidin enrichment retains biotin-run-on RNA; 5′ processing permits the second RNA-adapter ligation."),
        ("Reverse-transcribe and PCR-index", [Row(chunks=[("5′ PCR arm — cDNA of nascent RNA — index — 3′ PCR arm", None, True)])], "INFERRED — final run-arm bases are represented by the computed final-library panel below."),
    ]
