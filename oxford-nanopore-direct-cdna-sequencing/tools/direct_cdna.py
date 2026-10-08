"""Oxford Nanopore SQK-DCS109 direct-cDNA workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Row
from rna_special import nanopore_entry_rows, template_switch_scene

TITLE = "Oxford Nanopore Direct cDNA Sequencing — SQK-DCS109"
NOTES = "01_direct-cdna.html"
SOURCE = 'Commercial protocol: Oxford Nanopore <a href="https://nanoporetech.com/document/chemistry-technical-document">Chemistry Technical Document — Direct cDNA Sequencing Kit SQK-DCS109</a>.'
SUMMARY = "Poly(A) RNA is reverse-transcribed with strand switching, RNA is removed and second-strand DNA is synthesized. The PCR-free duplex cDNA is end-repaired, dA-tailed and ligated to a motor-loaded dT-overhang sequencing adapter."
CAVEAT = "SQK-DCS109 is a historical commercial protocol and its adapter bases are proprietary. The page distinguishes it from cDNA-PCR: no PCR or rapid-attachment tag is introduced."
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
SEQUENCING_ENDING = "There is no synthesis sequencing primer. A motor in the ligated sequencing adapter presents one strand of the PCR-free duplex cDNA to the pore. The adapter-ligated end enters first."

def sections():
    return [
        ("Make first-strand cDNA by strand switching", template_switch_scene(inferred=True).rows(), "INFERRED — primer bases are proprietary; their published functions are poly(A)-primed RT and capture of a completed cDNA end by strand switching."),
        ("Remove RNA and synthesize the second strand", [Row(chunks=[("RNA:cDNA hybrid → RNase cocktail → single-stranded cDNA", None, False)]), Row(chunks=[("second-strand synthesis → full-length duplex cDNA", "r1", False)])], "Unlike direct RNA sequencing, the RNA is removed and DNA is the pore substrate."),
        ("End-repair and dA-tail the duplex", [Row(chunks=[("5′ — full-length duplex cDNA — A 3′", None, False)]), Row(chunks=[("3′ — full-length duplex cDNA — T 5′", None, False)])], "Repair creates ligatable ends and dA tailing makes them compatible with the dT-overhang sequencing adapter."),
        ("Ligate the motor-loaded sequencing adapter", nanopore_entry_rows("PCR-free cDNA", amplified=False), "A dT-tailed sequencing adapter ligates to a dA-tailed cDNA end; ** marks the ligation boundary."),
    ]
