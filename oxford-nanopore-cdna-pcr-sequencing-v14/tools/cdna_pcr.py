"""Oxford Nanopore SQK-PCS114 cDNA-PCR workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Row
from rna_special import nanopore_entry_rows, template_switch_scene

TITLE = "Oxford Nanopore cDNA-PCR Sequencing V14 — SQK-PCS114"
NOTES = "01_cdna-pcr.html"
SOURCE = 'Commercial protocol: Oxford Nanopore <a href="https://nanoporetech.com/document/pcr-cdna-sequencing-v14-sqk-pcs114">cDNA-PCR Sequencing V14, SQK-PCS114</a>.'
SUMMARY = "Strand-switch reverse transcription adds a UMI and terminal handles, PCR selects and amplifies full-length cDNA while adding rapid-attachment tags, and a motor-loaded Rapid Adapter prepares one cDNA strand for nanopore sequencing."
CAVEAT = "SQK-PCS114 oligo and adapter bases are proprietary. Roles, UMI placement and reaction order are vendor-supported; unknown bases remain role placeholders."
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
SEQUENCING_ENDING = "There is no synthesis sequencing primer. The motor-loaded Rapid Adapter captures one strand of amplified cDNA at the pore and meters that strand through the nanopore; the adapter-bearing end enters first."

def sections():
    return [
        ("Reverse-transcribe and strand-switch", template_switch_scene(inferred=True).rows(), "INFERRED — proprietary CRTA/RTP/SSPII bases are not drawn. The current vendor protocol states that strand switching incorporates a UMI."),
        ("Select full-length molecules by PCR", [Row(chunks=[("5′ rapid-attachment tag — UMI — full-length cDNA — rapid-attachment tag 3′", None, True)]), Row(chunks=[("                             PCR amplification", "umi", False)])], "INFERRED — cDNA Primer amplification enriches molecules carrying both terminal handles and installs the rapid-adapter attachment tags."),
        ("Attach the motor-loaded Rapid Adapter", nanopore_entry_rows("amplified cDNA", amplified=True), "The Rapid Adapter assembles at the PCR-added terminal tag; ** marks the adapter attachment boundary."),
    ]
