"""SHAPE-MaP chemical probing and mutational profiling workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg, truseq_library
from chemdraw import Row
from rna_special import mutational_read_rows

TITLE = "SHAPE-MaP — RNA structure by mutational profiling"
NOTES = "01_shape-map.html"
SOURCE = 'Defining paper: <a href="https://doi.org/10.1038/nmeth.3029">Siegfried et al. (2014)</a>; implementation guidance: <a href="https://doi.org/10.1038/nprot.2018.044">Busan &amp; Weeks (2018)</a>.'
SUMMARY = "A SHAPE reagent acylates conformationally flexible RNA 2′-hydroxyls. MaP reverse transcription records those adducts as mutations in cDNA, which is converted to an indexed sequencing library."
CAVEAT = "SHAPE-MaP was selected over DMS-MaPseq because it is the defining general-purpose MaP workflow and has an authoritative full protocol. The chemical adduct is a mark, not an added sequencing adapter."
FINAL_LIBRARY, SEQ_PRIMERS = truseq_library([seg("MaP cDNA amplicon", "X"*44, placeholder=True)], "SHAPE-MaP library", inferred_adapters=True)
FINAL_CAPTION = "Indexed amplicon containing MaP-encoded substitutions at SHAPE-reactive positions."
SEQUENCING_INTRO = "Standard sequencing primers read the cDNA; structure is inferred from mutation frequencies rather than termination positions."

def sections():
    return [
        ("Modify flexible RNA nucleotides", [Row(chunks=[("folded RNA     paired stem ║ loop / flexible region ~~~~~", None, False)]), Row(chunks=[("SHAPE reagent                         ●  ●●   ●", "umi", False)])], "Electrophilic SHAPE reagent forms covalent 2′-O-adducts preferentially at flexible nucleotides; untreated and denatured controls accompany the modified sample."),
        ("Encode adducts as cDNA mutations", mutational_read_rows(), "MaP conditions encourage reverse transcriptase to traverse adducts while misincorporating, rather than stopping."),
        ("Amplify target and install sequencing arms", [Row(chunks=[("gene-specific PCR / library PCR → indexed double-stranded cDNA", None, False)]), Row(chunks=[("INFERRED — canonical flow-cell arms flank the assayed RNA-derived insert", None, True)])], "INFERRED — library-arm bases depend on the chosen amplicon implementation; the panel below states the modeled standard configuration."),
    ]
