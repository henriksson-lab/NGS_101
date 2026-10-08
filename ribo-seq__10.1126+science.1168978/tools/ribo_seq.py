"""Original ribosome-profiling molecular workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg, truseq_library
from chemdraw import Row
from rna_special import adapter_ligation_scene, circular_cdna_rows

TITLE = "Ribo-seq — ribosome profiling"
NOTES = "01_ribo-seq.html"
SOURCE = 'Defining paper: <a href="https://doi.org/10.1126/science.1168978">Ingolia et al. (2009)</a>; library details: <a href="https://doi.org/10.1038/nprot.2012.086">Ingolia et al. (2012)</a>.'
SUMMARY = "Nuclease leaves ribosome-protected RNA footprints. Size-selected footprints receive a 3′ adapter, are reverse-transcribed and circularized, then PCR adds the Illumina run structure."
CAVEAT = "The schematic preserves the published reaction order. Final flow-cell arms are shown as inferred canonical TruSeq roles because the historical protocol used experiment-specific index primers."
FINAL_LIBRARY, SEQ_PRIMERS = truseq_library([seg("ribosome footprint cDNA", "X"*30, placeholder=True)], "Ribo-seq library", inferred_adapters=True)
FINAL_CAPTION = "PCR-linearized library from a circular footprint cDNA; inferred canonical run arms are dotted."
SEQUENCING_INTRO = "Read 1 enters the reverse-transcribed ribosome footprint; index reads identify the sample."

def sections():
    return [
        ("Digest unprotected RNA", [Row(chunks=[("mRNA ——— [80S ribosome protects ~28–30 nt] ——— mRNA", None, False)]), Row(chunks=[("RNase I ↓          retained ribosome footprint", "umi", False)])], "RNase digestion and monosome isolation select ribosome-protected fragments."),
        ("Ligate the footprint 3′ adapter", adapter_ligation_scene(fragment_name="ribosome footprint").rows(), "A preadenylated adapter is joined to the footprint 3′ end before reverse transcription."),
        ("Reverse-transcribe and circularize cDNA", circular_cdna_rows([seg("footprint cDNA", "X"*30, placeholder=True), seg("adapter-derived PCR handle", "X"*18, "r2", placeholder=True)]), "The detailed protocol circularizes first-strand cDNA, providing the second PCR priming boundary without a second RNA ligation."),
    ]
