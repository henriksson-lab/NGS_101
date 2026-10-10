"""scGRO-seq clickable nascent-RNA chemistry."""
from chemdraw import Row
from multibranch import illumina_branch
from protocol_extensions import clicked_nascent_rna
TITLE="scGRO-seq — single-cell nascent RNA by click chemistry"
NOTES="01_scgro-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41586-024-07517-7">Mahat et al., <i>Nature</i> (2024)</a>.'
SUMMARY="Nuclear run-on incorporates a chain-terminating O-propargyl nucleotide; click ligation attaches the barcode/UMI/RT handle used to recover nascent RNA from single cells."
FINAL_LIBRARY,SEQ_PRIMERS=illumina_branch("scGRO-seq","nascent-RNA cDNA",cell=8,umi=8)
FINAL_CAPTION="Cell barcode and UMI precede cDNA copied from the run-on-labelled RNA 3-prime end."
SEQUENCING_INTRO="Declared Illumina primers bind the completed scGRO-seq library."
def sections():
    return [("Run on engaged polymerases with 3-prime O-propargyl NTPs",[Row(chunks=[("nascent RNA → RNA—O-propargyl-NMP (chain terminated)", "w1", False)])],"The terminal alkyne marks the active RNA 3-prime end."),("Click the azide barcode oligo to the alkyne",clicked_nascent_rna().rows(),"Copper-catalysed azide–alkyne cycloaddition creates a triazole linkage that RT can traverse."),("Reverse-transcribe and add sequencing adapters",[Row(chunks=[("clicked RNA → barcoded cDNA → PCR/index adapters", "r1", False)])],"The clicked oligo supplies molecule and cell identity before amplification."),]
