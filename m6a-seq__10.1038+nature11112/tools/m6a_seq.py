"""Original m6A-seq / MeRIP-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg,truseq_library
from chemdraw import Row
from terminal_rna import immunoprecipitation_rows

TITLE="m6A-seq / MeRIP-seq — antibody enrichment of methylated RNA"
NOTES="01_m6a-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/nature11112">Dominissini et al., <i>Nature</i> (2012)</a>.'
SUMMARY="Poly(A)-selected RNA is fragmented, enriched with anti-m6A antibody, competitively eluted with free m6A, and converted into an Illumina cDNA library alongside an input control."
CAVEAT="The paper defines the enrichment assay but does not establish a unique adapter sequence architecture; the final conventional Illumina arms are therefore marked inferred."
FINAL_LIBRARY,SEQ_PRIMERS=truseq_library([seg("m6A-enriched RNA-derived cDNA","X"*56,placeholder=True)],"m6A-seq IP library",inferred_adapters=True)
FINAL_CAPTION="The IP and matched input use the same cDNA-library chemistry; enrichment is biological selection, not a retained molecular barcode."
SEQUENCING_INTRO="Standard Illumina primers read the enriched cDNA insert and sample indexes."
READ_LENGTHS={"Read 1":50,"Read 2":50}

def sections(): return [
 ("Select polyadenylated RNA and fragment it",[Row(chunks=[("poly(A)+ RNA → approximately 100-nt RNA fragments",None,False)])],"Fragmentation limits immunoprecipitation enrichment to a local interval around a modified residue."),
 ("Immunoprecipitate m6A-containing fragments",immunoprecipitation_rows("m6A"),"Anti-m6A antibody retains methylated fragments while an aliquot is kept as input."),
 ("Competitively elute with free m6A",[Row(chunks=[("antibody-bound RNA + excess free m6A → released enriched RNA", "w1",False)])],"The supplementary validation identifies free-m6A competition as the specific elution."),
 ("Convert IP and input RNA into cDNA libraries",[Row(chunks=[("RNA fragment → reverse transcription → adapter-complete cDNA",None,False)])],"IP enrichment is quantified relative to the matched input library."),
]
