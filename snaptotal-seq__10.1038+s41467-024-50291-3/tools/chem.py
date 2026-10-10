"""snapTotal-seq MALBAC/template-switch total-RNA workflow."""
from batch_ngs import seg, truseq_library
from chemdraw import Row, feature
TITLE="snapTotal-seq — single-cell total-RNA sequencing"
NOTES="01_snaptotal-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-024-50291-3">Niu et al., <i>Nature Communications</i> (2024)</a>.'
SUMMARY="Low-temperature MALBAC priming and template switching amplify total RNA; a later double-strand conversion adds an 8-nt cell barcode to only one cDNA end."
FINAL_LIBRARY,SEQ_PRIMERS=truseq_library([seg("captured RNA sequence","X"*42,placeholder=True),seg("cell barcode","B"*8,"cbc",placeholder=True,feature=feature("cell_barcode","cell_barcode","whitelist",whitelist="Niu et al. Supplementary Data 12")),seg("barcode-to-UMI spacer","L"*35,placeholder=True),seg("UMI","U"*5,"umi",placeholder=True,feature=feature("umi","umi","random"))],"snapTotal-seq standard-primer library")
FINAL_CAPTION="The paper's standard-Illumina option reads RNA sequence in Read 1; Read 2 reaches the 8-nt cell barcode and 5-nt UMI."
SEQUENCING_INTRO="This page shows the paper's standard-primer alternative; the original run configuration also provides custom Read 2 and index primers."
def oligos(): return ["<b>Cell-barcode DSC primer</b> — 5′-GTGTGCTCTTCCGATCTNNNNNNNNAGGAGAGTGTGAGTGATGGTTGAGGATGTGTGGAG-3′"]
def sections(): return [
 ("Prime total RNA at low temperature",[Row(chunks=[("MALBAC primer → dispersed RNA sites","r1",False)])],"Random multiple annealing avoids dependence on poly(A)."),
 ("Template-switch and amplify cDNA",[Row(chunks=[("RNA → cDNA with common MALBAC-derived ends","tso",False)])],"Ten annealing/extension cycles precede cDNA PCR."),
 ("Add one cell barcode by double-strand conversion",[Row(chunks=[("one cDNA end ** 8-nt cell barcode ** UMI-bearing handle","cbc",False)])],"Adding the barcode after PCR avoids placing complementary barcode arms at both ends."),
 ("Pool, shear, enrich tagged fragments and dA-tail",[Row(chunks=[("pooled cDNA → 150–250-bp fragments → tagged subset → 3′ dA","me",False)])],"A second conversion enriches fragments carrying cell barcode and UMI."),
 ("Ligate Y adapter and complete sequencing ends",[Row(chunks=[("dA-tailed fragment ** Y adapter → indexed Illumina library","r2",False)])],"PCR produces either the custom-primer library or the standard-primer variant drawn below."),]
