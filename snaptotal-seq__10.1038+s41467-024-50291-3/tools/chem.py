"""snapTotal-seq MALBAC/template-switch total-RNA workflow."""
from batch_ngs import seg, truseq_library
from chemdraw import Construct, MolecularState, Row, Workflow, feature
from endprep import (ThreePrimeOverhangAdapter, cohesive_ligation_geometry,
                     dA_tailed_scene, ligation_product_scene, repair_and_dA_tail)
TITLE="snapTotal-seq — single-cell total-RNA sequencing"
NOTES="01_snaptotal-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-024-50291-3">Niu et al., <i>Nature Communications</i> (2024)</a>.'
SUMMARY="Low-temperature MALBAC priming and template switching amplify total RNA; a later double-strand conversion adds an 8-nt cell barcode to only one cDNA end."
FINAL_LIBRARY,SEQ_PRIMERS=truseq_library([seg("captured RNA sequence","X"*42,placeholder=True),seg("cell barcode","B"*8,"cbc",placeholder=True,feature=feature("cell_barcode","cell_barcode","whitelist",whitelist="Niu et al. Supplementary Data 12")),seg("barcode-to-UMI spacer","L"*35,placeholder=True),seg("UMI","U"*5,"umi",placeholder=True,feature=feature("umi","umi","random"))],"snapTotal-seq standard-primer library")
FINAL_CAPTION="The paper's standard-Illumina option reads RNA sequence in Read 1; Read 2 reaches the 8-nt cell barcode and 5-nt UMI."
SEQUENCING_INTRO="This page shows the paper's standard-primer alternative; the original run configuration also provides custom Read 2 and index primers."
Y_ADAPTER=ThreePrimeOverhangAdapter("snapTotal-seq Y-adapter ligating stem","GCTCTTCCGATCT","GATCGGAAGAGC")
def oligos(): return ["<b>Cell-barcode DSC primer</b> — 5′-GTGTGCTCTTCCGATCTNNNNNNNNAGGAGAGTGTGAGTGATGGTTGAGGATGTGTGGAG-3′"]
def workflow():
 w=Workflow(MolecularState("total RNA",(Row(chunks=[("cellular total RNA",None,False)]),)))
 w.react("Prime total RNA at low temperature",[Row(chunks=[("MALBAC primer → dispersed RNA sites","r1",False)])],name="multiply primed RNA",note="Random multiple annealing avoids dependence on poly(A).")
 w.react("Template-switch and amplify cDNA",[Row(chunks=[("RNA → cDNA with common MALBAC-derived ends","tso",False)])],name="amplified cDNA",note="Ten annealing/extension cycles precede cDNA PCR.")
 w.react("Add one cell barcode by double-strand conversion",[Row(chunks=[("one cDNA end ** 8-nt cell barcode ** UMI-bearing handle","cbc",False)])],name="barcoded cDNA",note="Adding the barcode after PCR avoids complementary barcode arms at both ends.")
 tailed=repair_and_dA_tail(Construct([seg("tagged cDNA fragment","X"*36,placeholder=True)]))
 w.react("Pool, shear, enrich tagged fragments and add 3′ dA",dA_tailed_scene(tailed).rows(),name="dA-tailed tagged fragments",note="A second conversion enriches fragments carrying cell barcode and UMI.")
 geom=cohesive_ligation_geometry(tailed,Y_ADAPTER,statement="Supplementary Data 12 prints both Y-adapter oligos; their ligating stem leaves a terminal 3′-T complementary to each insert 3′-dA.")
 product=ligation_product_scene(FINAL_LIBRARY,"captured RNA sequence","UMI",label="snapTotal-seq library")
 w.ligate_adapter("Ligate the Y adapter",geom,product.rows(),name="indexed Illumina library",note="PCR produces either the custom-primer or standard-primer variant drawn below.")
 return w
