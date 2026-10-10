"""CoBATCH molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct,Row,feature
from chromatin_epigenetics import antibody_enzyme_rows

TITLE="CoBATCH — combinatorial barcoding by targeted tagmentation"
NOTES="01_cobatch.html"
SOURCE='Defining source: <a href="https://doi.org/10.1016/j.molcel.2019.07.015">Wang et al., <i>Molecular Cell</i> (2019)</a>.'
SUMMARY="Antibody-tethered Protein-A–Tn5 transfers a pair of well barcodes at the chromatin target; pooled cells are redistributed and receive a second i5/i7 well identity during PCR."
CAVEAT="CoBATCH cell identity has four indexed parts: the paired T5/T7 transposome barcodes and the paired PCR-well indexes. They are modeled as cell barcodes, not sample indexes."
def cb(identifier,part): return feature(identifier,"cell_barcode","combinatorial",group="cobatch_cell",part=part,whitelist="published CoBATCH barcode plate")
FINAL_LIBRARY=Construct([
 seg("P5",il.P5,"p5"),seg("PCR-well i5","J"*8,"cbc",placeholder=True,feature=cb("cobatch_i5","PCR well i5")),seg("S5",nx.S5,"s5"),
 seg("T5 outer connector","TCCACGC"),seg("T5 well barcode","A"*8,"cbc",placeholder=True,feature=cb("cobatch_t5","PAT-T5 well")),seg("T5 inner connector","GCGATCGAGGACGGC"),seg("left mosaic end",nx.ME,"me"),
 seg("targeted chromatin fragment","X"*42,placeholder=True),seg("right mosaic end reverse complement",nx.ME_RC,"me"),
 seg("T7 inner connector reverse complement","GAGGCGGAGACGGTG"),seg("T7 well barcode reverse complement","B"*8,"cbc",placeholder=True,feature=cb("cobatch_t7","PAT-T7 well")),seg("T7 outer connector reverse complement","GGACAGGGACAG"),seg("S7 reverse complement",nx.S7_RC,"s7"),
 seg("PCR-well i7 reverse complement","I"*8,"cbc",placeholder=True,feature=cb("cobatch_i7","PCR well i7")),seg("P7 reverse complement",il.P7_RC,"p7")],name="CoBATCH library")
COBATCH_READ1=nx.S5
COBATCH_READ2=nx.S7
SEQ_PRIMERS=(
 sp.custom("Read 1","CoBATCH S5-side read primer",COBATCH_READ1,"Wang et al. supplementary oligo table"),
 sp.custom("Read 2","CoBATCH S7-side read primer",COBATCH_READ2,"Wang et al. supplementary oligo table"),
)
_problems=sp.verify(FINAL_LIBRARY,SEQ_PRIMERS,required_roles=("Read 1","Read 2"))
if _problems: raise ValueError("CoBATCH library: "+"; ".join(_problems))
FINAL_CAPTION="The first plate is encoded beside both mosaic ends; the second plate is read in i5 and i7. Their combination identifies a cell."
SEQUENCING_INTRO="Paired genomic reads map targeted chromatin fragments; index and inline cycles recover the two-round cell identity."
SEQUENCING_UNAVAILABLE="The fetched source establishes i5/i7 PCR indexes but does not state a standalone custom index-read primer sequence; the two genomic read-primer sites are shown exactly."
READ_LENGTHS={"Read 1":50,"Read 2":50}

def sections(): return [
 ("Bind antibody at the chromatin feature",antibody_enzyme_rows("Tn5","targeted adapter transfer"),"Protein A tethers Tn5 to the antibody-bound histone mark or chromatin protein."),
 ("Distribute cells into the first plate and target-tagment",[Row(chunks=[("PAT-T5 barcode → target fragment ← PAT-T7 barcode","cbc",False)])],"Each well uses a defined combination of barcoded T5 and T7 transposomes."),
 ("Pool cells and redistribute into PCR wells",[Row(chunks=[("T5/T7 well identity + i5/i7 PCR-well identity = cell identity","cbc",False)])],"Only the combinatorial pair of wells is treated as a single-cell label."),
 ("PCR-release and complete the sequencing library",[Row(chunks=[("P5 — i5 — T5 barcode — ME — insert — ME — T7 barcode — i7 — P7",None,False)])],"PCR both releases indexed fragments and installs complete flow-cell arms."),
]
