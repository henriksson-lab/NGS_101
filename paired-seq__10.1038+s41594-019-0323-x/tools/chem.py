"""Paired-seq joint RNA and chromatin accessibility."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import MolecularState,Row,Scene,Segment,Workflow,revcomp
from multibranch import illumina_branch
TITLE="Paired-seq — combinatorial RNA and chromatin accessibility"
NOTES="01_paired-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41594-019-0323-x">Zhu et al., <i>Nature Structural &amp; Molecular Biology</i> (2019)</a>.'
SUMMARY="Permeabilized nuclei undergo RNA reverse transcription and chromatin tagmentation before four rounds of split-pool cellular indexing."
CAVEAT="Barcode rounds and branch chemistry are source-established. The compact final Illumina shells are marked inferred because the full protocol-specific molecules have not yet been transcribed base by base."
BC=(("round 1 barcode",3),("round 2 barcode",7),("round 3 barcode",7),("round 4 barcode",7))
SPLINTS=("CGAATGCTCTGGCCTCTCAAGCACGTGGAT","GGTCTGAGTTCGCACACGCTCCTCGAATCC","AGTCGTACGCCGATGCGAAACATCGGCCAC")
RNA,RP=illumina_branch("Paired-seq RNA library","cDNA fragment",cell_parts=BC,umi=10,inferred_adapters=True)
ATAC,AP=illumina_branch("Paired-seq ATAC library","accessible chromatin fragment",architecture="nextera",cell_parts=BC,inferred_adapters=True)
FINAL_LIBRARIES=(("RNA library",RNA,RP,"cDNA reads carry the four-part cell barcode and UMI.","Illumina primers read the RNA branch."),("Chromatin-accessibility library",ATAC,AP,"Tagmented DNA reads carry the same four-part cell barcode.","Nextera-compatible primers read the ATAC branch."))
def splinted_round(round_no,splint):
 rc=revcomp(splint)
 sc=Scene(); sc.strand("splint",[Segment(f"round-{round_no} new-arm half",splint[15:]),Segment(f"round-{round_no} previous-arm half",splint[:15])],label=f"Linker-R0{round_no}")
 sc.anneal("ligating oligos",[Segment(f"new round-{round_no} 5-prime arm",rc[:15]),Segment("previous-round 3-prime arm",rc[15:])],to="splint",pair=(f"new round-{round_no} 5-prime arm",f"round-{round_no} new-arm half"),label="two adjacent barcode oligos",mod5="p")
 sc.ligation_end("ligating oligos",f"new round-{round_no} 5-prime arm",kind="splinted",polarity="5-prime",owner="adapter",text="5′-phosphorylated nick partner")
 sc.footer("splint holds the new 5′ end beside the previous 3′ end; ligase seals that nick","ligating oligos",f"new round-{round_no} 5-prime arm")
 return AdapterLigationGeometry.known("splinted",sc,statement=f"Supplementary oligos print Linker-R0{round_no} and both 15-base annealing arms; the barcode bodies remain unpaired outside this duplex.")

def sections():
 wf=Workflow(MolecularState("Permeabilized nuclei",(Row(chunks=[("RNA + accessible chromatin",None,False)]),)))
 wf.react("Reverse-transcribe RNA and tagment chromatin in nuclei",[Row(chunks=[("RNA → cDNA     accessible DNA ** Tn5 ends","me",False)])],note="Both analytes remain inside the same permeabilized nucleus; these reactions install BC1.")
 for n,s in enumerate(SPLINTS,2):
  wf.ligate_adapter(f"Ligate barcode round {n}",splinted_round(n,s),[Row(chunks=[(f"BC1 ** " + " ** ".join(f"BC{x}" for x in range(2,n+1)) + " ** RNA/ATAC molecule","cbc",False)])],note=f"The published Linker-R0{n} aligns the two covalent ends before ligation.")
 wf.react("Add a terminal homopolymer and anchor",[Row(chunks=[("barcoded molecule—poly(C) ← G-anchor extension","tso",False)])],note="Terminal transferase and anchored amplification make indexed material amplifiable.")
 wf.react("INFERRED — complete RNA and ATAC sequencing shells",[Row(chunks=[("restriction split → RNA library | ATAC library","r2",False)])],note="Branch-selective restriction is source-established; the compact outer adapter completion is inferred.")
 return [(r.action,r.after.rows,r.note) for r in wf.reactions]
