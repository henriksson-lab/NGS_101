"""Drop-BS molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from base_conversion import bisulfite_path
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import nextera_library,seg
from chemdraw import MolecularState,Row,Scene,Workflow,feature,oligo
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)

CELL_LEN=15
RANDOM_PRIMER="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"+"N"*9
FINAL_LIBRARY,SEQ_PRIMERS=nextera_library([
 seg("bisulfite-converted genomic insert","X"*38,placeholder=True),
 seg("bead-ligation handle","X"*20,placeholder=True),
 seg("15-nt bead cell barcode","D"*CELL_LEN,"cbc",placeholder=True,
     feature=feature("dropbs_cell","cell_barcode","whitelist",whitelist="Drop-BS bead barcode set")),
],"Drop-BS library")
TITLE="Drop-BS — droplet barcoding before bisulfite conversion"
NOTES="01_drop-bs.html"
SOURCE='<a href="https://doi.org/10.1038/s41467-023-40411-w">Fu et al., <i>Nature Communications</i> (2023)</a>.'
SUMMARY="MNase-fragmented DNA and a photocleavable barcode bead meet in one droplet. Ligation attaches cell identity before a second droplet conversion, and random priming plus indexed PCR completes the library."
FINAL_CAPTION="Nextera-ended Drop-BS library carrying the bead-derived 15-nt cell barcode on the converted genomic molecule."
SEQUENCING_INTRO="Paired-end genomic reads retain the bead barcode; i5 and i7 are added during the final indexing PCR."
BEAD_LIGATION_END="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
BLOCKED_SPLINT="CTGTCTCTTATACACATCT"  # 5′→3′; source figure prints this strand 3′→5′

def bead_ligation_geometry():
 sc=Scene();sc.strand("released bead oligo",[seg("bead barcode and PCR arm","CAAGCAGAAGACGGCATACGAGAT"+"D"*15,"cbc",placeholder=True),seg("ligation-end arm",BEAD_LIGATION_END[:-len(BLOCKED_SPLINT)],"me"),seg("splint-paired terminus",BEAD_LIGATION_END[-len(BLOCKED_SPLINT):],"me")],label="photocleaved bead oligo")
 sc.anneal("blocked complement",[seg("3′-blocked ligation splint",BLOCKED_SPLINT,"me")],to="released bead oligo",pair=("3′-blocked ligation splint","splint-paired terminus"),label="3′-SpC3; 5′-phosphorylated")
 sc.ligation_end("released bead oligo","splint-paired terminus",kind="splinted",polarity="3-prime",owner="adapter",text="published splinted bead-oligo 3′ end")
 sc.strand("insert DNA",[seg("prepared insert end","X"*20,placeholder=True)],label="MNase fragment after end repair")
 sc.ligation_end("insert DNA","prepared insert end",kind="splinted",polarity="5-prime",owner="substrate",text="insert end presented to splint")
 return AdapterLigationGeometry.known("splinted",sc,statement="Drop-BS Supplementary Figure 1 prints the bead ligation sequence and the 3′-SpC3-blocked, 5′-phosphorylated complementary oligo that presents the ligation site.")
def oligos(): yield oligo("random primer",[seg("S5+ME","TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG","me"),seg("N9","N"*9,placeholder=True)])
def sections(): return [
 ("Lyse and fragment one nucleus in a droplet",[Row(chunks=[("nucleus + MNase → genomic fragments",None,False)])],"MNase fragmentation is optimized independently of conversion."),
 ("Fuse with one barcode-bead droplet",bead_ligation_geometry().rows()+[Row(chunks=[("UV release: bead barcode ** genomic fragment", "cbc",False)])],"End repair and splint-supported ligation append the photocleaved 15-nt bead barcode."),
 ("Bisulfite-convert barcoded DNA in droplets",[Row(chunks=[(bisulfite_path(protected=False).text(),None,False)]),Row(chunks=[(bisulfite_path(protected=True).text(),"w1",False)])],"Cell identity is already covalent before droplets are broken."),
 ("Random-prime and index",[Row(chunks=[("S5—ME—N9 → extension; i5/i7 PCR completes P5/P7",None,False)])],"The printed N9 primer installs the second amplification handle."),]

def workflow():
 w=Workflow(MolecularState("MNase-fragmented nucleus and barcode bead",tuple(sections()[0][1])))
 w.ligate_adapter("Ligate the photocleaved bead oligo",bead_ligation_geometry(),[Row(chunks=[("bead barcode ** genomic fragment","cbc",False)])],note="The blocked complementary oligo aligns the published ligation end.")
 return w

def protocol_graph():
 source="Defining source: https://doi.org/10.1038/s41467-023-40411-w"
 graph=ProtocolGraph("Drop-BS",authority=source)
 nucleus=graph.input("nucleus","single nucleus",MaterialKind.BIOLOGICAL,primary=True)
 mnase=graph.input("mnase","micrococcal nuclease",MaterialKind.CATALYST)
 fragments=graph.product("mnase_fragments","MNase-fragmented genomic DNA",MaterialKind.MOLECULAR)
 graph.react("mnase_fragmentation","Lyse and fragment one nucleus with MNase",
             StepKind.FRAGMENTATION,
             ((nucleus,InputRole.TEMPLATE),(mnase,InputRole.CATALYST)),(fragments,))
 bead=graph.input("barcode_bead","photocleavable 15-nt barcode bead",MaterialKind.ADAPTER)
 fusion=graph.input("fusion_condition","bead-droplet fusion condition",MaterialKind.CONDITION)
 paired=graph.product("paired_droplet","fragment-and-barcode-bead droplet",MaterialKind.MOLECULAR)
 excluded=graph.product("fusion_loss","unpaired droplet material",MaterialKind.WASTE)
 graph.react("droplet_fusion","Fuse the fragment droplet with one barcode-bead droplet",
             StepKind.SAMPLING,
             ((fragments,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
              (fusion,InputRole.CONDITION)),(paired,excluded),changes_molecule=False)
 graph.discard(excluded,"unpaired or unusable fused droplets are not retained")
 repair=graph.input("repair_enzymes","DNA end-repair enzymes",MaterialKind.CATALYST)
 repaired=graph.product("repaired_fragments","end-repaired genomic fragments in droplets",
                        MaterialKind.MOLECULAR)
 graph.react("end_repair","End-repair the MNase fragments",StepKind.END_REPAIR,
             ((paired,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(repaired,))
 ligase=graph.input("ligase","DNA ligase",MaterialKind.CATALYST)
 barcoded=graph.product("barcoded_dna","bead-barcode-ligated genomic DNA",MaterialKind.MOLECULAR)
 graph.react("barcode_ligation","Photorelease and ligate the bead barcode",
             StepKind.LIGATION,
             ((repaired,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
              (ligase,InputRole.CATALYST)),(barcoded,),geometry=bead_ligation_geometry())
 bisulfite=graph.input("bisulfite","bisulfite-conversion condition",MaterialKind.CONDITION)
 converted=graph.product("converted_dna","bisulfite-converted barcoded DNA",MaterialKind.MOLECULAR)
 graph.react("bisulfite_conversion","Bisulfite-convert the barcoded DNA in droplets",
             StepKind.CHEMICAL_CONVERSION,
             ((barcoded,InputRole.TEMPLATE),(bisulfite,InputRole.CONDITION)),(converted,))
 random_primer=graph.input("random_primer","printed S5–ME–N9 random primer",MaterialKind.OLIGO)
 polymerase=graph.input("polymerase","random-extension and PCR polymerase",MaterialKind.CATALYST)
 extended=graph.product("extended_dna","random-primer-extended converted DNA",MaterialKind.MOLECULAR)
 graph.react("random_extension","Random-prime converted DNA to install the second handle",
             StepKind.EXTENSION,
             ((converted,InputRole.TEMPLATE),(random_primer,InputRole.STRUCTURAL),
              (polymerase,InputRole.CATALYST)),(extended,))
 i5=graph.input("i5_primer","indexed i5 PCR primer",MaterialKind.OLIGO)
 i7=graph.input("i7_primer","indexed i7 PCR primer",MaterialKind.OLIGO)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 graph.react("library_pcr","Add i5/i7 and PCR-complete the library",StepKind.PCR,
             ((extended,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
              (i7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
             geometry=PCRGeometry.unavailable(
                 FINAL_LIBRARY,reason="The compact final-library model does not expose separate PCR landing objects."))
 graph.sequence(library,SequencingPlan(
     library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
     acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
     read_uses={p.role:("demultiplex" if "Index" in p.role else "include") for p in SEQ_PRIMERS}))
 graph.declare_reviewed_scope((
     "nucleus","mnase","mnase_fragmentation","barcode_bead","fusion_condition","droplet_fusion",
     "repair_enzymes","end_repair","ligase","barcode_ligation","bisulfite",
     "bisulfite_conversion","random_primer","polymerase","random_extension",
     "i5_primer","i7_primer","library_pcr","library"),source)
 return graph.require_publishable()
