"""CoBATCH molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct,MolecularState,Row,Scene,Workflow,duplex_rows,feature,revcomp
from chromatin_epigenetics import antibody_enzyme_rows
from tagmentation import loaded_tn5_adapter
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)

TITLE="CoBATCH — combinatorial barcoding by targeted tagmentation"
NOTES="01_cobatch.html"
SOURCE='Defining source: <a href="https://doi.org/10.1016/j.molcel.2019.07.015">Wang et al., <i>Molecular Cell</i> (2019)</a>.'
SUMMARY="Antibody-tethered Protein-A–Tn5 transfers a pair of well barcodes at the chromatin target; pooled cells are redistributed and receive a second i5/i7 well identity during PCR."
CAVEAT="CoBATCH cell identity has four indexed parts: the paired T5/T7 transposome barcodes and the paired PCR-well indexes. They are modeled as cell barcodes, not sample indexes."
def cb(identifier,part): return feature(identifier,"cell_barcode","combinatorial",group="cobatch_cell",part=part,whitelist="published CoBATCH barcode plate")
T5_OUTER="TCCACGC"
T5_INNER="GCGATCGAGGACGGC"
T7_OUTER="CTGTCCCTGTCC"
T7_INNER="CACCGTCTCCGCCTC"
T5_BARCODES=("TATAGCCT","ATAGAGGC","CCTATCCT","GGCTCTGA")
T7_BARCODES=("CGAGTAAT","TCTCCGGA","AATGAGCG","GGAATCTC")
TN5_ME_REV=nx.ME_RC

def t5_adapter(barcode=T5_BARCODES[0]):
 if barcode not in T5_BARCODES: raise ValueError("representative T5 barcode must be published")
 return [seg("S5",nx.S5,"s5"),seg("T5 outer connector",T5_OUTER),
         seg("T5 well barcode",barcode,"cbc",feature=cb("cobatch_t5","PAT-T5 well")),
         seg("T5 inner connector",T5_INNER),seg("T5 mosaic end",nx.ME,"me")]

def t7_adapter(barcode=T7_BARCODES[1]):
 if barcode not in T7_BARCODES: raise ValueError("representative T7 barcode must be published")
 return [seg("S7",nx.S7,"s7"),seg("T7 outer connector",T7_OUTER),
         seg("T7 well barcode",barcode,"cbc",feature=cb("cobatch_t7","PAT-T7 well")),
         seg("T7 inner connector",T7_INNER),seg("T7 mosaic end",nx.ME,"me")]

def loaded_scene(parts,label):
 sc=loaded_tn5_adapter(parts,mosaic_end=parts[-1].name,label=label,
                       bottom_name="common annealing primer",
                       bottom_label="5′-phosphorylated common ME bottom",
                       bottom_sequence=TN5_ME_REV,mod5="p")
 bc=next(s.name for s in parts if s.feature)
 sc.mark("Tn5 top",bc,"first-plate well identity")
 return sc

def transferred_fragment():
 right=t7_adapter()
 right_rc=[seg(s.name+" reverse complement",revcomp(s.top),s.tag,feature=s.feature)
           for s in reversed(right)]
 return Construct([*t5_adapter(),seg("targeted chromatin fragment","X"*42,placeholder=True),
                   *right_rc],name="gap-filled target-coded fragment")

TRANSFERRED=transferred_fragment()
FINAL_LIBRARY=Construct([
 seg("P5",il.P5,"p5"),seg("PCR-well i5","J"*8,"cbc",placeholder=True,feature=cb("cobatch_i5","PCR well i5")),seg("S5",nx.S5,"s5"),
 seg("T5 outer connector",T5_OUTER),seg("T5 well barcode",T5_BARCODES[0],"cbc",feature=cb("cobatch_t5","PAT-T5 well")),seg("T5 inner connector",T5_INNER),seg("left mosaic end",nx.ME,"me"),
 seg("targeted chromatin fragment","X"*42,placeholder=True),seg("right mosaic end reverse complement",nx.ME_RC,"me"),
 seg("T7 inner connector reverse complement",revcomp(T7_INNER)),seg("T7 well barcode reverse complement",revcomp(T7_BARCODES[1]),"cbc",feature=cb("cobatch_t7","PAT-T7 well")),seg("T7 outer connector reverse complement",revcomp(T7_OUTER)),seg("S7 reverse complement",nx.S7_RC,"s7"),
 seg("PCR-well i7 reverse complement","I"*8,"cbc",placeholder=True,feature=cb("cobatch_i7","PCR well i7")),seg("P7 reverse complement",il.P7_RC,"p7")],name="CoBATCH library")
COBATCH_READ1=T5_INNER+nx.ME
COBATCH_READ2=T7_INNER+nx.ME
SEQ_PRIMERS=(
 sp.custom("Read 1","CoBATCH extended mosaic-end T5 primer",COBATCH_READ1,"Wang et al. supplementary oligo table"),
 sp.custom("Read 2","CoBATCH extended mosaic-end T7 primer",COBATCH_READ2,"Wang et al. supplementary oligo table"),
)
_problems=sp.verify(FINAL_LIBRARY,SEQ_PRIMERS,required_roles=("Read 1","Read 2"))
if _problems: raise ValueError("CoBATCH library: "+"; ".join(_problems))
FINAL_CAPTION="The first plate is encoded beside both mosaic ends; the second plate is read in i5 and i7. Their combination identifies a cell."
SEQUENCING_INTRO="Paired genomic reads map targeted chromatin fragments; index and inline cycles recover the two-round cell identity."
SEQUENCING_UNAVAILABLE="The fetched source establishes i5/i7 PCR indexes but does not state a standalone custom index-read primer sequence; the two genomic read-primer sites are shown exactly."
READ_LENGTHS={"Read 1":50,"Read 2":50}
SHOW_READ_CONTENT=True

def oligos(): return (
 f"<b>T5-1</b> — 5′-{''.join(s.top for s in t5_adapter())}-3′",
 f"<b>T7-2</b> — 5′-{''.join(s.top for s in t7_adapter())}-3′",
 f"<b>common annealing primer</b> — 5′-/5Phos/{TN5_ME_REV}-3′",
 f"<b>representative i5 PCR primer</b> — 5′-{il.P5}[8-nt PCR-well i5]{nx.S5}-3′",
 f"<b>representative i7 PCR primer</b> — 5′-{il.P7}[8-nt PCR-well i7]{nx.S7}-3′",
)

def workflow():
 start=Scene.duplex([seg("native chromatin","X"*42,placeholder=True)],label="antibody-bound chromatin")
 wf=Workflow(MolecularState("Antibody-bound chromatin in first-plate wells",tuple(start.rows())))
 wf.react("Anneal the T5 and T7 barcode oligos to common ME bottom and load pA–Tn5",
          [*loaded_scene(t5_adapter(),"representative PAT-T5").rows(),
           *loaded_scene(t7_adapter(),"representative PAT-T7").rows()],
          name="First-plate barcoded pA–Tn5",
          note="The representative T5-1 and T7-2 sequences are drawn; the source table supplies the full whitelist.")
 wf.react("Tether through antibody, activate tagmentation, and fill the staggered gaps",
          duplex_rows(TRANSFERRED,label="T5/T7 target-coded fragment"),
          name="First-plate-indexed chromatin fragments")
 wf.react("Pool cells, redistribute them, and PCR with second-well i5/i7 primers",
          duplex_rows(FINAL_LIBRARY,label="four-part CoBATCH library"),
          name="CoBATCH sequencing library",
          note="The two transposome barcodes plus the two PCR-well indexes jointly identify one cell.")
 return wf

def protocol_graph():
 source="Defining source: https://doi.org/10.1016/j.molcel.2019.07.015"
 graph=ProtocolGraph("CoBATCH",authority=source)
 chromatin=graph.input("chromatin","antibody-bound chromatin",MaterialKind.BIOLOGICAL,
                       primary=True,state=Scene.duplex([seg("native chromatin","X"*42,placeholder=True)]))
 transposomes=graph.input("barcoded_transposomes","T5/T7 first-plate pA-Tn5 complexes",
                          MaterialKind.ADAPTER)
 tn5=graph.input("pa_tn5","protein-A–Tn5 transposase",MaterialKind.CATALYST)
 tagged=graph.product("tagged_fragments","gapped first-plate-indexed chromatin fragments",
                      MaterialKind.MOLECULAR,state=TRANSFERRED)
 tagmentation_waste=graph.product("tagmentation_waste","non-library tagmentation products",
                                  MaterialKind.WASTE)
 graph.react("targeted_tagmentation","Tether pA-Tn5 and transfer the T5/T7 well adapters",
             StepKind.TAGMENTATION,
             ((chromatin,InputRole.TEMPLATE),(transposomes,InputRole.STRUCTURAL),
              (tn5,InputRole.CATALYST)),(tagged,tagmentation_waste))
 graph.discard(tagmentation_waste,"molecules without a usable paired T5/T7 insertion do not amplify")
 polymerase=graph.input("polymerase","gap-fill and PCR polymerase",MaterialKind.CATALYST)
 filled=graph.product("gap_filled_fragments","gap-filled first-plate-indexed fragments",
                      MaterialKind.MOLECULAR,state=TRANSFERRED)
 graph.react("gap_fill","Fill the staggered transposase gaps",StepKind.EXTENSION,
             ((tagged,InputRole.TEMPLATE),(polymerase,InputRole.CATALYST)),(filled,))
 allocation=graph.input("second_plate_allocation","second-plate redistribution rule",
                        MaterialKind.SELECTION)
 allocated=graph.product("allocated_fragments","fragments assigned to second-plate wells",
                         MaterialKind.MOLECULAR,state=TRANSFERRED)
 excluded=graph.product("allocation_loss","fragments excluded during redistribution",
                        MaterialKind.WASTE)
 graph.react("redistribute","Pool cells and redistribute them to second-plate wells",
             StepKind.SAMPLING,
             ((filled,InputRole.TEMPLATE),(allocation,InputRole.SELECTION)),(allocated,excluded),
             changes_molecule=False)
 graph.discard(excluded,"material not retained in an assayed second-plate well is not sequenced")
 i5=graph.input("i5_primer","PCR-well i5 primer",MaterialKind.OLIGO)
 i7=graph.input("i7_primer","PCR-well i7 primer",MaterialKind.OLIGO)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 graph.react("second_well_pcr","Amplify and add the second-well i5/i7 identities",StepKind.PCR,
             ((allocated,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
              (i7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
             geometry=PCRGeometry.unavailable(
                 FINAL_LIBRARY,reason="The source establishes the indexed PCR ends, but the page does not model exact primer annealing geometry."))
 graph.sequence(library,SequencingPlan(
     library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
     acquired_roles=tuple(p.role for p in SEQ_PRIMERS),read_lengths=READ_LENGTHS,
     read_uses={"Read 1":"include","Read 2":"include"},
     limitation=SEQUENCING_UNAVAILABLE))
 graph.declare_reviewed_scope((
     "chromatin","barcoded_transposomes","pa_tn5","targeted_tagmentation",
     "polymerase","gap_fill","second_plate_allocation","redistribute",
     "i5_primer","i7_primer","second_well_pcr","library"),source)
 return graph.require_publishable()
