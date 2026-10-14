"""EpiDamID / scDam&T-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct,MolecularState,Row,Scene,Segment,Workflow,annotation_rows,oligo
from adapter_ligation import AdapterLigationGeometry
from damid import DamMarkedGATC,damid_adapter
from single_cell_hic import unresolved_illumina_library
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)

ADRT="CTAATACGACTCACTATAGGGCAGCGTGGTCGCGGCCGAGGA"
ADRB="TCCTCGGCCGCG"
ADR_PCR="GGTCGCGGCCGAGGATC"
RANDOMHEX_RT="GCCTTGGCACCCGAGAATTCCANNNNNN"
RNA_PCR_1="AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA"
RNA_PCR_INDEX="CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
MARK=DamMarkedGATC("a reader-domain-bound histone modification")
ADAPTER=damid_adapter(ADRT,ADRB)

def mark_scene():
 sc=Scene.duplex(list(MARK.locus()),label="living-cell chromatin record");sc.note("top","targeting domain recruits Dam; m6A accumulates at nearby GATC sites");sc.labels("top");return sc
def adapter_scene():
 sc=Scene();sc.strand("fragment",[Segment("DpnI-cut genomic fragment","X"*34,placeholder=True)],label="methylated GATC fragment");sc.strand("adapter",list(ADAPTER),label="AdRt/AdRb adapter");sc.note("adapter","T4 ligase joins adapters to DpnI-selective fragment ends");sc.labels("adapter");return sc
def adapter_ligation_geometry():
 sc=Scene();sc.strand("DpnI fragment",[Segment("blunt DpnI end","X"*18,placeholder=True)],label="DpnI-cut methylated GATC")
 sc.strand("AdRt",[Segment("AdRt T7 arm",ADRT[:-len(ADRB)],"t7"),Segment("AdRt duplex end",ADRT[-len(ADRB):],"t7")],label="AdRt")
 sc.anneal("AdRb",[Segment("AdRb",ADRB,"t7")],to="AdRt",pair=("AdRb","AdRt duplex end"),label="AdRb")
 sc.ligation_end("DpnI fragment","blunt DpnI end",kind="blunt",owner="substrate",text="blunt DpnI end")
 sc.ligation_end("AdRt","AdRt duplex end",kind="blunt",owner="adapter",text="blunt DamID2 adapter end")
 return AdapterLigationGeometry.known("blunt",sc,statement="DpnI cleaves methylated GATC with blunt ends; the paper prints the annealed AdRt/AdRb DamID2 adapter strands.")

TAGGED=Construct([Segment("T7 DamID adapter",ADRT,"t7"),Segment("Dam-marked genomic insert","X"*38,placeholder=True),Segment("opposite DamID adapter","X"*18,placeholder=True,inferred=True)],name="adapter-ligated DamID fragment")
FINAL_LIBRARY=unresolved_illumina_library(TAGGED,"scDam&T / CEL-Seq2-derived EpiDamID library",indexed=True)
SEQ_PRIMERS=()

def workflow():
 wf=Workflow(MolecularState("Targeted m6A recording in living cells",tuple(mark_scene().rows())))
 wf.react("DpnI digestion of methylated GATC",mark_scene().rows(),note="DpnI makes recovery selective for Dam-marked sites.")
 wf.ligate_adapter("Ligate T7-bearing DamID2 adapters",adapter_ligation_geometry(),adapter_scene().rows(),note="The printed AdRt/AdRb blunt duplex installs the T7 promoter and PCR handle.")
 wf.react("T7 in-vitro transcription",[Row(chunks=[("one adapter-ligated DNA template → many amplified RNA copies", "t7", False)])],note="Pooled, cell-indexed templates are linearly amplified as RNA for 14 h.")
 wf.react("CEL-Seq2 fragmentation, reverse transcription and PCR",Scene.duplex(list(FINAL_LIBRARY),label="sequencing library",unpaired=tuple(s.name+"'" for s in FINAL_LIBRARY if s.is_role_token())).rows(),note="The DNA and transcriptome branches are completed using the published scDam&T/CEL-Seq2 workflow.")
 return wf

TITLE="EpiDamID — histone-targeted adenine recording with scDam&T-seq"
NOTES="01_epidamid.html"
SOURCE='<a href="https://doi.org/10.1016/j.molcel.2022.03.009">Rang and de Luca et al., <i>Molecular Cell</i> (2022)</a>.'
SUMMARY="A histone-binding domain or mintbody recruits Dam in living cells, leaving m6A at nearby GATC sites. DpnI-selective adapter ligation and T7 amplification recover that chromatin record; scDam&T-seq can pair it with transcriptomes."
CAVEAT="INFERRED — the source prints the DamID2 adapter, random-hexamer RT primer and terminal PCR primers, but not one complete strand-resolved post-CEL-Seq2 library. The final outer arms are therefore unresolved rather than silently assigned to a TruSeq generation."
FINAL_CAPTION="EpiDamID insert after the scDam&T/CEL-Seq2 branch. Exact DamID adapter sequence is retained; unresolved outer sequencing regions are visibly inferred."
SEQUENCING_ENDING="Libraries were sequenced on NextSeq 500 (75-base reads) or NextSeq 2000 (100-base reads)."
SEQUENCING_UNAVAILABLE="A complete final sequencing-primer binding construct cannot be verified from the paper's split DamID2 and CEL-Seq2 descriptions; the source-transcribed primers remain available below."
def oligos():
 yield oligo("AdRt",[Segment("adapter top",ADRT,"t7")])
 yield oligo("AdRb",[Segment("adapter bottom",ADRB,"t7")])
 yield oligo("AdR_PCR",[Segment("PCR primer",ADR_PCR,"r1")])
 yield oligo("RandomhexRT",[Segment("random-hexamer RT primer",RANDOMHEX_RT,"r1")])
 yield oligo("RNA PCR primer 1",[Segment("RNA PCR 1",RNA_PCR_1,"r1")],mods="terminal G*A phosphorothioate")
 yield oligo("RNA PCR index primer, example CGTGATGT",[Segment("indexed RNA PCR",RNA_PCR_INDEX,"r2")],mods="terminal C*A phosphorothioate")

def protocol_graph():
 source="Defining source: https://doi.org/10.1016/j.molcel.2022.03.009"
 unresolved=source+" (complete post-CEL-Seq2 strand architecture is not printed)"
 graph=ProtocolGraph("EpiDamID",authority=source)
 chromatin=graph.input("chromatin","living-cell chromatin",MaterialKind.BIOLOGICAL,
                       primary=True,state=mark_scene())
 dam=graph.input("targeted_dam","histone-targeted Dam methylase",MaterialKind.CATALYST)
 marked=graph.product("marked_chromatin","m6A-marked chromatin record",MaterialKind.MOLECULAR,
                      state=mark_scene())
 graph.react("m6a_recording","Record nearby GATC sites through targeted Dam methylation",
             StepKind.CHEMICAL_CONVERSION,
             ((chromatin,InputRole.TEMPLATE),(dam,InputRole.CATALYST)),(marked,))
 dpni=graph.input("dpni","DpnI",MaterialKind.CATALYST)
 fragments=graph.product("dpni_fragments","DpnI-selected methylated GATC fragments",
                         MaterialKind.MOLECULAR)
 graph.react("dpni_digest","Digest methylated GATC with DpnI",StepKind.DIGESTION,
             ((marked,InputRole.TEMPLATE),(dpni,InputRole.CATALYST)),(fragments,))
 adapter=graph.input("damid2_adapter","T7-bearing DamID2 AdRt/AdRb adapter",
                     MaterialKind.ADAPTER,state=adapter_scene())
 ligase=graph.input("ligase","T4 DNA ligase",MaterialKind.CATALYST)
 tagged=graph.product("adapter_ligated",TAGGED.name,MaterialKind.MOLECULAR,state=TAGGED)
 graph.react("adapter_ligation","Ligate the T7-bearing DamID2 adapter",StepKind.LIGATION,
             ((fragments,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
              (ligase,InputRole.CATALYST)),(tagged,),geometry=adapter_ligation_geometry())
 t7=graph.input("t7_polymerase","T7 RNA polymerase",MaterialKind.CATALYST)
 rna=graph.product("ivt_rna","linearly amplified DamID RNA copies",MaterialKind.MOLECULAR)
 graph.react("ivt","Amplify adapter-ligated templates by T7 transcription",
             StepKind.TRANSCRIPTION,
             ((tagged,InputRole.TEMPLATE),(t7,InputRole.CATALYST)),(rna,))
 fragmentation=graph.input("fragmentation_condition","CEL-Seq2 RNA fragmentation condition",
                           MaterialKind.CONDITION)
 fragmented=graph.product("fragmented_rna","fragmented DamID RNA",MaterialKind.MOLECULAR)
 graph.react("fragmentation","Fragment the T7-amplified RNA",StepKind.FRAGMENTATION,
             ((rna,InputRole.TEMPLATE),(fragmentation,InputRole.CONDITION)),(fragmented,))
 rt_primer=graph.input("rt_primer","source-printed RandomhexRT reverse-transcription primer",
                       MaterialKind.OLIGO,source=source)
 rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
 cdna=graph.product("cdna","reverse-transcribed DamID cDNA",MaterialKind.MOLECULAR,
                    source=unresolved,inferred=True)
 graph.react("reverse_transcription","Reverse-transcribe the fragmented DamID RNA",
             StepKind.REVERSE_TRANSCRIPTION,
             ((fragmented,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
              (rt,InputRole.CATALYST)),(cdna,),source=unresolved)
 p5=graph.input("p5_primer","source-printed RNA PCR primer 1",MaterialKind.OLIGO,
                source=source)
 p7=graph.input("p7_primer","source-printed RNA PCR index primer",MaterialKind.OLIGO,
                source=source)
 polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY,
                       source=unresolved,inferred=True)
 graph.react("library_pcr","PCR-complete the scDam&T/CEL-Seq2 library",StepKind.PCR,
             ((cdna,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
              (p7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
             geometry=PCRGeometry.unavailable(
                 FINAL_LIBRARY,reason="The paper does not print one complete strand-resolved post-CEL-Seq2 construct."),
             source=unresolved)
 graph.sequence(library,SequencingPlan(
     limitation=SEQUENCING_UNAVAILABLE,
     acquired_roles=("Read 1",),read_lengths={"Read 1":75},
     read_uses={"Read 1":"include"}))
 graph.declare_reviewed_scope((
     "chromatin","targeted_dam","m6a_recording","dpni","dpni_digest",
     "damid2_adapter","ligase","adapter_ligation","t7_polymerase","ivt",
     "fragmentation_condition","fragmentation","rt_primer","reverse_transcriptase",
     "reverse_transcription","p5_primer","p7_primer","pcr_polymerase",
     "library_pcr","library"),source)
 return graph.require_publishable()
