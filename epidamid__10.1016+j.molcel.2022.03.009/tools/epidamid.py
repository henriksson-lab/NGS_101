"""EpiDamID / scDam&T-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct,MolecularState,Row,Scene,Segment,Workflow,annotation_rows,oligo
from adapter_ligation import AdapterLigationGeometry
from damid import DamMarkedGATC,damid_adapter
from single_cell_hic import unresolved_illumina_library

ADRT="CTAATACGACTCACTATAGGGCAGCGTGGTCGCGGCCGAGGA"
ADRB="TCCTCGGCCGCG"
ADR_PCR="GGTCGCGGCCGAGGATC"
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
CAVEAT="INFERRED — the source prints the DamID2 adapter and terminal PCR primers, but not one complete strand-resolved post-CEL-Seq2 library. The final outer arms are therefore unresolved rather than silently assigned to a TruSeq generation."
FINAL_CAPTION="EpiDamID insert after the scDam&T/CEL-Seq2 branch. Exact DamID adapter sequence is retained; unresolved outer sequencing regions are visibly inferred."
SEQUENCING_ENDING="Libraries were sequenced on NextSeq 500 (75-base reads) or NextSeq 2000 (100-base reads)."
SEQUENCING_UNAVAILABLE="A complete final sequencing-primer binding construct cannot be verified from the paper's split DamID2 and CEL-Seq2 descriptions; the source-transcribed primers remain available below."
def oligos():
 yield oligo("AdRt",[Segment("adapter top",ADRT,"t7")])
 yield oligo("AdRb",[Segment("adapter bottom",ADRB,"t7")])
 yield oligo("AdR_PCR",[Segment("PCR primer",ADR_PCR,"r1")])
