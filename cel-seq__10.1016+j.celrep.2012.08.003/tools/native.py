"""Native original CEL-seq chemistry, source-bounded to the2012 originals."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import celseq as C
from chemdraw import Construct,Scene,Segment,complement_segments,revcomp,feature,oligo
from dataclasses import replace
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind,InputRole,StepKind,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
import seqprimers as sp

def rc_parts(parts):
 return [role(p.name+" complementary copy",molecular_feature=p.feature) if p.is_role_token() else replace(p,name=p.name+" complementary copy",top=revcomp(p.top)) for p in reversed(parts)]
def add_duplex(sc,key,parts,label):
 sc.strand(key,parts,label=label)
 bottom=complement_segments(parts)
 sc.anneal(key+" complement",bottom,to=key,pair=(parts[0].name+"'",parts[0].name),unpaired=tuple(p.name for p in bottom if p.is_role_token()),label=label)
def rna_start(population="endogenous cellular RNA"):
 sc=Scene();sc.strand("original RNA",[role(population+"; polyadenylated physical RNA U")],label=population+" population")
 return sc
def first_parts(): return [*rt_parts(),role("first-strand antisense transcript cDNA")]
def rna_parts():
 parts=rt_parts(); t7=next(i for i,p in enumerate(parts) if p.name=="T7 promoter")
 return [Segment("transcribed initial GGG","GGG","t7"),*[replace(p,top=p.top.replace("T","U")) for p in parts[t7+1:]],role("antisense RNA insert; physical U")]
def first_state(exo=False):
 sc=Scene();sc.strand("first cDNA",first_parts(),label="actual primer-extended single-stranded cDNA")
 sc.strand("original RNA",[role("independent original cellular/ERCC RNA; possible RNA:cDNA hybrid association")],label="original RNA retained; native fold not reconstructed")
 sc.strand("primer products",[role("ExoI hydrolysed free-primer DNA; possible unreacted primers" if exo else "possible unincorporated RT primers")],label="independent primer population; hydrolysis is not physical removal")
 return sc
def duplex_state():
 sc=Scene();add_duplex(sc,"T7 DNA template",first_parts(),"actual double-stranded cDNA including paired T7 promoter")
 sc.strand("RNA remnants",[role("original RNA/ERCC or RNA fragments; second-strand kit processing/recovery not fully assigned")],label="independent possible RNA lineage")
 sc.strand("primer remnants",[role("possible unincorporated primers/DNA hydrolysate; recovery unassigned")],label="independent possible remnants")
 return sc
def ivt_state(stage="intact",dnase=False):
 sc=Scene();sc.strand("aRNA",rna_parts(),label="new T7 antisense RNA copies; physical U; "+stage)
 if dnase:
  sc.strand("DNA products",[role("DNaseI-hydrolysed promoter-bearing DNA; possible undigested template")],label="independent DNA products; digestion is not physical discard")
 else: add_duplex(sc,"original promoter DNA",first_parts(),"original promoter-bearing DNA template; IVT does not consume it")
 sc.strand("prior RNA",[role("possible original RNA/ERCC and prior primer/RNA fragments; integrity/partition unassigned")],label="independent prior populations")
 if stage!="intact":
  sc.strand("other RNA fragments",[role("other aRNA fragments without the RT1 tag; not expected to enrich with both PCR handles")],label="nonamplification does not physically remove these fragments")
 return sc
def ligated_state():
 sc=Scene();sc.strand("tagged RNA",[*rna_parts(),*adapter_parts()],mod3=ADAPTER_MOD3,label="tag-retaining aRNA fragment joined to native adapter; RNA/DNA boundary explicit")
 sc.junction("tagged RNA","antisense RNA insert; physical U",adapter_parts()[0].name,"RNA3′–adapter5′ phosphodiester; detailed facing-end state bounded")
 sc.strand("other populations",[role("possible untagged/internal aRNA fragments, unligated adapters, original RNA/ERCC, promoter DNA/DNA fragments and free primers")],label="independent possible co-carried populations")
 return sc
def rt2_parts(): return [*rt2_head(),*rc_parts(rna_parts())]
def rt2_state():
 sc=Scene();sc.strand("second RT cDNA",rt2_parts(),label=RT2_LABEL)
 sc.strand("RNA template",[*rna_parts(),*adapter_parts()],mod3=ADAPTER_MOD3,label="independent original RNA–adapter template; not consumed by RT")
 sc.strand("other populations",[role("possible earlier original RNA/ERCC, promoter DNA/DNA hydrolysate, primers and non-tagged RNA fragments; recovery unresolved")],label="independent molecular populations")
 return sc
def amplified_state():
 sc=Scene();add_duplex(sc,"expected PCR daughter",list(FINAL_LIBRARY),"expected two-primer amplified library; bounded native precursor geometry")
 sc.strand("prior templates",[role("possible original RNA/ERCC, promoter templates/hydrolysate, RNA–adapter templates, RT2cDNA, free primers and nonamplifying fragments; partition unresolved")],label="independent possible carried original/copy populations")
 return sc

def graph_tools(g):
 def inp(i,label,kind=MaterialKind.CATALYST,state=None): return g.input(i,label,kind,state=state)
 def out(i,label,state,kind=MaterialKind.MOLECULAR): return g.product(i,label,kind,state=state)
 def op(i,label,kind,inputs,state,reason="",geometry=None):
  product=out(i+"_out",label,state)
  g.react(i,label,kind,inputs,(product,),**({"reason":reason} if reason else {}),geometry=geometry)
  return product
 def clean(i,label,x,state=None):
  selection=inp(i+"_selection",label,MaterialKind.SELECTION);product=out(i+"_out",label,state if state is not None else x.state,x.kind);waste=out(i+"_waste","physically removed cleanup flow-through/washes",None,MaterialKind.WASTE)
  g.react(i,label,StepKind.CLEANUP,((x,InputRole.TEMPLATE),(selection,InputRole.SELECTION)),(product,waste));g.discard(waste,"physically removed wash/flow-through fraction; exact molecular partition remains unassigned")
  return product
 return inp,out,op,clean
TITLE="CEL-seq — native barcoded T7 linear amplification"
NOTES="01_cel-seq.html"
SOURCE='<a href="https://doi.org/10.1016/j.celrep.2012.08.003">Hashimshony et al., Cell Reports2012</a>; original mmc4 CEL-Seq Protocol2012-07-16 and mmc5 article/S1.'
SUMMARY="Per-cell RT and second-strand synthesis precede pooling. A paired T7 promoter produces antisense RNA, whose tag-retaining fragments acquire the native RA3 adapter, a second RT copy and indexed PCR."
CAVEAT="This represents the original2012 protocol, including independent ERCC RNA. TruSeq SmallRNA RA3/RTP/RP1/RPIX are named in the detailed protocol; their ordered sequences, modifications, revision and full native end geometry remain unassigned within this capped pass. These open assignments do not establish source absence. IVT preserves its DNA template; RNaseH processing and physical column recoveries are separate. Internal RNA fragments are not physically discarded merely because they lack both PCR handles."
RT2_LABEL="actual single-stranded RTP-derived cDNA; native primer footprint unassigned"
ADAPTER_MOD3="native RA3 terminal modification unassigned"
INDEX=feature("cel_library_index","sample_index","unknown")
def rt_parts(): return C.rt_primer_segments("CEL-Seq")
def adapter_parts(): return [role("native RA3 RNA3′ adapter; sequence/modifications unassigned")]
def rt2_head(): return [role("native RTP primer-derived end; full strand unassigned")]
FINAL_LIBRARY=Construct([role("native RP1-derived end"),Segment("RA5 source handle",C.CEL1_RA5,"r1"),Segment("cell barcode","N"*8,"cbc",placeholder=True,feature=C.CELL_BARCODE),Segment("copied oligo(dT)","T"*24),Segment("source V anchor","V",placeholder=True),role("native tag insert"),role("RA3/RPIX-derived indexed end",molecular_feature=INDEX)],name="CEL-seq native library boundary")
SEQ_PRIMERS=()
FINAL_CAPTION="Source-positive barcode8 and tag insert; native outer arms/index extent unassigned."
SEQUENCING_ENDING="HiSeq2000 paired-end: Read1 at least15 bases, Read2 50; Illumina library barcode acquired when needed. Exact primer sites and selected index extent unassigned."
READ_SPECS=(sp.read_spec("Read 1","at least15","8-base cell barcode then poly(T)","extract"),sp.read_spec("Read 2",50,"transcript3′tag","include"),sp.read_spec("Index 1 (i7)","when needed; length unassigned","RPIX library identity","demultiplex"))
def oligos(): yield oligo("original CEL-seq RT primer",rt_parts())

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE);inp,out,op,clean=graph_tools(g)
 primary=g.input("mrna","single-cell poly(A) RNA",MaterialKind.BIOLOGICAL,state=rna_start(),primary=True)
 ercc=inp("ercc","ERCC poly(A) RNA spike-in; independent molecules",MaterialKind.MOLECULAR,rna_start("independent ERCC spike-in RNA"))
 primer=inp("rt_primer","2012 barcode8/T24V RT primer",MaterialKind.OLIGO,Construct(rt_parts()))
 rt=inp("array_script","MessageAmpII ArrayScript;42°C2h")
 first=op("reverse_transcription","Prime and copy endogenous/ERCC RNA in each well",StepKind.REVERSE_TRANSCRIPTION,((primary,InputRole.TEMPLATE),(ercc,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),first_state())
 mix=inp("second_mix","MessageAmpII DNA polymerase+RNaseH;16°C2h")
 ds=op("second_strand","Second-strand synthesis with RNaseH processing; actual paired T7 promoter",StepKind.EXTENSION,((first,InputRole.TEMPLATE),(mix,InputRole.CATALYST)),duplex_state())
 peers=inp("peer_wells","independently indexed ds-cDNA from other cells",MaterialKind.MOLECULAR,duplex_state())
 pooled=op("pool","Pool only after each well's second-strand synthesis",StepKind.POOL,((ds,InputRole.TEMPLATE),(peers,InputRole.TEMPLATE)),duplex_state())
 pooled=clean("cdna_column","Ambion cDNA column/washes/two warm elutions;speedvac6.4µl",pooled)
 t7=inp("t7_mix","MessageAmpII T7 polymerase+A/C/G/U rNTPs;37°C13h")
 arna=op("ivt","T7 IVT from actual ds-promoter; initialGGG and RNAU",StepKind.TRANSCRIPTION,((pooled,InputRole.TEMPLATE),(t7,InputRole.CATALYST)),ivt_state())
 frag=inp("fragment_conditions","MgOAc/KOAc/Tris-acetate94°C3min;EDTA stop",MaterialKind.CONDITION)
 arna=op("fragment_arna","Fragment RNA; retain tag-bearing and other fragment populations",StepKind.FRAGMENTATION,((arna,InputRole.TEMPLATE),(frag,InputRole.CONDITION)),ivt_state("fragmented"))
 arna=clean("arna_column","Ambion aRNA column/washes/two warm elutions",arna)
 phosphatase=inp("phosphatase","Antarctic phosphatase37°C30min;65°C5min")
 arna=op("dephosphorylate","Phosphatase treatment of fragmented RNA",StepKind.DEPHOSPHORYLATION,((arna,InputRole.TEMPLATE),(phosphatase,InputRole.CATALYST)),arna.state)
 pnk=inp("pnk","T4PNK+ATP37°C60min")
 arna=op("phosphorylate","PNK treatment of fragmented RNA",StepKind.PHOSPHORYLATION,((arna,InputRole.TEMPLATE),(pnk,InputRole.CATALYST)),arna.state)
 arna=clean("rna_rneasy","RNeasyMinElute washes/elution;speedvac5µl",arna)
 adapter=inp("ra3","TruSeqSmallRNA RA3 diluted5fold; native geometry unassigned",MaterialKind.ADAPTER)
 lig=inp("rna_ligase","T4RNA ligase2truncated28°C1h;StopSolution28°C15min")
 geom=AdapterLigationGeometry("unavailable","unavailable",(ligated_state(),),"Named native RA3 strand/end chemistry and complete RNA-facing geometry remain unassigned; known reaction order is explicit.")
 ligated=op("ligate_adapter","Ligate native RA3 to treated RNA3′ends",StepKind.LIGATION,((arna,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(lig,InputRole.CATALYST)),ligated_state(),geometry=geom)
 rtp=inp("rtp","TruSeq RNA RT primer RTP; full native sequence unassigned",MaterialKind.OLIGO)
 rt2=inp("ssii","SuperScriptII50°C1h after annealing70°C2min")
 copied=op("library_reverse_transcription","RT2 copies adapter-ligated RNA into actual ssDNA",StepKind.REVERSE_TRANSCRIPTION,((ligated,InputRole.TEMPLATE),(rtp,InputRole.STRUCTURAL),(rt2,InputRole.CATALYST)),rt2_state())
 f=inp("rp1","TruSeq RP1; ordered sequence unassigned",MaterialKind.OLIGO);r=inp("rpix","selected Truseq RPIX; member/extent unassigned",MaterialKind.OLIGO);pol=inp("pml","IlluminaPML;12cycles source example11–15range;30s extension")
 library=out("library","native CEL-seq amplified library",amplified_state(),MaterialKind.LIBRARY)
 g.react("library_pcr","Amplify actual RTP-derived ssDNA with RP1/RPIX",StepKind.PCR,((copied,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(library,),geometry=PCRGeometry.unavailable(copied.state,reason="Native ordered RP1/RPIX/RTP sites and early-cycle ssDNA geometry unassigned; no completed duplex used as input."))
 library=clean("final_cleanup1","FirstAMPure50µlbeads/50µlPCR;30µlrecovery",library)
 library=clean("final_cleanup2","SecondAMPure39µlbeads;10µlrecovery",library)
 g.sequence(library,SequencingPlan(limitation=SEQUENCING_ENDING,acquired_roles=("Read 1","Read 2"),read_lengths={"Read 2":50}))
 g.declare_reviewed_scope(("mrna","ercc","rt_primer","array_script","reverse_transcription","second_mix","second_strand","peer_wells","pool","cdna_column_selection","cdna_column","t7_mix","ivt","fragment_conditions","fragment_arna","arna_column_selection","arna_column","phosphatase","dephosphorylate","pnk","phosphorylate","rna_rneasy_selection","rna_rneasy","ra3","rna_ligase","ligate_adapter","rtp","ssii","library_reverse_transcription","rp1","rpix","pml","library_pcr","final_cleanup1_selection","final_cleanup1","final_cleanup2_selection","final_cleanup2","final_cleanup2_out"),SOURCE)
 return g.require_publishable()
