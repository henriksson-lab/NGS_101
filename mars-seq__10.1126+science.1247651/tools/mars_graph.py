"""Source-owned MARS-seq native reactions; polymer populations/fates explicit."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
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
IS_MARS2=False
PROTOCOL="MARS-seq"
WELL_NT,UMI_NT,RANDOM_NT=6,4,3
NOTES="01_mars-seq.html"
SOURCE="https://doi.org/10.1126/science.1247651"
RT2="TCTAGCCTTCTCGCAGCACATC"
READ_LIMIT="Original2014SupplementMethods:HiSeq2500 paired-end; cycle numbers/index acquisition unassigned. Canonical TruSeq sites provide only computed endpoint compatibility; loaded sequencing oligos remain unverified."
CAVEAT="INFERRED—Original SupplementTableS7 prints RT2 TCTAGCCTTCTCGCAGCACATC, whereas Methods calls it complementary to adapterAGATCGGAAGAGCGTCGTGTAG. Its reversal CTACACGACGCTCTTCCGATCT equals the computed adapter reverse complement; possible3′→5′listing is unresolved. The printed oligo is retained without choosing its landing; the expected Read1-side PCR daughter remains inferred. Later2019RT2 does not settle the2014order. Original2014fullMethods coauthor upload is independently distinguished from the two-page table extract. ExoI follows24wellintermediatepooling, then1.2×SPRI/192cellbatches; no2019reagent/order inheritance."
# Original2014TableS7 splits each printed PCR primer across sequence/role-label lines.
SOURCE_OLIGO_EXCEPTIONS={"pcr_forward":"Jaitin2014TableS7: role label interrupts two printed sequence lines; assembled directly from the original table.","pcr_reverse":"Jaitin2014TableS7: role label interrupts two printed sequence lines; assembled directly from original table."}

TITLE=PROTOCOL+" — native three-level indexed T7 amplification"
SUMMARY="Actual RT1-extended cDNA becomes a paired T7-promoter template. T7 makes new RNA copies, DNA-template digestion precedes physical recovery, and fragmented RNA acquires a phosphorylated/C3-blocked pool adapter before RT2 and PCR."
ADAPTER_MOD3="C3 spacer; no3′OH"
INDEX=feature("pool_barcode","sample_index","fixed",note="Representative ix1 source adapter; inline4-base identity, not an instrument index read")
CELL=feature("cell_barcode_well","cell_barcode","whitelist",whitelist=PROTOCOL+" RT1 set")
UMI=feature("umi","umi","random")
T7="TAATACGACTCACTATAGGG"
SPACER="CGATTGAGGCCGG"
HANDLE="GACGTGTGCTCTTCCGATCT"
LIG="AGATCGGAAGAGCGTCGTGTAG"
POOL="GACT"
def rt_parts():
 return [Segment("5′spacer",SPACER),Segment("T7 promoter",T7,"t7"),Segment("transcribed leader","GC"),Segment("RT1/PCR handle",HANDLE,"r2"),Segment("well barcode","N"*WELL_NT,"cbc",placeholder=True,feature=CELL),Segment("UMI","N"*UMI_NT,"umi",placeholder=True,feature=UMI),Segment("oligo(dT)","T"*20),Segment("N anchor","N",placeholder=True)]
def adapter_parts(): return [Segment("ordered pool barcode",POOL,"cbc",feature=INDEX),Segment("adapter random bases","N"*RANDOM_NT,placeholder=True),Segment("ligation/RT2 handle",LIG,"r1")]
def adapter_state():
 sc=Scene();sc.strand("adapter",adapter_parts(),mod5="p",mod3=ADAPTER_MOD3,label="DNA donor; source5′phosphate and3′C3block")
 return sc
def rt2_head():
 if IS_MARS2: return [Segment("native RT2",RT2,"r1"),*rc_parts(adapter_parts()[:2])]
 return [role("native printed-RT2-derived end; orientation/landing unresolved",inferred=True),*rc_parts(adapter_parts()[:2])]
RT2_LABEL="actual RT2-extended single-stranded DNA; no completed P5/P7 arms" if IS_MARS2 else "INFERRED—actual RT2copy class; printed2014 RT2 orientation/landing unresolved"
from illumina import P5,P7_RC,TRUSEQ_READ1,TRUSEQ_READ2,TRUSEQ_P5_FULL,P7
P5_PCR=TRUSEQ_P5_FULL
P7_PCR=P7+TRUSEQ_READ2
FINAL_LIBRARY=Construct([Segment("P5",P5,"p5"),Segment("TruSeqRead1 nonoverlap",TRUSEQ_READ1[4:],"r1",inferred=not IS_MARS2),*rc_parts(adapter_parts()[:2]),role("sense transcript3′tag"),Segment("copied N anchor","N",placeholder=True),Segment("copied oligo(dT)","A"*20),Segment("complementary UMI","N"*UMI_NT,"umi",placeholder=True,feature=UMI),Segment("complementary well barcode","N"*WELL_NT,"cbc",placeholder=True,feature=CELL),Segment("TruSeqRead2 complement",revcomp(TRUSEQ_READ2),"r2"),Segment("P7 complement",P7_RC,"p7")],name=PROTOCOL+" expected native PCR daughter")
# Reverse adapter order in rc_parts is random bases then pool identity, derived from its physical donor.
SEQ_PRIMERS=sp.TRUSEQ["R1"],sp.TRUSEQ["R2"]
FINAL_CAPTION="Expected native PCR daughter from source primer flanks; actual PCR input remains RT2-derived single-stranded DNA."
SEQUENCING_ENDING=READ_LIMIT
SEQUENCING_INTRO=READ_LIMIT
READ_SPECS=(sp.read_spec("Read 1",75 if IS_MARS2 else "unassigned",str(RANDOM_NT)+"adapter random bases+4-base pool barcode then transcript","extract_then_include"),sp.read_spec("Read 2",15 if IS_MARS2 else "unassigned",str(WELL_NT)+"well barcode+"+str(UMI_NT)+"UMI","extract"))
if IS_MARS2:
 READ_SPECS+= (sp.read_spec("Index 1 (i7)","not acquired","source procedure explicitly no indices","not_applicable"),sp.read_spec("Index 2 (i5)","not acquired","source procedure explicitly no indices","not_applicable"))
else:
 READ_SPECS+= (sp.read_spec("Index 1 (i7)","unassigned","acquisition not assigned; inlinepool does not establish instrument cycles","undetermined"),sp.read_spec("Index 2 (i5)","unassigned","acquisition not assigned","undetermined"))
def oligos():
 yield oligo("native RT1",rt_parts())
 yield oligo("representative ix1 pool adapter",adapter_parts(),mods="/5Phos/ … /3SpC3/")
 yield oligo("source-written RT2",[Segment("RT2",RT2)])
 yield oligo("source P5/Read1 PCR",[Segment("P5/Read1",P5_PCR)])
 yield oligo("source P7/Read2 PCR",[Segment("P7/Read2",P7_PCR)])

def protocol_graph():
 g=ProtocolGraph(PROTOCOL,authority=SOURCE);inp,out,op,clean=graph_tools(g)
 primary=g.input("mrna","single-cell poly(A) RNA in capture well",MaterialKind.BIOLOGICAL,state=rna_start(),primary=True)
 ercc=inp("ercc","independent ERCC poly(A) RNA spike-in",MaterialKind.MOLECULAR,rna_start("independent ERCC spike-in RNA"))
 primer=inp("rt_primer","source-native RT1",MaterialKind.OLIGO,Construct(rt_parts()))
 rt=inp("superscript_iii","SuperScriptIII;42°C2min/50°C50min/85°C5min")
 first=op("reverse_transcription","RT1 extension in each cell well",StepKind.REVERSE_TRANSCRIPTION,((primary,InputRole.TEMPLATE),(ercc,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),first_state())
 if not IS_MARS2:
  peers=inp("peer_wells","other independently barcoded RT wells",MaterialKind.MOLECULAR,first_state())
  first=op("pool_initial","Pool first-strand reactions into24intermediate wells",StepKind.POOL,((first,InputRole.TEMPLATE),(peers,InputRole.TEMPLATE)),first_state())
 exo=inp("exonuclease_i","ExonucleaseI37°C30min;80°C10min inactivation")
 first=op("remove_free_primers","Hydrolyse unincorporated RT1; fragments remain in mixture",StepKind.DIGESTION,((first,InputRole.TEMPLATE),(exo,InputRole.CATALYST)),first_state(True))
 peers=inp("peer_pools","other ExoI-treated wells/pools",MaterialKind.MOLECULAR,first_state(True))
 first=op("pool_final","Pool384well plate via12then6wells" if IS_MARS2 else "Pool192cell amplification batches during bead recovery",StepKind.POOL,((first,InputRole.TEMPLATE),(peers,InputRole.TEMPLATE)),first_state(True))
 first=clean("first_cleanup","0.9×AMPure recovery19µl" if IS_MARS2 else "1.2×SPRI recovery17µl",first)
 second=inp("second_strand_mix","NEB second-strand kit;16°C2h" if IS_MARS2 else "NEB second-strand kit;16°C2.5h")
 ds=op("second_strand","Convert actual first cDNA to dsDNA including paired T7 promoter",StepKind.EXTENSION,((first,InputRole.TEMPLATE),(second,InputRole.CATALYST)),duplex_state())
 ds=clean("ds_cleanup","1.4×SPRI recovery8µl; IVT on retained beads",ds)
 t7=inp("t7_mix","T7HighYield polymerase+A/C/G/U rNTPs;37°Covernight")
 arna=op("ivt","T7 transcription makes new initialGGG/GC-bearing RNA; DNA retained",StepKind.TRANSCRIPTION,((ds,InputRole.TEMPLATE),(t7,InputRole.CATALYST)),ivt_state())
 dnase=inp("turbo_dnase_i","TurboDNaseI37°C15min")
 arna=op("digest_dna","Hydrolyse promoter DNA template; no physical discard at digestion",StepKind.DIGESTION,((arna,InputRole.TEMPLATE),(dnase,InputRole.CATALYST)),ivt_state(dnase=True))
 arna=clean("rna_cleanup","1.2×SPRI aRNA recovery",arna)
 fragment=inp("fragment_conditions","Zn2+fragmentation70°C3min;stop solution",MaterialKind.CONDITION)
 arna=op("fragment_arna","Fragment RNA; tag-retaining and internal populations remain",StepKind.FRAGMENTATION,((arna,InputRole.TEMPLATE),(fragment,InputRole.CONDITION)),ivt_state("fragmented",dnase=True))
 arna=clean("fragment_cleanup","2×SPRI fragmented RNA recovery",arna)
 adapter=inp("pool_adapter","source5′phosphorylated/3′C3blocked DNA adapter",MaterialKind.ADAPTER,adapter_state())
 ligase=inp("rna_ligase","T4RNA ligase1+DMSO/ATP/PEG;22°C2h after70°C3min annealing")
 geom=AdapterLigationGeometry("unavailable","unavailable",(adapter_state(),),"Source-positive DNA donor5′phosphate/3′C3block; complete facing fragmented-RNA end geometry remains unassigned, without inferring a ligation-ready3′OH population.")
 arna=op("ligate_pool_adapter","Ligate native pool-adapter DNA to fragmented RNA",StepKind.LIGATION,((arna,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),ligated_state(),geometry=geom)
 if IS_MARS2: arna=clean("ligation_cleanup","1.5×SPRI ligation-product recovery",arna)
 rt2=inp("rt2","source-written second RT primer",MaterialKind.OLIGO,RT2)
 rt2mix=inp("affinityscript","AffinityScript42°C2min/45°C50min/85°C5min" if IS_MARS2 else "AffinityScript42°C2min/50°C45min/85°C5min")
 copied=op("second_reverse_transcription",RT2_LABEL,StepKind.REVERSE_TRANSCRIPTION,((arna,InputRole.TEMPLATE),(rt2,InputRole.STRUCTURAL),(rt2mix,InputRole.CATALYST)),rt2_state())
 copied=clean("rt2_cleanup","1.5×SPRI RT2cDNA recovery",copied)
 f=inp("pcr_forward","source P5/Read1 enrichment oligo",MaterialKind.OLIGO,P5_PCR);r=inp("pcr_reverse","source P7/Read2 enrichment oligo",MaterialKind.OLIGO,P7_PCR)
 pol=inp("pcr_mix","KAPAHiFi;example14cycles2at57°C+12at65°C;calibrate11–17" if IS_MARS2 else "KAPA PCR ready mix; source calls nested PCR, exact cycle count unassigned")
 library=out("library","native amplified "+PROTOCOL+" library",amplified_state(),MaterialKind.LIBRARY)
 g.react("library_pcr","Amplify actual recovered RT2ssDNA; full first-cycle geometry bounded",StepKind.PCR,((copied,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(library,),geometry=PCRGeometry.unavailable(copied.state,reason="Actual single-stranded RT2copy provided. Full sequential first-cycle geometry is not reconstructed;2014 printed RT2 orientation remains unresolved." if not IS_MARS2 else "Actual RT2single-stranded cDNA supplied; separate full-duplex PCR site view would not represent this input's first-cycle state."))
 library=clean("final_cleanup","0.7×SPRI library recovery",library)
 g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=SEQ_PRIMERS,acquired_roles=("Read 1","Read 2"),read_lengths={"Read 1":75,"Read 2":15} if IS_MARS2 else {},limitation=READ_LIMIT))
 scope=["mrna","ercc","rt_primer","superscript_iii","reverse_transcription","exonuclease_i","remove_free_primers","peer_pools","pool_final","first_cleanup_selection","first_cleanup","second_strand_mix","second_strand","ds_cleanup_selection","ds_cleanup","t7_mix","ivt","turbo_dnase_i","digest_dna","rna_cleanup_selection","rna_cleanup","fragment_conditions","fragment_arna","fragment_cleanup_selection","fragment_cleanup","pool_adapter","rna_ligase","ligate_pool_adapter","rt2","affinityscript","second_reverse_transcription","rt2_cleanup_selection","rt2_cleanup","pcr_forward","pcr_reverse","pcr_mix","library_pcr","final_cleanup_selection","final_cleanup","final_cleanup_out"]
 if IS_MARS2: scope.extend(("ligation_cleanup_selection","ligation_cleanup"))
 else: scope.extend(("peer_wells","pool_initial"))
 g.declare_reviewed_scope(scope,SOURCE)
 return g.require_publishable()
