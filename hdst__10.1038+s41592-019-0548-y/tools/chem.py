"""Protocol-owned native spatial chemistry; molecular populations and source limits explicit."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from dataclasses import replace
from chemdraw import Construct,Scene,Segment,feature,oligo,revcomp,complement_segments
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind,InputRole,StepKind,SequencingPlan,Terminal,Disposition,SourceClaim,SourceDiscrepancy
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
import seqprimers as sp

def add_duplex(sc,key,parts,label):
 sc.strand(key,parts,label=label);bottom=complement_segments(parts)
 sc.anneal(key+" complement",bottom,to=key,pair=(parts[0].name+"'",parts[0].name),unpaired=tuple(p.name for p in bottom if p.is_role_token()),label=label)
def rc_parts(parts):
 return [role("copied VNanchor complement;5′N then B(C/G/T)",length_bp=2) if p.name=="VNanchor" else role(p.name+" copied complement",molecular_feature=p.feature,inferred=p.inferred) if p.is_role_token() else replace(p,name=p.name+" copied complement",top=revcomp(p.top)) for p in reversed(parts)]
def tools(g):
 def inp(i,label,kind=MaterialKind.CATALYST,state=None):return g.input(i,label,kind,state=state)
 def out(i,label,state=None,kind=MaterialKind.MOLECULAR):return g.product(i,label,kind,state=state)
 def op(i,label,kind,inputs,state,geometry=None,reason=""):
  product=out(i+"_out",label,state);g.react(i,label,kind,inputs,(product,),geometry=geometry,**({"reason":reason} if reason else {}));return product
 def clean(i,label,x,state=None):
  selection=inp(i+"_selection",label,MaterialKind.SELECTION);kept=out(i+"_out",label,x.state if state is None else state,x.kind);waste=out(i+"_waste","physically removed wash/flowthrough fraction",None,MaterialKind.WASTE)
  g.react(i,label,StepKind.CLEANUP,((x,InputRole.TEMPLATE),(selection,InputRole.SELECTION)),(kept,waste));g.discard(waste,"Source physical recovery fraction; exact molecular partition/purity unassigned.");return kept
 return inp,out,op,clean

SOURCE_OLIGO_EXCEPTIONS={"native_rt_primer":"Assembled bead-bound capture primer, not a single ordered oligo: published Methods121 prints the precursor skeleton, donor roles, completed bridges and pool/UMI/tail widths separately. Local source tests guard those components; exact selected barcode members are unassigned."}
TITLE="HDST — native bead assembly and T7 spatial transcriptomics"
NOTES="01_hdst.html"
SOURCE='<a href="https://doi.org/10.1038/s41592-019-0548-y">Vickovic et al., Nature Methods2019</a>; directly cited Jemt et al.2016workflow.'
SUMMARY="Split-pool bead oligos capture RNA at decoded spatial coordinates. Actual bead-extended cDNA is USER-released, second-strand synthesized and processed through source-known T7/RNA-adapter/secondRT operations before native indexing PCR."
CAVEAT="Published synthesis15/14/15-base barcode pools and UMI5 differ from published preprocessing coordinates1–15,31–45,61–76 and77–82 (15/15/16 and UMI6). Physical molecules follow synthesis; the computational mask stays separate and unresolved. Superseded preprint UMI7 is not current authority. The paper describes vendor barcode/decoder pools as proprietary; exact native members remain unassigned in this capped pass, without a confirmedabsence claim. Jemt2016main supplies operation order/reaction cleanups; detailed enzymes, oligos and native ends delegated to Ståhl2016SI remain unread. No canonical index width, complete endpoint shell or USER end geometry is substituted."
B12="GACTTGTCTAGAGC";B34="TGATGCCACACTACTC"
T7="TAATACGACTCACTATAGGG";HANDLE="ACACGACGCTCTTCCGATCT"
PRECURSOR_WRITTEN="/AmC6/UUUUUGACTCGTAATACGACTCACTATAGGGACACGACGCTCTTCCGATCT-Spatial_barcode_Pool1-Bridge1"
BC=[feature("hdst.bc"+str(i),"spatial_barcode","whitelist",group="hdst.spatial",part="pool"+str(i),whitelist="HDST source pool"+str(i)) for i in (1,2,3)]
UMI=feature("hdst.umi","umi","random");INDEX=feature("hdst.index","sample_index","unknown")
def identity_parts():return [Segment("barcode1","N"*15,"cbc",placeholder=True,feature=BC[0]),Segment("bridge1+2",B12),Segment("barcode2","N"*14,"cbc",placeholder=True,feature=BC[1]),Segment("bridge3+4",B34),Segment("barcode3","N"*15,"cbc",placeholder=True,feature=BC[2]),Segment("UMI","N"*5,"umi",placeholder=True,feature=UMI),Segment("dT20","T"*20),Segment("VNanchor","VN",placeholder=True)]
def capture_parts(stage=3,released=False):
 prefix=([role("native USER-released5′boundary; end chemistry unassigned")] if released else [Segment("USER-sensitive dU5","U"*5)])+[Segment("source leader","GACTCG"),Segment("T7 promoter",T7,"t7"),Segment("source Illumina handle",HANDLE,"r1")]
 ids=identity_parts()
 if stage==1:return [*prefix,ids[0],Segment("bridge1",B12[:7])]
 if stage==2:return [*prefix,*ids[:3],role("bridge3; partition within source bridge3+4 unassigned")]
 return [*prefix,*ids]
def bead_scene(stage=3):
 sc=Scene();sc.strand("bead DNA",capture_parts(stage),mod5="AmC6 covalent bead linkage",label="one capture DNA; independently randomized UMI per molecule")
 return sc
def captured_scene():
 sc=bead_scene();sc.strand("RNA",[role("independent tissue poly(A)RNA; physical U")],label="hybridized RNA population; no covalent RNA–barcode fusion")
 return sc
def cdna_parts():return [*capture_parts(released=True),role("extended antisense transcript cDNA")]
def cdna_scene(duplex=False):
 sc=Scene()
 if duplex:add_duplex(sc,"promoter DNA",cdna_parts(),"actual paired T7-promoter dsDNA before IVT")
 else:sc.strand("cDNA",cdna_parts(),label="actual released primer-extended ssDNA; no completed library arms")
 sc.strand("original RNA",[role("independent original RNA/hydrolysates; association/recovery unassigned")],label="RNA lineage retained")
 return sc
def arna_parts():return [Segment("transcribed GGG","GGG"),Segment("transcribed handle",HANDLE.replace("T","U")),*[replace(p,top=p.top.replace("T","U")) for p in identity_parts()],role("antisense transcript RNA; physical U")]
def arna_scene(ligated=False):
 sc=Scene();parts=arna_parts()+([role("native ligated RNA3′adapter; strand/end chemistry unassigned")] if ligated else [])
 sc.strand("new aRNA",parts,label="new antisense RNA copies; RNA U distinct from DNA T")
 add_duplex(sc,"DNA template",cdna_parts(),"independent promoter DNA template; IVT does not consume it")
 sc.strand("prior populations",[role("possible original RNA, USER fragments, primers and other reaction products; recovery unresolved")],label="independent possible carried species")
 return sc
def second_cdna_scene():
 sc=arna_scene(True);sc.strand("second cDNA",[role("native secondRT primer/adapter-derived end; geometry unassigned"),*rc_parts(arna_parts())],label="actual secondRT-derived ssDNA entering native PCR")
 return sc
FINAL_LIBRARY=Construct([role("native PCR/Read1 end"),*identity_parts(),role("native transcript insert"),role("native opposite PCR/index end; width unassigned",molecular_feature=INDEX)],name="HDST source-limited library")
SEQ_PRIMERS=();FINAL_CAPTION="Physical15/14/15barcodes+UMI5; native outer ends/index scheme unassigned."
SEQUENCING_ENDING="Represented MOB run R1=125/R2=175; independent breast cancer cohort R1=150/R2=150. Index acquisition unassigned. The published computational R1 mask conflicts with physical synthesis lengths."
READ_SPECS=(sp.read_spec("Read 1","125MOB;150breast","reported barcode positions1–15,31–45,61–76; UMI77–82; unresolved versus synthesis","extract"),sp.read_spec("Read 2","175MOB;150breast","transcript","include"),sp.read_spec("Index 1 (i7)",None,"native index member/extent/acquisition unassigned","undetermined"),sp.read_spec("Index 2 (i5)",None,"native index scheme/acquisition unassigned","undetermined"))
def oligos():return ("Precursor as source-written: "+PRECURSOR_WRITTEN,"Ligated bridge1+2: "+B12,"Ligated bridge3+4: "+B34)
def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE);inp,out,op,clean=tools(g)
 bead=g.input("bead","AmC6-linked precursor-bearing silica beads",MaterialKind.BIOLOGICAL,state=bead_scene(1),primary=True)
 donor2=inp("donor2","5′phosphorylated Bridge2-pool2-Bridge3",MaterialKind.ADAPTER)
 helper2=inp("helper2","complementary14bp Bridge1′Bridge2′ splint; bases derived from B12",MaterialKind.OLIGO,revcomp(B12))
 lig=inp("ligase","T4DNA ligase")
 geo=AdapterLigationGeometry("unavailable","unavailable",(bead_scene(1),),"Source5′phosphorylated donors/complementary splints and completed bridges known; full ligation-facing end geometry/selected barcode sequences unassigned.")
 bead=op("ligate_round2","Split/pool and splint-ligate barcode2",StepKind.LIGATION,((bead,InputRole.TEMPLATE),(donor2,InputRole.STRUCTURAL),(helper2,InputRole.STRUCTURAL),(lig,InputRole.CATALYST)),bead_scene(2),geometry=geo)
 donor3=inp("donor3","5′phosphorylated Bridge4-pool3-UMI5-dT20VN",MaterialKind.ADAPTER);helper3=inp("helper3","complementary Bridge3′Bridge4′ helper; split unassigned",MaterialKind.OLIGO,revcomp(B34))
 bead=op("ligate_round3","Split/pool and splint-ligate barcode3/capture tail",StepKind.LIGATION,((bead,InputRole.TEMPLATE),(donor3,InputRole.STRUCTURAL),(helper3,InputRole.STRUCTURAL),(lig,InputRole.CATALYST)),bead_scene(),geometry=geo)
 decode=inp("decoders","14cycles complementaryFAM/Cy3/dark decoder pools and stripping",MaterialKind.OLIGO)
 bead=op("decode","Load array and decode native bead identities; strip decoder pools",StepKind.HYBRIDIZATION,((bead,InputRole.TEMPLATE),(decode,InputRole.STRUCTURAL)),bead_scene())
 rna=inp("rna","independent tissue poly(A)RNA after fixation/H&E/permeabilization",MaterialKind.MOLECULAR)
 captured=op("capture","Capture tissue RNA by poly(dT) hybridization",StepKind.HYBRIDIZATION,((bead,InputRole.TEMPLATE),(rna,InputRole.STRUCTURAL)),captured_scene())
 primer=inp("native_rt_primer","actual assembled bead capture primer; same molecular species",MaterialKind.OLIGO,Construct(capture_parts()))
 rt=inp("ssiii","SuperscriptIII+actinomycinD/RNaseOUT;42°Covernight")
 bound=bead_scene();bound.strands["bead DNA"].segs.append(role("extended antisense transcript cDNA"))
 bound.strand("RNA",[role("independent original RNA; physical U")],label="retained RNA template")
 cdna=op("rt","Extend native bead primer with first-strand DNA",StepKind.REVERSE_TRANSCRIPTION,((captured,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),bound)
 pk=inp("proteinase","ProteinaseK56°C1h")
 cdna=op("digest_tissue","Digest tissue; bead/cDNA not physically released by proteolysis",StepKind.DEPROTEINIZATION,((cdna,InputRole.TEMPLATE),(pk,InputRole.CATALYST)),cdna.state)
 user=inp("user","USER37°C2h")
 released=out("released","actual USER-released cDNA; exact5′boundary unassigned",cdna_scene());remnant=out("bead_remnants","bead-linked dU cleavage fragments; recovery unassigned",None)
 g.react("release","Cleave dU5release tract",StepKind.DIGESTION,((cdna,InputRole.TEMPLATE),(user,InputRole.CATALYST)),(released,remnant));g.mark_terminal(Terminal(remnant.id,Disposition.UNRESOLVED,"Chemical cleavage is not an asserted physical discard; downstream bead-remnant recovery unassigned."))
 second=inp("second_mix","Jemt second-strand synthesis reagents; exact enzymes deferred")
 ds=op("second_strand","Synthesize actual paired-promoter dsDNA",StepKind.EXTENSION,((released,InputRole.TEMPLATE),(second,InputRole.CATALYST)),cdna_scene(True));ds=clean("second_cleanup","Jemt reaction cleanup after second strand",ds)
 repair=inp("end_repair_mix","Jemt end-repair reagents; exact end chemistry deferred")
 ds=op("end_repair","End-repair promoter-bearing DNA",StepKind.END_REPAIR,((ds,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),ds.state);ds=clean("repair_cleanup","Jemt reaction cleanup after end repair",ds)
 t7=inp("t7","T7IVT reagents with RNA nucleotides; detailed conditions deferred")
 arna=op("ivt","Transcribe new antisense RNA from paired T7promoter",StepKind.TRANSCRIPTION,((ds,InputRole.TEMPLATE),(t7,InputRole.CATALYST)),arna_scene());arna=clean("ivt_cleanup","Jemt reaction cleanup after IVT",arna)
 adapter=inp("rna_adapter","native RNA3′adapter; ordered strand/end chemistry deferred",MaterialKind.ADAPTER);rl=inp("rna_ligase","Jemt RNA-adapter ligation mix; exact enzyme deferred")
 ageo=AdapterLigationGeometry("unavailable","unavailable",(arna_scene(),),"Native RNA/adapter ligation geometry deferred to uninspected Ståhl2016SI; actual RNA input retained.")
 arna=op("rna_ligate","Ligate native RNA3′adapter",StepKind.LIGATION,((arna,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(rl,InputRole.CATALYST)),arna_scene(True),geometry=ageo);arna=clean("rna_ligation_cleanup","Jemt reaction cleanup after RNAadapter ligation",arna)
 rtp=inp("second_rt_primer","native adapter-complementary secondRT primer; unassigned",MaterialKind.OLIGO);rt2=inp("second_rt","Jemt second cDNA-synthesis enzyme; detailed identity deferred")
 copied=op("second_rt_copy","Copy adapter-ligated RNA into actual ssDNA",StepKind.REVERSE_TRANSCRIPTION,((arna,InputRole.TEMPLATE),(rtp,InputRole.STRUCTURAL),(rt2,InputRole.CATALYST)),second_cdna_scene());copied=clean("second_rt_cleanup","Jemt reaction cleanup after second cDNA synthesis",copied)
 f=inp("forward","native library forward primer; unassigned",MaterialKind.OLIGO);r=inp("reverse","native library indexed opposite primer; member/width unassigned",MaterialKind.OLIGO);pol=inp("pcr_mix","native indexing PCR mix; detailed chemistry deferred")
 final=second_cdna_scene();add_duplex(final,"expected library",list(FINAL_LIBRARY),"native expected library plus possible carried originals/copies")
 library=out("library","HDST native amplified library",final,MaterialKind.LIBRARY)
 g.react("pcr","Amplify actual secondRTssDNA; native boundaries unassigned",StepKind.PCR,((copied,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(library,),geometry=PCRGeometry.unavailable(copied.state,reason="Actual recovered ssDNA supplied; native ordered primers and first-cycle geometry remain deferred, not a finished-duplex input."))
 g.sequence(library,SequencingPlan(acquired_roles=("Read 1","Read 2"),read_lengths={"Read 1":125,"Read 2":175},limitation=SEQUENCING_ENDING))
 g.declare_reviewed_scope(("bead","donor2","helper2","ligase","ligate_round2","donor3","helper3","ligate_round3","decoders","decode","rna","capture","native_rt_primer","ssiii","rt","proteinase","digest_tissue","user","release","second_mix","second_strand","second_cleanup_selection","second_cleanup","end_repair_mix","end_repair","repair_cleanup_selection","repair_cleanup","t7","ivt","ivt_cleanup_selection","ivt_cleanup","rna_adapter","rna_ligase","rna_ligate","rna_ligation_cleanup_selection","rna_ligation_cleanup","second_rt_primer","second_rt","second_rt_copy","second_rt_cleanup_selection","second_rt_cleanup","forward","reverse","pcr_mix","pcr","library"),SOURCE)
 return g.require_publishable()
