"""Hou2016: protocol-owned cytoplasmic/Tang and nuclear/scRRBS native routes."""
from chemdraw import Construct, Scene, Segment, complement_segments, revcomp
from protocol_model import role
from adapter_ligation import AdapterLigationGeometry
from endprep import EndPreparedDuplex, dA_tailed_scene
from pcr import PCRGeometry
from protocol_graph import Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
import seqprimers as sp
TITLE="scTrio-seq — genome copy number, methylome and transcriptome"
NOTES="01_sctrio-seq.html"
SOURCE='<a href="https://doi.org/10.1038/cr.2016.23">Hou etal., CellResearch2016</a>, Methods152–166; explicitly cited Guo2013 scRRBS Methods105 and Tang component scope named locally.'
SUMMARY="Soft lysis and centrifugation separate4µl cytoplasm/ERCC RNA from3µl nuclear fraction. Tang-lineage cDNA is size-selected, Covaris-sheared and processed with NEBNextUltraE7370; nuclear/lambda DNA enters native scRRBS and supplies both CNV and methylation reads."
CAVEAT="No canonical TruSeq shell, fixed index width or borrowed cell barcode is supplied. Native NEBNextE7370 oligo/strand assembly and premethylated scRRBS adapters remain unassigned. Tang2009 UP1/UP2 and5′amine variants are exact component references; the defining2016 paper cites the2010 Tang method, whose detailed order is retained from existing inherited scope, not newly reread in this three-original pass. Historical SOLiD adapters are not substituted for E7370."
UP1_HANDLE="ATATGGATCCGGCGCGCCGTCGAC"
UP2_HANDLE="ATATCTCGAGGGCGCGCCGGATCC"
UP1=UP1_HANDLE+"T"*24
UP2=UP2_HANDLE+"T"*24
RNA=Construct([role("nativeNEBNextUltraE7370 end A"),role("cytoplasmic/ERCC cDNA-derived fragment"),role("nativeE7370 end B; index family/member/width unassigned")],name="scTrio transcriptome native endpoint")
DNA=Construct([role("native premethylated adapter-derived end A"),role("nuclear or lambda DNA insert; converted-site PCR readout T/A"),role("native adapter-derived end B; index member/width unassigned")],name="scTrio genome/methylome native endpoint")
RP=DP=()
READ_LENGTHS={"Read 1":100,"Read 2":100}
FINAL_LIBRARIES=(("Transcriptome library",RNA,(),"NativeE7370 endpoint roles; Tang/SOLiD sequencing ends not inherited.","Published HiSeq2000/2500 PE100."),("Genome and methylome library",DNA,(),"One physical scRRBS library supplies CNV and methylation observations.","Published HiSeq2000/2500 PE100; index acquisition unassigned."))
SHOW_READ_CONTENT=False
READ_SPECS={name:(sp.read_spec("Read 1",100,"RNA-derived cDNA" if name=="Transcriptome library" else "nuclear DNA methylation/CNV insert","include"),sp.read_spec("Read 2",100,"opposite insert end","include"),sp.read_spec("Index 1 (i7)","unassigned","native index acquisition/width unassigned","undetermined"),sp.read_spec("Index 2 (i5)","unassigned","native index acquisition/width unassigned","undetermined")) for name,*_ in FINAL_LIBRARIES}

def populations(*labels):
 sc=Scene()
 for n,label in enumerate(labels):sc.strand(str(n),[role(label)],label=label)
 return sc

def rt_state(tailed=False,second=False,amine=False):
 sc=Scene()
 for name in ("cytoplasmic RNA","ERCC RNA"):
  sc.strand(name,[role(name)],label="independent RNA template; recovery/degradation unassigned")
  parts=[Segment("UP1",UP1),role("first-strand copy of "+name)]
  if tailed:parts.append(role("variable3′poly(dA) tail"))
  sc.strand(name+" first",parts,label="first-copy strand",mod5="amine" if amine else "")
  if second:sc.strand(name+" second",[Segment("UP2",UP2),role("UP2-initiated opposite copy; terminal extent unassigned")],label="separate complementary-copy role",mod5="amine" if amine else "")
 return sc

def oligo_state(name,seq,amine=False):
 sc=Scene();sc.strand(name,[Segment(name,seq)],mod5="NH2" if amine else "",label=name);return sc

def dna_population(stage,converted=False,carrier=False):
 sc=Scene()
 for name in ("nuclear genomic DNA","unmethylated lambda spike-in"):
  for n in range(2 if converted else 1):sc.strand(name+str(n),[role(name+": "+stage)],label="physical converted strand; U before PCR" if converted else name)
 if stage.startswith("PCR"):
  for name in ("nuclear genomic DNA","unmethylated lambda spike-in"):
   sc.strand(name+" carried originals",[role("carried converted U-bearing original "+name+"; native survival/recovery unassigned")],label="original templates are not assumed chemically erased by PCR")
 if carrier:sc.strand("tRNA",[role("10ng tRNA carrier; source added during conversion-column recovery")],label="independent RNA carrier; later fate unassigned")
 return sc

def restriction_state():
 rows=[]
 for name in ("nuclear DNA","lambda DNA"):
  core=[role(name+" paired MspI fragment body")]
  sc=Scene();sc.strand("top",[Segment("5′CG overhang","CG"),*core],label=name,mod5="p")
  sc.anneal("bottom",[Segment("opposite5′CG overhang","CG"),*complement_segments(core)],"top",pair=(core[0].name+"'",core[0].name),label=name,mod5="p",unpaired=(core[0].name+"'",))
  rows.extend(sc.rows())
 return rows

def prepared_state():
 rows=[]
 for name in ("nuclear DNA","lambda DNA"):
  prepared=EndPreparedDuplex(Construct([role(name+" filled fragment body")],name=name),five_prime_phosphate=True,three_prime_overhang="A",overhang_inferred=False)
  rows.extend(dA_tailed_scene(prepared,label=name+" opposed3′A ends").rows())
 return rows

def cleanup(g,key,template,label):
 select=g.input(key+"_selection",label,MaterialKind.SELECTION)
 kept=g.product(key+"_retained","physically recovered "+template.label,MaterialKind.MOLECULAR,state=template.state)
 waste=g.product(key+"_waste","removed physical fraction; molecular composition unassigned",MaterialKind.WASTE)
 g.react(key,label,StepKind.CLEANUP,((template,InputRole.TEMPLATE),(select,InputRole.SELECTION)),(kept,waste),changes_molecule=False)
 g.discard(waste,"Physical non-retained fraction; species-specific recovery not inferred.")
 return kept

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cell=g.input("cell","one intact cell",MaterialKind.BIOLOGICAL,state=populations("RNA","nuclear genomic DNA"),primary=True)
 lysis=g.input("soft_lysis","7µl soft buffer;4°C30min, vortex1min,1000g5min",MaterialKind.CONDITION)
 cyto=g.product("cytoplasm","upper4µl cytoplasm transferred to RT tube",MaterialKind.MOLECULAR,state=populations("cytoplasmic poly(A) RNA"))
 nuc=g.product("nucleus","remaining3µl containing intact nucleus",MaterialKind.BIOLOGICAL,state=populations("nuclear genomic DNA"))
 g.react("separate","Soft-lyse, centrifuge and transfer upper cytoplasm",StepKind.SPLIT,((cell,InputRole.TEMPLATE),(lysis,InputRole.CONDITION)),(cyto,nuc),changes_molecule=False)
 ercc=g.input("ercc","ERCC spike-in RNA in RT tube",MaterialKind.MOLECULAR)
 p1=g.input("up1","Tang2009 UP1 component reference",MaterialKind.OLIGO,state=oligo_state("UP1",UP1))
 rt=g.input("rt","Tang-lineage RT; selected enzyme unassigned",MaterialKind.CATALYST)
 first=g.product("first","UP1-primed ss-cDNA/RNA mixture",MaterialKind.MOLECULAR,state=rt_state())
 g.react("rt_copy","Reverse-transcribe cytoplasmic RNA and ERCC",StepKind.REVERSE_TRANSCRIPTION,((cyto,InputRole.TEMPLATE),(ercc,InputRole.SPIKE_IN),(p1,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,))
 exo=g.input("exo","inherited Tang ExonucleaseI",MaterialKind.CATALYST)
 dig=g.product("digested","ExoI-treated RT mixture",MaterialKind.MOLECULAR,state=first.state)
 debris=g.product("primer_hydrolysate","UP1 digestion products; physical removal unassigned",MaterialKind.MOLECULAR)
 g.react("digest_up1","Digest unused UP1",StepKind.DIGESTION,((first,InputRole.TEMPLATE),(exo,InputRole.CATALYST)),(dig,debris))
 g.mark_terminal(Terminal(debris.id,Disposition.UNRESOLVED,"Chemical hydrolysis is not physical removal; later column/gel recovery unassigned."))
 tdt=g.input("tdt","inherited Tang TdT+dATP",MaterialKind.CATALYST)
 tailed=g.product("tailed","variable3′dA-tailed ss-cDNA mixture",MaterialKind.MOLECULAR,state=rt_state(tailed=True))
 g.react("tail","Poly(dA)-tail first copies",StepKind.TAILING,((dig,InputRole.TEMPLATE),(tdt,InputRole.CATALYST)),(tailed,))
 p2=g.input("up2","Tang2009 UP2 component reference",MaterialKind.OLIGO,state=oligo_state("UP2",UP2))
 copying=g.input("copy_polymerase","inherited Tang second-strand polymerase; identity unassigned",MaterialKind.CATALYST)
 copies=g.product("copies","UP2-initiated complementary copies",MaterialKind.MOLECULAR,state=rt_state(tailed=True,second=True))
 g.react("second_copy","Copy tailed first strand with UP2",StepKind.EXTENSION,((tailed,InputRole.TEMPLATE),(p2,InputRole.STRUCTURAL),(copying,InputRole.CATALYST)),(copies,))
 pol=g.input("cdna_polymerase","Tang cDNA PCR polymerase; selected identity unassigned",MaterialKind.CATALYST)
 amp=g.product("cdna_amplified","UP1/UP2 amplified cDNA population",MaterialKind.MOLECULAR,state=copies.state)
 g.react("cdna_pcr","Amplify copied cDNA with UP1/UP2",StepKind.PCR,((copies,InputRole.TEMPLATE),(p1,InputRole.STRUCTURAL),(p2,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(copies.state,reason="Actual UP2 copy mixture; native first-cycle/end derivation not supplied by an expected completed amplicon."))
 a1=g.input("aup1","Tang2009 NH2-UP1 reference; applicability inherited not newly resolved",MaterialKind.OLIGO,state=oligo_state("NH2-UP1",UP1,True))
 a2=g.input("aup2","Tang2009 NH2-UP2 reference",MaterialKind.OLIGO,state=oligo_state("NH2-UP2",UP2,True))
 blocked=g.product("blocked","amine-primer-derived cDNA reference population",MaterialKind.MOLECULAR,state=rt_state(tailed=True,second=True,amine=True))
 retained_unmodified=g.product("unmodified_cdna_precursors","unmodified UP1/UP2 precursors co-carried through amine-primer PCR",MaterialKind.MOLECULAR,state=amp.state)
 g.react("amine_pcr","INFERRED inherited Tang amine-primer amplification stage",StepKind.PCR,((amp,InputRole.TEMPLATE),(a1,InputRole.STRUCTURAL),(a2,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(blocked,retained_unmodified),geometry=PCRGeometry.unavailable(amp.state,reason="Existing inherited Tang stage; table verifies primer bases/5′NH2 but2010 detailed method and2016 end retention remain deferred."))
 g.mark_terminal(Terminal(retained_unmodified.id,Disposition.UNRESOLVED,"Amine-primer PCR creates modified copies without chemically replacing all unmodified originals. Their recovery through the column/gel/native-kit route is unassigned; no physical isolation is implied."))
 selected=cleanup(g,"cdna_column",blocked,"DNA Clean&Concentrator5 column")
 selected=cleanup(g,"cdna_gel",selected,"2%agarose500–3000bp selection;VisClean recovery; source removes amplification primers")
 shear=g.input("covaris","CovarisS2 fragmentation to150–250bp",MaterialKind.CONDITION)
 frag=g.product("rna_fragments","150–250bp cDNA fragments; terminal amine classes unassigned",MaterialKind.MOLECULAR,state=populations("fragmented cytoplasmic cDNA","fragmented ERCC cDNA","residual RNA/native amine-terminal fragments; class recovery unassigned"))
 g.react("shear","Shear selected cDNA",StepKind.FRAGMENTATION,((selected,InputRole.TEMPLATE),(shear,InputRole.CONDITION)),(frag,))
 kit=g.input("e7370","NEBNextUltraDNA E7370 native preparation; historical oligos/order unassigned",MaterialKind.CONDITION)
 native=g.product("rna_native","native barcodedE7370 prePCR material",MaterialKind.MOLECULAR,state=populations("nativeE7370-prepared cDNA; end/strand geometry unassigned"))
 g.react("e7370_prepare","NativeE7370 library-preparation boundary before amplification",StepKind.MAGIC,((frag,InputRole.TEMPLATE),(kit,InputRole.CONDITION)),(native,),reason="Defining2016 source names barcoding/amplification usingE7370. Native historical endprep/adapter oligos and internal preparation order are outside inspected cap; a2023 manual is not substituted.")
 rf=g.input("rna_forward","nativeE7370 PCR primerA",MaterialKind.OLIGO)
 rr=g.input("rna_reverse","nativeE7370 PCR primerB",MaterialKind.OLIGO)
 rp=g.input("rna_polymerase","nativeE7370 PCR polymerase",MaterialKind.CATALYST)
 rna=g.product("rna_library",RNA.name,MaterialKind.LIBRARY,state=RNA)
 g.react("rna_index_pcr","Amplify source barcoded cDNA library",StepKind.PCR,((native,InputRole.TEMPLATE),(rf,InputRole.STRUCTURAL),(rr,InputRole.STRUCTURAL),(rp,InputRole.CATALYST)),(rna,),geometry=PCRGeometry.unavailable(native.state,reason="Actual nativeE7370 prepared input; exact kit primer/junction assembly unassigned."))
 protease=g.input("protease","protease nuclear release",MaterialKind.CATALYST)
 lam=g.input("lambda","60fg unmethylated lambda DNA",MaterialKind.MOLECULAR)
 naked=g.product("naked","released nuclear DNA and independent lambda DNA",MaterialKind.MOLECULAR,state=dna_population("released DNA"))
 g.react("release_dna","Release nuclear DNA and add lambda control",StepKind.DEPROTEINIZATION,((nuc,InputRole.TEMPLATE),(lam,InputRole.SPIKE_IN),(protease,InputRole.CATALYST)),(naked,))
 mspi=g.input("mspi","inherited scRRBS9U MspI,37°C3h",MaterialKind.CATALYST)
 restriction=g.product("restriction","MspI fragments with opposed5′CG ends",MaterialKind.MOLECULAR,state=restriction_state())
 g.react("digest_mspi","Digest nuclear/lambda DNA with MspI",StepKind.DIGESTION,((naked,InputRole.TEMPLATE),(mspi,InputRole.CATALYST)),(restriction,))
 klenow=g.input("klenow","5U Klenowexo−,1mM dATP/0.1mM dGTP/0.1mM dCTP",MaterialKind.CATALYST)
 prep=g.product("dna_prepared","filled and3′dA-tailed fragments",MaterialKind.MOLECULAR,state=prepared_state())
 g.react("endprep","Fill recessed ends and add3′A with source Klenow mixture",StepKind.END_REPAIR,((restriction,InputRole.TEMPLATE),(klenow,InputRole.CATALYST)),(prep,))
 adapter=g.input("dna_adapter","standard premethylated indexed Illumina adapters; native strands/member unassigned",MaterialKind.ADAPTER)
 ligase=g.input("ligase","30U concentratedT4DNA ligase",MaterialKind.CATALYST)
 ligated=g.product("ligated","adapter-ligated nuclear/lambda DNA",MaterialKind.MOLECULAR,state=dna_population("native premethylated adapter-ligated DNA"))
 g.react("ligate","Ligate source premethylated indexed adapters",StepKind.LIGATION,((prep,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(ligated,),geometry=AdapterLigationGeometry.unavailable("Source names premethylated indexed adapters, not a chosen canonical native strand assembly."))
 bs=g.input("bisulfite","MethyCode:98°C10min then64°C2.5h",MaterialKind.CONDITION)
 converted=g.product("converted","physical converted single strands",MaterialKind.MOLECULAR,state=dna_population("unprotectedC→U; protected methyl-C retained",converted=True))
 g.react("convert","Bisulfite-convert adapter-ligated DNA",StepKind.CHEMICAL_CONVERSION,((ligated,InputRole.TEMPLATE),(bs,InputRole.CONDITION)),(converted,))
 trna=g.input("trna","10ng tRNA carrier at conversion-column recovery",MaterialKind.MOLECULAR)
 recovery=g.input("bs_column","on-column desulfonation/ZymoSpin recovery;30µl elution",MaterialKind.SELECTION)
 recovered=g.product("recovered","column-recovered converted DNA plus carrier; exact carrier recovery unassigned",MaterialKind.MOLECULAR,state=dna_population("converted U-bearing originals",converted=True,carrier=True))
 waste=g.product("conversion_waste","physical column fractions",MaterialKind.WASTE)
 g.react("recover_conversion","Desulfonate/purify converted DNA with tRNA carrier",StepKind.CLEANUP,((converted,InputRole.TEMPLATE),(trna,InputRole.SPIKE_IN),(recovery,InputRole.SELECTION)),(recovered,waste))
 g.discard(waste,"Source column flowthrough/washes; individual carrier/DNA recovery unassigned.")
 current=recovered
 for key,cycles,enzyme in (("dna_round1",25,"PfuTurboCx uracil-stalling-free"),("dna_round2",22,"PhusionHighFidelity")):
  f=g.input(key+"_forward","native scRRBS forward primer; sequence unassigned",MaterialKind.OLIGO)
  r=g.input(key+"_reverse","native scRRBS reverse primer; relation between rounds unassigned",MaterialKind.OLIGO)
  enzyme_port=g.input(key+"_polymerase",enzyme,MaterialKind.CATALYST)
  product=g.product(key,"amplified nuclear/lambda DNA; copiedU yieldsA/T",MaterialKind.MOLECULAR,state=dna_population("PCR-copied DNA; converted-site A/T",carrier=True))
  g.react(key+"_pcr",str(cycles)+"cycle "+enzyme+" amplification",StepKind.PCR,((current,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(enzyme_port,InputRole.CATALYST)),(product,),geometry=PCRGeometry.unavailable(current.state,reason="Actual preceding converted/nativePCR mixture; primer strands and inter-round transfer/cleanup remain unassigned."))
  current=product
 page=g.input("page","12%nativeTBE PAGE200–500bp",MaterialKind.SELECTION)
 dna=g.product("dna_library",DNA.name,MaterialKind.LIBRARY,state=DNA)
 rejected=g.product("page_waste","nonrecovered PAGE size fraction",MaterialKind.WASTE)
 g.react("size_select","Retain200–500bp DNA library after twoPCRrounds",StepKind.SELECTION,((current,InputRole.TEMPLATE),(page,InputRole.SELECTION)),(dna,rejected))
 g.discard(rejected,"Physically nonrecovered PAGE size fraction; tRNA/other species class recovery unassigned.")
 for lib,endpoint in ((rna,RNA),(dna,DNA)):g.sequence(lib,SequencingPlan(library=endpoint,acquired_roles=("Read 1","Read 2"),read_lengths=READ_LENGTHS,limitation="PublishedPE100; native read-primer/index acquisition unassigned."))
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()

REVIEWED=("cell","soft_lysis","separate","ercc","up1","rt","rt_copy","exo","digest_up1","tdt","tail","up2","copy_polymerase","second_copy","cdna_polymerase","cdna_pcr","aup1","aup2","amine_pcr","unmodified_cdna_precursors","cdna_column_selection","cdna_column","cdna_gel_selection","cdna_gel","covaris","shear","e7370","e7370_prepare","rna_forward","rna_reverse","rna_polymerase","rna_index_pcr","rna_library","protease","lambda","release_dna","mspi","digest_mspi","klenow","endprep","dna_adapter","ligase","ligate","bisulfite","convert","trna","bs_column","recover_conversion","dna_round1_forward","dna_round1_reverse","dna_round1_polymerase","dna_round1_pcr","dna_round2_forward","dna_round2_reverse","dna_round2_polymerase","dna_round2_pcr","page","size_select","dna_library")
