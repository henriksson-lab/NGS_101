"""Source-owned ZipSeq tumor-v3 chemistry; unknown native kit boundaries explicit."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct, Scene, Segment, feature, oligo, revcomp, complement_segments
from protocol_model import role
from endprep import EndPreparedDuplex, dA_tailed_scene
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal, Disposition
import seqprimers as sp
ZIP1="GTAGCAACCACAGATCGCACCCGAGAATTCCATGATGC"+"A"*28

TITLE="ZipSeq — sequential photo-uncaging of spatial zipcodes"
NOTES="01_zipseq.html"
SOURCE='<a href="https://doi.org/10.1038/s41592-020-0880-2">Hu et al., Nature Methods (2020)</a>'
SUMMARY="Independent RNA and hybridized DNA zipcodes acquire the same GEM cell identifier. Physical SPRI fractionation separates transcript and zipcode products before native library construction and a 1:10 zipcode:RNA molar sequencing pool."
CAVEAT="This is the tumor experiment's 10x 3-prime v3 branch, not the paper's v2 experiments. CG000183 Rev A supplies the represented 28/8/91 acquisition and 16-base cell barcode/12-base UMI; the actual experimental guide revision is unpinned. Native kit oligo bases, adapter geometry, additive-primer sequence and CITE RPI member/index extent remain unassigned in this capped pass. The published typeset oligo table was inspected: its text describes red NPOM-modified T positions, but those marks were not resolved; no preprint positions are substituted. Chemical uncaging and blocking do not themselves remove DNA. Computational RNA-cell whitelisting of zipcode counts does not physically discard GEM molecules."
CELL=feature("cell_barcode","cell_barcode","whitelist",whitelist="10x Chromium v3")
UMI=feature("umi","umi","random")
INDEX=feature("sample_index_i7","sample_index","unknown")
REGION=feature("zipseq_region","spatial_barcode","whitelist",whitelist="Hu named zipcodes")
O1_BLOCK="CGATCTGTGGTTGCTAC"
O1_CAGED=O1_BLOCK+"CCGTCAACATCAAAAAAAAAAAA"
O1_PRIME=revcomp(O1_BLOCK)
ZIPCODE1=ZIP1[-36:-28]
ZIP_HANDLE=ZIP1[len(O1_PRIME):-36]
CITE_RPI1="CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
LIMIT="Native end/primer geometry and exact zipcode sequencing entry remain unassigned; 28/8/91 are guide acquisition cycles, not proof of the native RPI molecular index width."
def cell_parts():
 return [role("native bead-derived handle"),Segment("cell barcode","B"*16,"cbc",placeholder=True,feature=CELL),Segment("UMI","N"*12,"umi",placeholder=True,feature=UMI),role("oligo(dT) priming tract; length unassigned")]
def copied_parts(z=False):
 return [*cell_parts(),*([Segment("copied zipcode",revcomp(ZIPCODE1),"cbc",feature=REGION),Segment("copied handle",revcomp(ZIP_HANDLE)),Segment("copied O1 prime",revcomp(O1_PRIME))] if z else [role("transcript cDNA; physical RNA U copied as DNA A")]),role("native opposite end; RT/switch/PCR footprint unassigned")]
def final(z=False):
 return Construct([role("native P5/Read1 end"),*copied_parts(z),role("native sample-index-bearing end",molecular_feature=INDEX),role("native P7 end")],"ZipSeq zipcode" if z else "ZipSeq transcript")
RNA_LIB=final(); ZIP_LIB=final(True)
HEADINGS=("Gene-expression library (tumor v3)","Zipcode 1 library (tumor v3)")
READ_LENGTHS={h:{"Read 1":28,"Index 1 (i7)":8,"Read 2":91} for h in HEADINGS}
READ_SPECS={h:(sp.read_spec("Read 1",28,"16-base cell barcode +12-base UMI","extract"),sp.read_spec("Index 1 (i7)",8,"sample index; native member/extent unassigned","demultiplex"),sp.read_spec("Read 2",91,"zipcode identity; exact entry unassigned" if n else "transcript insert","extract" if n else "include")) for n,h in enumerate(HEADINGS)}
SEQUENCING_UNAVAILABLE={h:LIMIT for h in HEADINGS}
FINAL_LIBRARIES=((HEADINGS[0],RNA_LIB,(),"Native GEX end roles; actual kit bases unassigned.",LIMIT),(HEADINGS[1],ZIP_LIB,(),"Native zipcode end roles; no assumed RPI member or index width.",LIMIT))
def anchor_scene(tagged=False,with_rna=False):
 sc=Scene(); a=[Segment("anchor prefix","ATCCAG"),Segment("anchor docking","T"*12+"GATGTTGACGG")]
 sc.strand("anchor",a,mod5="anti-CD45–Amine-C6",label="antibody-linked DNA anchor")
 sc.junction("anchor","anchor prefix","anchor docking","internal Cy5")
 c=[Segment("O1 overhang",O1_BLOCK),Segment("anchor complement",O1_CAGED[len(O1_BLOCK):])]
 sc.anneal("caged",c,to="anchor",pair=("anchor complement","anchor docking"),unpaired=("O1 overhang",),label="DNA; NPOM T positions unassigned before illumination")
 if tagged:
  sc.anneal("zipcode",[Segment("O1 prime",O1_PRIME),Segment("zipcode handle",ZIP_HANDLE),Segment("region zipcode",ZIPCODE1,"cbc",feature=REGION),Segment("poly(A)28","A"*28)],to="caged",pair=("O1 prime","O1 overhang"),unpaired=("zipcode handle","region zipcode","poly(A)28"),label="independent hybridized DNA zipcode")
 if with_rna: sc.strand("RNA",[role("independent cellular poly(A) RNA population; physical U")],label="RNA, not covalent to anchor or zipcode")
 return sc

def add_duplex(sc,key,parts,label):
 sc.strand(key,parts,label=label)
 bottom=complement_segments(parts)
 sc.anneal(key+" complement",bottom,to=key,pair=(bottom[0].name,parts[0].name),unpaired=tuple(p.name for p in bottom if p.is_role_token()),label=label)

def copies(which="both",amp=False):
 sc=Scene()
 for z,key in ((False,"transcript"),(True,"zipcode")):
  if which not in ("both",key): continue
  parts=copied_parts(z)
  if amp: add_duplex(sc,key,parts,"amplified "+key)
  else: sc.strand(key,parts,label="actual single-stranded "+key+" copy; native switched end bounded")
 sc.strand("residual",[role("possible original RNA, zipcode templates, anchor/blocker DNA and amplification precursors; partition/integrity unresolved")],label="independent possible co-carried lineages; no purity assumption")
 return sc

def gex_intermediate(ligated=False):
 parts=[*cell_parts(),role("fragmented transcript cDNA")]
 if ligated:
  parts.extend([Segment("source dA junction","A"),role("ligated native v3 adapter; strand/ligation footprint unassigned")])
  sc=Scene();add_duplex(sc,"adapted",parts,"adapter-ligated actual PCR precursor; native geometry bounded")
 else:
  sc=dA_tailed_scene(EndPreparedDuplex(Construct(parts),five_prime_phosphate=False,three_prime_overhang="A"),label="source guide3.1 fragmented/end-prepared dA-tailed cDNA;5′P marking unassigned")
 sc.strand("residual",[role("possible original RNA/zipcode/anchor/blocker/copied products; partition unassigned")],label="independent unresolved residual lineages")
 return sc

def library_scene(z=False):
 sc=Scene(); lib=ZIP_LIB if z else RNA_LIB
 add_duplex(sc,"library",list(lib),"expected amplified native library")
 sc.strand("residual",[role("possible co-carried original/copy/template/anchor/blocker products; unresolved recovery")],label="independent possible molecular residuals")
 return sc

def oligos():
 yield oligo("tumor antibody anchor",[Segment("prefix","ATCCAG"),Segment("docking","T"*12+"GATGTTGACGG")],mods="Amine-C6; internal Cy5 between prefix and docking")
 yield oligo("O1 caged",[Segment("source bases",O1_CAGED)],mods="NPOM-modified T; exact published marks unassigned in this pass")
 yield oligo("O1 blocker",[Segment("source bases",O1_BLOCK)])
 yield oligo("Zipcode 1",[Segment("O1 prime",O1_PRIME),Segment("handle",ZIP_HANDLE),Segment("region",ZIPCODE1,"cbc",feature=REGION),Segment("poly(A)28","A"*28)])

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 def inp(i,l,k=MaterialKind.CONDITION,state=None): return g.input(i,l,k,state=state)
 def prod(i,l,state,k=MaterialKind.MOLECULAR): return g.product(i,l,k,state=state)
 def op(i,l,k,inputs,state):
  out=prod(i+"_out",l,state); g.react(i,l,k,inputs,(out,),geometry=AdapterLigationGeometry.unavailable("Native v3 adapter strand/end geometry not assigned in this capped pass.") if k is StepKind.LIGATION else None); return out
 def clean(i,l,x,state=None):
  sel=inp(i+"_selection",l,MaterialKind.SELECTION); out=prod(i+"_out",l,state or x.state); waste=prod(i+"_waste","physically removed wash/supernatant fraction",None,MaterialKind.WASTE)
  g.react(i,l,StepKind.CLEANUP,((x,InputRole.TEMPLATE),(sel,InputRole.SELECTION)),(out,waste));g.discard(waste,"physically removed fraction; exact molecular partition unassigned");return out
 starting=Scene();starting.strand("RNA",[role("cellular poly(A) RNA population; physical U")],label="RNA within live cells")
 cells=g.input("cells","150 µm live tumor sections; independent cellular RNA",MaterialKind.BIOLOGICAL,state=starting,primary=True)
 antibody=inp("antibody_anchor","anti-CD45 internal-Cy5 DNA conjugate",MaterialKind.MOLECULAR)
 cage=inp("caged_oligo","O1 NPOM strand",MaterialKind.OLIGO,O1_CAGED)
 pre=op("prehybridize","Prehybridize anchor and caged strand 1:1,37°C15min; cool",StepKind.HYBRIDIZATION,((antibody,InputRole.TEMPLATE),(cage,InputRole.STRUCTURAL)),anchor_scene())
 coat=op("coat","ssDNA/BSA block; antibody-anchor stain 1h4°C",StepKind.BINDING,((cells,InputRole.TEMPLATE),(pre,InputRole.ANALYTE)),anchor_scene(with_rna=True))
 washed=clean("coat_wash","Wash unbound conjugates",coat)
 light=inp("light","365 nm patterned1000ms pulses")
 unc=op("uncage","Photolyze NPOM cages; DNA sequence retained",StepKind.CHEMICAL_CONVERSION,((washed,InputRole.TEMPLATE),(light,InputRole.CONDITION)),anchor_scene(with_rna=True))
 zipin=inp("zipcode","Zipcode1 poly(A) DNA",MaterialKind.OLIGO,ZIP1)
 tagged=op("hybridize","Hybridize zipcode to exposed O1 overhang",StepKind.HYBRIDIZATION,((unc,InputRole.TEMPLATE),(zipin,InputRole.STRUCTURAL)),anchor_scene(True,True))
 washed=clean("tag_wash","Wash excess zipcode",tagged)
 blocker=inp("blocker","O1 blocking strand",MaterialKind.OLIGO,O1_BLOCK)
 blocked=op("block","Hybridize blocker to carried zipcode; repeat spatial rounds as required",StepKind.HYBRIDIZATION,((washed,InputRole.TEMPLATE),(blocker,InputRole.STRUCTURAL)),anchor_scene(True,True))
 blocked.state.strand("blocker products",[role("independent blocked excess zipcode; not chemically destroyed")],label="possible carried blocking complexes")
 washed=clean("block_wash","Wash following blocking/labeling rounds",blocked)
 sort=inp("sort_gate","collagenase I/IV +mechanical dissociation; live Cy5+ FACS",MaterialKind.SELECTION)
 selected=prod("sorted","live Cy5+ cells",washed.state); rejected=prod("sort_waste","gate-excluded cells/debris",None,MaterialKind.WASTE)
 g.react("sort","Dissociate, strain and sort live Cy5+ cells",StepKind.SELECTION,((washed,InputRole.TEMPLATE),(sort,InputRole.SELECTION)),(selected,rejected));g.discard(rejected,"physical FACS gate")
 selected=clean("cell_wash","PBS0.04%BSA wash",selected)
 gem=inp("gem_conditions","10x3′v3 GEM partitioning;8k cell target")
 occupied=prod("occupied","represented cell-containing GEM population",selected.state)
 other=prod("other_gems","other GEM populations; unresolved molecular fates",None)
 g.react("partition","Partition cells into GEMs",StepKind.PARTITION,((selected,InputRole.TEMPLATE),(gem,InputRole.CONDITION)),(occupied,other))
 g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,"Empty/multiple/unassigned GEM populations; no physical discard inferred from barcode QC."))
 bead=inp("bead_primer","v3 native cell16/UMI12 oligo(dT) primer; full strand unassigned",MaterialKind.OLIGO)
 rt=inp("rt_mix","native v3 RT mix",MaterialKind.CATALYST)
 copied=op("rt","Copy independent RNA and zipcode-DNA templates with bead primers",StepKind.REVERSE_TRANSCRIPTION,((occupied,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),copies())
 recovery=inp("gem_recovery_conditions","CG000183 native post-GEM recovery/cleanup boundary; detailed operations not assigned in inspected sections")
 recovered=prod("recovered_copies","native recovered mixed RT copies; physical partition/integrity unresolved",copies())
 g.react("gem_recovery","Native post-GEM recovery before cDNA amplification",StepKind.MAGIC,((copied,InputRole.TEMPLATE),(recovery,InputRole.CONDITION)),(recovered,),reason="Guide process has post-GEM recovery before cDNA amplification; detailed native recovery operations/partition were not inspected within the frozen source cap.")
 copied=recovered
 f=inp("cdna_primer","native10x cDNA amplification primer species",MaterialKind.OLIGO);r=inp("additive","partial SmallRNA Read2 additive; native sequence/revision unassigned",MaterialKind.OLIGO);pol=inp("amp_mix","native cDNA amplification mix",MaterialKind.CATALYST)
 amp=prod("amplified","amplified transcript and short zipcode molecules",copies(amp=True))
 g.react("cdna_pcr","Amplify cDNA with partial SmallRNA Read2 additive",StepKind.PCR,((copied,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(copied.state,reason="Native primer species/end geometry not assigned in inspected sources."))
 spri=inp("split_selection","post-amplification SPRI",MaterialKind.SELECTION)
 rna=prod("rna_fraction","bead-bound transcript cDNA",copies("transcript",True));z=prod("zip_fraction","decanted <200bp zipcode supernatant",copies("zipcode",True))
 g.react("split","SPRI physical bead/supernatant split; purity unassigned",StepKind.SELECTION,((amp,InputRole.TEMPLATE),(spri,InputRole.SELECTION)),(rna,z))
 z=clean("zip_clean1","First3×SPRI zipcode recovery",z);z=clean("zip_clean2","Second3×SPRI zipcode recovery",z)
 endmix=inp("fragment_endprep_mix","CG0001833.1 fragmentation/endrepair/A-tailing mix",MaterialKind.CATALYST)
 rna=op("fragment_endprep","Fragment/end-repair/A-tail cDNA (CG0001833.1)",StepKind.END_REPAIR,((rna,InputRole.TEMPLATE),(endmix,InputRole.CATALYST)),gex_intermediate())
 rna=clean("rna_size","Double-sided SPRI3.2",rna)
 adapter=inp("rna_adapter","native v3 adapter; ligation geometry unassigned",MaterialKind.OLIGO);lig=inp("ligase","native adapter-ligation mix",MaterialKind.CATALYST)
 rna=op("rna_ligate","Adapter ligation3.3; native geometry remains bounded",StepKind.LIGATION,((rna,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(lig,InputRole.CATALYST)),gex_intermediate(True))
 rna=clean("rna_lig_cleanup","Post-ligation cleanup3.4",rna)
 outputs=[]
 for prefix,x,iszip in (("rna",rna,False),("zip",z,True)):
  f=inp(prefix+"_forward","CITE10xP5 class" if iszip else "native10x forward index-PCR primer",MaterialKind.OLIGO)
  r=inp(prefix+"_reverse","CITE smallRNA RPIx class; member/extent unassigned" if iszip else "native10x sample-index primer",MaterialKind.OLIGO)
  pol=inp(prefix+"_polymerase","KAPA HiFi" if iszip else "native10x PCR mix",MaterialKind.CATALYST)
  out=prod(prefix+"_library","native amplified "+prefix+" library",library_scene(iszip),MaterialKind.LIBRARY)
  g.react(prefix+"_index_pcr","Amplify actual native "+prefix+" precursor",StepKind.PCR,((x,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(out,),geometry=PCRGeometry.unavailable(x.state,reason="Native ordered primer/member geometry unassigned; no completed library used as input."))
  if not iszip: out=clean("rna_final_size","Final double-sided SPRI3.6",out)
  outputs.append(out)
 poolstate=Scene()
 for name,lib in (("RNA",RNA_LIB),("zipcode",ZIP_LIB)):
  add_duplex(poolstate,name,list(lib),"independent "+name+" library molecules")
 poolstate.strand("residual",[role("possible co-carried original/templates/anchor/blocker fragments; fate unresolved")],label="independent residual populations")
 pooled=prod("sequencing_pool","1:10 zipcode:RNA molar pool; Bioanalyzer/qPCR",poolstate,MaterialKind.LIBRARY)
 g.react("pool","Mix zipcode and associated RNA libraries1:10molar",StepKind.POOL,tuple((x,InputRole.TEMPLATE) for x in outputs),(pooled,),changes_molecule=False)
 g.sequence(pooled,SequencingPlan(primers=(),acquired_roles=("Read 1","Index 1 (i7)","Read 2"),read_lengths=READ_LENGTHS[HEADINGS[0]],limitation=LIMIT))
 # Manual source scope: local known operations only; no absence/full-review claim.
 g.declare_reviewed_scope(("gem_recovery","gem_recovery_conditions","prehybridize","coat","coat_wash","uncage","hybridize","tag_wash","block","block_wash","sort","cell_wash","partition","rt","cdna_pcr","split","zip_clean1","zip_clean2","fragment_endprep","rna_size","rna_ligate","rna_lig_cleanup","rna_index_pcr","zip_index_pcr","rna_final_size","pool","cells","antibody_anchor","caged_oligo","coat_wash_selection","light","zipcode","tag_wash_selection","blocker","block_wash_selection","sort_gate","cell_wash_selection","gem_conditions","bead_primer","rt_mix","cdna_primer","additive","amp_mix","split_selection","zip_clean1_selection","zip_clean2_selection","fragment_endprep_mix","rna_size_selection","rna_adapter","ligase","rna_lig_cleanup_selection","rna_forward","rna_reverse","rna_polymerase","zip_forward","zip_reverse","zip_polymerase","rna_final_size_selection","sequencing_pool"),SOURCE)
 return g.require_publishable()
