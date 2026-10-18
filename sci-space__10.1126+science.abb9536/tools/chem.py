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
 return [role(p.name+" copied complement",molecular_feature=p.feature,inferred=p.inferred) if p.is_role_token() else replace(p,name=p.name+" copied complement",top=revcomp(p.top)) for p in reversed(parts)]
def tools(g):
 def inp(i,label,kind=MaterialKind.CATALYST,state=None):return g.input(i,label,kind,state=state)
 def out(i,label,state=None,kind=MaterialKind.MOLECULAR):return g.product(i,label,kind,state=state)
 def op(i,label,kind,inputs,state,geometry=None,reason=""):
  product=out(i+"_out",label,state);g.react(i,label,kind,inputs,(product,),geometry=geometry,**({"reason":reason} if reason else {}));return product
 def clean(i,label,x,state=None):
  selection=inp(i+"_selection",label,MaterialKind.SELECTION);kept=out(i+"_out",label,x.state if state is None else state,x.kind);waste=out(i+"_waste","physically removed wash/flowthrough fraction",None,MaterialKind.WASTE)
  g.react(i,label,StepKind.CLEANUP,((x,InputRole.TEMPLATE),(selection,InputRole.SELECTION)),(kept,waste));g.discard(waste,"Source physical recovery fraction; exact molecular partition/purity unassigned.");return kept
 return inp,out,op,clean
HISTORICAL_SOURCE_NOTICES = ('The printed hash, section hash and transcript are separate molecules that are associated through shared nucleus indices; the final construct shows the spatial-hash branch, not a fictitious covalent fusion to a transcript. The readable sources do not establish the section-hash tail sequence, the barcoded oligo(dT) RT primer or the index members. '
          '<b>Source differences.</b> The published Supplementary Materials (Materials and Methods, pp. 3-5) differ from the linked protocols.io guides on several conditions; the model follows the study-specific Supplementary Materials and the differences remain visible here: '
          'fixation of the scraped embryo section in 4% paraformaldehyde (SM p. 3, "Spatial transfer") versus 5% PFA in the transfer guide (steps 8-10) and in the SM hash-concentration test (p. 2); '
          'indexed oligo(dT) at 25 uM (SM p. 4) versus 10 uM barcoded RT primer (two-level guide step 15); an RT ramp starting at 10 C (SM) versus 4 C (guide); '
          'index PCR 75 C 3 min then 18 cycles with 1-min extension (SM p. 5) versus 72 C 5 min then 18-22 cycles (guide step 53); '
          'pooled Zymo DNA Clean & Concentrator plus 0.8x AMPure XP (SM p. 5) versus pooled 1x AMPure (guide step 55); and only a custom n7-loaded Tn5 in the SM, whereas the guide also allows dual-loaded Illumina Tn5. '
          'These are differences between a general guide and the study methods, not established publication errors.')

TITLE="sci-Space — native independent transcript and DNA-hash populations"
NOTES="01_sci_space.html"
SOURCE='<a href="https://doi.org/10.1126/science.abb9536">Srivatsan et al., Science2021</a>; published SupplementaryMethodspp1–5 and cited gridprotocol bm64k9gw,13Feb2021.'
SUMMARY="Independent spot, sector and section DNA hashes are associated with endogenous RNA through nucleus RT/PCR indices. This graph selects Slides1–6: pooledRT nuclei are sorted200perPCRwell. The independently prepared Slides7–14 dilution250perwell cohort is described separately, without an invented split."
CAVEAT="INFERRED identity-window placement — RTprimer sequence/native layout is unassigned; displayed UMI8/RTbarcode10 follow reported R1content and do not reconstruct a verified full primer. Native customn7adapter, ordered indexedPCRprimers and exact fragment/end geometry remain unassigned in this capped pass. Printed spot/sector DNAhasher sequences follow the exact cited2021gridguide; section-hash tail is unassigned. Physical RNAU, originalDNAhashes and DNAcopies remain distinct. Three DNAhashes do not become one molecule with a transcript. Source18/10/10/52cycles and explicitly10bpPCRindices are separate statements; no generic8nt index is used. "
CAVEAT+=" Historical previously reviewed guide comparisons (not newly expanded in this pass): "+HISTORICAL_SOURCE_NOTICES.replace("The readable sources do not establish the section-hash tail sequence, the barcoded oligo(dT) RT primer or the index members.","Those native strand/member assignments remain unassigned in this model.")
CELL=feature("sci_space.rt","cell_barcode","combinatorial",group="sci_space.cell",part="RTwell")
UMI=feature("sci_space.umi","umi","random")
IDX5=feature("sci_space.pcr5","sample_index","unknown",group="sci_space.cell",part="PCRwellP5index",note="Secondcellindex allocation, not spatial/section identity")
IDX7=feature("sci_space.pcr7","sample_index","unknown",group="sci_space.cell",part="PCRwellP7index",note="PairedP5/P7members identify PCRwell")
SPOT=feature("sci_space.spot","spatial_barcode","whitelist",group="sci_space.position",part="spot",whitelist="cited2021grid oligo panel")
SECTOR=feature("sci_space.sector","spatial_barcode","whitelist",group="sci_space.position",part="sector",whitelist="cited2021grid oligo panel")
SECTION=feature("sci_space.section","sample_index","unknown",note="Slide/section hash on independent DNA molecule; not PCRindex")
HASH_PREFIX="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
KINDS=("transcript","spot","sector","section")
def hash_parts(kind):
 if kind=="section":return [role("native section-hash strand/tail geometry unassigned",molecular_feature=SECTION),role("source poly(A)capture tail; extent unassigned")]
 return [Segment("native grid5′prefix",HASH_PREFIX),Segment(kind+"barcode","N"*10,"cbc",placeholder=True,feature=SPOT if kind=="spot" else SECTOR),Segment("sourceBanchor(C/G/T)","B",placeholder=True),Segment("sourcepoly(A)32","A"*32)]
def grid_scene():
 sc=Scene()
 for k in ("spot","sector"):sc.strand(k,hash_parts(k),label="independent "+k+"DNAhash")
 return sc

def native_populations(*,section=False,grid=False):
 sc=Scene();sc.strand("RNA",[role("endogenous polyadenylated RNA; physical U")],label="independent RNA population")
 if section:sc.strand("sectionDNA",hash_parts("section"),label="independent slide/sectionDNAhash")
 if grid:
  for k in ("spot","sector"):sc.strand(k+"DNA",hash_parts(k),label="independent "+k+"DNAhash; no covalent transcript fusion")
 return sc
def rt_end():
 return [role("native RTprimer handle; full species unassigned",inferred=True),Segment("UMIwindow","N"*8,"umi",placeholder=True,feature=UMI,inferred=True),Segment("RTwellbarcodewindow","N"*10,"cbc",placeholder=True,feature=CELL,inferred=True),role("native oligo(dT)tract; extent/placement unassigned",inferred=True)]
def copy_parts(kind,tagged=False):
 if kind=="transcript":payload=[role("antisense transcriptcDNA; copied RNA U represented in DNA alphabet")]
 elif kind=="section":payload=[role("copied section-hash DNA; native tail/primer junction unassigned",molecular_feature=SECTION)]
 else:payload=[role("complement of sourceBanchor(C/G/T); A/C/Gbase unassigned",length_bp=1),*rc_parts(hash_parts(kind)[:2])]
 return [*rt_end(),*payload,*([role("native n7-tagged or retained hash-prefix end; gap/nick/fragment boundaries unassigned")] if tagged else [])]
def copies(duplex=False,tagged=False):
 sc=Scene()
 for k in KINDS:
  parts=copy_parts(k,tagged)
  if duplex:add_duplex(sc,k,parts,"native copied "+k+" dsDNA; INFERRED identity-window placement")
  else:sc.strand(k,parts,label="actual ssDNA "+k+"copy; INFERRED identity-window placement")
 sc.strand("prior molecules",[role("independent possible originalRNA, originalspot/sector/sectionDNAhashes, RTprimers, copies and nonamplifying fragments; integrity/partition unresolved")],label="prior lineages are not erased by copying or computational filters")
 return sc
FINAL_LIBRARY=Construct([role("native P5-derived end"),Segment("PCRwellP5index","N"*10,"cbc",placeholder=True,feature=IDX5),*rt_end(),role("representative copied spatialhash DNA",molecular_feature=SPOT),role("native n7/opposite end; exactgeometry unassigned"),Segment("PCRwellP7index","N"*10,"cbc",placeholder=True,feature=IDX7),role("native P7-derived end")],name="representative sci-Space native spatialhash library")
SEQ_PRIMERS=();FINAL_CAPTION="Independent transcript/spot/sector/sectionlibrary populations; spatialhash representative only. Native full primer/adapter layout unassigned."
SEQUENCING_ENDING="PublishedSMp5:NextSeq500 onboard primers,R1=18(UMI8+RTbarcode10),I1=10,I2=10,R2=52(transcript or DNAhash). Native primer entry/layout and ordered index members unassigned."
READ_SPECS=(sp.read_spec("Read 1",18,"UMI8+RTbarcode10; nucleus association includes PCRwellpair","extract"),sp.read_spec("Index 1 (i7)",10,"PCRwellindex1; sourcealso explicitly says10basepair","demultiplex"),sp.read_spec("Index 2 (i5)",10,"PCRwellindex2; sourcealso explicitly says10basepair","demultiplex"),sp.read_spec("Read 2",52,"transcript or independent hash molecule; processing differs by class","extract_then_include"))
def library_state():
 sc=copies(True,True)
 for k in KINDS:
  feature_obj=SPOT if k=="spot" else SECTOR if k=="sector" else SECTION if k=="section" else None
  parts=[role("native P5PCRend"),Segment("P5PCRwellindex","N"*10,"cbc",placeholder=True,feature=IDX5),*copy_parts(k,True),Segment("P7PCRwellindex","N"*10,"cbc",placeholder=True,feature=IDX7),role("native P7PCRend")]
  add_duplex(sc,"daughter "+k,parts,"expected indexed "+k+"daughter; original/copy residuals retained")
 return sc
def oligos():
 for k in ("spot","sector"):yield oligo("cited2021native "+k+"DNAhash",hash_parts(k))
 yield "Native RTprimer/indexPCRprimer full species remain unassigned; physical UMI8/RTbarcode10 and PCRindex10+10 are source-declared independently of run lengths."
def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE);inp,out,op,clean=tools(g)
 tissue=g.input("tissue","freshfrozen embryosection selectedSlides1–6",MaterialKind.BIOLOGICAL,state=native_populations(),primary=True)
 section=inp("section_hash","independent slide-specific poly(A)DNAhash",MaterialKind.OLIGO)
 buffer=inp("permeabilization","cold IGEPAL/RNaseinhibitor permeabilization buffer+DAPI",MaterialKind.CONDITION)
 labelled=op("section_label","Permeabilize and introduce independent sectionDNAhash",StepKind.BINDING,((tissue,InputRole.TEMPLATE),(section,InputRole.STRUCTURAL),(buffer,InputRole.CONDITION)),native_populations(section=True))
 grid=inp("grid","printed agarose grid with independent spot/sectorDNAhashes and SYBRfiducials",MaterialKind.OLIGO,grid_scene())
 labelled=op("grid_transfer","Diffuse gridDNAhashes into nuclei; image tissue/gridregistration",StepKind.BINDING,((labelled,InputRole.TEMPLATE),(grid,InputRole.STRUCTURAL)),native_populations(section=True,grid=True))
 pfa=inp("pfa","scrape sourceSlides1–6 into4%PFA;15min onice",MaterialKind.CONDITION)
 fixed=op("fix","Scrape and fix labeled nuclei; no Slides7–14sonication branch imposed",StepKind.CHEMICAL_CONVERSION,((labelled,InputRole.TEMPLATE),(pfa,InputRole.CONDITION)),labelled.state)
 fixed=clean("nuclei_recovery","800g10min pellet;NSBpool;flashfreeze",fixed)
 perm=inp("thaw_perm","thaw/pellet;NSB0.25%Triton3min;pellet/resuspend;12secsonication foraggregates",MaterialKind.CONDITION)
 nuclei=op("permeabilize","Recover and permeabilize nuclei before indexedRTallocation",StepKind.PERMEABILIZATION,((fixed,InputRole.TEMPLATE),(perm,InputRole.CONDITION)),fixed.state)
 primer=inp("rt_primer","uniquely indexed25µM oligo(dT); native full strand unassigned",MaterialKind.OLIGO)
 rt=inp("ssiv","SuperscriptIV ramp10/20/30/40/50°C2mineach,55°C15min")
 first=op("rt","INFERRED primer-layout—copy RNA and each separate DNAhash in indexedRTwells",StepKind.REVERSE_TRANSCRIPTION,((nuclei,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),copies())
 peers=inp("peer_rt_wells","other independent indexedRTnuclear populations",MaterialKind.MOLECULAR,copies())
 pooled=op("pool_rt","Pool indexedRTwells,filtercap,DAPIstain",StepKind.POOL,((first,InputRole.TEMPLATE),(peers,InputRole.TEMPLATE)),copies())
 sort=inp("sort_rule","Slides1–6:FACS200nuclei per96PCRwell;8RTplates. IndependentSlides7–14dilution250notthisbranch",MaterialKind.SELECTION)
 seeded=out("seeded","200sorted nuclei per downstream well",copies());other=out("unallocated","other/unallocated poolednuclei; fate unassigned",copies())
 g.react("redistribute","Sort independently sampled200nucleus wellallocations",StepKind.SELECTION,((pooled,InputRole.TEMPLATE),(sort,InputRole.SELECTION)),(seeded,other));g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,"Unallocated nuclear population is not a fabricated alternative-cohort branch or computational QC discard."))
 second=inp("second_mix","NEBmRNAsecondstrand mix16°C150min")
 ds=op("second_strand","Second-strand synthesis on actual separate first-copy molecules",StepKind.EXTENSION,((seeded,InputRole.TEMPLATE),(second,InputRole.CATALYST)),copies(True))
 transposome=inp("n7_transposome","customn7loadedTn5; native full strand/endgeometry unassigned",MaterialKind.ADAPTER);tn5=inp("tn5","source0.01µlTn5in2×TD;55°C5min")
 tagged=op("tagment","React copied populations with customn7Tn5; retainedhash/nativefragment classes unresolved",StepKind.TAGMENTATION,((ds,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),copies(False,True))
 stop=inp("zymo_stop","12µlZymoDNAbindingbuffer5minRT",MaterialKind.CONDITION)
 tagged=op("stop_tagment","Stop native tagmentation before physical recovery",StepKind.QUENCHING,((tagged,InputRole.TEMPLATE),(stop,InputRole.CONDITION)),tagged.state)
 tagged=clean("tag_cleanup","36µlAMPure;17µlEBeluate;16µltransferred toPCR",tagged)
 f=inp("p5","indexedP5primer; ordered full strand/member unassigned",MaterialKind.OLIGO);r=inp("p7","indexedP7primer; ordered full strand/member unassigned",MaterialKind.OLIGO);pol=inp("pcr_mix","NEBNextHiFi;75°C3min then18cycles98/66/72°C1minextension")
 tagged=op("initial_extension","Native75°C3min pre-cycle extension; completion/gaps bounded",StepKind.EXTENSION,((tagged,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),copies(False,True))
 amp=out("amplified","perPCRwell indexed transcript/hash libraries",library_state(),MaterialKind.LIBRARY)
 g.react("pcr","INFERRED native-end layout—PCR actual recovered native copied/tagged populations",StepKind.PCR,((tagged,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(tagged.state,reason="Actual native copied/tagged mixture supplied, not final duplexes. RTprimer/customn7/fullindexprimer geometry and firstcycle strand boundaries remain unassigned."))
 peerlibs=inp("peer_pcr_wells","other independently indexedPCRwelllibraries",MaterialKind.LIBRARY,library_state())
 pool=out("pooled_libraries","physically pooledPCRwell libraries",library_state(),MaterialKind.LIBRARY)
 g.react("pool_pcr","Pool actual PCRwell libraries",StepKind.POOL,((amp,InputRole.TEMPLATE),(peerlibs,InputRole.TEMPLATE)),(pool,),changes_molecule=False)
 pool=clean("concentrate","ZymoDNACleanandConcentrator recovery",pool);pool=clean("final_cleanup","0.8×AMPure recovery;Qubit/TapeStation/2pMloading",pool)
 g.sequence(pool,SequencingPlan(acquired_roles=("Read 1","Index 1 (i7)","Index 2 (i5)","Read 2"),read_lengths={"Read 1":18,"Index 1 (i7)":10,"Index 2 (i5)":10,"Read 2":52},limitation=SEQUENCING_ENDING))
 g.declare_reviewed_scope(("tissue","section_hash","permeabilization","section_label","grid","grid_transfer","pfa","fix","nuclei_recovery_selection","nuclei_recovery","thaw_perm","permeabilize","rt_primer","ssiv","rt","peer_rt_wells","pool_rt","sort_rule","redistribute","second_mix","second_strand","n7_transposome","tn5","tagment","zymo_stop","stop_tagment","tag_cleanup_selection","tag_cleanup","p5","p7","pcr_mix","initial_extension","pcr","peer_pcr_wells","pool_pcr","concentrate_selection","concentrate","final_cleanup_selection","final_cleanup","final_cleanup_out"),SOURCE)
 return g.require_publishable()
