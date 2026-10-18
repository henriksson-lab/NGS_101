"""Native published two-level sciPlex-ATAC graph; independent hash/DNA classes."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature,revcomp
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind as M,InputRole as I,StepKind as K,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry
import nextera as nx
import seqprimers as sp
TITLE='sciPlex-ATAC-seq — native two-level route'
NOTES='01_sciplex_atac_seq.html'
SOURCE='Booth et al., <a href="https://doi.org/10.1186/s12864-023-09832-1">BMCGenomics2023 Methods169–197/247–271 and published SupplementaryTable4</a>; linked author two-level protocol, printedID87966, modified4Oct2023.'
SUMMARY='Sample hashes are independent DNA oligos retained in nuclei; extended capture molecules and tagmented genomic DNA receive independent first-well components before common PCR-well indexing and separate gel recovery.'
CAVEAT=('Selected published two-level route; the three-level experiment is distinct and outside this graph. <b>Published/linked-source discrepancies:</b> published Methods179 prints TTT+TTT…TTTVN and defines + as LNA on the following base (nine modified Ts after the initial unmodified T3); linked author protocol ID87966 pp3–4 prints unmodified T30VN. The graph uses the published formulation; sample-specific choice is unassigned. Published Methods197 uses hash199bp/ATAC~200–3000bp; linked protocol step8 uses~200bp/300–1200bp. Published windows are modelled. Published sequencing Methods261 assigns hashRead1 toUMI+PCRP5index andIndex2 to capture-well index; the printed capture oligo instead places its own wellindex afterUMI at a compatibleR1 site, with PCRP5index atIndex2. The source run lengths are retained without assigning these disputed cycles to a molecular barcode. Native complete PCR-primer strands and first-cycle/gap/anchor geometry remain unassigned; canonical final-library shells are not supplied.')
FINAL_LIBRARY=None
SEQ_PRIMERS=()
SEQUENCING_UNAVAILABLE=''
SEQUENCING_ENDING='Separate native ATAC and hash libraries: ATACcustomR1/R2=51, I1=8+27dark+10, I2=10+21dark+8 (18 acquired bases each); hashstandardR1=18,R2≥16,I1/I2=10. Physical barcode identity and acquired-cycle assignment remain separate; the published hashR1/I2 assignment is disputed as described above.'
HASH=feature('sample_hash','sample_index','whitelist',whitelist='Booth2023 publishedTable3 sample-hash identities',note='independent hashDNA, not a genomic insert segment')
UMI=feature('hash_umi','umi','random')
CAPTURE_WELL=feature('capture_well','cell_barcode','whitelist',group='first_well',part='capture',whitelist='source96firstwell capture oligos; exactmemberunassigned')
N5=feature('tn5_n5','cell_barcode','whitelist',group='first_well',part='n5',whitelist='sourceTable4 N5 plate8members')
N7=feature('tn5_n7','cell_barcode','whitelist',group='first_well',part='n7',whitelist='sourceTable4 N7 plate12members')
P5=feature('pcr_p5','cell_barcode','unknown',group='pcr_well',part='p5',note='sourcePCR-well coordinate; full orderedprimerunassigned')
P7=feature('pcr_p7','cell_barcode','unknown',group='pcr_well',part='p7',note='sourcePCR-well coordinate; full orderedprimerunassigned')
HASH_PREFIX=nx.ADAPTOR_S7
CAPTURE_PREFIX=nx.ADAPTOR_S5
CAPTURE_T30_WRITTEN='+'.join(['TTT']*10)
P5_I5_1_TRANSFERRED='TCGTCGGCAGCGTCTCCACGCTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG'
P7_I7_1_TRANSFERRED='GTCTCGTGGGCTCGGCTGTCCCTGTCCCGAGTAATCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG'
TN5ME_REV=nx.ME_RC
CUSTOM_R1='GCGATCGAGGACGGCAGATGTGTATAAGAGACAG'
CUSTOM_R2='CACCGTCTCCGCCTCAGATGTGTATAAGAGACAG'
CUSTOM_I1='CTGTCTCTTATACACATCTGAGGCGGAGACGGTG'
CUSTOM_I2='CTGTCTCTTATACACATCTGCCGTCCTCGATCGC'
CUSTOM_SEQ_PRIMERS=tuple(sp.custom(role_,name,seq,SOURCE) for role_,name,seq in [('Read 1','sourcecustomR1',CUSTOM_R1),('Read 2','sourcecustomR2',CUSTOM_R2),('Index 1 (i7)','sourcecustomI1',CUSTOM_I1),('Index 2 (i5)','sourcecustomI2',CUSTOM_I2)])
# Declared as a source oligo inventory: no exact native final molecule is claimed
# on which to assign these sequencing-primer landing sites.

def hash_parts():return [Segment('hash prefix',HASH_PREFIX),Segment('sample hash','N'*10,placeholder=True,feature=HASH),Segment('sourceBanchor','B',placeholder=True),Segment('hash polyA30','A'*30)]
def capture_parts():return [Segment('capture prefix',CAPTURE_PREFIX),Segment('hash UMI','N'*8,placeholder=True,feature=UMI),Segment('capture-well index','N'*10,placeholder=True,feature=CAPTURE_WELL),Segment('capture T30','T'*30,note='publishedphysicalLNA atT4,7,10,13,16,19,22,25,28'),Segment('capture VNanchor','VN',placeholder=True)]
def oligo_state(label,parts):
 sc=Scene();sc.strand('oligo',parts,label=label);return sc

def loading_scene():
 sc=Scene()
 for key,seq,offset,f in [('N5',P5_I5_1_TRANSFERRED,21,N5),('N7',P7_I7_1_TRANSFERRED,27,N7)]:
  me_start=len(seq)-len(nx.ME)
  sc.strand(key,[Segment(key+'outer',seq[:offset]),Segment(key+'index',seq[offset:offset+8],feature=f),Segment(key+'extendedME',seq[offset+8:me_start]),Segment(key+'ME',seq[me_start:])],label='representative transferred'+key)
  sc.anneal(key+'nontransfer',[Segment(key+'MErev',TN5ME_REV)],to=key,pair=(key+'MErev',key+'ME'),label='5Pnontransfer',mod5='p')
 return sc

def state(*,hashed=False,captured=False,extended=False,tagged=False,copied=False,atac_only=False,hash_only=False):
 sc=Scene()
 if copied:
  if not hash_only:
   sc.strand('ATAC daughterclass',[role('indexed ATAC daughterDNA; native complete primer ends/gap completion unassigned')],label='independent genomic daughterDNA')
   sc.strand('ATAC identifiers',[role('N5 barcode inventory',molecular_feature=N5,length_bp=8),role('N7 barcode inventory',molecular_feature=N7,length_bp=8)],label='unplaced ATAC firstwell identifiers; not a joined free molecule')
  if not atac_only:
   sc.strand('hash daughterclass',[role('hash-derived ordinary DNA daughters; native ends/anchor junction unassigned')],label='independent hash-derived daughterDNA')
   sc.strand('hash identifiers',[role('samplehash inventory',molecular_feature=HASH,length_bp=10),role('UMI inventory',molecular_feature=UMI,length_bp=8),role('capture-well barcode inventory',molecular_feature=CAPTURE_WELL,length_bp=10)],label='unplaced hash-class features; not an extra freeDNA chain')
   sc.note('hash daughterclass','Physical LNA residues belong to original capture molecules; DNA copying does not reproduce LNA chemistry. Native full primer-end placement is unassigned.')
  sc.strand('PCRwell indexinventory',[role('PCRP5 coordinate',molecular_feature=P5,length_bp=10),role('PCRP7 coordinate',molecular_feature=P7,length_bp=10)],label='unplaced PCRwell coordinates of daughterDNA; sample-specific members unassigned')
  sc.strand('carried originals',[role('possible freeRNA/originalhash/capture/genomicDNA/primers and non-targetclasses; physical recovery unassigned')],label='independent possible co-carried molecular populations')
  return sc
 if not hash_only:
  sc.strand('genome',[role('native gapped indexed genomicDNA' if tagged else 'native genomic DNA')],label='independent genomic DNA population')
  if tagged:
   sc.strand('N5 inventory',[role('transferredN5 end; nativegap/strand placement unassigned',molecular_feature=N5,length_bp=8)],label='N5 feature inventory, not free genomic hash')
   sc.strand('N7 inventory',[role('transferredN7 end; nativegap/strand placement unassigned',molecular_feature=N7,length_bp=8)],label='N7 feature inventory')
 if hashed and not atac_only:sc.strand('hash',hash_parts(),label='independent samplehash DNA originals')
 if captured and not atac_only:
  sc.strand('capture',capture_parts()+([role('new hash-template-derived DNA; nativeVNanchor junction unassigned',molecular_feature=HASH)] if extended else []),label='extendedcapture DNA' if extended else 'physical capture primer')
 sc.strand('possible carriedclasses',[role('possible originalRNA/protein/primers/nonproductiveDNA; source-specificphysical recovery unassigned')],label='independent possible co-carried molecules')
 return sc

def oligos():
 yield 'Sample hash:5′-'+HASH_PREFIX+'[10ntsamplehash]B'+'A'*30+'-3′'
 yield 'Published capture:5′-'+CAPTURE_PREFIX+'NNNNNNNN[10ntfirstwellindex]'+CAPTURE_T30_WRITTEN+'VN-3′; + denotes followingLNA-T.'
 for label,seq in [('N5-1 transferred',P5_I5_1_TRANSFERRED),('N7-1 transferred',P7_I7_1_TRANSFERRED),('5Pnontransfer','/5Phos/'+TN5ME_REV),('customR1',CUSTOM_R1),('customR2',CUSTOM_R2),('customI1',CUSTOM_I1),('customI2',CUSTOM_I2)]:yield label+':5′-'+seq+'-3′'

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 def inp(key,label,kind=M.CONDITION,state_=None):return g.input(key,label,kind,state=state_,source=SOURCE)
 def op(key,label,kind,x,state_,extra=(),**kw):
  y=g.product(key+'_out',label,M.MOLECULAR,state=state_);g.react(key,label,kind,((x,I.TEMPLATE),*extra),(y,),**kw);return y
 def unresolved(x,why):g.mark_terminal(Terminal(x.id,Disposition.UNRESOLVED,why))
 def recover(key,label,x):
  selection=inp(key+'_selection',label,M.SELECTION);keep=g.product(key+'_out','retained '+x.label,M.MOLECULAR,state=x.state);waste=g.product(key+'_waste','removed physical wash/flowthrough fraction',M.WASTE)
  g.react(key,label,K.CLEANUP,((x,I.TEMPLATE),(selection,I.SELECTION)),(keep,waste),changes_molecule=False);g.discard(waste,'Physical non-retained fraction; species-complete depletion not inferred.');return keep
 nuclei=g.input('nuclei','permeabilized nuclei separated bysample;selected two-level route',M.BIOLOGICAL,state=state(),primary=True)
 hash_=inp('hash_oligo','source ssDNA samplehash;source0.5pmol/1000cells',M.OLIGO,oligo_state('sourcehash',hash_parts()))
 x=op('hash','Absorb hash DNA into nuclei;ice5min',K.BINDING,nuclei,state(hashed=True),((hash_,I.STRUCTURAL),))
 fix=inp('fixative','1%finalformaldehyde15minice')
 x=op('fix','Fix hash-labelled nuclei',K.CROSS_LINKING,x,x.state,((fix,I.CONDITION),))
 peers=inp('sample_peers','independently hashed/fixed other sample populations',M.MOLECULAR,state(hashed=True))
 x=op('pool_samples','Pool independently fixed sample nuclei',K.POOL,x,x.state,((peers,I.TEMPLATE),),changes_molecule=False)
 x=recover('fix_wash','500gpellet/removefixative;resuspend/tweenwash/repellet',x)
 freeze=inp('freeze','sourcefreezingbuffer/flashfreeze−80C')
 x=op('freeze_thaw','Freeze then thaw poolednuclei',K.THERMAL_TREATMENT,x,x.state,((freeze,I.CONDITION),))
 allocation=inp('firstwell_rule','source5000nuclei/well;independent96wellallocations')
 first=g.product('firstwell_out','selectedfirst-roundwellnuclei',M.MOLECULAR,state=x.state)
 remainder=g.product('unallocated_nuclei','poolednuclei outside selectedallocation',M.MOLECULAR,state=x.state)
 g.react('firstwell','Allocate5000nuclei toselectedfirstwell',K.SAMPLING,((x,I.TEMPLATE),(allocation,I.CONDITION)),(first,remainder),changes_molecule=False)
 unresolved(remainder,'Unallocatedpoolmaterial outside selectedfirstwell lineage; physicaldispositionunassigned.')
 x=first
 capture=inp('capture_oligo','publishedLNAcapture/UMI8/firstwell10',M.OLIGO,oligo_state('captureprimer',capture_parts()))
 anneal=inp('anneal_conditions','55C5min thenice5min')
 x=op('anneal','Anneal capture to independent hashDNA',K.HYBRIDIZATION,x,state(hashed=True,captured=True),((capture,I.STRUCTURAL),(anneal,I.CONDITION)))
 pol=inp('capture_pol','NEBNextHighFidelity2X;55C10min',M.CATALYST)
 x=op('extend_hash','Extend captured hash DNA',K.EXTENSION,x,state(hashed=True,captured=True,extended=True),((pol,I.CATALYST),))
 adapters=inp('indexed_tn5_adapters','source8N5×12N7 indexedloadingoligos;representative1/1',M.ADAPTER,loading_scene())
 tn5=inp('tn5','sourceTn5/DMFtagbuffer+CLB55C15min',M.CATALYST)
 tagged=g.product('tagged','tagmentedgenome plus independent capturedhashes',M.MOLECULAR,state=state(hashed=True,captured=True,extended=True,tagged=True));other=g.product('tag_other','other/unconverted genomicDNA',M.MOLECULAR,state=x.state)
 g.react('tagment','Indexedtransposition ofgenomicDNA',K.TAGMENTATION,((x,I.TEMPLATE),(adapters,I.STRUCTURAL),(tn5,I.CATALYST)),(tagged,other));unresolved(other,'Transposition is not physical disposal; unconverted/offproduct recovery unassigned.')
 stop=inp('stop_mix','EDTA/spermidine;37C15min')
 x=op('stop','Stop transposition',K.QUENCHING,tagged,tagged.state,((stop,I.CONDITION),))
 otherwells=inp('firstwell_peers','independently indexed/captured first-roundwell populations',M.MOLECULAR,x.state)
 x=op('pool_firstwell','Pool first-roundwells andaddDAPI',K.POOL,x,x.state,((otherwells,I.TEMPLATE),),changes_molecule=False)
 sort=inp('sort_rule','DAPI FACS;limitednuclei persecondwell, experiment-specific')
 selected=g.product('sorted','selectednuclei in secondPCRwell',M.MOLECULAR,state=x.state);unused=g.product('unsorted','nuclei notsortedinto selectedPCRwell',M.MOLECULAR,state=x.state)
 g.react('sort','Physically sort limitednuclei into reversecrosslinkbuffer',K.PARTITION,((x,I.TEMPLATE),(sort,I.CONDITION)),(selected,unused),changes_molecule=False);unresolved(unused,'Outside selectedPCRwell;actualcollection/disposal notassigned.')
 proteinase=inp('proteinase','proteinaseK/SDS65C13.5h',M.CATALYST)
 x=op('reverse_crosslinks','Reverse crosslinks',K.DEPROTEINIZATION,selected,selected.state,((proteinase,I.CATALYST),))
 f=inp('pcr_p5','indexedNexteraP5PCRprimer;nativeorderedstrandunassigned',M.OLIGO);r=inp('pcr_p7','indexedNexteraP7PCRprimer;nativeorderedstrandunassigned',M.OLIGO);pcrmix=inp('pcr_mix','NEBNextHighFidelity2X+sourceTween20',M.CATALYST)
 initial=inp('initial_hold','72C5min beforecycling')
 x=op('initial_extension','Source72Cinitialextension withPCRprimers/mix',K.EXTENSION,x,x.state,((f,I.STRUCTURAL),(r,I.STRUCTURAL),(pcrmix,I.CATALYST),(initial,I.CONDITION)))
 product=g.product('amplified','coamplifiedDNA/hashdaughterclasses',M.MOLECULAR,state=state(hashed=True,captured=True,extended=True,tagged=True,copied=True));original=g.product('pcr_originals','originalmodifiedcapture/hash/genomicDNA andnonselectedclasses',M.MOLECULAR,state=x.state)
 g.react('pcr','PCRsource22cycles plusfinal72C5min',K.PCR,((x,I.TEMPLATE),(f,I.STRUCTURAL),(r,I.STRUCTURAL),(pcrmix,I.CATALYST)),(product,original),geometry=PCRGeometry.unavailable(x.state,reason='Actualmixednative gappedgenomicDNA andcapture-extendedhash DNA enterPCR. Fullorderedindexedprimers andnativefirst-cycle/anchorcompletionunassigned.'))
 unresolved(original,'PCR doesnot erase originalhash/captureLNA/genomicclasses; component-specificlaterrecoveryunassigned.')
 pcrpeers=inp('pcrwell_peers','independently PCR-indexed otherwell libraries',M.MOLECULAR,product.state)
 x=op('pool_pcr','Pool independentPCRwells',K.POOL,product,product.state,((pcrpeers,I.TEMPLATE),),changes_molecule=False)
 x=recover('column','Zymoclean/concentrate5Xbindingbuffer;elute100µlEB',x)
 gel=inp('gel_selection','1%agarose;publishedhash199bp/ATAC~200–3000bp',M.SELECTION)
 atac=g.product('atac_library','separate genomicATAC library;withoutcovalentlyjoinedhash',M.LIBRARY,state=state(tagged=True,copied=True,atac_only=True))
 hashes=g.product('hash_library','separate hash-derivedlibrary;LNAoriginalrecoveryunassigned',M.LIBRARY,state=state(hashed=True,captured=True,extended=True,copied=True,hash_only=True))
 rest=g.product('gel_rest','gelregionsoutside sourcecutwindows',M.WASTE)
 g.react('gel','Gel recover both libraryfractions;Nucleospin extraction',K.SELECTION,((x,I.TEMPLATE),(gel,I.SELECTION)),(atac,hashes,rest));g.discard(rest,'Physically notcut/retained;co-elutingmolecularclasseswithinwindows remainunassigned.')
 g.sequence(atac,SequencingPlan(acquired_roles=sp.ROLES,read_lengths={'Read 1':51,'Read 2':51,'Index 1 (i7)':18,'Index 2 (i5)':18},limitation='Sourcecustomsequencingprimerinventory shownabove. I1=8+27dark+10 andI2=10+21dark+8;18acquiredbases each. Complete nativeprimerlanding/endstructureunassigned; no samplehash isagenomicinsertsegment.'))
 g.sequence(hashes,SequencingPlan(acquired_roles=sp.ROLES,read_lengths={'Read 1':18,'Index 1 (i7)':10,'Index 2 (i5)':10},limitation='SourcehashRead2≥16,Read1=18,index1/index2=10. PublishedUMI/P5-versuscapture-well assignmentconflictswithorderedcapturegeometry; cycles arenotassignedtobarcodeuntilresolved. Complete nativeprimerlandingunassigned.'))
 g.declare_reviewed_scope(REVIEWED,SOURCE);return g.require_publishable()
REVIEWED=('nuclei','hash_oligo','hash','fixative','fix','sample_peers','pool_samples','fix_wash_selection','fix_wash','freeze','freeze_thaw','firstwell_rule','firstwell','capture_oligo','anneal_conditions','anneal','capture_pol','extend_hash','indexed_tn5_adapters','tn5','tagment','stop_mix','stop','firstwell_peers','pool_firstwell','sort_rule','sort','proteinase','reverse_crosslinks','pcr_p5','pcr_p7','pcr_mix','initial_hold','initial_extension','pcr','pcrwell_peers','pool_pcr','column_selection','column','gel_selection','gel','atac_library','hash_library')
