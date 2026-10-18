"""Protocol-owned native molecules and source-bounded graph."""
from chemdraw import Scene, Segment, revcomp, feature, oligo
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind as M,InputRole as I,StepKind as K,SequencingPlan,Terminal,Disposition
FINAL_LIBRARY=None
SEQ_PRIMERS=()
SEQUENCING_UNAVAILABLE=''
SHOW_READ_CONTENT=False
import seqprimers as sp
def population(*rows):
 sc=Scene()
 for key,parts,label in rows:sc.strand(key,tuple(parts),label=label)
 return sc
def recover(g,key,x,label):
 reagent=g.input(key+'_selection',label,M.SELECTION)
 kept=g.product(key+'_retained','recovered '+x.label,M.MOLECULAR,state=x.state)
 waste=g.product(key+'_waste','physical non-retained supernatant/washes',M.WASTE)
 g.react(key,label,K.CLEANUP,((x,I.TEMPLATE),(reagent,I.SELECTION)),(kept,waste),changes_molecule=False)
 g.discard(waste,'Actual non-retained fraction; complete species-specific depletion is not asserted.')
 return kept
def unresolved(g,x,reason):g.mark_terminal(Terminal(x.id,Disposition.UNRESOLVED,reason))
def op(g,key,x,label,kind,state,reagent,port=I.CATALYST):
 out=g.product(key+'_out',label,M.MOLECULAR,state=state)
 g.react(key,label,kind,((x,I.TEMPLATE),(reagent,port)),(out,),**({'changes_molecule':False} if kind in (K.GENERIC,K.QUENCHING) else {}))
 return out
import nextera as nx
from topology import HairpinTopology
TITLE='SAMOSA-Tag — native marked DNA and source hairpin transposition'
NOTES='01_samosa-tag.html'
SOURCE='Dennis et al. NatureGenetics2024 doi:10.1038/s41588-024-01748-0; publishedPMC11176058 Methods/Results/Fig1, originalSupplementaryTable5, MOESM1 repair/geometry notes.'
SUMMARY='EcoGII marks original native DNA before barcoded hairpin-Tn5 transposition. Explicit9ntgap filling/nick sealing produces candidate closed templates; ExoIII hydrolysis and three successive source bead recoveries do not imply complete chemical purity.'
CAVEAT='Selected in-situ cellline methylated route. Reserved no-EcoGII control is a distinct allocation retaining any native marks; ex-situ andPDX experiments andoptional high/low-size selections remain separate. Native methylation positions/efficiency, exact original-to-repair junction coordinates and polymerase/primer binding polarity are unassigned. Full repair uses short unmodified daughter patches, not a replacement amplified genomic duplex. Published oligo does not establish proprietary sequencing-primer absence.'
BC='CGGAAGAA'
LOOP='ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT'
BARCODE=feature('samosa_sample_barcode','sample_index','whitelist',whitelist='publishedSupplementTable5',note='representativeSMRT-A_bc001; not short-read acquired index cycles')
ADAPTER=(revcomp(nx.ME)+revcomp(BC)+LOOP+BC+nx.ME)
STEM=len(nx.ME)+len(BC)
def adapter_parts():
 return (Segment('mosaic-end reversecomplement',revcomp(nx.ME)),Segment('samplebarcode reversecomplement',revcomp(BC),placeholder=False,feature=BARCODE),Segment('hairpinloop',LOOP),Segment('samplebarcode',BC,feature=BARCODE),Segment('mosaic-end',nx.ME))
def adapter_scene():
 sc=population(('hairpin',adapter_parts(),'actual source singleDNA oligo; intramolecularstem shown separately'))
 sc.strands['hairpin'].mod5='p'
 return sc
TOPOLOGY_PANELS=(('Published native hairpin substrate',HairpinTopology('SMRT-A_bc001',ADAPTER,0,len(ADAPTER)-STEM,STEM).folded_panel('Source5′phosphorylated oligo; stem/barcode folding derived from Table5. This substrate is not a completed closed genomic template.')),)
LIMIT='Source closedDNA topology requires gapfill/nickseal after transfer. SourceSequelII8MSMRTcell30hmovies/2hpreextension/4himmobilization; polymerase2.1or2.2 library-size dependent. Native sequencing-primer/fullbinding geometry and which strand enters first are unassigned, not confirmed absent.'
SEQUENCING_ENDING=LIMIT
READ_SPECS=(sp.read_spec('Sequencing output',None,'Repeated native template observations and methylation kinetics; samplebarcode physical8bases is not8 acquired index cycles','include'),)
def native(stage='unmarked native',tagged=False,closed=False,digested=False):
 sc=population(('originalA',(role('original genomicDNAstrandA; '+stage+'; native/accessibility methylation positions unassigned'),),'original sideA of native duplex/template'),('originalB',(role('original genomicDNAstrandB; '+stage+'; independent native/accessibility methylation pattern'),),'original sideB; not a copied replacement'))
 if tagged:
  for side in ('left','right'):
   sc.strand(side+'hairpin',adapter_parts(),label='covalently transferred '+side+'hairpin component; native orientedjunction unassigned')
  sc.note('originalA','Source9nt gaps remain at transfer boundaries' if not closed else 'Bothsource9ntgaps repaired/nickssealed; candidate one covalentlyclosed template includes both originalDNA strands and both hairpins')
 if closed and not digested:sc.strand('openclasses',(role('possible incompletely filled/nick-bearing/open nativeDNA co-carried with candidatecircles; exact recovery unassigned'),),label='actual possible openDNA in mixed repair/recovery and ExoIII input; not a purified closed fraction')
 if closed:sc.strand('repairpatches',(role('short normalDNA gap-repair patches;9ntgap target per source; exact junction sequence/patch extent unassigned'),),label='newly synthesized repairDNA, distinct from marked original strands')
 if digested:sc.strand('digestclasses',(role('ExoIII hydrolysates and possible surviving incompletelyclosed DNA; subsequent bead partition unassigned'),),label='independent digestion/off-product classes')
 sc.strand('other',(role('RNA/protein/native fragments or their digestion products; physical recovery unassigned'),),label='independent possible carried constituents')
 return sc

def oligos():
 yield oligo('SMRT-A_bc001 sourceoligo',adapter_parts(),mods='/5Phos/')

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 nuclei=g.input('nuclei','source permeabilizedOS152/E14 nuclei inBufferM; selectedOS152 route',M.BIOLOGICAL,state=native(),primary=True)
 mix=g.input('sam','1mMSAM; replenish1.16mM at15min inmarked/control',M.CONDITION)
 marked_input=g.product('marking_aliquot','200µl source methylation allocation',M.MOLECULAR,state=nuclei.state)
 control=g.product('control','200µl reserved no-EcoGIIcontrol allocation',M.MOLECULAR,state=nuclei.state)
 g.react('reserve_control','Reserve separate no-enzyme control before EcoGII',K.SPLIT,((nuclei,I.TEMPLATE),(mix,I.CONDITION)),(marked_input,control));unresolved(g,control,'Separate sourcecontrol experiment, not physically discarded or falsely assumed fullyunmethylated. Its later preparation is outside selectedmarked-route graph.')
 ecogii=g.input('ecogii','source250UEcoGII;37C30min',M.CATALYST)
 x=op(g,'mark',marked_input,'Mark accessible adenines on originalDNA',K.CHEMICAL_CONVERSION,native('EcoGII m6dAmarked native'),ecogii)
 x=recover(g,'pellet',x,'Sourcepellet600g10min25C; resuspendOmniATACbuffer;40µmfilter')
 aliquot=g.product('tagmentation_aliquot','10000–50000 markednuclei sourcealiquot',M.MOLECULAR,state=x.state)
 remainder=g.product('remaining_nuclei','remaining methylatednuclei allocation',M.MOLECULAR,state=x.state)
 rule=g.input('aliquot_rule','source10000–50000nuclei allocation',M.CONDITION)
 g.react('take_nuclei','Take physical markednuclei aliquot',K.SAMPLING,((x,I.TEMPLATE),(rule,I.CONDITION)),(aliquot,remainder),changes_molecule=False);unresolved(g,remainder,'Unallocatednuclei not destroyed by sampling.')
 adapter=g.input('hairpin','source5Pbarcodedhairpin annealed95C5min/25C30min',M.ADAPTER,state=adapter_scene())
 tn5=g.input('tn5','sourceTn5R27S,E54K,L372P;loaded23C55min;9.4–18.8pmol selectedtagmentation',M.CATALYST)
 tagged=g.product('tagged','native hairpin-transferredDNA withsource9ntgaps',M.MOLECULAR,state=native('EcoGII marked',True))
 off=g.product('transfer_offproducts','other/unproductive tagmentation classes',M.MOLECULAR,state=native('EcoGII marked'))
 g.react('tagment','Hairpin transfer/tagmentation55C45–60min inOmniATAC',K.TAGMENTATION,((aliquot,I.TEMPLATE),(adapter,I.STRUCTURAL),(tn5,I.CATALYST)),(tagged,off));unresolved(g,off,'Nonproductive end identity is not physical discard; sourcebead recovery unassigned.')
 rnase=g.input('rnase','10µgRNaseA37C15min',M.CATALYST)
 x=op(g,'rna_digest',tagged,'SourceRNaseA digestion beforeproteinase stop',K.DIGESTION,tagged.state,rnase)
 stop=g.input('stop','50µgProteinaseK/SDS/EDTA;60C1–2h',M.CATALYST)
 x=op(g,'deproteinize',x,'Stop transposition and deproteinize',K.DEPROTEINIZATION,native('EcoGII marked',True),stop)
 x=recover(g,'first_spri',x,'Source2XSPRI;two80%EtOHwashes;20µlEB37C15min')
 x=recover(g,'size_spri',x,'Separate0.6XSPRI enriches>500bp; no complete species purity inferred')
 phusion=g.input('phusion','PhusionHF2U/0.8mMdNTP;37C1h',M.CATALYST)
 ligase=g.input('taq_ligase','80UTaqDNAligase/1XTaqbuffer;sourcecoupledrepair reaction',M.CATALYST)
 filled=g.product('filled','native gapfilled candidateDNA; nicksunsealed class',M.MOLECULAR,state=native('EcoGII original marks retained;normalDNA fill patches',True))
 g.react('fill_gaps','Fill source9nt transposition gaps inactual nativeDNA',K.EXTENSION,((x,I.TEMPLATE),(phusion,I.CATALYST)),(filled,))
 closed=g.product('closed','candidate closedtemplates and possible co-carried incomplete/open originalDNA',M.MOLECULAR,state=native('EcoGII marked',True,True))
 unclosed=g.product('unclosed','incomplete gapfill/sealing classes',M.MOLECULAR,state=x.state)
 g.react('seal_nicks','Seal repairednicks in same source37C1h coupledreaction',K.NICK_SEALING,((filled,I.TEMPLATE),(ligase,I.CATALYST)),(closed,unclosed));unresolved(g,unclosed,'SourceExoIII attacks susceptible unclosedDNA; hydrolysates and remaining intactspecies are not physically removed until recovery, exact partition unassigned.')
 x=recover(g,'repair_spri',closed,'Separate source2XSPRI aftergaprepair;12µlEB')
 exo=g.input('exo','ExoIII100Uper160ng/NEBuffer2;37C1h',M.CATALYST)
 x=op(g,'exo_digest',x,'Digest susceptible openDNA; hydrolysates remain molecular',K.DIGESTION,native('EcoGII marked',True,True,True),exo)
 x=recover(g,'final_spri',x,'Separate final2XSPRI afterExoIII;12µlEB')
 lib=g.product('library','nativeclosed SAMOSA-Tag library withunassigned residualpartition',M.LIBRARY,state=x.state)
 g.react('retain_library','Retain recovered native circular-template population',K.GENERIC,((x,I.TEMPLATE),),(lib,),changes_molecule=False)
 g.sequence(lib,SequencingPlan(limitation=LIMIT))
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()
REVIEWED=('nuclei','sam','reserve_control','ecogii','mark','pellet_selection','pellet','aliquot_rule','take_nuclei','hairpin','tn5','tagment','rnase','rna_digest','stop','deproteinize','first_spri_selection','first_spri','size_spri_selection','size_spri','phusion','taq_ligase','fill_gaps','seal_nicks','repair_spri_selection','repair_spri','exo','exo_digest','final_spri_selection','final_spri','retain_library','library')
