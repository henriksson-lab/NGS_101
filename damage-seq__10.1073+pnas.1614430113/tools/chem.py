"""Main-original Damage-seq partial graph; relevant SI unavailable, no guessed oligos."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
from chemdraw import Construct,Scene,Segment,feature
from protocol_graph import ProtocolGraph,InputRole as R,MaterialKind as K,StepKind as S,SequencingPlan,Terminal,Disposition
from protocol_model import molecule,role
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
import seqprimers as sp
PAGE_FILE='damage-seq.html';TITLE='Damage-seq — copy the boundary before a bulky adduct';NOTES='01_damage-seq.html'
SOURCE='<a href="https://doi.org/10.1073/pnas.1614430113">Hu et al.(2016)</a>,publishedmainPMC5068337,Results/Figure1/Methods. RelevantSIAppendixremainsunavailableinthecurrentsourcebundle.'
SUMMARY='Fragment genomic DNA, ligate an adapter and select Pt-adduct-bearing strands, then extend a biotinylated primer with Q5. Stalled and unblocked products both proceed to library preparation; full-length background reads are filtered after sequencing.'
CAVEAT='This is a source-limited partial model from full main text,not a complete chemistry review. SIAppendix download is blocked by the source gate; exact oligos,ligatingends/conditions,strand-separationdetails and runmask remain unassembled,not claimed unpublished. New primer-derived DNA does not contain the template Pt adduct. Antibody enrichment and polymerase stalling are not modeled as perfect or quantitative. Both stalled and unblocked classes remain sequenced; computational filtering is not molecular disposal.'
@dataclass(frozen=True)
class DamageState:
 stage:str
 def __post_init__(self):
  if self.stage not in {'template','stalled','unblocked','associated','library'}:raise ValueError('unknownDamage-seqstate')
 def rows(self):
  sc=Scene()
  if self.stage in {'template','associated'}:
   sc.strand('template',[role('originalPt-adduct-bearingDNA;lesionpositionunassigned')],label='originaltemplate;Ptadductnotcopied')
  if self.stage in {'stalled','associated'}:
   sc.strand('stalled',[role('biotinprimer-derivedDNAendingadjacenttoblockingadduct')],label='stallednewstrand;doesnotincludetemplateadduct')
   sc.note('stalled','Biotinpositionandexactnascent3-primeendgeometrynotinspected;noassignedgenomiccoordinate.')
  if self.stage in {'unblocked','associated'}:sc.strand('unblocked',[role('biotinprimer-derivedfull-lengthextensionproduct')],label='unblockedbackgroundnewstrand')
  if self.stage=='library':
   sc.strand('stalled_library',[role('adapter2/index-PCRlibraryfromstalledextension')],label='source-sequencedstalledclass')
   sc.strand('unblocked_library',[role('adapter2/index-PCRlibraryfromunblockedextension')],label='source-sequencedbackgroundclass;filterreadsaftersequencing')
  return sc.rows()
INDEX=feature('damage_sample_index','sample_index','unknown')
FINAL_LIBRARY=Construct((role('unassembledlibraryend'),role('stalledorunblockedprimer-derivedsequence'),Segment('sourceindex;member,length,sideunassigned','[sampleindex]',placeholder=True,feature=INDEX),role('unassembledoppositelibraryend')),name='Damage-seqlibrary-roleoverview')
FINAL_CAPTION='Role-only overview;not an exact adapter order or index-position assignment. Both stalled/background molecular classes are retained in graph library state.'
SEQ_PRIMERS=();SEQUENCING_UNAVAILABLE='Platform, read lengths, acquired read roles, index mask and loaded primers await the chemistry-relevant SI Appendix. Main-text lesion-boundary analysis does not establish exact run geometry.'
SEQUENCING_ENDING='The source filters adapter-1-starting undamaged-derived reads after sequencing. Lesions map immediately 5′ upstream of aligned boundary reads; the lesion is not physically inserted in the new strand.'
READ_SPECS=(sp.read_spec('Read 1','unassigned','damage-boundaryreportingpath;exactacquiredroleassignmentunresolved','undetermined'),)

def unknown_ligation(label):
 sc=Scene();sc.strand('end',[role('adapterendgeometrynotinspected;relevantSIunavailable')],label=label)
 return AdapterLigationGeometry('unavailable','unavailable',(sc,),label+':sourcepositiveligation,exactoligos/endconditionsawaitSI;notanabsenceclaim.')

def protocol_graph():
 g=ProtocolGraph(TITLE,authority='PublishedHu2016mainResults/Figure1;relevantSIunavailable,scopeexplicitlypartial')
 dna=g.input('dna','isolatedgenomicDNAfromPt-treatedGM12878cells;preciseexperimentassignmentunreviewed',K.MOLECULAR,state=DamageState('template'),primary=True)
 sonication=g.input('sonication','sourcesonication;settingsunreadSI',K.CONDITION)
 fragments=g.product('fragments','sonicatedDNAfragments',K.MOLECULAR,state=dna.state)
 g.react('fragment','Sonicate isolated genomic DNA',S.FRAGMENTATION,((dna,R.TEMPLATE),(sonication,R.CONDITION)),(fragments,))
 enzyme=g.input('endrepair','sourceendrepairreagents;exactidentityunreadSI',K.CATALYST)
 repaired=g.product('endrepaired','end-repairedDNA',K.MOLECULAR,state=fragments.state)
 g.react('repair_ends','Repair fragment ends',S.END_REPAIR,((fragments,R.TEMPLATE),(enzyme,R.CATALYST)),(repaired,))
 a1=g.input('adapter1','sourcefirstadapter;exactpolymer/endsunreadSI',K.ADAPTER,state=molecule('adapter1',(role('unassembledfirstadapter'),)))
 ligated=g.product('adapter1_ligated','first-adapter-ligatedDNA',K.MOLECULAR,state=repaired.state)
 g.react('ligate_first','Ligate source first adapter',S.LIGATION,((repaired,R.TEMPLATE),(a1,R.STRUCTURAL)),(ligated,),geometry=unknown_ligation('firstadapter'))
 heat=g.input('denaturation','sourcedenaturation;conditionsunreadSI',K.CONDITION)
 denatured=g.product('denatured','denaturedfirst-adapter-bearingDNA',K.MOLECULAR,state=ligated.state)
 g.react('denature','Denature before anti-Pt selection',S.THERMAL_TREATMENT,((ligated,R.TEMPLATE),(heat,R.CONDITION)),(denatured,))
 antibody=g.input('anti_pt','antiPtantibodyselection;exactIPconditionsunreadSI',K.SELECTION)
 selected=g.product('selected','Pt-enrichedDNApool;purityunassigned',K.MOLECULAR,state=DamageState('template'))
 flow=g.product('ip_flowthrough','nonselectedIPmaterial',K.MOLECULAR,state=molecule('nonselectedmaterial',(role('DNAandothermaterialnotretainedinIP'),)))
 g.react('immunoprecipitate','Select Pt-adduct-bearing DNA strands',S.AFFINITY_CAPTURE,((denatured,R.TEMPLATE),(antibody,R.SELECTION)),(selected,flow));g.mark_terminal(Terminal(flow.id,Disposition.UNRESOLVED,'Nonselectedfractionnotcarriedintorepresentedlibrary;physicaldispositionnotreconstructedwithoutSI.'))
 primer=g.input('biotin_primer','biotinylatedadapter1-bindingprimer;exactpolymer/tagpositionunreadSI',K.OLIGO,state=molecule('biotinprimer',(role('sourcebiotinylatedprimer;positionsunassigned'),)))
 q5=g.input('q5','NEBNextQ5DNApolymerase',K.CATALYST)
 stalled=g.product('stalled','biotinylatedstalledextensionclass',K.MOLECULAR,state=DamageState('stalled'))
 unblocked=g.product('unblocked','biotinylatedunblockedextensionclass',K.MOLECULAR,state=DamageState('unblocked'))
 templates=g.product('remaining_templates','originaltemplatesandpossiblefreeprimerremaininpreparation',K.MOLECULAR,state=DamageState('template'))
 g.react('extend','Q5 primer extension; damaged template can stall the new strand',S.EXTENSION,((selected,R.TEMPLATE),(primer,R.STRUCTURAL),(q5,R.CATALYST)),(stalled,unblocked,templates))
 beads=g.input('streptavidin','streptavidin-coatedbeads;exactseparationconditionsunreadSI',K.SELECTION)
 recovered=g.product('extension_products','selectedbiotin-extensionproductpool;associatedtemplate/freeprimerfatesunassigned',K.MOLECULAR,state=DamageState('associated'))
 nonselected=g.product('extension_flowthrough','nonretainedcapturefraction',K.MOLECULAR,state=molecule('unretainedfraction',(role('nonretainedmaterial;cannotassigneveryoriginaltemplatehere'),)))
 g.react('purify_extension','Purify biotin-bearing extension products',S.AFFINITY_CAPTURE,((stalled,R.TEMPLATE),(unblocked,R.TEMPLATE),(templates,R.TEMPLATE),(beads,R.SELECTION)),(recovered,nonselected));g.mark_terminal(Terminal(nonselected.id,Disposition.UNRESOLVED,'Sourceaffinityrecoverypositive;exactwash/template-separationfateawaitsSI.'))
 a2=g.input('adapter2','sourcesecondadapter;exactpolymer/endsunreadSI',K.ADAPTER,state=molecule('adapter2',(role('unassembledsecondadapter'),)))
 tagged=g.product('adapter2_ligated','second-adapter-ligatedpoolwithstalled/backgroundclasses',K.MOLECULAR,state=DamageState('associated'))
 g.react('ligate_second','Ligate source second adapter; retain both product classes',S.LIGATION,((recovered,R.TEMPLATE),(a2,R.STRUCTURAL)),(tagged,),geometry=unknown_ligation('secondadapter'))
 qc_condition=g.input('qc_condition','reserve1%ligationproductsfortwo-primer-setQC',K.CONDITION)
 qc=g.product('qc','1%ligation-productQCaliquot',K.QC,state=tagged.state)
 retained=g.product('retained','remaininglibrarypreparationmaterial',K.MOLECULAR,state=tagged.state)
 g.react('qc_split','Reserve diagnostic PCR QC aliquot',S.SPLIT,((tagged,R.TEMPLATE),(qc_condition,R.CONDITION)),(qc,retained));g.mark_terminal(Terminal(qc.id,Disposition.QC,'Sourcecomparestwoprimersetsforunblocked/allproducts;exactoligos/reactionconditionsunassembled,notmodeledasquantitativeproof.'))
 forward=g.input('forward_index_primer','sourceindexed-PCRforwardpolymer;exactidentityunreadSI',K.OLIGO,state=molecule('forwardPCRprimer',(role('unassembledsourceforwardprimer'),)))
 reverse=g.input('reverse_index_primer','sourceindexed-PCRreversepolymer;exactidentityunreadSI',K.OLIGO,state=molecule('reversePCRprimer',(Segment('index-bearingprimer;positionunassigned','[unassembledindexedprimer]',placeholder=True,feature=INDEX),)))
 mix=g.input('pcr_mix','sourcelibraryPCRmixture;exactconditionsunreadSI',K.CATALYST)
 library=g.product('library','indexedsequencinglibrarywithstalledandunblockedclasses',K.LIBRARY,state=DamageState('library'))
 g.react('pcr','Source indexed PCR; both molecular classes proceed to sequencing',S.PCR,((retained,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(mix,R.CATALYST)),(library,),geometry=PCRGeometry.unavailable(molecule('sourcePCRpool',(role('bothprimer-derivedproductclasses;exactendsunassembled'),)),reason='SIAppendixrelevantbutunavailable;exactsourceprimerpolymers/geometrynotinspected. Noexactendpointinvented.'))
 g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_UNAVAILABLE))
 g.declare_reviewed_scope(('dna','sonication','fragments','fragment','endrepair','endrepaired','repair_ends','adapter1','adapter1_ligated','ligate_first','denaturation','denatured','denature','anti_pt','selected','ip_flowthrough','immunoprecipitate','biotin_primer','q5','stalled','unblocked','remaining_templates','extend','streptavidin','extension_products','extension_flowthrough','purify_extension','adapter2','adapter2_ligated','ligate_second','qc_condition','qc','retained','qc_split','forward_index_primer','reverse_index_primer','pcr_mix','library','pcr'),g.authority)
 return g.require_publishable()
def render_page():
 from batch_page import render
 return render(sys.modules[__name__])
def selftest():
 protocol_graph();render_page();print('PASS source-limitedDamagegraphandrendering;noexactoligoclaims')
