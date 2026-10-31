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
from page import head, info, caveat, table
from protocol_views import WorkflowSection, PlanLimitationSection, assemble_page
PAGE_FILE='damage-seq.html';TITLE='Damage-seq — copy the boundary before a bulky adduct';NOTES='01_damage-seq.html'
SOURCE='<a href="https://doi.org/10.1073/pnas.1614430113">Hu et al. (2016)</a>, published main PMC5068337, Results, Figure 1 and Methods. The chemistry-relevant SI Appendix remains unavailable in the current source bundle.'
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
   sc.strand('template',[role('Original Pt-adduct-bearing DNA; lesion position unassigned')],label='Original template; Pt adduct is not copied')
  if self.stage in {'stalled','associated'}:
   sc.strand('stalled',[role('Biotin-primer-derived DNA ending near a blocking adduct')],label='Stalled new strand; does not contain the template adduct')
   sc.note('stalled','Biotin position and exact nascent 3-prime geometry uninspected; no genomic coordinate assigned.')
  if self.stage in {'unblocked','associated'}:sc.strand('unblocked',[role('Biotin-primer-derived full-length extension product')],label='Unblocked background new strand')
  if self.stage=='library':
   sc.strand('stalled_library',[role('Adapter-2/index-PCR library from stalled extension')],label='Source-sequenced stalled class')
   sc.strand('unblocked_library',[role('Adapter-2/index-PCR library from unblocked extension')],label='Source-sequenced background class; filter reads after sequencing')
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
 g=ProtocolGraph(TITLE, authority='Published Hu2016 main Results/Figure1; relevant SI unavailable, scope explicitly partial', representation='published_instance')
 dna=g.input('dna','Genomic DNA from Pt-treated GM12878 cells; precise experiment assignment unreviewed',K.MOLECULAR,state=DamageState('template'),primary=True)
 sonication=g.input('sonication','sourcesonication;settingsunreadSI',K.CONDITION)
 fragments=g.product('fragments','Sonicated DNA fragments',K.MOLECULAR,state=dna.state)
 g.react('fragment','Sonicate isolated genomic DNA',S.FRAGMENTATION,((dna,R.TEMPLATE),(sonication,R.CONDITION)),(fragments,))
 enzyme=g.input('endrepair','Source end-repair reagents; exact identity awaits SI',K.CATALYST)
 repaired=g.product('endrepaired','End-repaired DNA',K.MOLECULAR,state=fragments.state)
 g.react('repair_ends','Repair fragment ends',S.END_REPAIR,((fragments,R.TEMPLATE),(enzyme,R.CATALYST)),(repaired,))
 a1=g.input('adapter1','Source first adapter; exact polymer/ends await SI',K.ADAPTER,state=molecule('adapter1',(role('unassembledfirstadapter'),)))
 ligated=g.product('adapter1_ligated','First-adapter-ligated DNA',K.MOLECULAR,state=repaired.state)
 g.react('ligate_first','Ligate source first adapter',S.LIGATION,((repaired,R.TEMPLATE),(a1,R.STRUCTURAL)),(ligated,),geometry=unknown_ligation('firstadapter'))
 heat=g.input('denaturation','Source denaturation; conditions await SI',K.CONDITION)
 denatured=g.product('denatured','Denatured first-adapter-bearing DNA',K.MOLECULAR,state=ligated.state)
 g.react('denature','Denature before anti-Pt selection',S.THERMAL_TREATMENT,((ligated,R.TEMPLATE),(heat,R.CONDITION)),(denatured,))
 antibody=g.input('anti_pt','Anti-Pt antibody selection; exact IP conditions await SI',K.SELECTION)
 selected=g.product('selected','Pt-enriched DNA pool; purity unassigned',K.MOLECULAR,state=DamageState('template'))
 flow=g.product('ip_flowthrough','Nonselected IP material',K.MOLECULAR,state=molecule('nonselectedmaterial',(role('DNA and other material not retained in IP'),)))
 g.react('immunoprecipitate','Select Pt-adduct-bearing DNA strands',S.AFFINITY_CAPTURE,((denatured,R.TEMPLATE),(antibody,R.SELECTION)),(selected,flow));g.mark_terminal(Terminal(flow.id,Disposition.UNRESOLVED,'Nonselectedfractionnotcarriedintorepresentedlibrary;physicaldispositionnotreconstructedwithoutSI.'))
 primer=g.input('biotin_primer','Biotinylated adapter-1-binding primer; exact polymer/tag position awaits SI',K.OLIGO,state=molecule('biotinprimer',(role('Source biotinylated primer; positions unassigned'),)))
 q5=g.input('q5','NEBNext Q5 DNA polymerase; identity reported in main Results',K.CATALYST)
 # Molecular contribution is required by extension, but its source formulation is not assigned.
 extension_supply=g.input('extension_supply','Unassigned nucleotide-bearing molecular supply for source primer extension; identity and dose await SI',K.MOLECULAR)
 stalled=g.product('stalled','Biotinylated stalled-extension class',K.MOLECULAR,state=DamageState('stalled'))
 unblocked=g.product('unblocked','Biotinylated unblocked-extension class',K.MOLECULAR,state=DamageState('unblocked'))
 templates=g.product('remaining_templates','Original templates and possible free primer remain in preparation',K.MOLECULAR,state=DamageState('template'))
 g.react('extend','Q5 primer extension; damaged template can stall the new strand',S.EXTENSION,((selected,R.TEMPLATE),(primer,R.STRUCTURAL),(q5,R.CATALYST),(extension_supply,R.STRUCTURAL)),(stalled,unblocked,templates))
 beads=g.input('streptavidin','Streptavidin-coated beads; exact separation conditions await SI',K.SELECTION)
 recovered=g.product('extension_products','Selected biotin-extension product pool; associated template/free-primer fates unassigned',K.MOLECULAR,state=DamageState('associated'))
 nonselected=g.product('extension_flowthrough','Nonretained capture fraction',K.MOLECULAR,state=molecule('unretainedfraction',(role('nonretainedmaterial;cannotassigneveryoriginaltemplatehere'),)))
 g.react('purify_extension','Purify biotin-bearing extension products',S.AFFINITY_CAPTURE,((stalled,R.TEMPLATE),(unblocked,R.TEMPLATE),(templates,R.TEMPLATE),(beads,R.SELECTION)),(recovered,nonselected));g.mark_terminal(Terminal(nonselected.id,Disposition.UNRESOLVED,'Sourceaffinityrecoverypositive;exactwash/template-separationfateawaitsSI.'))
 a2=g.input('adapter2','Source second adapter; exact polymer/ends await SI',K.ADAPTER,state=molecule('adapter2',(role('unassembledsecondadapter'),)))
 tagged=g.product('adapter2_ligated','Second-adapter-ligated pool with stalled/background classes',K.MOLECULAR,state=DamageState('associated'))
 g.react('ligate_second','Ligate source second adapter; retain both product classes',S.LIGATION,((recovered,R.TEMPLATE),(a2,R.STRUCTURAL)),(tagged,),geometry=unknown_ligation('secondadapter'))
 qc_condition=g.input('qc_condition','Reserve 1% ligation products for two-primer-set QC',K.CONDITION)
 qc=g.product('qc','1% ligation-product QC aliquot',K.QC,state=tagged.state)
 retained=g.product('retained','Remaining library preparation material',K.MOLECULAR,state=tagged.state)
 g.react('qc_split','Reserve diagnostic PCR QC aliquot',S.SPLIT,((tagged,R.TEMPLATE),(qc_condition,R.CONDITION)),(qc,retained));g.mark_terminal(Terminal(qc.id,Disposition.QC,'Sourcecomparestwoprimersetsforunblocked/allproducts;exactoligos/reactionconditionsunassembled,notmodeledasquantitativeproof.'))
 forward=g.input('forward_index_primer','Source indexed-PCR forward polymer; exact identity awaits SI',K.OLIGO,state=molecule('forwardPCRprimer',(role('unassembledsourceforwardprimer'),)))
 reverse=g.input('reverse_index_primer','Source indexed-PCR reverse polymer; exact identity awaits SI',K.OLIGO,state=molecule('reversePCRprimer',(Segment('index-bearingprimer;positionunassigned','[unassembledindexedprimer]',placeholder=True,feature=INDEX),)))
 mix=g.input('pcr_mix','Unassigned molecular reagent supply for source library PCR; no kit recipe inferred',K.MOLECULAR)
 polymerase=g.product('pcr_enzyme','Unassigned catalytic portion of that same source PCR supply',K.CATALYST)
 nts=g.product('pcr_nucleotides','Unassigned nucleotide-bearing portion of that same source PCR supply',K.MOLECULAR)
 conditions=g.product('pcr_conditions','Unassigned nonincorporated portion of that same source PCR supply',K.CONDITION)
 g.react('pcr_supply_components','Inventory view of one unresolved PCR supply; no physical separation or extra additions',S.MAGIC,((mix,R.TEMPLATE),),(polymerase,nts,conditions),changes_molecule=False,reason='Conceptual inventory partition required by PCR, not a disclosed formulation. Exact ingredients, identities and conditions await the SI Appendix.')
 library=g.product('library','Indexed sequencing library with stalled and unblocked classes',K.LIBRARY,state=DamageState('library'))
 g.react('pcr','Source indexed PCR; both molecular classes proceed to sequencing',S.PCR,((retained,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(polymerase,R.CATALYST),(nts,R.STRUCTURAL),(conditions,R.CONDITION)),(library,),geometry=PCRGeometry.unavailable(molecule('sourcePCRpool',(role('bothprimer-derivedproductclasses;exactendsunassembled'),)),reason='SIAppendixrelevantbutunavailable;exactsourceprimerpolymers/geometrynotinspected. Noexactendpointinvented.'))
 g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_UNAVAILABLE))
 g.declare_reviewed_scope(('dna','sonication','fragments','fragment','endrepair','endrepaired','repair_ends','adapter1','adapter1_ligated','ligate_first','denaturation','denatured','denature','anti_pt','selected','ip_flowthrough','immunoprecipitate','biotin_primer','q5','extension_supply','stalled','unblocked','remaining_templates','extend','streptavidin','extension_products','extension_flowthrough','purify_extension','adapter2','adapter2_ligated','ligate_second','qc_condition','qc','retained','qc_split','forward_index_primer','reverse_index_primer','pcr_mix','pcr_supply_components','pcr_enzyme','pcr_nucleotides','pcr_conditions','library','pcr'),g.authority)
 return g.require_publishable()
def render_page():
 g=protocol_graph()
 return assemble_page(g,(
  head(TITLE),'<div class="wrap" data-protocol-complete="true"><h1>'+TITLE+'</h1>',
  '<p class="research-notes"><a href="'+NOTES+'">Research notes and source checklist</a></p>',info(SOURCE),
  '<h2>Copy the lesion boundary</h2>',
  '<p>Damage-seq maps bulky cisplatin adducts by extending a primer until the polymerase stalls. The platinum adduct stays on the original genomic template. The newly synthesized strand reports the nearby boundary without containing the adduct itself.</p>',
  '<p>Genomic DNA is fragmented, end-repaired and given a first adapter. Anti-Pt enrichment selects damaged strands after denaturation. A biotinylated adapter-binding primer is extended with Q5; streptavidin recovery is followed by a second adapter and indexed PCR.</p>',
  caveat('Source-limited published preparation: the main article and Figure 1 establish this order, while the detailed SI Appendix remains behind the recorded download gate. Exact oligo bases, end geometry, separation conditions, sequencing platform and run mask remain unassembled. They are not claimed absent from publication. Unassigned molecular reagent supplies expose model gaps; they do not specify a buffer or kit formulation.'),
  '<h2>Keep molecular background separate from read filtering</h2>',
  table(('Material or output','Interpretation'),(
   ('Stalled extension','Primer-derived DNA ends near a blocking template lesion; exact molecular geometry remains unresolved.'),
   ('Unblocked extension','Background products remain in the preparation and sequencing library; no perfect antibody selection or stalling efficiency is assumed.'),
   ('One-percent QC aliquot','The source compares two primer sets on a reserved ligation-product aliquot. The model records the analytical branch without inventing its oligos or treating it as a sequencing library.'),
   ('Read filtering','Adapter-1-starting background reads are discarded after sequencing. This is a computational filter, not a physical discard of every background molecule.'),
   ('Lesion coordinate','The source places the lesion immediately upstream of the aligned boundary read. This readout interpretation does not locate a Pt group in the new DNA strand.'))),
  '<h2>Sequencing scope</h2>',info(SEQUENCING_ENDING),PlanLimitationSection('library'),
  WorkflowSection('Source-qualified preparation and unresolved molecular fates'),
  '<p>The distinct Pt-XR-seq branch in the same paper uses different chemistry and is outside this graph.</p></div>'))


def selftest():
 from checks import Check, Source, have
 from storage import sources_root
 c=Check()
 source=Source(sources_root()/HERE.parent.name/'pnas.1614430113_PMC5068337.html.txt','https://pmc.ncbi.nlm.nih.gov/articles/PMC5068337/')
 if have(c,source,label='published Damage-seq main'):
  text=' '.join(source.path.read_text().split())
  c('named extension catalyst retains published identity','extended by NEBNext Q5 DNA polymerase' in text)
  c('QC fraction comes from ligation products','1% of the ligation products' in text)
  c('background discard is post-sequencing','After sequencing, reads from undamaged DNA' in text and 'discarded before alignment and data analysis' in text)
  c('reported lesion boundary is upstream of reads','immediately 5′ (upstream) of our reads' in text)
 protocol_graph().require_publishable()
 c.report()
