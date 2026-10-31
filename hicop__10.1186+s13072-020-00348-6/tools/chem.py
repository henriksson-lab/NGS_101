"""Source-qualified public preparation with partial molecular detail; no full source/absence audit."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem,SourceCoverage
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from chemdraw import Scene,Segment

def unknown_scene():
 s=Scene();s.strand('unassembled',[Segment('unknown_template','[unassembled]',placeholder=True,inferred=True)]);return s
class T:
 def __init__(self,g,x):self.g=g;self.x=x
 def inp(self,k,label,kind=K.CONDITION,role=R.CONDITION,**kw):return(self.g.input(k,label,kind,**kw),role)
 def op(self,k,label,ports=(),kind=S.MAGIC,chemical_state=None,**kw):
  y=self.g.product(k+'_product','Source preparation: '+label+'; whole occupied molecules/yield/sideproducts unassembled',K.MOLECULAR,chemical_state=chemical_state)
  if kind==S.MAGIC:kw.setdefault('reason','Positive source operation; detailed molecular products/occupancy unassembled, not source-absence evidence.')
  if kind==S.GENERIC:kw['changes_molecule']=False
  self.g.react(k,label,kind,((self.x,R.TEMPLATE),)+tuple(ports),(y,),**kw);self.x=y;return y
 def clean(self,k,label,ports=()):
  y=self.g.product(k+'_kept','Actual retained source workup fraction',K.MOLECULAR);w=self.g.product(k+'_excluded','Actual excluded source fraction; purity not asserted',K.WASTE)
  self.g.react(k,label,S.CLEANUP,((self.x,R.TEMPLATE),)+tuple(ports)+(self.inp(k+'_condition',label),),(y,w));self.g.discard(w,'Actual excluded workup fraction; no invented molecular partition or complete purification');self.x=y;return y
 def aliquot(self,k,label,qc=False):
  y=self.g.product(k+'_selected','Actual allocation: '+label,K.QC if qc else K.MOLECULAR);r=self.g.product(k+'_remainder','Untransferred source stock; actual amount/later fate unassigned',K.MOLECULAR)
  self.g.react(k,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r))
  self.g.mark_terminal(Terminal(y.id if qc else r.id,Disposition.QC if qc else Disposition.UNRESOLVED,label if qc else 'Untransferred stock not silently consumed/discarded'))
  self.x=r if qc else y;return y
 def library(self,k,label,platform,reads):
  y=self.g.product(k+'_library',label+'; whole library ends/primerlanding unassembled',K.LIBRARY)
  self.g.react(k+'_loading','Actual sequencing loading contribution',S.GENERIC,((self.x,R.TEMPLATE),),(y,),changes_molecule=False)
  self.g.sequence(y,SequencingPlan(read_lengths=reads,limitation=platform+'; source run lengths only; whole adapter/index/primerlanding geometry unassembled.'));return y

def intended_atail_scene():
 s=Scene();s.strand('substrate',[Segment('unknown_source_fragment','[genomic fragment]',placeholder=True,inferred=True),Segment('intended_source_3A','A',inferred=True,note='Source dA-tailing class; actual occupied end/whole fragment context unassembled')])
 s.ligation_end('substrate','intended_source_3A',kind='cohesive',polarity='3-prime',owner='substrate',text='INFERRED occupied substrate view: source intended3A class only; adapter body/overhang and actual occupancy unassembled')
 return s

def contact_geometry():
 return AdapterLigationGeometry('unavailable','unavailable',(unknown_scene(),),'Source physical contact joining; actual genomic partners/occupied end geometry unassembled in this selected declaration. No source-absence claim or borrowed Hi-C junction.')

def adapter_geometry():
 return AdapterLigationGeometry('unavailable','unavailable',(intended_atail_scene(),),'Source dA-tailing establishes intended substrate3A class; actual occupied fragment/adapter bodies and junction geometry unassembled in this selected declaration. No source-absence claim.')

AUTHORITY='Published2020 PMC7329539 mainMethods141–143/Results67–71 and actualFig1PDFp2/Fig3p5; publishedMOESM1DOCXcaptions/roles. Sourcequalified HiCoP contact/column route only, not CoP-only Tn5/HiSeq or unpublished/preprint defaults.'

def protocol_graph():
 g=ProtocolGraph('HiCoP source-qualified K562 accessible-contact library',authority=AUTHORITY,representation='published_instance')
 cells=g.input('crosslinked_cells','Actual source-prepared10^7 crosslinkedK562cells; upstream fixation/lysis allocation unassembled for selectedHiCoP',K.BIOLOGICAL,primary=True,chemical_state={'source':'HiCoPMethods143crosslinkedcells','upstream':'No silentCoP2percent ormouse1percentfixationrecipe inheritance'})
 b=T(g,cells)
 g.add_source_item(SourceItem('selection_boundary','Main67–71/Fig1/3 distinguishes OCEAN-C FAIRE enrichment from HiCoP silica-column purification of protein-free sonicated contact fragments.','Physical selection boundary retained; no exclusiveaccessiblecontact specificity claimed. OCEAN-C originalunread; comparison qualified to definingHiCoP account, no separateCoP-only count.'))
 g.add_source_item(SourceItem('branch_scope','SelectedHiCoPMethods143: NovaSeqPE150/13–19PCRcycles/0.55–0.75xbeads.','PrecedingCoP-seqNGSsection Tn5/5cycleqPCR/HiSeq alternatives not merged; no actualHiCoPqPCRaliquot orTn5addition invented.'))
 b.op('mboi_digest','ActualMboIdigestion of crosslinkedcells; dose/endoccupancy unassembled',(b.inp('mboi','ActualMboI enzyme',K.CATALYST,R.CATALYST),),kind=S.DIGESTION)
 supplies=g.input('fill_supplies','Aggregate actual source fill-in supplies; number/stockpackaging/doses unassembled',K.MOLECULAR)
 pol=g.product('fill_polymerase','Actual polymerase catalyst class performing sourcefill-in; identity/dose unassembled',K.CATALYST)
 bio=g.product('biotin_dctp','Actual source biotin-labelled dCTP nucleotide contributor; attachment atom/linker unassigned',K.MOLECULAR,chemical_state={'nucleotide':'biotin-labelled dCTP','meaning':'Affinity handle; not sampleindex/UMI','attachment':'Source name alone does not define exact modified atom/linker'})
 other=g.product('fill_other','Other actual fill-in molecular supplies; exactNTmix/buffers unassembled',K.MOLECULAR)
 g.react('fill_inventory','Component inventory view, not documented premix or extra fresh doses',S.GENERIC,((supplies,R.TEMPLATE),),(pol,bio,other),changes_molecule=False)
 b.op('fill','ActualbiotindCTPfill-in; localcovalentincorporationclass only',((pol,R.CATALYST),(bio,R.STRUCTURAL),(other,R.STRUCTURAL)),kind=S.EXTENSION,chemical_state={'DNA':'Intendedsourcebiotinylatedfilledcontact-fragmentclass; occupiedends/wholeproducts/unincorporatedNTclassesunassembled'})
 current=g.product('joining_template','Actualcontactfragmentend inventoryview within same filledpreparation',K.MOLECULAR);partners=g.product('joining_partners','Othercopresentcontactfragmentends from same preparation; not fresh addedDNA',K.MOLECULAR)
 g.react('joining_inventory','Bookkeeping copresentendclasses; no physicalsplit/purification',S.GENERIC,((b.x,R.TEMPLATE),),(current,partners),changes_molecule=False);b.x=current
 b.op('contact_join','ActualT4ligasecontactjoining',((partners,R.STRUCTURAL),b.inp('t4_ligase','ActualT4DNA ligase;doseunassigned',K.CATALYST,R.CATALYST)),kind=S.LIGATION,geometry=contact_geometry())
 b.clean('contact_wash','Actualsourcewashbefore400ulTEsonication; detailedworkupunassembled')
 b.op('sonicate','400ulcoldTE;20min15son25soff/Q800R2/60percentamplitude',(b.inp('sonication_conditions','Sourceactualmechanicalprogramme'),),chemical_state={'DNA':'Sonicatedcontactfragments withbiotin class; exactfragmentends/breakpositionsunassembled'})
 selected=g.product('column_input','Actual360ulremainingchromatin forCoPcolumn',K.MOLECULAR);saved=g.product('saved_input','Actual40ulsparedinputfrom400ulsonicate; separatecontrolstock',K.MOLECULAR)
 g.react('input_split','Actual40ulinput/360ulcolumnallocation',S.SPLIT,((b.x,R.TEMPLATE),),(selected,saved));g.mark_terminal(Terminal(saved.id,Disposition.UNRESOLVED,'Savedsourceinputcontrolprocessedseparately; notpooled intoselectedHiCoPlibrary'));b.x=selected
 column=g.input('cop_column','ActualPCRpurificationcolumn silica selection;sourceDongShengN1093, not generic DNAcleanupsubstitution',K.SELECTION)
 captured=g.product('column_bound','Actualcolumn-retainedprotein-freecontactfragmentpopulation; specificitynotabsolute',K.MOLECULAR);flow=g.product('column_nonselected','Actualnonselectedprotein-bound/othercolumnfraction',K.WASTE)
 g.react('cop_select','SourceCoPsilica-columnselection BEFORE decrosslinking',S.SELECTION,((b.x,R.TEMPLATE),(column,R.SELECTION)),(captured,flow));g.discard(flow,'Actualnonselectedcolumnfraction; noperfectaccessibilityclaim');b.x=captured
 b.clean('column_elute','Actual40ulEBelution; retainedcolumnmaterial excluded, not magicallydecrosslinked')
 b.op('rnase','Actual1ul100mgmlRNase A37C30min',(b.inp('rnase_a','ActualRNase A1ul100mgml',K.CATALYST,R.CATALYST),),kind=S.DIGESTION)
 b.op('pk','Actual5ul10mgmlProteinaseK55C30min',(b.inp('proteinase_k','ActualProteinaseK5ul10mgml',K.CATALYST,R.CATALYST),),kind=S.DEPROTEINIZATION)
 b.op('reverse_crosslinks','Actual65Covernight afterPK',(b.inp('reverse_conditions','Actualsourceovernight65C'),))
 b.clean('dna_purify','ActualPCRcolumnDNApurification/elute50ulEB')
 b.aliquot('dna_qc','ActualQubitDNAquantificationallocation;volumeunassigned',True)
 raw=g.input('c1_raw','Actual5ul10mgmlC1streptavidinbeads',K.MOLECULAR)
 washed=g.product('c1_washed','ActualC1beadsafter1x200ulTWBwash',K.MOLECULAR);waste=g.product('c1_prewash_waste','ActualTWBprewashsupernatant',K.WASTE)
 g.react('c1_prewash','Sourcebeadprewash200ulTWB',S.CLEANUP,((raw,R.TEMPLATE),b.inp('c1_prewash_conditions','TWB5mMTris7.5/.5mMEDTA/1MNaCl/.05percentTween20')),(washed,waste));g.discard(waste,'Actualprewashsupernatant')
 support=g.product('c1_prepared','ActualwashedC1beads resuspended50ul2xbindingbuffer',K.SELECTION,chemical_state={'support':'SourceC1streptavidin','binding':'Biotin-affinity, not covalentbiotinbondcleavage'})
 g.react('c1_resuspend','Actual50ul2xbindingbuffer10Tris7.5/1EDTA/2MNaCl',S.GENERIC,((washed,R.TEMPLATE),b.inp('c1_binding_buffer','Actual50ul2xbindingbuffer')),(support,),changes_molecule=False)
 bound=g.product('biotin_bound','ActualC1boundbiotinDNApreparation;supportremainscopresent',K.MOLECULAR);unbound=g.product('biotin_nonselected','ActualunboundC1capturesupernatant',K.WASTE)
 g.react('biotin_capture','AddpreparedbeadstoactualDNAstock;RT15minrotation',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),(support,R.SELECTION)),(bound,unbound));g.discard(unbound,'Actualnonselectedcapturesupernatant, notperfectbiotinpurity');b.x=bound
 b.clean('capture_wash','Actual200ulTWB55C2minmix/magnetdiscard/repeat;retainedbeadDNA')
 finish=g.input('finishing_supplies','Actualaggregateend-repair/dA-tail/adapterfinishinginventory;packaging/order/dosesunassembled, notprovenpremix',K.MOLECULAR)
 enz=g.product('end_repair_enzymes','Actualsourceend-repair/dA-tailcatalystclasses;identitiesunassembled',K.CATALYST);nts=g.product('end_repair_nt','Actualend-repair/dA-tailnucleotideclasses;exactmixunassembled',K.MOLECULAR);ad=g.product('adapters','Actualsourceadapteroligoinventory;wholebodies/modificationsunassembled',K.OLIGO);rest=g.product('end_repair_other','Otheractualfinishingreagents;stockpackagingunassembled',K.MOLECULAR)
 g.react('finishing_inventory','Bookkeepactualfinishingcomponents;nofreshduplicatedoses',S.GENERIC,((finish,R.TEMPLATE),),(enz,nts,ad,rest),changes_molecule=False)
 b.op('end_repair_atail','Sourceend-repair/dA-taildelegatedin-situHiCrecipe;beadDNA retained',((enz,R.CATALYST),(nts,R.STRUCTURAL),(rest,R.STRUCTURAL)),chemical_state={'DNA':'Intendedrepaired/dA-tailedfragmentclass;actualendoccupancy/wholeproductunassembled'})
 b.op('adapter_ligation','SourceadapteradditiontoDNAfragments;endgeometryunassembled',((ad,R.STRUCTURAL),b.inp('adapter_ligation_supplies','Actualadapterjoiningmolecularsupplies;enzymestockdetailsunassembled',K.MOLECULAR,R.STRUCTURAL)),kind=S.LIGATION,geometry=adapter_geometry())
 pcr=g.input('pcr_supplies','ActualPCRaggregateinventory;stocknumber/packagingunassigned',K.MOLECULAR)
 f=g.product('pcr_forward','ActualIlluminaprimerforwardfunctionalclass;body/role/indexlayoutunassembled',K.OLIGO);r=g.product('pcr_reverse','ActualIlluminaprimerreversefunctionalclass;body/role/indexlayoutunassembled',K.OLIGO);cat=g.product('pcr_enzyme','ActualPCRpolymeraseclassfromsourcesupplies;identityunassigned',K.CATALYST);nt=g.product('pcr_nt','ActualPCRnucleotideclassfromsourcesupplies;dosesunassigned',K.MOLECULAR);buff=g.product('pcr_other','OtheractualPCRsupplies;packagingunassembled',K.MOLECULAR)
 g.react('pcr_inventory','ActualPCRcomponents inventory;notnewindependentadditions',S.GENERIC,((pcr,R.TEMPLATE),),(f,r,cat,nt,buff),changes_molecule=False)
 b.op('pcr','Source13–19cyclesPCRwithIlluminaprimers',((f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST),(nt,R.STRUCTURAL),(buff,R.STRUCTURAL)),kind=S.PCR,geometry=PCRGeometry.unavailable(unknown_scene(),reason='Wholeprimers/landing/index/moleculesunassembled; no canonicalkitbodyinheritance.'))
 b.clean('size_select','Actual0.55–0.75xVAHTSDNACleanbeadselection;keptsizefraction/unselectedfractionsexplicit',ports=(b.inp('size_beads','ActualVAHTSDNACleanbeadcontribution',K.MOLECULAR,R.STRUCTURAL),))
 b.aliquot('library_qc','ActualsourceQubitlibraryQCallocation;amountunassigned',True)
 b.aliquot('loading_allocation','Actualsequencingloadingallocation;amountunassigned')
 b.library('hicop','ActualHiCoPlibrary','IlluminaNovaSeq',{'Read1':150,'Read2':150})
 return _publication_coverage(g)

def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();g=protocol_graph()
 c('BiotindCTPfromactualsupplies','biotin_dctp'in g.steps['fill_inventory'].outputs,True)
 c('NoextraexternalbiotindCTP',not any('biotin_dctp'in s.outputs for s in g.steps.values()),False)
 c('Joiningpartnerssamepreparation','joining_partners'in g.steps['joining_inventory'].outputs,True)
 c('Columnuses360ulselectedallocation',g.steps['cop_select'].inputs[0].material,'column_input')
 c('Inputcontrolnotdiscarded',g.terminals['saved_input'].disposition,Disposition.UNRESOLVED)
 c('RNaseafterphysicalcolumneluate',g.steps['rnase'].inputs[0].material,'column_elute_kept')
 c('PKbeforecrosslinkreversal',g.steps['reverse_crosslinks'].inputs[0].material,'pk_product')
 c('CaptureusespostQCstock',g.steps['biotin_capture'].inputs[0].material,'dna_qc_remainder')
 c('C1supportpreparedfromactualwashedbeads',g.steps['c1_resuspend'].inputs[0].material,'c1_washed')
 c('FinishingusesretainedbeadDNA',g.steps['end_repair_atail'].inputs[0].material,'capture_wash_kept')
 c('Knownadapterinventoryconsumed',any(i.material=='adapters'for i in g.steps['adapter_ligation'].inputs),True)
 c('PCRbothfunctionalprimersconsumed',all(any(i.material==p for i in g.steps['pcr'].inputs)for p in ('pcr_forward','pcr_reverse')),True)
 c('Bothprimersproducedsameinventory',all(p in g.steps['pcr_inventory'].outputs for p in ('pcr_forward','pcr_reverse')),True)
 c('NoCoPqPCRaddition',any('qpcr'in k for k in g.steps),False)
 c('NoTn5branch',any('tn5'in k for k in g.materials),False)
 c('Noimaginedphysicalpool',any(s.kind==S.POOL for s in g.steps.values()),False)
 c('Loadingremainderexplicit','loading_allocation_remainder'in g.terminals,True)
 c('SourcePE150structured',g.sequencing['hicop_library'].read_lengths,{'Read1':150,'Read2':150})
 p=sources_root()/HERE.parent.name/'s13072-020-00348-6_PMC7329539.xml.txt'
 if have(c,Source(p,'https://doi.org/10.1186/s13072-020-00348-6')):
  t=p.read_text();c('Actualsourceinputallocation','40 μl of chromatin was spared as input'in ' '.join(t.split()),True);c('ActualsourcePCRrange','13 to 19 cycles'in t,True);c('ActualsourceNovaSeq','Novaseq-PE150'in t,True)
 c('SourceintendedAendretainedwithoutknownadapter',[(e.owner,e.polarity)for e in intended_atail_scene().ligation_ends],[('substrate','3-prime')])
 c('Occupiedsubstrateviewinferredbydefault',all(seg.inferred for st in intended_atail_scene().strands.values() for seg in st.segs),True)
 c('Publication graph validates with declared limitations',g.publication_problems(),[]);c.report()


def _cover_reviewed_targets(g, groups):
 """Only the literal curator-reviewed targets below receive these source assignments."""
 for source,targets in groups:
  for target in targets:
   label=g.materials[target].label if target in g.materials else g.steps[target].label
   key='reviewed_'+target
   g.add_source_item(SourceItem(key,label,source))
   g.cover_source_item(SourceCoverage(key,'modelled',target))

PAGE_FILE='hicop.html'
SOURCE_TARGETS=(
 ('Published HiCoP2020 mainMethods143/Fig3PDFp5: preparedcrosslinkedcells/MboI/biotindCTPfill/T4joining/wash/400ulTEsonication and40ulinput/360ulcolumnallocation. Ingredient inventories are bookkeeping; full ends unassembled.',
  ('crosslinked_cells','mboi','mboi_digest','fill_supplies','fill_inventory','fill','t4_ligase','joining_inventory','contact_join','contact_wash_condition','contact_wash','sonication_conditions','sonicate','input_split')),
 ('Same publishedMethods143 and Results67–71/Fig1PDFp2: actual silica selection before de-crosslinking,40ulEB/RNase A/PK/65Cover night/DNAcleanup50ulEB/Qubit.',
  ('cop_column','cop_select','column_elute_condition','column_elute','rnase_a','rnase','proteinase_k','pk','reverse_conditions','reverse_crosslinks','dna_purify_condition','dna_purify','dna_qc')),
 ('Same publishedMethods143:5ulC1/200ulprewash/50ulbindingbuffer/RT15mincapture/repeatedTWB55C2minwash. Whole capture occupancy unassembled.',
  ('c1_raw','c1_prewash_conditions','c1_prewash','c1_binding_buffer','c1_resuspend','biotin_capture','capture_wash_condition','capture_wash')),
 ('Same publishedMethods143: endrepair/dAtail/adapteraddition/13–19cyclesPCR/0.55–0.75xbeads/Qubit/NovaSeqPE150. Actual supply classes and inferred intended3Aprimeend only; ref17kitrecipes/wholeoligos not inherited.',
  ('finishing_supplies','finishing_inventory','end_repair_atail','adapter_ligation_supplies','adapter_ligation','pcr_supplies','pcr_inventory','pcr','size_beads','size_select_condition','size_select','library_qc','loading_allocation','hicop_loading','hicop_library')),
)

def _publication_coverage(g):
 g.cover_source_item(SourceCoverage('selection_boundary','modelled','cop_select'))
 g.cover_source_item(SourceCoverage('branch_scope','modelled','hicop_library'))
 _cover_reviewed_targets(g,SOURCE_TARGETS)
 return g.require_publishable()

def render_page():
 from page import head,caveat,table
 from chemdraw import panel
 from protocol_views import assemble_page,WorkflowSection
 g=protocol_graph()
 return assemble_page(g,[head('HiCoP — column-selected accessible chromatin contacts'),
 '<main class="wrap" data-protocol-complete="true"><h1>HiCoP</h1><p>HiCoP combines contact ligation with silica-column enrichment of protein-free chromatin fragments. Biotin capture then retains contact-library DNA for adapter addition, PCR and paired-end sequencing.</p><p><a href="https://doi.org/10.1186/s13072-020-00348-6">Zhang et al., Epigenetics &amp; Chromatin (2020)</a> · <a href="01_hicop.html">Research note and evidence scope</a></p>',
 caveat('<strong>Partial molecular declaration from the published HiCoP procedure.</strong> Complete genomic contact products, occupied ends, exact adapter/primer bodies and read-primer landing sites remain unassembled. Column and affinity selection do not establish exclusive regulatory contacts. The intended substrate A-tail view below is inferred in an unresolved genomic context.'),
 '<h2>Preparation</h2><ol><li>Digest 10⁷ source-prepared crosslinked cells with MboI, fill with biotin-labelled dCTP and join contacts with T4 ligase.</li><li>Sonicate in 400 µl TE. Save 40 µl as input and column-select the remaining 360 µl before de-crosslinking.</li><li>Elute 40 µl, treat with RNase A and proteinase K, reverse crosslinks and purify DNA into 50 µl EB; retain an actual Qubit allocation.</li><li>Capture DNA on prepared C1 streptavidin beads and wash. Retain bead-associated DNA for end repair/dA-tailing and adapter addition.</li><li>Amplify 13–19 cycles, select with 0.55–0.75× VAHTS beads, quantify and allocate the library for sequencing.</li></ol>',
 '<h2>Known molecular parts</h2>',table(('Part','Role and limit'),(('Biotin-labelled dCTP','Source fill-in nucleotide and affinity handle; exact modified atom/linker and occupied sites unassembled.'),('Silica column','Physical enrichment occurs before crosslink reversal. Excluded material and saved input remain separate.'),('C1 streptavidin beads','5 µl at 10 mg/ml, prewashed and resuspended in 50 µl binding buffer; retained DNA proceeds to library finishing.'),('Illumina primer pair','Both PCR functional classes are explicit; whole bodies, sample-index architecture and occupied sites remain unresolved.'))),
 panel(intended_atail_scene().rows(),caption='INFERRED — intended source dA-tailed substrate class; occupied genomic fragment and adapter geometry unassembled'),
 '<h2>Read setup</h2><p>Illumina NovaSeq, paired-end 150 bases. Exact read-primer sites, index cycles and complete library molecules remain unresolved. No physical sample pool is invented.</p>',
 '<h2>Source boundary</h2><p>The defining article contrasts OCEAN-C’s FAIRE selection with HiCoP’s silica-column selection. This comparison is qualified to the HiCoP article; OCEAN-C’s full original was not reviewed. CoP-only accessibility, its Tn5/HiSeq alternatives and its five-cycle/qPCR workflow are not merged into the selected HiCoP preparation. The upstream fixation allocation and the detailed cited in-situ Hi-C finishing recipe remain unresolved; no canonical kit oligos are supplied by name.</p>',
 '<details><summary>Full molecular workflow, physical fractions and unresolved chemistry</summary>',WorkflowSection('Source-qualified molecular workflow'),'</details></main>'])

if __name__=='__main__':selftest()
