"""Private source-qualified partial preparation; no full source/absence audit."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem
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
  y=self.g.product(k+'_selected','Actual allocation: '+label,K.MOLECULAR);r=self.g.product(k+'_remainder','Untransferred source stock; actual amount/later fate unassigned',K.MOLECULAR)
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
 bio=g.product('biotin_dctp','Actual source biotin-labelleddCTP nucleotide contributor; attachment atom/linker unassigned',K.MOLECULAR,chemical_state={'nucleotide':'biotin-labelleddCTP','meaning':'Affinity handle; not sampleindex/UMI','attachment':'Source name alone does not define exact modified atom/linker'})
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
 b.op('rnase','Actual1ul100mgmlRNaseA37C30min',(b.inp('rnase_a','ActualRNaseA1ul100mgml',K.CATALYST,R.CATALYST),),kind=S.DIGESTION)
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
 return g

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
 c('Privateincompletegraph',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
