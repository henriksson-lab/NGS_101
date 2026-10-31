from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, feature, workflow_from_protocol_graph, workflow_panel
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan
from page import head,info
PAGE_FILE='chem.html'

def operations(g):
    def op(key,x,label,reagent=None,state=None,kind=S.GENERIC,reagent_kind=K.CONDITION,reagent_role=R.CONDITION):
        inputs=[(x,R.TEMPLATE)]
        if reagent:
            z=g.input(key+'_reagent',reagent,reagent_kind);inputs.append((z,reagent_role))
        y=g.product(key+'_product',label+'; preparation includes unremoved reaction components',K.MOLECULAR,state=state)
        kw={'reason':'Actualsourceoperation; wholeoccupiedgeometry/stoichiometryunassembled.'} if kind==S.MAGIC else {}
        if kind==S.GENERIC:kw['changes_molecule']=False
        g.react(key,label,kind,inputs,(y,),**kw);return y
    def clean(key,x,label):
        z=g.input(key+'_selection',label,K.SELECTION)
        y=g.product(key+'_kept','Physically retained source fraction: '+label,K.MOLECULAR)
        w=g.product(key+'_removed','Physically removed fractions: '+label,K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(z,R.SELECTION)),(y,w));g.discard(w,'Source physical workup only; no perfect specificity or purity assertion.');return y
    return op,clean

def end(g,x,label):g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=label))
def incomplete(g):
    g.declare_reviewed_scope(('targets','foundation','forward_cololinker','reverse_cololinker','hybridize','chip','foundation_mix_supplies','foundation_mix','foundation_allocate','foundation_parts','seed','seed_wash_selection','seed_wash','universal_bridge','bridge_allocation','cycle1_bsa','cycle2_bsa','cycle3_bsa','cycle4_bsa','cycle5_bsa','cycle6_bsa','cycle1_digest','cycle2_digest','cycle3_digest','cycle4_digest','cycle5_digest','cycle6_digest','ngs_pair_supplies','ngs_pair_parts','ngs_adapter','ngs_bridge','ngs_mix_supplies','ngs_prepare','ngs_allocate','ngs_parts','ngs_ligate','adapter_wash_selection','adapter_wash','phix','ht1','cartridge_control','ready'),SOURCE+'; literalreviewedselectedpreparation/cyclefunctions only, notwholeassignment/junction/sourceaudit.')
    g.add_source_item(SourceItem('unassembled','Unreviewed linked originals and full molecular geometry',LIMIT));return g


from dataclasses import dataclass,replace
from protocol_graph import Terminal,Disposition

def join(g,key,inputs,label):
    y=g.product(key+'_product',label+';unremoved companion species remain',K.MOLECULAR)
    g.react(key,label,S.GENERIC,tuple((x,R.TEMPLATE if i==0 else R.STRUCTURAL) for i,x in enumerate(inputs)),(y,),changes_molecule=False);return y

def render_page():
    return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_chem.html">Research notes</a></p>'+workflow_panel(workflow_from_protocol_graph(protocol_graph(),ENDPOINT))+'</main>'

from illumina import P5_RC
TITLE='ProtSeq — demonstrated chip barcode-cycle sequencing'
SOURCE='Published PMC8717419 Results262–278,STARMethods588–609/643–657/669–673;publishedSI TableS2/S4 andFigS3'
LIMIT='Private incomplete model of demonstrated DNA-DNA barcode-cycle NGS. This is not complete amino-acid/protein sequencing. Peptide–nanobody capture and TargetSwitchSELEX are separate demonstrations. Source-selected populations,not perfect cognate binding or six successful ligations per molecule. Complete foundation/cololinker/binder junction geometry remains unassembled.'
ENDPOINT='ready_product'
FOUNDATION='CGACTGCGAGCTGATGGCCTTGATGATAACG'
NGS_PREFIX='AGATCGGAAGAGCGTCGTGTAGGGAAAGA'
NGS_ADAPTER=NGS_PREFIX+P5_RC
BRIDGE_PREFIX='CTGCGCCTATAGGAATTCGTTATC';BRIDGE_SUFFIX='GGACACGGCCGTTATC'
@dataclass(frozen=True)
class UniversalBridge:
    prefix:str=BRIDGE_PREFIX
    universal_residues:tuple=('5-Nitroindole',)*12
    suffix:str=BRIDGE_SUFFIX
    def notation(self):return self.prefix+'/i5NitInd/'*len(self.universal_residues)+self.suffix
    def __post_init__(self):
        if self.universal_residues!=('5-Nitroindole',)*12:raise ValueError('Source12actual5NitIndresidues,not12randomDNA or inosines')
@dataclass(frozen=True)
class ObservedChain:
    recorded:tuple=()
    def record_observed_event(self,binder,cycle):
        if not binder or type(cycle)is not int or cycle<1 or (self.recorded and cycle<=self.recorded[-1][1]):raise ValueError('Explicit increasing source observation needed')
        return replace(self,recorded=self.recorded+((binder,cycle),))
def foundation():
    sc=Scene();sc.strand('foundation',[Segment('common source foundationbody','CGACTGCGAGCTGATG'),Segment('sourceFd7identity','GCCTTGAT',feature=feature('protseq.fd7','feature_barcode','fixed',note='Construct/target foundation identity,not cell or UMI.')),Segment('source foundationtail','GATAACG')],mod5='phosphate',label='Source TableS2 Fd7 example;exact experimental target assignment unassembled');return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance');op,clean=operations(g)
    def aliquots(k,stock,labels,remainder):
        parts=[g.product(k+'_'+str(i+1),label,stock.kind,state=stock.state,chemical_state=stock.chemical_state) for i,label in enumerate(labels)]
        rest=g.product(k+'_rest',remainder,stock.kind,state=stock.state,chemical_state=stock.chemical_state)
        g.react(k,'Actualsourceallocation toindependentportions, notcopyingstock',S.SPLIT,((stock,R.TEMPLATE),),parts+[rest]);g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Untransferredstock, possiblyzero/handlingloss; quantityunassigned.'));return parts
    def transformed(k,x,label,ports):
        y=g.product(k+'_product',label+'; allcleavage/unligated/reagentcompanions retaineduntilactualwash',K.MOLECULAR)
        g.react(k,label,S.MAGIC,((x,R.TEMPLATE),)+tuple(ports),(y,),reason='ActualsourceDNA-changingoperation; exactoccupiedsubstrate/junction/cleavagegeometryunassembled, no perfectevent/yieldclaim.');return y
    targets=g.input('targets','Actualsix-cycleDNA-target populations; individualmemberassignmentunexpanded',K.OLIGO,primary=True)
    fd=g.input('foundation','Actualsix-cyclefoundationpool; wholemember/targetassignmentunassembled. Fd7helper isONLYsourceinventoryreference, notassertedusedinmeasuredsixcyclerun.',K.OLIGO,chemical_state={'source_identity':feature('protseq.foundation','feature_barcode','unknown',note='Actualfoundationidentity; memberassignmentunassembled, noFd7 experimentalassignment')})
    fc=g.input('forward_cololinker','SourceFC wholepolymerunassembled',K.OLIGO);rc=g.input('reverse_cololinker','SourceRC wholepolymerunassembled',K.OLIGO)
    x=join(g,'hybridize',(targets,fd,fc,rc),'Compoundactualone-target-per-wellhybridizations BEFOREmixing,95→25Cstagedcool;10nMFC/othersinexcess. Independentmembercontributionsunexpanded, notpostreactioncopies.')
    chip=g.input('chip','ActualMiSeqflowcell surfaceoligos; historicalwholegeometryunreviewed',K.MOLECULAR)
    supplies=g.input('foundation_mix_supplies','Actual10ul2xBluntTAMM+HybridizationBuffer to100ulFoundationMix; sourceequalworkingstocks500pM→120pMfinal, exactmembershipunexpanded',K.MOLECULAR)
    prepared=g.product('foundation_mix_stock','Actual100ulFoundationMix fromcompletedsourcehybridizationarray andactualsupplies',K.MOLECULAR)
    g.react('foundation_mix','Sourceequalhybridizedtargetfoundationcontributions mixedwithligationsupplies',S.GENERIC,((x,R.TEMPLATE),(supplies,R.STRUCTURAL)),(prepared,),changes_molecule=False)
    portions=aliquots('foundation_allocate',prepared,['Actualfirst30ulFoundationMixload','Actualsecond30ulFoundationMixload'],'UntransferredFoundationMix; nominal40ulbeforeunquantifiedhandlingloss')
    dna=g.product('foundation_loaded_dna','CopresentDNAcomponents ofsameactualtwo30ulloads',K.MOLECULAR);cat=g.product('foundation_loaded_ligase','CopresentactualT4ligaseofsameFoundationMixloads',K.CATALYST);buf=g.product('foundation_loaded_buffer','OthercopresentFoundationMixsupplies/28C15min',K.CONDITION)
    g.react('foundation_parts','Components oftwoactualsequentialloads, notanextra enzymeaddition',S.GENERIC,((portions[0],R.TEMPLATE),(portions[1],R.STRUCTURAL)),(dna,cat,buf),changes_molecule=False)
    x=transformed('seed',chip,'ActualT4ligation oftetheredfoundationclasses,28C15min',((dna,R.STRUCTURAL),(cat,R.CATALYST),(buf,R.CONDITION)))
    x=clean('seed_wash',x,'Source100percentformamidewash40C90s; retainchipclasses')
    blocking=g.input('block_supplies','Actualsourceblockingbuffer/ChipBlockingSolution inventory; exactcomponentbodies/packagingunassembled',K.MOLECULAR)
    blockers=g.product('block_oligos','ActualsourceblockingDNAfunctions; fullbodiesunassembled',K.OLIGO);blockother=g.product('block_other','Otheractualblockingmaterial; fullrecipeunassembled',K.MOLECULAR)
    g.react('block_parts','Actualblockinginventoryfunctions, notextra freesupplies',S.GENERIC,((blocking,R.TEMPLATE),),(blockers,blockother),changes_molecule=False)
    x=join(g,'block',(x,blockers,blockother),'Sourceblocksurface/exposedDNA,37C15min; occupancyunassigned')
    bridge=g.input('universal_bridge','Sourceordered12actual5NitInd bridge stock, notN12/inosine; wholevisualpairingunassembled',K.OLIGO,state=UniversalBridge(),chemical_state=UniversalBridge())
    bridge_parts=aliquots('bridge_allocation',bridge,['Actualsource bridgeportion forround'+str(i) for i in range(1,7)],'Untransferredsamebridge stock; amountsunassigned')
    for cycle,portion in enumerate(bridge_parts,1):
        name='cycle'+str(cycle)
        binder=g.input(name+'_binder','Actualround-specificDNA-binder/barcodepopulation; wholemember/junctionunassembled',K.OLIGO)
        mix=join(g,name+'_bridge',(binder,portion),'Source1:2binder/bridgeannealing95C5minthenRT1h; separateactualroundportion, notwhole-stockreuse')
        bsa=g.input(name+'_bsa','ActualsourceBSA AFTERannealing, BEFOREchipincubation; final100ugml',K.MOLECULAR)
        mix=join(g,name+'_bsa_add',(mix,bsa),'ActualpostannealingBSAaddition tosourcebinder solution')
        x=join(g,name+'_bind',(x,mix),'ActualDNA-DNAbinding RT30min; specificity/occupancyunassigned')
        x=clean(name+'_wash',x,'Sourcewashout ofunboundbinders, retainchip-associatedpopulation')
        ligstock=g.input(name+'_ligase_stock','Actualsource2xT4BluntTAligasemix addition; exactDNA-sixcycledosesunassigned, noSpotTagdoseinheritance',K.MOLECULAR)
        lig=g.product(name+'_ligase','CopresentactualT4ligasecomponent',K.CATALYST);other=g.product(name+'_ligase_other','OtheractualsourceMMsupplies; noATPidentity/doseinferredfromunreadguide',K.MOLECULAR)
        g.react(name+'_ligase_parts','Sameactualroundligasemixfunctions, noextrafreshdoses',S.GENERIC,((ligstock,R.TEMPLATE),),(lig,other),changes_molecule=False)
        x=transformed(name+'_ligate',x,'ActualsourceT4ligation retainsoriginalDNAbodies, no perfectincorporation',((lig,R.CATALYST),(other,R.STRUCTURAL)))
        eco=g.input(name+'_ecori','ActualsourceEcoRI catalyst; selectedDNA-sixcycledoseunassigned',K.CATALYST)
        restriction=g.input(name+'_restriction_bridge','Knownsourcecomplementaryrestrictionbridge, fulloccupiedbindinggeometryunassembled',K.OLIGO,state=Construct([Segment('source restrictionbridge','CTGCGCCTATACGAATTCGTTATC')]))
        conditions=g.input(name+'_digest_buffer','ActualsourceCutSmart/digestionconditions; noSpotTagspecificvolumeimport',K.CONDITION)
        digested=g.product(name+'_digest_product','CleavedbinderDNA plusretainedchipbarcodeclasses andsourcecompanions; no aminoacidremoval',K.MOLECULAR)
        g.react(name+'_digest','ActualEcoRI DNAcleavage, releasedDNA remainsuntilwash',S.DIGESTION,((x,R.TEMPLATE),(eco,R.CATALYST),(restriction,R.STRUCTURAL),(conditions,R.CONDITION)),(digested,));x=digested
        x=clean(name+'_release_wash',x,'Sourceformamide/bufferwash physicallyremovesreleasedDNAclasses; retainedchainyieldunassigned')
    pairs=g.input('ngs_pair_supplies','Actualsource1.5ul2:1 1uMUniversalNGSadapter/bridgepreparation; exactconcentrationconventionunassigned',K.MOLECULAR)
    ad=g.product('ngs_adapter','Copresentknown5P universalNGSadaptor',K.OLIGO,state=adapter_scene());abr=g.product('ngs_bridge','CopresentknownNGSadapterbridge9/5',K.OLIGO,state=Construct([Segment('source bridge','TTCCGATCTCGTTA')]))
    g.react('ngs_pair_parts','Actualsameadapter/bridge inventory, notnewfreeadditions',S.GENERIC,((pairs,R.TEMPLATE),),(ad,abr),changes_molecule=False)
    supplied=g.input('ngs_mix_supplies','Actual25ul2xBluntTAligaseMM+10ul10xCutSmart+63.5ulwater',K.MOLECULAR)
    stock=g.product('ngs_mix_stock','Actual100ulNGSligationpreparation',K.MOLECULAR)
    g.react('ngs_prepare','Actualsourceadapter/bridge+ligationsupplies preparation',S.GENERIC,((ad,R.TEMPLATE),(abr,R.STRUCTURAL),(supplied,R.STRUCTURAL)),(stock,),changes_molecule=False)
    loads=aliquots('ngs_allocate',stock,['Actualfirst30ulNGSligationload','Actualsecond30ulNGSligationload'],'UntransferredNGSligationstock, nominal40ulbeforehandlingloss')
    loadedDNA=g.product('ngs_loaded_dna','Sameactualsourceadapter/bridgecomponents oftwo30ulloads',K.MOLECULAR);cat=g.product('ngs_loaded_ligase','SameactualT4ligasecomponent ofloadedmixture',K.CATALYST);buf=g.product('ngs_loaded_conditions','Sameactualsourcebuffer;40C2min45s',K.CONDITION)
    g.react('ngs_parts','Loadedmixturefunctions, noextraenzyme/adapterdoses',S.GENERIC,((loads[0],R.TEMPLATE),(loads[1],R.STRUCTURAL)),(loadedDNA,cat,buf),changes_molecule=False)
    x=transformed('ngs_ligate',x,'Actualsourcebridge-assistedNGSadaptorligation40C2min45s',((loadedDNA,R.STRUCTURAL),(cat,R.CATALYST),(buf,R.CONDITION)))
    x=clean('adapter_wash',x,'Actualtwo500ulwaterwashes90sapart; sourcechipDNAretained')
    phix=g.input('phix','Actual20ul20pMdenaturedPhiX DNAcontrol',K.MOLECULAR);ht1=g.input('ht1','Actual580ulHT1buffer',K.CONDITION)
    control=g.product('cartridge_phix','Actualprepared600ulPhiXcartridgecontrolsample; notligatedtobarcodechains',K.QC)
    g.react('cartridge_control','SourcePhiXcontrolloadedincartridgesamplewell, notcovalentchainincorporation',S.GENERIC,((phix,R.TEMPLATE),(ht1,R.CONDITION)),(control,),changes_molecule=False);g.mark_terminal(Terminal(control.id,Disposition.QC,'Separateactualcartridgecontrol; mixing/acquisitionbeyondsourceinstructionsnotinferred.'))
    x=op('ready',x,'Sourcechip directlytoadaptedMiSeqrun; omitsinitialwashing/libraryloading',kind=S.GENERIC)
    end(g,x,'ActualMiSeqV2,45–600cycleexperiment-specific read. Primerlanding/exactcycles/unassembledwholechain unknown. DemonstratedDNA-DNABCScapture, notcompleteaminoacidsequencing/SELEX-PCR.')
    return incomplete(g)

def adapter_scene():
    sc=Scene();sc.strand('adapter',[Segment('sourcebody',NGS_PREFIX),Segment('sourcecanonicalP5reversecomplement',P5_RC)],mod5='phosphate');return sc

def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph();n=ObservedChain().record_observed_event('sourcebinder',1).record_observed_event('sourceother',3)
    c('Observed missingcycle retainedwithout invented incorporation',tuple(z[1] for z in n.recorded),(1,3))
    c('5NitInd is actualmodifiedresidue,not DNA N/I',set(UniversalBridge().universal_residues),{'5-Nitroindole'})
    c('Bridge stocksplit feeds independentroundportions',g.steps['cycle2_bridge'].inputs[1].material,'bridge_allocation_2')
    c('BSAadded afteractualbinder/bridgeannealing',g.steps['cycle1_bsa_add'].inputs[0].material,'cycle1_bridge_product')
    c('Ligationisactualmoleculartransformation',g.steps['cycle1_ligate'].kind,S.MAGIC)
    c('Restrictionbridge isactualDNAcontributor',g.steps['cycle1_digest'].inputs[2].material,'cycle1_restriction_bridge')
    c('PhiXcartridgecontrol separatefrombarcodechain',g.terminals['cartridge_phix'].disposition,Disposition.QC)
    c('Foundation retains5phosphate',foundation().strands['foundation'].mod5,'phosphate')
    c('Foundation identity is not UMI',foundation().strands['foundation'].segs[1].feature.role,'feature_barcode')
    c('Nextcycle receiveswashedretainedchip material',g.steps['cycle2_bind'].inputs[0].material,'cycle1_release_wash_kept')
    c('NGSadapter aftersixsourcecycles',g.steps['ngs_ligate'].inputs[0].material,'cycle6_release_wash_kept')
    c('NGSadapter physical5phosphate',adapter_scene().strands['adapter'].mod5,'phosphate')
    c('No aminoaciddegradation step in demonstrated DNA branch','aminoacid_degrade' in g.steps,False)
    c('Incomplete sourcegeometry staysprivate',g.publishable,False)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'j.isci.2021.103586_PMC8717419.xml.txt'
    if p.exists():
        t=p.read_text();c('Wholemodifieduniversalbridge source transcription',UniversalBridge().notation() in t,True);c('WholeproducedNGSadapter verbatim',''.join(s.top for s in g.materials['ngs_adapter'].state.strands['adapter'].segs) in t,True)
        c('WholeproducedNGSbridge verbatim',''.join(s.top for s in g.materials['ngs_bridge'].state.segments) in t,True)
    else:c.skip('Published originaltranscription','Externaloriginalabsent')
    c.report()
