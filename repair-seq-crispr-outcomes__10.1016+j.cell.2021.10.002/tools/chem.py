from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition

def unknown_scene(label='Unassembled whole molecular product'):
    sc=Scene();sc.strand('unassembled',[Segment('unassembled','[unassembled]',placeholder=True)]);return sc

def add(g,key,label,kind=K.MOLECULAR,role=R.STRUCTURAL,**kw):
    return (g.input(key,label,kind,**kw),role)
def condition(g,key,label):return add(g,key,label,K.CONDITION,R.CONDITION)
def opaque(g,key,label,x,ins=(),kind=S.MAGIC):
    y=g.product(key+'_product',label+'; whole molecular identities/junctions unassembled',K.MOLECULAR)
    kw={'reason':'Source-positive operation; whole kit composition, substrate context or product geometry unassembled, not publication absence.'} if kind==S.MAGIC else {}
    g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(ins),(y,),**kw);return y
def cleanup(g,key,label,x):
    y=g.product(key+'_kept','Source retained preparation; recovery/carryover unassigned',K.MOLECULAR);w=g.product(key+'_removed','Source workup fraction',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),condition(g,key+'_condition',label)),(y,w));g.discard(w,'Physical workup, not perfect molecular selectivity.');return y
def aliquot(g,key,label,x,remainder_label='Untransferred source stock; later fate unassigned'):
    y=g.product(key+'_chosen',label,K.MOLECULAR);r=g.product(key+'_remainder',remainder_label,K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,remainder_label));return y
def qc(g,key,label,x):
    a=g.product(key+'_analytical',label,K.MOLECULAR);k=g.product(key+'_kept','Retained production stock',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(a,k));g.mark_terminal(Terminal(a.id,Disposition.QC,'Analytical aliquot not returned; amount unassigned.'));return k
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
TITLE='Repair-seq: associated-manuscript pooled-screen outcome-library preparation'
AUTHORITY='AssociatedNIHMS1750972/PMC9074467.1 Methods264–270/292 and source TableS2 identity. One original; final-body equivalence unverified; no supersededpreprint or arrayed/Perturb-seq route merged.'
def partial_adapter_end(inferred=True):
    sc=Scene();sc.strand('partial_adapter',[Segment('cohesive_end','NNNN',placeholder=True,inferred=inferred),Segment('unassembled_body','[adapter body/UMI position unassembled]',placeholder=True,inferred=inferred)])
    sc.ligation_end('partial_adapter','cohesive_end',kind='cohesive',polarity='5-prime',owner='adapter',text='Source remaining4ntNotI-compatible interface; literal adapter strands/TableS2 unassembled;5prime assignment derives from NotI cleavage class, not whole-end proof')
    return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('screen_cells','Source cells collected at end of pooled Repair-seq screen; prior integrated perturbation/outcome constructs and editing treatments outside selected library-prep boundary',K.BIOLOGICAL,primary=True)
    x=opaque(g,'extract','Source NucleoSpinBloodXLMaxi gDNA extraction,50–100millioncells percolumn',x,(add(g,'extraction_preparation','Source actual DNAextraction preparation; full formulation/physical subfractions unassembled'),))
    x=opaque(g,'noti','Source NotI-HF digestion ofgDNA',x,(add(g,'noti_enzyme','Source NotI-HF',K.CATALYST,R.CATALYST),condition(g,'noti_conditions','Source enzymebuffer/conditions unassembled')),S.DIGESTION)
    x=cleanup(g,'target_gel','Source large0.8percentagarosegel, multiple wells as needed; wide target-containing band excision/NucleoSpingelcleanup; exact boundaries/recovery unassigned',x)
    x=opaque(g,'hindiii','Source HindIII-HF removes one NotI end; retained opposite4ntNotI-compatible interface, exact targetfragment sequences unassembled',x,(add(g,'hindiii_enzyme','Source HindIII-HF',K.CATALYST,R.CATALYST),condition(g,'hindiii_conditions','Source digestion condition unassembled')),S.DIGESTION)
    x=cleanup(g,'pre_ligation_clean','Source0.8xSPRIselect',x)
    x=aliquot(g,'ligation_input','Source1ugHindIIIdigestedDNA selected for110ulligation',x)
    adapters=g.input('adapter_preparation','Sourceequimolar oBA676/oBA677 HPLC oligos annealed95C5min→70x1mingradual95to25C→1Chold; whole strands unassembled because actualTableS2gated',K.ADAPTER,chemical_state={'source_umi_length':12,'identifier_meaning':'sourceunique molecular identifier beforePCR','whole_orderings':'unassembled, no borrowedcanonicalbody','annealing':'sourceequimolar two oligos; whole duplex unassembled'})
    adapter=g.product('annealed_adapters','Source annealedcustom12ntUMIadapterstock, wholeduplex unassembled',K.ADAPTER,state=partial_adapter_end());g.react('adapter_anneal','Sourceactualtwo-oligoannealing preparation',S.MAGIC,((adapters,R.TEMPLATE),),(adapter,),reason='Whole source strand orderings/UMI window position unassembled; known source annealing procedure retained, not source absence.')
    product=g.product('ligated_mixture','SourceUMIadapter-bearing fragments plus unligated/adaptor/offproduct classes; complete junctions unassembled',K.MOLECULAR)
    geometry=AdapterLigationGeometry('known','cohesive',(partial_adapter_end(),),'Knownsource4ntremainingNotIcompatibleend; full source oBA676/oBA677 ordering/strand chemistry/pairedinterface unassembled, not unreported.')
    g.react('adapter_ligate','Source110ul:30ulligationbuffer/10ulligase/KAPAHyperPrepPCR-free,200:1adapter:insert,4Covernight/lid30C',S.LIGATION,((x,R.TEMPLATE),(adapter,R.STRUCTURAL),add(g,'ligation_enzyme','Source10ulKAPAHyperPrep ligase',K.CATALYST,R.CATALYST),condition(g,'ligation_buffer','Source30ulKAPAHyperPrepbuffer/water110ul; containsPEG')),(product,),geometry=geometry);x=product
    x=cleanup(g,'post_lig1','Sourcecalculated0.65xSPRIselect; effective ratio higher becausePEG, not normalized',x);x=cleanup(g,'post_lig2','Source0.8xSPRIselect second cleanup',x)
    one=g.product('one_pcr_template','Source30–50ngselectedaliquot forone50ulPCR',K.MOLECULAR);others=g.product('other_pcr_templates','Actualother aliquots ofsamecleanedstockforparallelPCRs, count/amounts unassembled',K.MOLECULAR);rem=g.product('post_lig_remainder','Untransferred postligationstock; nearlyentire—not complete—consumption sourcewording',K.MOLECULAR)
    g.react('pcr_allocation','SourceparallelPCRs use nearlyentire cleanedstock; selected representative plus other source-derived aliquots',S.SPLIT,((x,R.TEMPLATE),),(one,others,rem));g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,'Actualuntransferredstocknotsequencingreads.'))
    bulk=g.input('pcr_preparation','Source50ulPCR:KAPAHiFiHotStartReadyMix1x,eachprimer0.6uM,3percentDMSO; sourceP7+indexedScreenindexing primers wholeorderingsunassembled',K.MOLECULAR)
    cat=g.product('pcr_enzyme','SamepreparedPCR enzymecomponent',K.CATALYST);nt=g.product('pcr_nt','SamepreparedPCR nucleotidecomponent',K.MOLECULAR);f=g.product('pcr_f','SourceP7primercomponent; wholeorderingunassembled',K.OLIGO);r=g.product('pcr_r','SourceindexedScreenindexingprimercomponent; wholeorderingunassembled',K.OLIGO);buf=g.product('pcr_buffer','SamePCRbuffer/DMSO',K.CONDITION)
    g.react('pcr_components','Onepreparedstock, copresentcomponents—notseparatefreshadditions/physicalfractionation',S.GENERIC,((bulk,R.TEMPLATE),),(cat,nt,f,r,buf),changes_molecule=False)
    amplified=g.product('one_pcr_product','SourceDNAamplicon/remainingstockmixture; wholeamplicongeometryunassembled',K.MOLECULAR)
    g.react('pcr','Source16or20cycles98C15s/70C15s; selectedcyclecountunassigned',S.PCR,((one,R.TEMPLATE),(cat,R.CATALYST),(nt,R.STRUCTURAL),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(buf,R.CONDITION)),(amplified,),geometry=PCRGeometry.unavailable(unknown_scene(),reason='SourceactualTableS2wholeprimers/completeoutcomefragmentunassembled; notsourceabsence.'))
    other_amp=opaque(g,'parallel_pcrs','OtheractualparallelPCRs fromsameallocatedDNAstock; number/cycles/wholemixtureunassembled',others,(add(g,'parallel_pcr_preparations','SourceactualparallelPCRpreparations; totalstockallocationunassembled'),))
    kept=g.product('pcr_cleanup_kept','Source purified products across actual parallelPCRs; per-PCR versus combined workup allocation unassembled',K.MOLECULAR)
    removed=g.product('pcr_cleanup_removed','Actual source0.8xSPRI workup fractions; molecular recovery/carryover unassigned',K.WASTE)
    g.react('pcr_cleanup','Source0.8xSPRIselect cleanup of actual amplifiedDNA; timing of parallel-product combination unassembled',S.MAGIC,((amplified,R.TEMPLATE),(other_amp,R.TEMPLATE),condition(g,'pcr_cleanup_condition','Source0.8xSPRIselect workup')),(kept,removed),reason='Associatedmain270 specifies parallelPCRs and subsequent cleanup, but not whether those PCRs are combined before or aftercleanup. Actual produced inputs retained without inventing a precleanup physical pool.')
    g.discard(removed,'Actual cleanup excludes workup fractions; no claim of perfect specificity or yield.');x=kept
    x=qc(g,'qc','Sourcequality/yieldaliquots NanoDrop/Bioanalyser/orNovexgel; selectedassay/amountunassigned',x)
    chosen=aliquot(g,'index_pool_input','Selectedindexed sample library forinter-samplepool, amountunassigned',x)
    others=g.input('other_indexed_libraries','Actualotherindexed samplelibraryaliquots forsourceindexmix; upstream preparations/volumesunassembled',K.MOLECULAR)
    pool=g.product('index_pool','Actualindexed librarymixture',K.MOLECULAR);g.react('pool_indexes','Sourceindexed samples mixed, ratiosunassigned',S.POOL,((chosen,R.TEMPLATE),(others,R.TEMPLATE)),(pool,))
    lib=g.product('library','Sourcepooled-screen Repair-seq outcome library',K.LIBRARY);g.react('loading','SourceNovaSeq6000sequencing',S.GENERIC,((pool,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='AssociatedAMpooledscreen:I1=12ntUMI,I2=8ntsampleindex,R1=45ntCRISPRisgRNA,R2=258ntoutcome. Wholetarget/adapters/primers/readlandingsunassembled; NOTarrayedMiSeq/Perturbseq branches.'))
    g.add_source_item(SourceItem('unassembled','ActualTableS2oligoworkbook gated; fullfinalarticlebody/version equivalence and linked historicalgel guide ordinaryunreviewed.','PartialknownNotIinterface drawnwithout inventing fulladapter/sealedproduct. Noabsence/fullvalidationclaim.'));return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph();c('NotI precedes targetgel',list(g.steps).index('noti')<list(g.steps).index('target_gel'),True);c('HindIIIaftergel',g.steps['hindiii'].inputs[0].material,'target_gel_kept');c('Ligationonlyselected1ug',g.steps['adapter_ligate'].inputs[0].material,'ligation_input_chosen');c('UMIadapteraddedbeforePCR',list(g.steps).index('adapter_ligate')<list(g.steps).index('pcr'),True);c('SamePCRstockstructuralNT',g.producer_of('pcr_nt').id,'pcr_components');c('OnePCRusesactualaliquot',g.steps['pcr'].inputs[0].material,'one_pcr_template');c('Nearlyentirestockremainderretained','post_lig_remainder'in g.terminals,True);c('Indexpoolhasactualotherlibraries',len(g.steps['pool_indexes'].inputs),2);c('Privateincomplete',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
