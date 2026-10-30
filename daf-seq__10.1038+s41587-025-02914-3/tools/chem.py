"""Source-owned private partial DAF-seq, published SI cell-input benchmark route.

Associated AM is separately scoped; no claim it equals inaccessible final Methods.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry
TITLE='DAF-seq selected published SI targeted cell-input route'
AUTHORITY='PublishedSI MOESM1 SupplementaryNote pp9–15, actualprimerrecipeprintedp14/physicalp15 andcontinuationphysicalp16; finalNaturemainpreview/maps only; associatedAM PMC13063407 main150–168/244–250 separatelyread, no final-equivalence claim. Three originals closed.'
FORWARD='AATCACTATATTTCCATACTACTCA'
REVERSE='TACACCTAATRCAATAACARCCT'
TARGET='hg38 chr3:179229955-179234791'
def primer_scene(seq):
    sc=Scene();parts=[];start=0
    for i,b in enumerate(seq):
        if b=='R':
            if i>start:parts.append(Segment('fixed_'+str(start),seq[start:i]))
            parts.append(Segment('source_purine_'+str(i),'R',placeholder=True,note='Source constrained primer pool R=A/G, not N or a molecular identifier'));start=i+1
    if start<len(seq):parts.append(Segment('fixed_'+str(start),seq[start:]))
    sc.strand('primer',parts);return sc
def unknown_scene():
    sc=Scene();sc.strand('unassembled',[Segment('unassembled','[source target/product geometry unassembled]',placeholder=True)]);return sc
def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('nuclei','Selected source100000 permeabilized cells/nuclei; particularbenchmarkspecimen unassigned; upstreampermeabilization outside selected publishedSIrecipe',K.BIOLOGICAL,primary=True)
    def op(key,label,kind,inputs=(),product='Source reaction mixture including incompletelyprocessed substrate andfree reagents',geometry=None):
        nonlocal x
        kw={}
        if kind==S.GENERIC:kw['changes_molecule']=False
        if kind==S.MAGIC:kw['reason']='Historicalsource-specific sequencer-service finish and whole library geometry unassembled.'
        if geometry is not None:kw['geometry']=geometry
        y=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);x=y
    def input(key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(g.input(key,label,kind,**kw),role)
    def clean(key,label):
        nonlocal x
        c=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','Source physically recovered DNA; exact recovery/carryover unassigned',K.MOLECULAR);w=g.product(key+'_waste','Source separated extraction/workup fractions',K.WASTE);g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(y,w));g.discard(w,'Source physical extraction/workup, not an assumed sequence-specific exclusion.');x=y
    op('deaminate','Source4uMSsDddA25C10min; selected100k input, higherinputs separatebenchmarkbranches',S.CHEMICAL_CONVERSION,(input('ssddda','SourceSsDddA family(namedinpublishedSInote),4uM;WTvsSsDddA5variantunassigned inselectedcell-inputparagraph'),input('deam_condition','Source stateddeaminationconditions; reactionbufferrecipe notassembled from inaccessiblefinalMethods',K.CONDITION,R.CONDITION)),product='Source chromatin-stenciled DNA carryingC-to-U lesions at accessiblecytidines plusprotected/unconvertedDNA; no100percentconversion or directbaseTclaim')
    op('quench','Source20percentSDSstockaddedto5percentfinal',S.QUENCHING,(input('sds','SourceSDS',K.CONDITION,R.CONDITION),))
    clean('extract','SourceMonarchSpinGenomicDNAExtractionT3010 afterquench; historicalinternalfractionrecipeunassembled')
    production=g.product('qubit_kept','ProductionDNAretained',K.MOLECULAR);qc=g.product('qubit_sample','SourceQubit1xdsDNAHS analyticalaliquot;allocationunassigned',K.QC);g.react('qubit','SourceDNAquantificationseparatefromproduction',S.SPLIT,((x,R.TEMPLATE),),(production,qc));g.mark_terminal(Terminal(qc.id,Disposition.QC,'SourceQubitQC'));x=production
    chosen=g.product('pcr_input','Source15–70ngDddA-treatedgDNAselectedforPCR',K.MOLECULAR);rest=g.product('dna_remainder','UntransferredsourceDNAstock;amount/laterfateunassigned',K.MOLECULAR);g.react('pcr_aliquot','SourceDNAaliquotforselectedtargetedPCR',S.SPLIT,((x,R.TEMPLATE),),(chosen,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'SourceDNAstocknotentirelyconsumedbyPCR.'));x=chosen
    f=input('forward','Source1.5ul10uMforwardprimer',K.OLIGO,R.STRUCTURAL,state=primer_scene(FORWARD),chemical_state={'polymer':'DNA','source_ordering':FORWARD})
    r=input('reverse','Source1.5ul10uMreverseprimer;R=A/G',K.OLIGO,R.STRUCTURAL,state=primer_scene(REVERSE),chemical_state={'polymer':'DNA','source_ordering':REVERSE,'degenerate_base_classes':{'R':'A/G'},'purpose':'priming convertedDNA, not UMI'})
    op('pcr','PublishedSI50ul:25ulrepliQaHiFiToughMix95200,1.5uleach10uMprimer,15–70ngDNA;98C30s,30x98C10s53C5s68C25s,68C1min',S.PCR,(f,r,input('repliqa','SourcefreshrepliQaHiFiToughMix')),product='Source targeted '+TARGET+' ampliconpopulations plus off-target/truncated/chimericproducts and residualoriginalDNA/freeprimers; UreadasT aftercopying, not perfectchemicalconversion',geometry=PCRGeometry.unavailable(unknown_scene(),reason='Wholeconvertedtargetpopulation and bothstrandprimerlanding unassembled; no exactreferenceamplicon or absenceclaim.'))
    clean('pcr_clean','SourceMonarchSpinPCRandDNACleanupT1130')
    op('service','Source standardPremiumPCRservicePlasmidsaurus; historicallibraryadapters/indexes/runconditionsunassembled',S.MAGIC)
    lib=g.product('library','Source long-readservice sequencing preparation',K.LIBRARY);g.react('loading','Source service loading',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='Selected publishedSIPlasmidsaurusPremiumPCRroute; no currentservicekit substitution. Other demonstrated routes include targetedONT/PacBio andsingle-cellcustomPTA; notimposedhere. C->U marking vsPCR-reportedT and computationalY/Rrecoding/deduplication are distinct; no barcode/UMI invented from deaminationpatterns.'))
    g.add_source_item(SourceItem('pending','FinalmainMethods subscriptionpreviewonly; associatedAM notclaimedfullyequivalent. PublishedSupplementaryTablesworkbook remainsunreadwithincap, historicalservicefinish/wholeproductsunassembled.','Source-positive publishedSIpreparation scope only; no fullsource/absenceaudit.'));return g
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('SDSquenchesactualdeamination',g.steps['quench'].inputs[0].material,'deaminate_product')
    c('DNAextractionafterSDS',g.steps['extract'].inputs[0].material,'quench_product')
    c('SourcePCRusesaliquot',g.steps['pcr'].inputs[0].material,'pcr_input')
    c('UntransferredDNAexplicit',g.terminals['dna_remainder'].disposition,Disposition.UNRESOLVED)
    c('ReverseRnotreplacedwithN',g.materials['reverse'].chemical_state['source_ordering'],REVERSE)
    c('EveryRinstatepreserved',''.join(s.top for s in primer_scene(REVERSE).strands['primer'].segs),REVERSE)
    c('Qubitphysicallyseparate',g.terminals['qubit_sample'].disposition,Disposition.QC)
    c('NounsupportedPTAinselectedroute',any('pta'in s.id for s in g.steps.values()),False)
    c('HistoricalserviceafterphysicalPCRcleanup',g.steps['service'].inputs[0].material,'pcr_clean_kept')
    p=sources_root()/Path(__file__).resolve().parents[1].name/'s41587-025-02914-3_41587_2025_2914_MOESM1_ESM.pdf.txt'
    if have(c,Source(p,'https://doi.org/10.1038/s41587-025-02914-3')):
        compact=''.join(p.read_text().split())
        for name,seq in [('forward',FORWARD),('reverse',REVERSE)]:c('WholepublishedSIprimer '+name,seq in compact,True)
    c('Privatepartialsourcereview',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
