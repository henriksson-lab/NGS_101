"""Private CAT-ATAC: source-owned modified primer/splint and ordered capture routes."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
class B:
    def __init__(self,title,authority,start,kind=K.BIOLOGICAL):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',start,kind,primary=True)
    def stock(self,key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind,inputs=(),geometry=None,product='Source mixture including incomplete/unreacted/off-target classes; occupancy/recovery unassigned',reason=None):
        y=self.g.product(key+'_product',product,K.MOLECULAR);kw={}
        if kind==S.MAGIC:kw['reason']=reason or 'Positive source operation; complete occupied molecular inputs/products and historical formulation unassembled, not a missing-information conclusion.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label):
        y=self.g.product(key+'_kept','Actual source retained physical fraction; yield/carryover unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual separated source workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label)),(y,w));self.g.discard(w,'Actual workup fraction excluded; no perfect specificity or molecular yield claim.');self.x=y
    def split(self,key,label,selected,remainder):
        y=self.g.product(key+'_selected',selected,K.MOLECULAR);r=self.g.product(key+'_remainder',remainder,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r));self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,remainder));self.x=y;return y

def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('source products unassembled','[unassembled]',placeholder=True)]);return s

from nextera import ADAPTOR_S5
from illumina import P5,P7,TRUSEQ_READ2
AUTHORITY='PublishedPMC12859477main122/379–425/486 andDocumentS1MethodS1originalp16–20; actualp16/p19 tables read. Two originals only.'
CS1='TTGCTAGGACCGGCCTTAAAGC'
RT_ORDERING='/5Phos/'+ADAPTOR_S5+'N'*12+CS1
SPLINT_BASES='TGACGCTGCCGACGACAGACGCG'
SPLINT_ORDERING='TGACGCTGCC+G+A+C+G+ACAGACGCG/3Phos/'
R2_TSO=TRUSEQ_READ2+'AAGCAGTGGTATCAACGCAGAG'
PARTIAL_P5=P5[:21]
INDEXES={'6a':'TTCGCAGT','6b':'CGAGACTA','6c':'ACAGCTCA','6d':'AAGTGTCG','6e':'AGTAAACC','6f':'CCGTTTAG'}
SAMPLE_SHEET_I7={k:revcomp(v)for k,v in INDEXES.items()}
def index_primer(which):return P7+INDEXES[which]+TRUSEQ_READ2

def rt_parts():
    return [Segment('RT first5',ADAPTOR_S5[:5]),Segment('RT next10',ADAPTOR_S5[5:15]),Segment('rest NexteraR1',ADAPTOR_S5[15:]),Segment('source UMI','N'*12,placeholder=True,feature=feature('cat_atac.umi','umi','random')),Segment('source CS1 complement',CS1)]
def splint_parts():
    return [Segment('DNA pair10',SPLINT_BASES[:10]),Segment('LNA pair5',SPLINT_BASES[10:15],note='Five sourceLNA residues; unmodifiedDNA thermodynamics not asserted'),Segment('bead-complementary arm',SPLINT_BASES[15:],note='Source sequence retained; whole occupied bead arm not assembled')]
def rt_scene():
    s=Scene();s.strand('rt',rt_parts(),mod5='P');return s

def splint_scene():
    s=Scene();s.strand('splint',splint_parts(),mod3='P');return s

def duplex():
    # Exact source local15bp interface; full bead/gRNA binding and copiedproduct omitted.
    s=Scene();s.strand('rt',rt_parts(),mod5='P',label='Source original RT-primer component; copiedgRNA extension unassembled')
    s.anneal('splint',splint_parts(),to='rt',pair=('DNA pair10','RT next10'),unpaired=('bead-complementary arm',),mod3='P')
    s.ligation_end('rt','RT first5',kind='splinted',polarity='5-prime',owner='adapter',text='SourceRT5P localcapturecomponent; completebead3end/occupiedgRNAcDNAjunctionunassembled')
    return s

def primer(seq,index_family=False):
    s=Scene()
    if index_family:segs=[Segment('P7',P7),Segment('selected source i7','N'*8,placeholder=True,feature=feature('cat_atac.i7','sample_index','unknown',note='One of six printed index choices; actual per-sample choice unassigned')),Segment('R2 tail',TRUSEQ_READ2)]
    else:segs=[Segment('source whole primer',seq)]
    s.strand('primer',segs);return s

def protocol_graph():
    b=B('CAT-ATAC published three-modality capture preparation',AUTHORITY,'Source preparedperturbedcells; constructdelivery/editingtreatments outside graph');g=b.g
    def qc(key,label,amount=None):
        a=g.product(key+'_analytical',label,K.MOLECULAR);r=g.product(key+'_kept','Actual retained productionstock afterQC; recoveryunassigned',K.MOLECULAR)
        g.react(key,label,S.SPLIT,((b.x,R.TEMPLATE),),(a,r));g.mark_terminal(Terminal(a.id,Disposition.QC,'Actualanalyticalaliquot notreturnedtoproduction.'));b.x=r
    def pcr(key,label,reverse,index_family=False):
        bulk=g.input(key+'_stock','Actual source100ulPCRstock:50ulKapaHiFiHotStartReadyMix+4ul10uMforward+4ul10uMreverse+37ulwater; template5ulseparate',K.MOLECULAR)
        f=g.product(key+'_f','Sourceoligo8partialP5 component ofsamepreparedstock',K.OLIGO,state=primer(PARTIAL_P5))
        r=g.product(key+'_r','Sourceoligo'+('6a-fselectedindexfamily'if index_family else'5R2TSO')+'component ofsamepreparedstock',K.OLIGO,state=primer(reverse,index_family),chemical_state={'source_choices':{k:index_primer(k)for k in INDEXES}}if index_family else{'source_ordering':reverse})
        e=g.product(key+'_enzyme','SamepreparedKAPAHiFi enzymecomponent',K.CATALYST);nt=g.product(key+'_nt','SamepreparedReadyMix nucleotidecomponent',K.MOLECULAR);buf=g.product(key+'_buffer','SamepreparedPCRbuffer/waterconditions',K.CONDITION)
        g.react(key+'_components','OnepreparedPCRstockcomponents, not freshduplicateadditions',S.GENERIC,((bulk,R.TEMPLATE),),(f,r,e,nt,buf),changes_molecule=False)
        b.op(key,label,S.PCR,((f,R.STRUCTURAL),(r,R.STRUCTURAL),(e,R.CATALYST),(nt,R.STRUCTURAL),(buf,R.CONDITION)),geometry=PCRGeometry.unavailable(unknown_scene(),reason='Wholebead/gRNA/TSO/PCRproduct and actualsampleindexchoice/primerlandingunassembled. Whole sourceprimers/choice inventory retained.'))
    b.op('extract_nuclei','Sourceextractnuclei accordingto linked10xcelllineprotocol',S.MAGIC,(b.stock('nuclear_preparation','Actualhistorical10xnucleuspreparationstocks, internalsunassembled',K.MOLECULAR,R.STRUCTURAL),))
    b.op('tagment','Sourcetagmentnuclei accordingto historical10xMultiomeprotocol',S.MAGIC,(b.stock('tagment_preparation','ActualhistoricalTn5/loadedoligo/molecularpreparationstock, wholecompositionunassembled',K.MOLECULAR,R.STRUCTURAL),))
    nuclei=b.x
    rt=g.input('rt_primer','Sourceoligo1'+RT_ORDERING,K.OLIGO,state=rt_scene(),chemical_state={'source_ordering':RT_ORDERING,'five_prime':'P','UMI_length':12})
    sp=g.input('splint','Sourceoligo2'+SPLINT_ORDERING,K.OLIGO,state=splint_scene(),chemical_state={'source_ordering':SPLINT_ORDERING,'LNA_one_based_positions':(11,12,13,14,15),'three_prime':'P'})
    paired=g.product('rt_splint','SourcepreannealedRT/splintduplex, exactlocal15bpinterface only',K.MOLECULAR,state=duplex(),chemical_state={'source_orderings':(RT_ORDERING,SPLINT_ORDERING),'fullbindingcontexts':'unassembled'})
    g.react('anneal','Sourcepreannealoligos1/2 perlinkedIDTmethod; exactthermalprogrammeunreviewed',S.HYBRIDIZATION,((rt,R.TEMPLATE),(sp,R.STRUCTURAL),b.cond('anneal_conditions','SourceIDTannealinginstructions ordinaryunread')),(paired,))
    b.x=nuclei
    b.op('gem','SourceGEMformationwithpreannealedduplex/HiT4 andsource64ulbarcodingmix',S.MAGIC,
       ((paired,R.STRUCTURAL),b.stock('barcoding_mix','Source49.5ulbarcodingreagentmix, molecularNTs/componentsunassembled',K.MOLECULAR,R.STRUCTURAL),b.stock('gem_tso','Source1.1ultemplateswitcholigo, wholeorderingunassembled; not automaticallyPCRoligo5',K.OLIGO,R.STRUCTURAL),b.stock('barcoding_enzymes','Source7.5ulbarcodingenzymemix',K.CATALYST,R.CATALYST),b.stock('hi_t4','Source3ulHiT4ligaseNEBM2622S, sameadditionthroughthermalprogramme',K.CATALYST,R.CATALYST),b.stock('bead_arms','ActualMultiomegelbeads/capturearms; wholeoligos/physicalallocationunassembled',K.MOLECULAR,R.STRUCTURAL,chemical_state={'cell_barcode_feature':feature('cat_atac.cell','cell_barcode','unknown'),'source_diagram_BC_length':16}),b.cond('gem_conditions','Source1.9ulreducingagent;1ulRTsplintdescribed0.5uM; total64ul; stock/finalconventionunassigned')),
       reason='Positive source coencapsulation; exact capturebead/occupiedgRNA/ATACfragments/internal molecularstocks unassembled.')
    b.op('atac_preligation','Source25C30min ATACcaptureligation precedes highertemperatureRT',S.MAGIC,(b.cond('first25','Source25C30min'),),reason='Source operationorder is explicit. Full ATACend/beadligation substrates/products and competingRTfill-inclasses unassembled, not assumed frompaperhypothesis.')
    b.op('rt53','Source53C45min gRNART/templateswitching andRNAcapture',S.MAGIC,(b.cond('rt_conditions','Source53C45min; catalysts/primers/NTs retainedfromsameGEMmix, nofreshaddition'),),reason='KnownsourceRT/TSO positive; whole gRNA/RNA/bead substrates/products andquantitativeconversion unassembled.')
    b.op('post25','Source25C30min afterRT, then4Chold; HiT4availablefromsameGEMstock',S.MAGIC,(b.cond('post25_conditions','SourcepostRT25C30min then4C'),),reason='SourcepostRTcaptureligation retained; knownRT5P/LNAsplintlocalinterface represented inGEMmaterial, fullbead/occupiedcopiedgRNAjunctions and sealed products unassembled.')
    b.op('gem_recover','SourcefollowhistoricalMultiomepostGEMworkflowuntilcleanup',S.MAGIC,(b.stock('gem_recovery_preparation','ActualhistoricalGEMrecovery/workupmolecularstocks, allocationsunassembled',K.MOLECULAR,R.STRUCTURAL),))
    b.clean('modified_cleanup','Source3.2xSPRI160ul+1.8xisopropanol90ul replacing1.8xbeads; intendedssgRNAcDNAretention, recoveryunassigned')
    b.op('preamp','SourcehistoricalMultiomepreamp stage; exactstocks/cyclesandwholeprimersunassembled',S.MAGIC,(b.stock('preamp_preparation','Actualsourcehistoricalpreampstock withmolecularnucleotides/primers/catalysts, wholeformulationunassembled',K.MOLECULAR,R.STRUCTURAL),))
    atac=g.product('atac_preparation','Actual sourceATAClibrarypreparationinputfromsameGEMpopulation, exactallocationunassembled',K.MOLECULAR)
    gex=g.product('gex_preparation','Actual sourceGEXlibrarypreparationinputfromsameGEMpopulation, exactallocationunassembled',K.MOLECULAR)
    cap=g.product('crispr_preamp','Actual sourceCRISPRpreampstock, exactrelation/amountallocatedfromstandardworkflowunassembled',K.MOLECULAR)
    rem=g.product('preamp_remainder','Untransferredstandardworkflowstock/fractions, fateunassembled',K.MOLECULAR)
    g.react('library_allocations','SourceATAC/GEXstandardroutes andCRISPRpreampaliquots; wholephysicalallocationsequenceunassembled',S.MAGIC,((b.x,R.TEMPLATE),),(atac,gex,cap,rem),reason='Sourcepositive3libraryroutes; exact historicalfractionation/aliquot/retention scheduleunassembled. RemovedlargeCRISPRfractionnotinventedasGEXinput.');g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,'Otheractualworkflowstockfatesunassembled.'))
    b.x=gex;b.op('gex_library_prep','SourceGEXlibrary construction perhistorical10xprotocol',S.MAGIC,(b.stock('gex_stock','ActualhistoricalGEX molecularconstructionstocks, end/oligo/programmeunassembled',K.MOLECULAR,R.STRUCTURAL),))
    b.split('gex_loading','SourceGEXsequencingallocation, amountunassigned','ActualGEXloadingaliquot','UntransferredGEXstock, fateunassigned')
    glib=g.product('gex_library','SourceCATATACGEXlibrary readout',K.LIBRARY);g.react('gex_load','SourceNovaSeqXGEX28x10x10x90',S.GENERIC,((b.x,R.TEMPLATE),),(glib,),changes_molecule=False);g.sequence(glib,SequencingPlan(limitation='SourceNovaSeqX28/10/10/90, >=20000readpairscell. WholeRNAcapturebeadarms/library/primerlandingunassembled; separatefromATACCRISPRpool.'))
    b.x=atac;atac_finished=b.op('atac_library_prep','SourceATAClibrary construction perhistorical10xprotocol',S.MAGIC,(b.stock('atac_stock','ActualhistoricalATAC molecularconstructionstocks, wholeends/primers/indexchoicesunassembled',K.MOLECULAR,R.STRUCTURAL),))
    b.x=cap;b.split('capture_input','Source40ulpreampselectedforCRISPRbranch','Actual40ulpreampaliquot','UntransferredCRISPRpreampstock, amount/fateunassigned')
    b.op('dilute_capture','Sourceadd60ulEBto40ulpreamp, total100ul',S.GENERIC,(b.cond('eb60','Source60ulEB'),))
    sn=g.product('capture_first_sn','ActualfirstSPRI supernatant retainedforsecondaddition',K.MOLECULAR);large=g.product('capture_large','Actualfirst60ulSPRIlargecDNAboundfractionremovedfromCRISPRroute',K.MOLECULAR)
    g.react('capture_first_sp','Sourcefirst60ulSPRI on100ulstock(.6x), exclude largecDNAfraction',S.SPLIT,((b.x,R.TEMPLATE),b.cond('first60','Actualfirst60ulSPRI')),(sn,large));g.mark_terminal(Terminal(large.id,Disposition.UNRESOLVED,'ExcludedlargefractionnotautomaticallyGEXproduction; exactfateunassigned.'));b.x=sn
    b.clean('capture_second_sp','Sourceactualsecond40ulSPRIaddition(total1.0xoriginalbasis), retainedselectedmaterial/elute30ul')
    b.split('feature_input','Source5ulselectedfor100ulFeaturePCR','Actual5ulFeaturePCRtemplate','OtherelutedCRISPRcapturestock, fateunassigned')
    pcr('feature_pcr','Source10cycles98C20s/58C5s/72C5s, initial98C45s/final72C1min',R2_TSO)
    qc('feature_qc','Actual2ulTapeStationD1000HSaliquot',2)
    b.clean('feature_cleanup','Source.8x–.9xSPRIcleanup/elute30ul; exactstepallocationunassembled')
    qc('feature_post_qc','ActualpostcleanupTapeStationconfirmationaliquot, amountunassigned')
    b.split('index_input','Source5ulselectedfor100ulIndexPCR','Actual5ulIndexPCRtemplate','UntransferredpostFeaturePCRstock, fateunassigned')
    pcr('index_pcr','Source9cycles98C20s/60C30s/72C20s, initial98C45s/final72C1min',None,index_family=True)
    qc('index_qc','Actual1ulprecleanupTapeStationaliquot',1)
    b.clean('index_cleanup','Source.7x–.8xdouble-sidedSPRI/elute30ul')
    qc('index_post_qc','Actual1ulpostcleanupTapeStationaliquot',1)
    crispr_chosen=b.split('crispr_pool_input','ActualCRISPRlibraryaliquotforATACCRISPRpool; amountunassigned','SelectedCRISPRpoolaliquot','OtherCRISPRlibrarystock, fateunassigned')
    b.x=atac_finished;atac_chosen=b.split('atac_pool_input','ActualATAClibraryaliquotforATACCRISPRpool; amountunassigned','SelectedATACpoolaliquot','OtherATAClibrarystock, fateunassigned')
    pool=g.product('atac_crispr_pool','ActualsourceATAC:CRISPR5:1pool; ratio basisunassigned',K.MOLECULAR);g.react('pool_atac_crispr','Sourcephysical5:1ATAC/CRISPRlibrarypool',S.POOL,((atac_chosen,R.TEMPLATE),(crispr_chosen,R.TEMPLATE)),(pool,))
    lib=g.product('atac_crispr_library','SourceATACCRISPRpool readout',K.LIBRARY);g.react('atac_crispr_load','SourceNovaSeqX151x8x24x151',S.GENERIC,((pool,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='SourceactualNovaSeqX151/8/24/151, ATAC25000/CRISPR5000readpairscell; sourceSI recommends100/8/24/100and>=75bpPE. Sourcei5diagramBC16 isnot24cyclephysicalreadorstandalonesampleindex. Actualindexchoice/wholemolecules/primerlandingunassembled.'))
    g.add_source_item(SourceItem('local_partial','Exactsource5P/N12UMI/LNA5/3P/partial15bpinterface preserved; wholecapturebead/gRNA/template-switch/PCRproducts andLNATm unassembled.','PublishedmainandSIonly; IDT/historical10xlinkedoriginals ordinaryunread; no fullsource/reaction/absenceaudit.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Exactlocal15bpinterface',revcomp(SPLINT_BASES[:15]),ADAPTOR_S5[:15])
    c('WholemodifiedsplintnotordinaryDNA',len(g.materials['splint'].chemical_state['LNA_one_based_positions']),5)
    c('ATACcapture precedeshigherRT',list(g.steps).index('atac_preligation')<list(g.steps).index('rt53'),True)
    c('PostRTcaptureusesactualRToutput',g.steps['post25'].inputs[0].material,'rt53_product')
    c('HiT4notfreshlyaddedperthermalstage',sum(any(x.material=='hi_t4'for x in s.inputs)for s in g.steps.values()),1)
    c('Secondsizefractionusesactualfirstsupernatant',g.steps['capture_second_sp'].inputs[0].material,'capture_first_sn')
    c('FeaturePCRusesactual5ulaliquot',g.steps['feature_pcr'].inputs[0].material,'feature_input_selected')
    c('IndexPCRusesactual5ulpostQCalc',g.steps['index_pcr'].inputs[0].material,'index_input_selected')
    c('Feature2ulQCprecedescleanup',g.steps['feature_cleanup'].inputs[0].material,'feature_qc_kept')
    c('SourceGEXreadoutnotmergedwithATACCRISPR',set(g.sequencing),{'gex_library','atac_crispr_library'})
    c('PhysicalATACCRISPRpooltwoactualfractions',len(g.steps['pool_atac_crispr'].inputs),2)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'j.crmeth.2025.101222_supp_mmc1.pdf.txt'
    if have(c,Source(p,'https://doi.org/10.1016/j.crmeth.2025.101222')):
        t=''.join(p.read_text().split()).upper()
        for name,seq in [('RT',RT_ORDERING),('LNAsplint',SPLINT_ORDERING),('R2TSO',R2_TSO),('partialP5',PARTIAL_P5)]+[(k,index_primer(k))for k in INDEXES]:c('Sourcewholeordering '+name,seq.upper()in t,True)
        for k,seq in SAMPLE_SHEET_I7.items():c('SourceprintedRCsamplesheet '+k,seq in t,True)
    c('Privatepartialremainsunready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
