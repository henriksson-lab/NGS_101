"""Protocol-owned selected chemical classes; incomplete source-qualified graph."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, revcomp, feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry

def scene(name,text):
    sc=Scene();sc.strand(name,[Segment(name,text,placeholder=True)]);return sc

def clean(g,key,x,label):
    condition=g.input(key+'_condition',label,K.CONDITION)
    kept=g.product(key+'_kept','Source recovered material; no perfect purity or quantitative partition assumed',K.MOLECULAR)
    waste=g.product(key+'_waste','Actually separated cleanup fractions with unquantified molecular carryover',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(condition,R.CONDITION)),(kept,waste));g.discard(waste,'Source physical separation, not perfect molecular classification.');return kept

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source mixture',catalyst=False,product=None,state=None,chemical_state=None):
    stock=g.input(key+'_reagent',reagent,K.CATALYST if catalyst else K.CONDITION)
    y=g.product(key+'_product',product or label+'; all physical reaction companions retained until separation',K.MOLECULAR,state=state,chemical_state=chemical_state)
    kwargs={}
    if kind==S.MAGIC:kwargs['reason']='Source operation established; exact internal chemistry/geometry not assembled from unreviewed historical originals.'
    if kind==S.PCR:kwargs['geometry']=PCRGeometry.unavailable(Scene(),reason='Exact source template/primer endpoint geometry unassembled, no ordinary complete library asserted.')
    if kind in (S.LIGATION,S.PROXIMITY_LIGATION):kwargs['geometry']=AdapterLigationGeometry('unavailable','unavailable',(scene('unassembled','[source-known ligation; full substrate/end geometry unassembled]'),),'Unreviewed geometry; not a claim of unavailable disclosure.')
    inputs=[(x,R.TEMPLATE),(stock,R.CATALYST if catalyst else R.CONDITION)]
    if kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.PROXIMITY_LIGATION} and not catalyst:
        enzyme=g.input(key+'_enzyme','Explicit catalytic component of source mixture: '+reagent,K.CATALYST);inputs.append((enzyme,R.CATALYST))
    if kind==S.LIGATION:
        adapter=g.input(key+'_adapter','Source ligating molecular adapter; exact historical oligo/end inventory unassembled',K.ADAPTER,state=scene(key+'_adapter','[source adapter geometry unassembled]'));inputs.append((adapter,R.STRUCTURAL))
    if kind==S.PCR:
        for arm in ('f','r'):
            primer=g.input(key+'_primer_'+arm,'Source PCR primer family; exact oligo/endpoint unassembled',K.OLIGO,state=scene(key+'_primer_'+arm,'[source primer body/index geometry unassembled]'));inputs.append((primer,R.STRUCTURAL))
    g.react(key,label,kind,inputs,(y,),**kwargs);return y

def qc_split(g,key,x,label):
    assay=g.product(key+'_assay','Physically separate source analytical aliquot',K.QC)
    kept=g.product(key+'_retained','Unassayed retained preparation for downstream library work',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(assay,kept));g.mark_terminal(Terminal(assay.id,Disposition.QC,'Source analytical aliquot never feeds retained library.'));return kept
from dataclasses import replace
from illumina import TRUSEQ_SMALL_RNA_RA3,TRUSEQ_SMALL_RNA_RA5,P5,P7
SOURCE='PublishedVORPMC8782182 Fig1/overview87–119, Methods185–247; publishedSupplementaryDOCX TableS5 originalXMLcells independentlyread inclmodifications. DOI10.1080/15476286.2021.1999105.'
@dataclass(frozen=True)
class FirstAdapter:
    umi:str='NNNNNNNN'
    activation:str='P'
    def __post_init__(self):
        if len(self.umi)!=8 or set(self.umi)-set('ACGTN') or self.activation not in {'P','rApp','ligated'}:raise ValueError('Source8ntUMIandexplicitadapteractivationrequired')
    @property
    def sequence(self):return 'TCAGATCG'+self.umi+TRUSEQ_SMALL_RNA_RA3
    def adenylate(self):
        if self.activation!='P':raise ValueError('Sourceadenylationrequires5phosphateadapter')
        return replace(self,activation='rApp')
    def scene(self):
        sc=Scene();sc.strand('first_adapter',[Segment('prefix','TCAGATCG'),Segment('umi',self.umi,placeholder='N' in self.umi,feature=feature('cyphyrna.umi','umi','unknown')),Segment('ra3',TRUSEQ_SMALL_RNA_RA3)],mod5=self.activation,mod3='Amino-C7');return sc

RT_PRIMER='G'+revcomp(TRUSEQ_SMALL_RNA_RA3)
SECOND_ADAPTER=revcomp(TRUSEQ_SMALL_RNA_RA5)
HAIRPIN='CGACACTGTCGGTACCGACGGGAGAAGTCGGTACCGACAGTGTCGTGGN'
def rt_binding():
    sc=FirstAdapter(activation='rApp').scene();sc.anneal('rt_primer',[Segment('source_extra_g','G'),Segment('rt_core',revcomp(TRUSEQ_SMALL_RNA_RA3))],to='first_adapter',pair=('rt_core','ra3'));return sc
@dataclass(frozen=True)
class LigatedRNA:
    rna:str
    source_end:str
    adapter:FirstAdapter
    junction_two_p:bool
    def __post_init__(self):
        if not self.rna or set(self.rna)-set('ACGU') or self.source_end not in {'cP','OH'} or self.adapter.activation!='ligated':raise ValueError('ExplicitRNA/sourceendandcovalentlyligatedDNAadapterrequired')
    def tpt1_selected(self):
        if not self.junction_two_p:raise ValueError('SelectedTpt1substratemusthavejunction2phosphate')
        return replace(self,junction_two_p=False)
    def selected_cdna_extension(self):
        # Extension upstream of annealedRA3core; doesnotcopyphysicalAmino-C7block.
        return Construct([Segment('umi_copy',revcomp(self.adapter.umi),placeholder='N' in self.adapter.umi,feature=feature('cyphyrna.umi','umi','unknown')),Segment('prefix_copy',revcomp('TCAGATCG')),Segment('rna_copy',revcomp(self.rna.replace('U','T')))],name='Conditional copiedextension, actualRTprimerbodyseparate')

def ligate_first(rna,end,adapter,enzyme):
    if adapter.activation!='rApp':raise ValueError('Sourcefirstadaptermustbepreadenylated')
    if (end,enzyme) not in {('cP','AtRNL_AA'),('OH','T4RNL2_TKQ')}:raise ValueError('Selectedsourcemolecularend/ligasecontractnotmet')
    return LigatedRNA(rna,end,replace(adapter,activation='ligated'),end=='cP')
@dataclass(frozen=True)
class SpacerPrimer:
    nodes:tuple=(('DNA',SECOND_ADAPTER),('HEG','iSp18'),('DNA','CACTCA'),('HEG','iSp18'),('DNA',RT_PRIMER))
    circular:bool=False
    extended:bool=False
    def __post_init__(self):
        if type(self.circular)is not bool or type(self.extended)is not bool or self.circular and not self.extended or sum(k=='HEG' for k,v in self.nodes)!=2 or any(k not in {'DNA','HEG'} or (k=='DNA' and (not v or set(v)-set('ACGT'))) or (k=='HEG' and v!='iSp18') for k,v in self.nodes):raise ValueError('OrderedsourceDNA/twoactualHEGnodesandconsistentextension/closurestate required')
    def extend(self):
        if self.circular or self.extended:raise ValueError('LinearunextendedsourceprimermustprecedeRT')
        return replace(self,extended=True)
    def close(self):
        if self.circular or not self.extended:raise ValueError('Sourcecircleclosurefollowsextended5Pprimer')
        return replace(self,circular=True)
    def flatten(self):raise ValueError('ActualHEGspacerscannotbeerasedtoordinaryDNA')
    def scene(self):
        sc=Scene();segs=[Segment('node'+str(i),value if kind=='DNA' else '[sourcehexaethyleneglycolspacer]',placeholder=kind!='DNA') for i,(kind,value) in enumerate(self.nodes)]
        if self.extended:segs.append(Segment('cdna','[sourceRNA/adaptor-derivedDNAextension; exactfragmentunassigned]',placeholder=True))
        sc.strand('spacer_primer',segs,mod5='' if self.circular else 'phosphate',label='Sourcecircleopenedfordisplay; actualspacersretained' if self.circular else 'Sourcespacer-containing5phosphateRTprimer');return sc

def oligo(name,seq,mod5='',mod3=''):
    sc=Scene();sc.strand(name,[Segment(name,seq)],mod5=mod5,mod3=mod3);return sc

def primer(which):
    if which=='f':return oligo('forward',P5+TRUSEQ_SMALL_RNA_RA5[:-5])
    if which!='r':raise ValueError('SourcePCRf/rrequired')
    sc=Scene();sc.strand('reverse',[Segment('p7',P7),Segment('index','N'*6,placeholder=True,feature=feature('cyphyrna.i7','sample_index','unknown')),Segment('source_arm','GTGACTGGAGTT'+revcomp(TRUSEQ_SMALL_RNA_RA3))]);return sc

def graph_for_scope(end='cP',route='linear'):
    if end not in {'cP','OH'} or route not in {'linear','circular_spikein'}:raise ValueError('Actualsourceend/librarybranchrequired')
    g=ProtocolGraph('cyPhyRNA-seq — '+end+'/'+route,authority=SOURCE,representation='published_instance')
    x=g.input('total_rna','SourceTRIzolisolatedtotalRNA; D.vulgarislinear branch or1ugE.colitotalRNAcircularspike-inbranch',K.MOLECULAR,primary=True)
    if route=='circular_spikein':
        spike=g.input('ribozyme_spike','Source1ngself-cleavedPAGEpurified '+('HHR5primefragmentwithterminal2prime3primecyclicphosphate' if end=='cP' else 'twister3primefragmentwith5primeOH')+'; actualfragmentsequenceunassigned',K.MOLECULAR)
        mixture=g.product('spiked_total_rna','SourceE.colitotalRNAplusindependentpurifiedribozymeRNAfragmentbeforeend-enrichment',K.MOLECULAR)
        g.react('spike_pool','Source1ngribozymefragmentaddedto1ugE.colitotalRNA',S.POOL,((x,R.TEMPLATE),(spike,R.TEMPLATE)),(mixture,));x=mixture
    if end=='cP':
        cond=g.input('peg_salt','Source5%PEG8000/0.5MNaCl,minus20C30min,10000g4C30min',K.SELECTION);small=g.product('short_rna','SourcecollectedsupernatantenrichedforshortRNA, notperfectrRNAremoval',K.MOLECULAR);large=g.product('large_rna','PhysicallyprecipitatedlargeRNAcompanion, exactfateunassigned',K.MOLECULAR)
        g.react('peg_select','SourcecPbranchlargeRNAprecipitation; retainsupernatant',S.SELECTION,((x,R.TEMPLATE),(cond,R.SELECTION)),(small,large));g.mark_terminal(Terminal(large.id,Disposition.UNRESOLVED,'ActuallargeRNAfractionnotsequencedbyselectedbranch; nolaterfateassumed.'))
        x=clean(g,'short_recover',small,'Source3volEtOHminus20Covernight/17000grecovery/70%EtOHwash')
    else:
        x=op(g,'rpph_xrn1',x,'SourceRppH/Xrn137C15min enrich5OHclasses',S.DIGESTION,reagent='SourceRppHconverts5PPPto5P; Xrn1digests5PRNA',product='Enriched5OH-RNAclasses plusdigestedRNAproducts/remainingclasses untilactualcleanup; noorigin-specificityasserted')
        x=clean(g,'end_cleanup',x,'SourceGeneJETcleanupafterRppH/Xrn1')
    h=g.input('blocking_hairpin','SourceCEW369tRNAblockingDNAhairpinwithterminaldegenerateNand5P/3P; notUMI',K.OLIGO,state=oligo('hairpin',HAIRPIN,mod5='phosphate',mod3='phosphate'))
    cat=g.input('blocking_ligase','SourceT4DNAligase/RNaseinhibitor32C5h',K.CATALYST);b=g.product('blocked_mix','Sourcehairpin-ligatedtRNAclasses plusunblockedtRNA/otherRNA/freehairpin/reagents retained',K.MOLECULAR)
    g.react('trna_block','SourcehairpinligationtoblocktRNA3ends,65C10mininactivation',S.LIGATION,((x,R.TEMPLATE),(h,R.STRUCTURAL),(cat,R.CATALYST)),(b,),geometry=AdapterLigationGeometry('unavailable','unavailable',(oligo('hairpin',HAIRPIN,'phosphate','phosphate'),),'Sourcehairpinsequence/blocksretained; exacttRNA/ligationsplintgeometryunassembled.'))
    x=clean(g,'block_cleanup',b,'Sourcephenol/chloroform+EtOHaftertRNAblocking')
    adapter=FirstAdapter();a=g.input('first_adapter','PublishedCEW3165phosphate/Amino-C7DNAadapterwith8N-UMI',K.ADAPTER,state=adapter.scene(),chemical_state=adapter);atp=g.input('adapter_atp','SourceATPforTS2126RNL1preadenylation',K.MOLECULAR);cat=g.input('adapter_adenylase','SourceTS2126RNL1',K.CATALYST);aden=g.product('adenylated_adapter_mix','SourcepreadenylatedrAppDNAadapter plusATP/enzyme/byproducts untilcleanup',K.MOLECULAR,state=adapter.adenylate().scene(),chemical_state=adapter.adenylate())
    g.react('adenylate','Sourceadapterpreadenylation',S.CHEMICAL_CONVERSION,((a,R.TEMPLATE),(atp,R.STRUCTURAL),(cat,R.CATALYST)),(aden,));aden=clean(g,'adenylation_cleanup',aden,'Sourcephenol/chloroformandEtOHadenylatedadapterrecovery')
    cat=g.input('first_ligase','SourceAtRNLK152A/D726AvariantAA' if end=='cP' else 'SourceT4RNL2TKQ',K.CATALYST);y=g.product('first_ligated','RNA-DNAadapterclasses withsourceUMI/Amino-C7retained; '+('RNAjunction2phosphateretainedbeforeTpt1' if end=='cP' else '3OHligationof5OH-enrichedRNA, notdirectligationtoits5OH')+';unligatedRNA/adapter/remnantsretained',K.MOLECULAR)
    g.react('first_ligation','Source40pmolrAppUMIadapterfirstligation',S.LIGATION,((x,R.TEMPLATE),(aden,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=AdapterLigationGeometry('unavailable','unavailable',(adapter.adenylate().scene(),),'SourceRNAend/ligaseknown; populationfragmentends/junctiongeometryunassembled.'))
    x=y
    if end=='cP':
        nad=g.input('nad','Source10mMNADphosphotransfercofactor; notRNA/DNAsequence',K.MOLECULAR);cat=g.input('tpt1_enzyme','SourceS.cerevisiaeTpt1',K.CATALYST);y=g.product('junction_dephosphorylated','Selectedjunction2phosphateremoved; RNA/DNAadapterbodies/UMI/Amino-C7retained, NADtransferproducts/remnantsstillpresent',K.MOLECULAR)
        g.react('tpt1','SourceTpt1/NAD30C30min afterAtRNLAA',S.CHEMICAL_CONVERSION,((x,R.TEMPLATE),(nad,R.STRUCTURAL),(cat,R.CATALYST)),(y,));x=y
    if route=='linear':p=g.input('rt_primer','PublishedCEW422RTprimer with5Goverhang relativeRA3core',K.OLIGO,state=oligo('rt_primer',RT_PRIMER))
    else:
        pp=SpacerPrimer();p=g.input('spacer_rt_primer','SourcecircularizationRTprimer; mainexplicit5P, SIrowomits5Pnotation, notestablishedcontradiction',K.OLIGO,state=pp.scene(),chemical_state=pp)
    rt=g.input('ssiv','SourceSuperScriptIV55C30min',K.CATALYST);mix=g.product('rt_mix','SourcecDNAplusoriginalRNA/DNAadaptertemplates/unextendedprimers/remnants beforeNaOH; fullRTextentunassigned',K.MOLECULAR)
    g.react('rt','SourcesourceRNAadapterRT',S.REVERSE_TRANSCRIPTION,((x,R.TEMPLATE),(p,R.STRUCTURAL),(rt,R.CATALYST)),(mix,))
    x=op(g,'rna_hydrolyze',mix,'Source200mMNaOH95C3minthen200mMHClneutralization',S.CHEMICAL_CONVERSION,reagent='SourceNaOH/HCl',product='cDNAplusRNAhydrolysisproducts, DNAadapter/Amino-C7/fragments/reagentsuntilgel;RNAhydrolysisdoesnoteraseDNAtemplateadapter')
    x=clean(g,'cdna_gel',x,'Source15%denaturingcDNAgelandEtOH' if route=='linear' else 'Source8%denaturingcDNAgelandEtOH')
    if route=='linear':
        a=g.input('second_adapter','SourceCEW4325phosphate/Amino-C6DNAadapter',K.ADAPTER,state=oligo('second_adapter',SECOND_ADAPTER,'phosphate','Amino-C6'));cat=g.input('rnl1','SourceT4RNL1/ATP/hexamminecobalt/PEG16C16h',K.CATALYST);y=g.product('second_ligated','cDNA/secondadapterclasses retainoldtemplates/unligatedadapter/Amino-C6remnants; fulljunctionunassembled',K.MOLECULAR)
        g.react('second_ligation','SourcesecondDNAadapterligationtocDNA3OH then65C10min',S.LIGATION,((x,R.TEMPLATE),(a,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=AdapterLigationGeometry('unavailable','unavailable',(oligo('second_adapter',SECOND_ADAPTER,'phosphate','Amino-C6'),),'SourcecDNA3OH/adapter5P; fullfragment/cDNAgeometryunassembled.'));x=y
    else:
        x=op(g,'circularize',x,'SourceTS2126RNL160C2h and80C5mininactivation',reagent='Sourcecircularizationbuffer/ATP/Mn/DTT; exactspacerpolymerasepathunassembled',product='SourcecircularcDNAclasses retainingtwoHEGnodes plusunclosedtemplates/remnants, notordinaryspacer-freeDNAcircle',state=SpacerPrimer().extend().close().scene(),chemical_state=SpacerPrimer().extend().close())
    # Source PCR physically samples preparation; keep the actual unallocated remainder explicit.
    taken=g.product('pcr_input','SourceallocatedreactionaliquotforPCR',K.MOLECULAR);remaining=g.product('reaction_remainder','UnallocatedsourceRT/ligation/circlepreparation; laterfateunassigned',K.MOLECULAR);g.react('allocate_pcr','Source5ulsecond-ligationreaction or2.5–5ulcircularspikeinreaction allocatedtoPCR',S.SPLIT,((x,R.TEMPLATE),),(taken,remaining));g.mark_terminal(Terminal(remaining.id,Disposition.UNRESOLVED,'No source storage/discardasserted.'))
    if route=='linear':
        f=g.input('pcr_forward','PublishedCEW116Illuminaforwardprimer',K.OLIGO,state=primer('f'));r=g.input('pcr_reverse','PublishedIlluminai7reverseprimerfamily;6basememberunassigned',K.OLIGO,state=primer('r'));cat=g.input('phusion','SourcePhusionGCbuffer16cycles',K.CATALYST);amp=g.product('amplified','PCRDNAcopies plusoriginalblockedtemplates/cDNA/primers/reagents; originalAmino-C6/C7notmovedintonewunmodifiedcopies',K.MOLECULAR)
        g.react('pcr','SourceTruSeqcompatibleprimerPCR',S.PCR,((taken,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(Scene(),reason='Exactfragment/product/fullprimerlanding/readgeometryunassembled; declaredsourcefamiliesonly.'));x=amp
    else:x=op(g,'pcr',taken,'SourcespikeinAmpliconEZpartialadapterPCR',S.PCR,reagent='SourceCEW317/318; fullspacertraversal/endgeometryunassembled')
    x=qc_split(g,'agarose_qc',x,'Source3%agarosePCRproductanalysis; analyticalfractionseparatefromproduction')
    x=clean(g,'library_cleanup',x,'SourceQIAquickPCRcleanuporVAHTSbeadslinearlibrary;QIAquickspikeinlibrary')
    g.sequence(x,SequencingPlan(acquisition='unspecified',limitation='SourceD.vulgarisHiSeq150bpPE~30millionreads; circularspikeinAmpliconEZ~10000–25000reads, sourcegeometryunassembled. Endcapturedoesnotproveribozymeorigin;UMIdeduplicationandcleavagesignalscomputational.'))
    g.add_source_item(SourceItem('pending','Citedhairpinblockingmethod52, ligaseTpt1/TS2126 predecessororiginals andhistoricalPCR/readguides','Selectedpublishedarticle/SITableS5supportlocaloperators; completefragment/end/pair/readgeometry,circularHEGtraversal anddirectlinkedoriginals remainunassembled, noabsenceclaim.'))
    return g

SOURCE_OLIGO_EXCEPTIONS={'pcr_reverse': 'Bare source token Index is not expanded by the generic source-text guard. Coordinator independently inspected original published DOCX word/document.xml complete row73 globally (TableS5/table5row46 within table),2026-10-29: exact P7 + Index + GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA; adjacent cell explicitly says Index containing six bases. One primer-family text-match exemption only; actual N6/fixed flanks retained, no absence claim. Local Source/have test checks the whole original row and width.'}

def protocol_graph():return graph_for_scope()
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();a=FirstAdapter('ACGTACGT').adenylate();p=ligate_first('ACGU','cP',a,'AtRNL_AA');q=p.tpt1_selected()
    c('AdenylationdoesnotchangeadapterDNAbody',a.sequence,FirstAdapter('ACGTACGT').sequence)
    c('AtRNLAAligationretainsjunction2PuntilTpt1',(p.junction_two_p,q.junction_two_p),(True,False))
    c('UMIidentitysurvivescopying',q.selected_cdna_extension().segments[0].feature.id,FirstAdapter().scene().strands['first_adapter'].segs[1].feature.id)
    c('Tpt1doesnotdeleteRNAorUMI',(q.rna,q.adapter.umi),('ACGU','ACGTACGT'))
    c('RTcopycontainscomplementofRNAandadapterUMIprefix',q.selected_cdna_extension().top(),revcomp('ACGTTCAGATCGACGTACGT'))
    try:ligate_first('ACGU','cP',a,'T4RNL2_TKQ');bad=False
    except ValueError:bad=True
    c('Selected3OHligasecontractdoesnotpretendcPligation',bad,True)
    binding=rt_binding();c('RTprimerbindingconstructedwithrealRA3core','rt_primer' in binding.strands,True)
    circle=SpacerPrimer().extend().close();c('CircularspikeinretainsbothactualHEGnodes',sum(k=='HEG' for k,v in circle.nodes),2)
    try:circle.flatten();bad=False
    except ValueError:bad=True
    c('CannotflattenHEGcircletoordinaryDNA',bad,True)
    g=protocol_graph();h=graph_for_scope('OH');c('5OHcaptureusesend-enrichmentnotPEG','rpph_xrn1' in h.steps and 'peg_select' not in h.steps,True)
    for end in ('cP','OH'):
        spiked=graph_for_scope(end,'circular_spikein')
        c('Independentphysicalspikepoolbefore'+end+'enrichment',(spiked.steps['spike_pool'].kind,{i.material for i in spiked.steps['spike_pool'].inputs},spiked.steps['peg_select' if end=='cP' else 'rpph_xrn1'].inputs[0].material),(S.POOL,{'total_rna','ribozyme_spike'},'spiked_total_rna'))
    c('Linearendogenouslibrarydoesnotinventspikeinput','ribozyme_spike' not in g.materials,True)
    c('SourcepostRTmixtureprecedesRNAhydrolysis',g.steps['rna_hydrolyze'].inputs[0].material,'rt_mix')
    c('SourceunusedPCRreactionremainderunresolved',g.terminals['reaction_remainder'].disposition,Disposition.UNRESOLVED)
    src=Source(sources_root()/'cyphyrna-seq__10.1080+15476286.2021.1999105'/'15476286.2021.1999105_supp_KRNB_A_1999105_SM4973_20211015_Supplementary_Data_RNA_Biology_final.docx.txt','get_sources --slugs cyphyrna-seq__10.1080+15476286.2021.1999105; TableS5originalXMLcellsverified')
    if have(c,src,label='publishedoligos'):
        t=''.join(src.path.read_text().split());c('Sourcefirstadapter8NandbothAminoblocks','TCAGATCGNNNNNNNNTGGAATTCTCGGGTGCCAAGG' in t and 'Amino-C7' in t and 'Amino-C6' in t,True)
    original=Source(src.path.with_suffix(''),'get_sources --slugs cyphyrna-seq__10.1080+15476286.2021.1999105; actualpublishedDOCXTableS5/table5row46')
    if have(c,original,label='wholeoriginalprimerfamily'):
        import zipfile,xml.etree.ElementTree as ET
        ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
        with zipfile.ZipFile(original.path) as archive:
            document=ET.fromstring(archive.read('word/document.xml'))
        rows=[]
        for row in document.findall('.//w:tr',ns):
            cells=[''.join(t.text or '' for t in cell.findall('.//w:t',ns)) for cell in row.findall('w:tc',ns)]
            if cells and cells[0]=="Illumina 3'-PCR primer":rows.append(cells)
        c('Exactlyonesourceprimerfamilyrow',len(rows),1)
        if rows:
            written=''.join(rows[0][1].replace("5'-",'').replace("-3'",'').split())
            c('WholeoriginalIndexprimerfixedflanksmatchcanonicalsource',written,P7+'Index'+'GTGACTGGAGTT'+revcomp(TRUSEQ_SMALL_RNA_RA3))
            c('Originaladjacentcellestablishessixbaseindex','Index containing six bases' in rows[0][2],True)
    c.report()
if __name__=='__main__':selftest()
