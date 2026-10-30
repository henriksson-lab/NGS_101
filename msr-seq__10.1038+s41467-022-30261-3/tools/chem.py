"""Private source-qualified partial preparation; no complete molecular/absence claim."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature,revcomp
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
class B:
    def __init__(self,title,authority,label,kind=K.MOLECULAR):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',label,kind,primary=True)
    def stock(self,key,label,kind=K.MOLECULAR,role=R.STRUCTURAL,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cat(self,key,label):return self.stock(key,label,K.CATALYST,R.CATALYST)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind=S.MAGIC,inputs=(),product='Actual source preparation; whole occupied substrates/products/side products and conversion unassembled',**kw):
        y=self.g.product(key+'_product',product,K.MOLECULAR)
        if kind==S.MAGIC:kw.setdefault('reason','Positive source operation; whole occupied molecules, end distributions, detailed molecular products and conversion unassembled, not a source-absence claim.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label,inputs=()):
        y=self.g.product(key+'_kept','Actual retained physical source fraction, carryover/recovery unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual excluded workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label))+tuple(inputs),(y,w));self.g.discard(w,'Actual source workup fraction excluded; no quantitative purification claim.');self.x=y;return y
    def aliquot(self,key,label,kept='Selected source aliquot',rest='Untransferred stock, amount/fate unassigned',qc=False):
        y=self.g.product(key+'_selected',kept,K.MOLECULAR);r=self.g.product(key+'_remainder',rest,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r))
        if qc:self.g.mark_terminal(Terminal(y.id,Disposition.QC,kept));self.x=r;return r
        self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,rest));self.x=y;return y



from adapter_ligation import AdapterLigationGeometry
from illumina import P5,P7,TRUSEQ_READ1,TRUSEQ_READ2
from hairpin import Hairpin
from dataclasses import dataclass
import re
SOURCE_ORDERS={'cho': '/5Phos/rACT GGAA AGA TCG Gaa gag cac acg at /iBiodT/agaCGT GTG CTC TTC CGA TCT TTCC AG rU/3Phos/', 'l2':'/5Phos/NN NNN NGA TCG TCG GAC TGT AGA A/3ddC/','pcr_f':'AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC','pcr_r':'CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'}
SOURCE_OLIGO_EXCEPTIONS={'cho':'Published Data5 Read2-4nt CHO contains /iBiodT/, an actual biotinylated deoxythymidine nucleotide. Model retains its physical T at zero-based position27 plus terminal rU. Generic text parsing drops the internal modified nucleotide; whole original-order fixture and model-polymer/token comparison verify the complete ordering. Narrow modified-nucleotide transcription exception, not chemistry or source absence.'}
AUTHORITY='PublishedPMC9072684 main49–65/165–203/301; actualData5MOESM8 Read2-4ntCHO/L2/customPCR wholeorders; SI1Fig1controls. Standarddeacylation/shortRT branch; noCMC/AlkB/charging/overnightRT/mRNA merge.'
STEM5='AGATCGGAAGAGCACACG';STEM3=revcomp(STEM5)
BARCODE=feature('msr.capture_sample','sample_index','whitelist',whitelist='PublishedData5Read2four-baseCHO12-memberfamily',note='Source selectedreferenceTTCC; actualexperimentalCHO/indexpairingunassigned')
INDEX=feature('msr.pcr_sample','sample_index','whitelist',whitelist='PublishedData5customPCR24indexes',note='SelectedsourceI1CGTGAT reference, notcommercialRPLprimerbody')
@dataclass(frozen=True)
class SelectedHairpin(Hairpin):
    inferred:bool=True

def cho_parts():
    return [Segment('5prime rA','A'),Segment('5prime DNA CT','CT'),Segment('barcode complement','GGAA',feature=BARCODE,note='Reversecomplementcopyofsource4ntsamplecode'),Segment('left stem',STEM5),Segment('loop prefix','AT'),Segment('internal biotin dT','T',note='Actual nucleotide with internal biotin modification'),Segment('loop suffix','AGA'),Segment('right stem',STEM3),Segment('read2 sample code','TTCC',feature=BARCODE),Segment('terminal DNA AG','AG'),Segment('terminal rU','U')]
def tokens(name):
    raw=''.join(SOURCE_ORDERS[name].split());raw=raw.replace('/iBiodT/','{BIOTIN_T}');raw=re.sub('/[^/]+/','',raw)
    out=re.findall(r'\{BIOTIN_T\}|r[ACGUTN]|[ACGUTN]',raw,flags=re.I)
    seq=''.join('T'if x=='{BIOTIN_T}'else x[-1].upper() for x in out)
    if name=='l2':seq+='C'
    return seq

def fold_candidate():
    seq=''.join(p.top for p in cho_parts());return SelectedHairpin('INFERRED source-compatible CHO core; known5prime7ntoverhangoutsidecore',seq[7:],18,overhang_3_len=7)

def oligo_scene(name):
    s=Scene()
    if name=='cho':
        fold_candidate();s.strand('oligo',cho_parts(),mod5='P',mod3='P');s.ligation_end('oligo','5prime rA',kind='single-strand',polarity='5-prime',owner='adapter',text='Sourceaccessible5P/rA CHO ligationend; fullRNA/hairpinoccupiedproductunassembled')
    elif name=='l2':
        seq=tokens(name);s.strand('oligo',[Segment('N6complexitywindow','N'*6,placeholder=True,note='BaselineData5diversitywindow; CMCspecificUMImeaningnotautomaticallyinherited'),Segment('sourceL2fixedDNA',seq[6:-1]),Segment('terminalddCnucleotide','C')],mod5='P',mod3='ddC');s.ligation_end('oligo','N6complexitywindow',kind='single-strand',polarity='5-prime',owner='adapter',text='Source5P second-adapter end; occupiedcDNA/productunassembled')
    elif name=='pcr_f':
        seq=SOURCE_ORDERS[name];parts=[Segment('canonical P5 body',P5),Segment('source multiplex suffix',seq[len(P5):])]if seq.startswith(P5)else[Segment('source multiplex fullbody',seq)];s.strand('oligo',parts)
    else:
        seq=SOURCE_ORDERS[name];prefix=seq[:24];s.strand('oligo',[Segment('source P7 prefix',P7[:24]),Segment('source I1sampleindex','CGTGAT',feature=INDEX),Segment('sourcecustomPCRbindingbody',seq[30:])])
    return s

def source_oligo(b,key,name):
    chem={'source_ordering':SOURCE_ORDERS[name]}
    if name=='cho':chem.update({'sugar_blocks':[{'start':0,'end':1,'sugar':'RNA'},{'start':1,'end':55,'sugar':'DNA'},{'start':55,'end':56,'sugar':'RNA'}],'modified_bases':[{'index':27,'base':'T','modification':'biotin'}],'five_prime':'P','three_prime':'P','hairpin_core_candidate':{'inferred':True,'five_prime_overhang':7,'stem':18,'loop':6,'three_prime_overhang':7}})
    if name=='l2':chem.update({'five_prime':'P','terminal_nucleotide':'additionalphysicalddC','N6meaning':'baselinecomplexity; CMCcontextUMIseparate'})
    return b.g.input(key,'Source modified '+name,K.OLIGO,state=oligo_scene(name),chemical_state=chem)

def protocol_graph():
    b=B('MSR-seq selectedcapturehairpin/beadretainedsmallRNApreparation',AUTHORITY,'Source up-to1ugtotalRNA, standarddeacylationbranch');g=b.g
    g.add_source_item(SourceItem('selected_reference','SourceData5 offers multiple3/4ntread1/read2CHO and24customPCRindexes. Model selects sourceprovidedRead2-4ntbc1 andI1 aslocalreferenceorders.','ActualstudyRNA/CHO/PCRindexpairingunassigned; not a claimthesetwoorderswereusedtogetherinaspecificrun. No commercialRPLbody orCMC100bp readlength substituted.'))
    g.add_source_item(SourceItem('degenerate_meaning','Data5 describesL2N6forreadcomplexity; main315 callsrandom6aUMI specificallyinCMCprocessing.','BaselineN6purposekeptscoped; sampleCHO4nt andPCR6nt aredistinctsampleindexes. CMCmechanism/run notimportedinto baseline.'))
    b.op('deacylate','Source100mMTris9/37C30min',inputs=(b.cond('deacyl_condition','Actualsourcebasicdeacylation'),))
    b.op('neutralize','SourceNaOAc4.8to180mM',S.QUENCHING,(b.cond('sodium_acetate','Actualneutralizationconditions'),));b.clean('deacyl_clean','Selectsourceethanolprecipitation/water, spincolumnalternativeoutsidetrace')
    cho=source_oligo(b,'cho','cho')
    b.op('first_ligate','Source50ulRnl1/CHO1uM/16C12+h RNA3primebarcoding',S.LIGATION,((cho,R.STRUCTURAL),b.cat('first_rnl1','Actual1UulT4RNAligaseI'),b.stock('first_atp','Actual50uMATP'),b.cond('first_lig_conditions','Sourcebuffer/15percentPEG/1mMcobalt/5percentDMSO')),geometry=AdapterLigationGeometry.known('single-strand',oligo_scene('cho'),statement='Knownsource5P/rA localCHOend; internalbiotindT/3rU3P retained. FullRNA/occupiedhairpin/sealedjunctionunassembled.'),product='SourceRNA-CHO/freeCHO/offtargetpreparation; fullcovalentproducts/occupancyunassembled')
    b.op('dilute','Actualequalvolumewater beforecapture',S.GENERIC,(b.cond('dilute_water','Source50ulwater for50ulreaction'),))
    bound=g.product('capture_bound','ActualCHO/RNAassociatedC1beads, residualfreeCHOmayalsobecaptured',K.MOLECULAR);free=g.product('capture_supernatant','Actualremovedcapture supernatant',K.WASTE)
    g.react('capture','SourceC1RT15min/1.2capacityratio;12ulper50pmolexample',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),b.stock('c1','ActualsourceC1beadallocation',K.SELECTION,R.SELECTION)),(bound,free));g.discard(free,'Actualsourceunboundfractionremoved; no complete RNAcaptureclaim.');b.x=bound;b.clean('capture_wash','Sourcehighsalt1MNaCl/20mMTris7.4 thenlowsalt100mMNaCl/20mMTris7.4')
    b.aliquot('pool_aliquot','Selectedsourcebarcodedbeadcontribution, amountunassigned')
    others=g.input('other_barcoded_bead_aliquots','Otherindividuallybarcoded/captured/washedsourcebeadcontributions; actualsamplecount/CHOorders/amountsunassigned',K.MOLECULAR)
    pool=g.product('barcoded_pool','Actualselectedbarcodedbeadpool',K.MOLECULAR);g.react('pool','Selectsourceoptionalbarcodedbeadpoolingroute',S.POOL,((b.x,R.TEMPLATE),(others,R.TEMPLATE)),(pool,));b.x=pool
    b.op('cip','Source50ul.04UulCIP/10mMMg/.5mMZn/20mMHEPES7.3/37C30min',S.DEPHOSPHORYLATION,(b.cat('cip_enzyme','Actualcalf-intestinephosphatase'),b.cond('cip_conditions','Actualsourceconditions')),product='Sourcebead-associatedCHO3OHclasses readyforRT; fullsealedRNA/hairpinproducts/unreactedclassesunassembled');b.clean('cip_wash','Sourcehighsalt/lowsaltwash,resuspend20ulwater')
    master=g.input('vilo_master','Actual5ulSSIVVILO5xmasteraddedto20ulbeads; formulationunassembled',K.MOLECULAR);enzyme=g.product('vilo_ssiv','SameactualSSIVVILOenzymecomponent',K.CATALYST);nt=g.product('vilo_nt','SameactualVILOnucleotideinventory; exactdoseunassembled',K.MOLECULAR);other=g.product('vilo_other','OtheractualVILOsupplies/conditions, fullcompositionunassembled',K.MOLECULAR)
    g.react('vilo_components','Copresentactualmastercomponents, notnewindependentadditions',S.GENERIC,((master,R.TEMPLATE),),(enzyme,nt,other),changes_molecule=False)
    rtproduct=g.product('rt_product','Sourcehairpin-extendedcDNA/RTstop/fullproducts+unextendedCHO mixture',K.MOLECULAR,chemical_state={'primer_provenance':{'kind':'intramolecularCHO3primeend','whole_covalent_occupied_geometry':'unassembled'},'terminal_classes':'ExtendedDNA3OH versus unextendedCHO3rU3OH; unreacted/TdTclassesretained'})
    g.react('rt','Selectedsource25ulSSIVVILO55C10min; covalentlyretainedCHOselfprimingcontext',S.MAGIC,((b.x,R.TEMPLATE),(enzyme,R.CATALYST),(nt,R.STRUCTURAL),(other,R.STRUCTURAL)),(rtproduct,),reason='SourceCIP/CHOselfprimingpositive; fullRNA/hairpin/cDNA/bondgeometry andVILOinternalsunassembled. No new externalRTprimeror overnightbranchinvented.');b.x=rtproduct;b.clean('rt_wash','Sourcehighsaltthenlowsaltwash')
    b.op('rnase_h','Source50ul.4UulRNaseH37C15min',S.DIGESTION,(b.cat('rnase_h_enzyme','ActualRNaseH'),b.cond('rnase_h_buffer','Actual1xRNaseHbuffer')));b.clean('rnase_h_wash','Sourcehighsalt/lowsaltwash/resuspend40ulwater')
    b.op('periodate','Source10ul250mMperiodate/.5MNaOAc5 added/RT30min',inputs=(b.stock('periodate_reagent','Actualperiodatechemicalcontributor'),b.cond('periodate_conditions','Actualsourceconditions')),product='SourceunextendedterminalrUoxidizedclasses and DNA-endedcDNAremainingligationclasses; fullconversion/offproductsunassembled')
    b.op('quench_periodate','Source167mMriboseRT5min',S.QUENCHING,(b.stock('ribose','Actualribosequenchreagent'),b.cond('ribose_conditions','Source167mMfinal/RT5min')));b.clean('periodate_wash','Sourcehighsalt/lowsaltwash')
    l2=source_oligo(b,'l2','l2')
    b.op('second_ligate','Source50ul2UulRnl1/L22uM/RT12+h',S.LIGATION,((l2,R.STRUCTURAL),b.cat('second_rnl1','ActualfreshT4RNAligaseI2Uul'),b.stock('second_atp','Actual50uMATP'),b.cond('second_lig_conditions','Sourcebuffer/25percentPEG/7.5percentDMSO/1mMcobalt')),geometry=AdapterLigationGeometry.known('single-strand',oligo_scene('l2'),statement='KnownsourceL25P andterminalddC; fullextendedDNA3OH/RNAhairpinoccupiedproductunassembled.'),product='Sourceon-beadsecond-adaptercDNA/freeCHO/TdT/unreactedpreparation; no completeoxidation/rejection/yieldclaim')
    b.op('second_dilute','Sourceadd50ulwater',S.GENERIC,(b.cond('second_water','Actualsourcewateraddition'),));b.clean('second_wash','Sourcehighsalt/lowsaltwash,resuspend6ulperoriginalsample; samplecountunassigned')
    b.aliquot('pcr_template','Actual5–10percentbeadslurryPCRallocation','Actualselected5–10percentbeadslurry','Actualuntransferredbarcodedbeadstock')
    pf=source_oligo(b,'pcr_f','pcr_f');pr=source_oligo(b,'pcr_r','pcr_r')
    b.op('pcr','Source50ulQ5PCR9–15cycles98C10s55C15s72C20s',inputs=((pf,R.STRUCTURAL),(pr,R.STRUCTURAL),b.cat('q5','Actual.02UulQ5polymerase'),b.stock('pcr_nt','Actual.2mMdNTP'),b.cond('pcr_buffer','Actual1xQ5buffer/.5uMeachprimer')),product='Sourceoff-beadPCRamplicons/beads/unreactedclasses; wholePCRproductsunassembled')
    libstock=g.product('off_bead_dna','SourceZymoDNAcleanedoff-beadPCRproduct',K.MOLECULAR);beadremain=g.product('post_pcr_beads','Residualbeads/otherclasses duringoff-beadDNAworkup; exactseparationscheduleunassembled',K.MOLECULAR);waste=g.product('column_waste','Actualexcludedcolumnworkupfractions',K.WASTE)
    g.react('pcr_cleanup','SourceDNAcleanconcentratorafterPCR, productoffbead; exactbeadseparationunassembled',S.MAGIC,((b.x,R.TEMPLATE),b.cond('pcr_cleanup_condition','ActualZymoDNAcleancolumn')),(libstock,beadremain,waste),reason='Sourceoff-beadPCRendpoint/columncleanup, no extra magneticrelease timing invented.');g.discard(waste,'Actualsourcecolumnworkupfractions.');g.mark_terminal(Terminal(beadremain.id,Disposition.UNRESOLVED,'Actualresidualbeadfate outsidepartialtrace.'));b.x=libstock
    b.op('gel','Source10percentnativeTBEdesiredproductsizecut/crush/500mMNaOAc5/overnightextract',inputs=(b.cond('gel_conditions','Actualsourcegel/extractionconditions'),),product='Sourcegel-selectedPCRproducts; exactsizewindowunassigned, TdT/offproductrejectionnotperfectlyvalidated');b.clean('gel_precip','Sourceethanolprecipitationofgel-extractedDNA')
    b.aliquot('loading','Actualsequencingloadingallocation, pool/massvolumeunassigned')
    lib=g.product('library','SourceMSRsmallRNA-seqbarcodedlibraryclass',K.LIBRARY);g.react('load','SourceHiSeq/NEXTseqpaired-endsequencing',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='Publishedmain301baselineHiSeq/NextSeqpaired-end; exactbaselinecycles/indexmask/primerlandingunassembled. CMCspecific100bpandN6UMIcontextnotautomaticallyinherited.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    import zipfile,xml.etree.ElementTree as E
    c=Check();g=protocol_graph();h=fold_candidate()
    c('WholeCHOphysicalorderingwithbiotindT',''.join(p.top for p in cho_parts()),tokens('cho'))
    c('Uncertainfoldcandidatedefault',h.inferred,True)
    c('BiotinisphysicalT',tokens('cho')[27],'T')
    c('CIPafterpool',g.steps['cip'].inputs[0].material,'barcoded_pool')
    c('NoexternalRTprimerstock',all('rt_primer'!=m.id for m in g.materials.values()),True)
    c('RTselfprimerprovenance',g.materials['rt_product'].chemical_state['primer_provenance']['kind'],'intramolecularCHO3primeend')
    c('PeriodateafterRNaseHwash',g.steps['periodate'].inputs[0].material,'rnase_h_wash_kept')
    c('PCRusesactualbeadaliquot',g.steps['pcr'].inputs[0].material,'pcr_template_selected')
    c('Remainingbeadstocknotconsumed','pcr_template_remainder'in g.terminals,True)
    c('OffbeadPCRproductnotresidualbeads',g.steps['gel'].inputs[0].material,'off_bead_dna')
    p=sources_root()/Path(__file__).resolve().parents[1].name/'s41467-022-30261-3_supp_41467_2022_30261_MOESM8_ESM.xlsx'
    if have(c,Source(p,'https://doi.org/10.1038/s41467-022-30261-3')):
        z=zipfile.ZipFile(p);ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'};ss=[''.join(x.itertext())for x in E.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)];actual={''.join(x.split()).upper()for x in ss}
        for k,v in SOURCE_ORDERS.items():c('Wholemodifiedsourceorder '+k,''.join(v.split()).upper()in actual,True)
    c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
