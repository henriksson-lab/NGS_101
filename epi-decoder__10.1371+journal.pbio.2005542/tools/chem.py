from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry
class B:
    def __init__(self,title,authority,start,kind=K.BIOLOGICAL,state=None):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',start,kind,primary=True,state=state)
    def stock(self,key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind,inputs=(),geometry=None,product='Source mixture including incomplete/unreacted/off-target classes; occupancy/recovery unassigned',reason=None):
        y=self.g.product(key+'_product',product,K.MOLECULAR);kw={}
        if kind==S.MAGIC:kw['reason']=reason or 'Positive source operation; complete occupied molecular inputs/products and historical formulation unassembled, not a missing-information conclusion.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label,extra=()):
        y=self.g.product(key+'_kept','Actual source retained physical fraction; yield/carryover unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual separated source workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label))+tuple(extra),(y,w));self.g.discard(w,'Actual workup fraction excluded; no perfect specificity or molecular yield claim.');self.x=y
    def split(self,key,label,selected,remainder):
        y=self.g.product(key+'_selected',selected,K.MOLECULAR);r=self.g.product(key+'_remainder',remainder,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r));self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,remainder));self.x=y;return y

def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('source products unassembled','[unassembled]',placeholder=True)]);return s
AUTHORITY='Published PLOS2018 PMC6059479 Methods174/178/182, Fig1 original PDFp3; actual S1Fig TIFF and S7Table XLSX. HO/KanMX route only, not ADE2.'
ORDERS={
'2B_UP_Fw':'AATGATACGGCGACCACCGAGATCTACACCTTCCGATCTAGATGTCCACGAGGTCTCT',
'2B_DN_Fw':'AATGATACGGCGACCACCGAGATCTACACTCTTCCGATCTACGGTGTCGGTCTCGTAG',
'2B_UP_Rv_all':'CAAGCAGAAGACGGCATACGAGATNNNNNNGTCGACCTGCAGCGTACG',
'2B_DN_Rv_all':'CAAGCAGAAGACGGCATACGAGATNNNNNNAACGAGCTCGAATTCATCGA',
'2B_UP_Seq':'ACACCTTCCGATCTAGATGTCCACGAGGTCTCT',
'2B_DN_Seq':'ACACTCTTCCGATCTACGGTGTCGGTCTCGTAG'}
def oligo(name):
    sc=Scene();s=ORDERS[name]
    if 'NNNNNN' in s:
        a,z=s.split('NNNNNN');parts=[Segment(name+' P7-side',a),Segment('Six-base sample index','NNNNNN',feature=feature('epidecoder.sample','sample_index','unknown',note='Input/IP plate index; selected values unassigned')),Segment(name+' genomic foot',z)]
    else:parts=[Segment(name,s)]
    sc.strand(name,parts);return sc

def qc(b,key,label):
    y=b.g.product(key+'_qc',label,K.MOLECULAR);r=b.g.product(key+'_stock','Stock remaining after analytical allocation; amount unassigned',K.MOLECULAR);b.g.react(key,label,S.SPLIT,((b.x,R.TEMPLATE),),(y,r));b.g.mark_terminal(Terminal(y.id,Disposition.QC,label));b.x=r

def genomic_reporter():
    sc=Scene();sc.strand('HO reporter partial context',[Segment('BC_UP20','[20-base-BC_UP]',placeholder=True,feature=feature('epidecoder.ho.up','feature_barcode','unknown',note='Clone-linked tagged-protein reporter identity, not cell/UMI')),Segment('Intervening KanMX genomic context','[unassembled HO/KanMX interval]',placeholder=True),Segment('BC_DN20','[20-base-BC_DN]',placeholder=True,feature=feature('epidecoder.ho.dn','feature_barcode','unknown',note='Paired clone-linked locus reporter; exact barcode lookup deferred'))]);return sc

def protocol_graph():
    b=B('Epi-Decoder HO reporter: input and TAP-ChIP barcode libraries',AUTHORITY,'Actual pool of yeast clones: paired HO/KanMX genomic barcodes and one endogenous TAP-tagged protein per clone',state=genomic_reporter());g=b.g
    b.op('crosslink','One-tenth culture volume FixSolution, 20min; source stock1% formaldehyde/50mMHEPES/100mMNaCl/1mMEDTA',S.MAGIC,(b.stock('fix','Actual freshly prepared FixSolution allocation containing formaldehyde and buffer; source final-concentration convention unassigned',K.MOLECULAR,R.STRUCTURAL),b.cond('fix_conditions','Log culture OD660~0.4; actual volume ratio retained without silently interpreting source1% as finalculture concentration')))
    b.op('quench','125mM final glycine,5min',S.MAGIC,(b.stock('glycine','Actual glycine125mM final',K.MOLECULAR,R.STRUCTURAL),))
    b.clean('cells','ColdTBS/.2mMPMSF wash; pellet frozen -80C')
    b.op('lyse','~1.5e9 cells bead-beaten in200ul breakingbuffer',S.MAGIC,(b.stock('lysis_supplies','Actual zirconia/silica beads and100mMTris7.9/20%glycerol/EDTA-freePIC supplies; no presumed premix',K.MOLECULAR,R.STRUCTURAL),))
    b.clean('fa','Two1mlFA washes: HEPES/NaCl/EDTA/Triton/deoxycholate/PIC')
    b.op('sonicate','BioruptorPICO10min,30s intervals',S.MAGIC,(b.cond('sonication','Source chromatin preparation; full fragments and crosslink occupancy unassembled'),))
    b.clean('clarify','Clarify4000rpm5min4C, retain chromatin lysate')
    inp=g.product('input100','Actual100ul INPUT chromatin allocation',K.MOLECULAR);ip=g.product('ip1000','Actual1ml chromatin for HO IgG-Sepharose IP',K.MOLECULAR);rem=g.product('chromatin_remainder','Untransferred clarified chromatin; amount unassigned',K.MOLECULAR)
    g.react('chromatin_allocate','Separate source INPUT and IP allocations',S.SPLIT,((b.x,R.TEMPLATE),),(inp,ip,rem));g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,'Untransferred chromatin'))
    rawbeads=g.input('igg_beads_raw','Actual IgGSepharose6FastFlow allocation before source prewash',K.MOLECULAR)
    pbsb=g.input('pbsb_stock','Actual aggregate PBSB stock allocations for three bead prewashes, PBS+5mg/ml BSA; total volume unassigned',K.MOLECULAR)
    bsa=g.product('pbsb_bsa','BSA molecular component of same actual PBSB wash allocations',K.MOLECULAR);pbs=g.product('pbsb_buffer','PBS component of same actual PBSB wash allocations',K.CONDITION)
    g.react('pbsb_parts','Expose copresent BSA and PBS from same source wash stock',S.GENERIC,((pbsb,R.TEMPLATE),),(bsa,pbs),changes_molecule=False)
    washed=g.product('igg_beads_washed','Actual three-times PBSB-prewashed IgGSepharose; residual BSA carryover unassigned',K.MOLECULAR);washout=g.product('igg_prewash_removed','Actual removed PBSB prewash fractions',K.WASTE)
    g.react('igg_prewash','Three source PBSB prewashes BEFORE chromatin capture',S.CLEANUP,((rawbeads,R.TEMPLATE),(bsa,R.STRUCTURAL),(pbs,R.CONDITION)),(washed,washout));g.discard(washout,'Actual removed prewash fractions; no assertion zero BSA carryover on washed beads')
    b.x=ip
    b.op('tap_capture','IgGSepharose6FastFlow;6h4C turning',S.MAGIC,((washed,R.STRUCTURAL),),reason='Actual IgG capture of TAP-associated chromatin. Protein/IgG are not DNA; occupied crosslinks and recovery remain unassembled.')
    b.clean('ip_wash','TwoFA/two500mMNaClFA/twoRIPA/oneTE washes; retain Sepharose')
    b.clean('ip_elute','Actual100ul50mMTris8/10mMEDTA/1%SDS,65C10min; retained eluate versus beads')
    purified={}
    for fate,stock in [('input',inp),('ip',b.x)]:
        b.x=stock
        b.op(fate+'_deproteinize','Separate INPUT or IP:0.5ulRNaseA10mg/ml+10ulProteinaseK10mg/ml in70ulTE;50C1h then65Covernight',S.MAGIC,(b.stock(fate+'_rnase','Actual0.5ulRNaseA10mg/ml'),b.stock(fate+'_pk','Actual10ulProteinaseK10mg/ml'),b.cond(fate+'_te','Actual70ulTE allocation; reverse-crosslink programme')))
        b.clean(fate+'_column','QIAquickPCRpurification');purified[fate]=b.x
    products=[]
    for fate,stock in purified.items():
        up=g.product(fate+'_up_template','Selected BC_UP PCR template allocation',K.MOLECULAR);dn=g.product(fate+'_dn_template','Selected BC_DN PCR template allocation',K.MOLECULAR);r=g.product(fate+'_dna_remaining','Untransferred purified '+fate+' DNA',K.MOLECULAR)
        g.react(fate+'_pcr_allocate','Separate UP/DN PCR allocations; quantities unassigned',S.SPLIT,((stock,R.TEMPLATE),),(up,dn,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,'Untransferred DNA'))
        for locus,x in [('UP',up),('DN',dn)]:
            k=fate+'_'+locus.lower();b.x=x
            primers=[b.stock(k+'_'+end,'Actual '+name,K.OLIGO,R.STRUCTURAL,state=oligo(name)) for end,name in [('fw','2B_'+locus+'_Fw'),('rv','2B_'+locus+'_Rv_all')]]
            supplies=(b.stock(k+'_pol','Actual PCR polymerase functional allocation; identity/dose/packaging unassigned'),b.stock(k+'_supplies','Actual PCR nucleotide/buffer supplies aggregate; identities/packaging unassigned',K.MOLECULAR,R.STRUCTURAL))
            b.op(k+'_pcr','Separate '+fate+' BC_'+locus+' PCR; plate-specific six-base index',S.PCR,tuple(primers)+supplies,geometry=PCRGeometry.unavailable(unknown_scene(),reason='Whole source primers declared, but selected genomic barcode/flanking allele and occupied binding/junction states unassembled; no canonical genome or actual index value borrowed.'))
            b.split(k+'_pool_select','Actual equimolar contribution selected from this PCR','Selected PCR pool contribution','Untransferred PCR stock');products.append(b.x)
    pool=g.product('amplicon_pool','Actual equimolar pool of four input/IP UP/DN PCR contributions',K.MOLECULAR);g.react('equimolar_pool','Actual source equimolar PCR pooling',S.POOL,tuple((x,R.TEMPLATE)for x in products),(pool,));b.x=pool
    b.clean('gel','Source100–150bp agarose fraction; recovery/other molecular classes unassigned')
    b.split('master_select','Selected share for indexed sequencing masterpool','Selected sample pool contribution','Untransferred size-selected stock')
    b.op('masterpool','Source indexed multiple-pool sequencing masterpool',S.MAGIC,(b.stock('other_pools','Other actual indexed sample-pool contributions; count/amounts unassigned',K.MOLECULAR,R.STRUCTURAL),))
    b.split('load_select','Selected final pooled loading allocation','Actual loading aliquot','Untransferred sequencing masterpool')
    seq=[b.stock('read_'+loc,'Actual custom sequencing primer '+name,K.OLIGO,R.STRUCTURAL,state=oligo(name))for loc,name in [('up','2B_UP_Seq'),('dn','2B_DN_Seq')]]
    lib=g.product('library','Selected barcode library loading mixture; actual complete allele/index geometry unassembled',K.LIBRARY);g.react('load','Custom UP/DN primer mixture; single-end HiSeq2500/MiSeq',S.MAGIC,((b.x,R.TEMPLATE),)+tuple(seq),(lib,),reason='Source custom primers one or mixture; selected mixed UP/DN library uses both. Complete selected genomic insert and landing/read-index allocation unresolved.')
    g.sequence(lib,SequencingPlan(read_lengths={},limitation='Source single-end >50bp; exact cycle count unassigned; index read-cycle allocation and whole allele remain unassembled. Computational six-index/20barcode extraction is not a physical trim.'))
    g.add_source_item(SourceItem('source_scope',AUTHORITY,'Private source-qualified partial graph; barcode identifies clone-linked protein/locus, not a cell or molecular deduplication UMI.'))
    g.add_source_item(SourceItem('run_scope','Source single-end >50bp; exact cycle count unassigned.','No invented integer cycle count or sequencing construct/landing proof.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    import zipfile,xml.etree.ElementTree as ET
    c=Check();g=protocol_graph()
    c('Actual input and IP allocated before capture',g.producer_of('input100').id,g.producer_of('ip1000').id)
    c('Capture consumes prewashed beads, not fresh PBSB',g.steps['tap_capture'].inputs[1].material,'igg_beads_washed')
    c('BSA wash component from actual PBSB stock',g.producer_of('pbsb_bsa').id,'pbsb_parts')
    c('Removed PBSB wash not carried as fresh capture addition',g.terminals['igg_prewash_removed'].disposition,Disposition.DISCARD)
    c('Only IP stock captured',g.steps['tap_capture'].inputs[0].material,'ip1000')
    c('Four source PCR categories',len([s for s in g.steps.values()if s.kind==S.PCR]),4)
    c('No early input/IP pooling',g.steps['input_deproteinize'].inputs[0].material,'input100')
    c('Only selected PCR contributions pooled',[x.material for x in g.steps['equimolar_pool'].inputs],[f+'_'+l+'_pool_select_selected'for f in ('input','ip')for l in ('up','dn')])
    c('Final loading has stock remainder','load_select_remainder'in g.terminals,True)
    c('Custom primers are suffixes of actual forward orders',all(ORDERS['2B_'+x+'_Fw'].endswith(ORDERS['2B_'+x+'_Seq'])for x in ('UP','DN')),True)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'journal.pbio.2005542_supp_pbio.2005542.s012.xlsx'
    if have(c,Source(p,'https://doi.org/10.1371/journal.pbio.2005542.s012')):
        with zipfile.ZipFile(p)as z:
            ss=ET.fromstring(z.read('xl/sharedStrings.xml'));ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'};strings=[''.join(t.itertext())for t in ss.findall('m:si',ns)]
            for name,seq in ORDERS.items():c('Whole original order '+name,seq in strings,True)
    c.report()
if __name__=='__main__':selftest()
