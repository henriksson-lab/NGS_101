"""Private source-qualified partial preparation; no complete molecular/absence claim."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry
from illumina import P5,P7,TRUSEQ_READ1,TRUSEQ_READ2
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
    def clean(self,key,label):
        y=self.g.product(key+'_kept','Actual retained physical source fraction, carryover/recovery unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual excluded workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label)),(y,w));self.g.discard(w,'Actual source workup fraction excluded; no quantitative purification claim.');self.x=y;return y
    def aliquot(self,key,label,kept='Selected source aliquot',rest='Untransferred stock, amount/fate unassigned',qc=False):
        y=self.g.product(key+'_selected',kept,K.MOLECULAR);r=self.g.product(key+'_remainder',rest,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r))
        if qc:self.g.mark_terminal(Terminal(y.id,Disposition.QC,kept));self.x=r;return r
        self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,rest));self.x=y;return y
    def pcr(self,key,label,forward,reverse,**kw):
        # Source one bulk preparation: produced components are copresent, not physical fractions.
        chemicals=kw.pop('primer_chemicals',{})
        stock=self.g.input(key+'_stock','Actual source prepared PCR reagent stock; quantitative internal recipe specified in step, components same tube',K.MOLECULAR)
        specs=[('f',K.OLIGO,forward,R.STRUCTURAL),('r',K.OLIGO,reverse,R.STRUCTURAL),('enzyme',K.CATALYST,None,R.CATALYST),('nt',K.MOLECULAR,None,R.STRUCTURAL),('buffer',K.CONDITION,None,R.CONDITION)]
        out=[];ins=[]
        for suffix,kind,state,role in specs:
            z=self.g.product(key+'_'+suffix,'Same prepared PCR stock '+suffix,kind,state=state,chemical_state=chemicals.get(suffix));out.append(z);ins.append((z,role))
        self.g.react(key+'_components','Copresent components of ONE source PCR stock, not extra additions',S.GENERIC,((stock,R.TEMPLATE),),tuple(out),changes_molecule=False)
        return self.op(key,label,S.PCR,ins,geometry=PCRGeometry.unavailable(unknown_scene(),reason='Whole occupied primer landing and heterogeneous amplicon products unassembled; source exact primers retained where assigned.'),**kw)
    def finish(self,key,label,plan):
        y=self.g.product(key,label,K.LIBRARY);self.g.react(key+'_load',label,S.GENERIC,((self.x,R.TEMPLATE),),(y,),changes_molecule=False);self.g.sequence(y,SequencingPlan(limitation=plan));return y

def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('whole source preparation unassembled','[unassembled]',placeholder=True)]);return s

def oligo(seq,name='whole source oligo',mod5=None,mod3=None):
    import re
    parts=[Segment(name+' '+str(i),part,placeholder=bool(re.search('[NRYKMSWBDHV]',part))) for i,part in enumerate(re.split('([NRYKMSWBDHV]+)',seq)) if part]
    s=Scene();s.strand('oligo',parts,mod5=mod5,mod3=mod3);return s

FIRST=TRUSEQ_READ2[:-1]+'T'*30+'V'
SECOND_PREFIX='C'+TRUSEQ_READ1[-21:]
SECOND=SECOND_PREFIX+'N'*5+'GGG'
PCR_P5=P5+TRUSEQ_READ1[4:]
PCR_P7=P7+'N'*6+TRUSEQ_READ2
SOURCE_ORDERINGS={'first':'/5Biosg/'+FIRST,'second':'/5Biosg/'+SECOND_PREFIX+'N'*5+'rGrGrG','P5':PCR_P5[:-1]+'*'+PCR_P5[-1],'P7':P7+'[i7]'+TRUSEQ_READ2[:-1]+'*'+TRUSEQ_READ2[-1]}
AUTHORITY='Published PMC6836736main136–146/154 plus SupplementalFile2 printedJan8,2019 standardintactRNA/option2/actualoligopage10. No FFPE/LCM merge.'
def first_scene():return oligo(FIRST,mod5='biotin')
def second_scene():
    s=Scene();s.strand('tso',[Segment('source DNA handle',SECOND_PREFIX),Segment('source random molecularUMI','N'*5,placeholder=True,feature=feature('smart3seq.umi','umi','random')),Segment('source terminalrG3','GGG',note='Three source riboguanosines, not moleculeUMI')],mod5='biotin');return s

def p7_scene():
    s=Scene();s.strand('primer',[Segment('source P7',P7),Segment('source i7 six class','N'*6,placeholder=True,feature=feature('smart3seq.i7','sample_index','unknown')),Segment('source reverse readtail',TRUSEQ_READ2)]);return s

def protocol_graph():
    b=B('SMART-3SEQ standard RNA andpre-SPRIpool preparation',AUTHORITY,'Source1ul purifiednondegradedtotalRNA; concentration andactualbatchN unassigned');g=b.g
    first=g.input('first_primer','Source.5ul20uM1Swholebiotinprimer',K.OLIGO,state=first_scene(),chemical_state={'source_ordering':SOURCE_ORDERINGS['first'],'five_prime':'biotin','terminal_T_run':30})
    b.op('fragment','SourceRNAfragment80C5min in5ultotal beforeRT',S.FRAGMENTATION,((first,R.STRUCTURAL),b.stock('fragment_nt','Source1ul10mMeachdNTPs; retainedforRT'),b.cond('frag_buffer','Source2ul5xSMARTScribebuffer+.5ul80mMMgCl2;4ulmix added to1ulRNA')))
    tso=g.input('second_primer','Source.5ul20uM2SwholebiotinRNA-endedoligo',K.OLIGO,state=second_scene(),chemical_state={'source_ordering':SOURCE_ORDERINGS['second'],'terminal_sugar':'three riboguanosines','five_prime':'biotin'})
    b.op('rt_ts','Sourceadd5ulTSRTmix at42C then42C30min70C10min',inputs=((tso,R.STRUCTURAL),b.cat('smartscribe','Source1ul100UulSMARTScribe addedonce'),b.cat('rt_inhib','Source.5ul20xRNaseinhibitor'),b.cond('tsrt_buffer','Source2ul5MTMG/1ul20mMDTT; prior1S/dNTP retainedfromfragmentmixture')),product='SourcefragmentedRNA/cDNA/template-switch populations, terminalbiotin/UMI-bearing classes; whole junctions and copying products unassembled')
    b.pcr('pcr','Sourceadd2.5ul2uMeachP5/P7pair+12.5ul2xKAPA;25ultotal;cyclesdependonRNAandN;98C15s/60C30s/72C10s final72C60s',oligo(PCR_P5),p7_scene(),primer_chemicals={'f':{'source_ordering':SOURCE_ORDERINGS['P5'],'PS_bond_before_last_base':True},'r':{'source_ordering':SOURCE_ORDERINGS['P7'],'PS_bond_before_last_base':True}})
    # The source selects entire prepared PCR tubes for pre-SPRI pooling.
    other=g.input('other_indexed_pcr','ActualotherdistinctlyindexedPCRproducts prepared bysamepublishedbranch; number/inputamounts/indexchoices unassembled',K.MOLECULAR,primary=True)
    pool=g.product('pre_spri_pool','ActualcombinedindexedPCRsamples beforeSPRI, N unassigned; not equimolarityassertion',K.MOLECULAR)
    g.react('pool','Sourceoption2combinecooledPCRsamples, eachdifferentindex',S.POOL,((b.x,R.TEMPLATE),(other,R.TEMPLATE)),(pool,));b.x=pool
    b.clean('spri1','Source15xNulSPRIbeads/5minRT/magnet; discardactualsupernatant')
    b.clean('respri','Source96ulReSPRI60ulwater+36ulbeadsorbindingbuffer/5minRT; discardactualsupernatant;2x1ml80%ethanol/airdry2min/elute10ul')
    b.aliquot('tapestation','SourceTapeStationHSsizeQC analyticalamountunassigned',qc=True)
    b.aliquot('qpcr','SourceqPCRdualprobequantification analyticalamountunassigned',qc=True)
    b.aliquot('loading','Sourcequantifiedsequencingpoolallocation; equimolnormalizationamounts/membershipunassembled','Actualsequencingpoolaliquot','Untransferredsourcepurifiedindexedpool')
    phix=g.input('phix','Actualsource1%PhiX174DNAcontrolsequencinglibrary, concentration/volumeunassigned',K.MOLECULAR)
    mix=g.product('spiked_pool','ActualDNAlibrarysequencingmixture withPhiXcontrol',K.MOLECULAR)
    g.react('spike','Sourcephysical1%PhiXcontroladdition tosequencingpreparation',S.GENERIC,((b.x,R.TEMPLATE),(phix,R.SPIKE_IN),b.cond('loading_buffer','Actualinstrumentloadingdilutionbuffer')),(mix,),changes_molecule=False);b.x=mix
    b.finish('library','SMART-3SEQselectedNextSeqreferenceRNApool','PublishedreferenceRNA NextSeq500R1=86/i7=6; first5R1basesUMI,next3G. Whole occupied primerlandings/TSjunctions andactualindices unassembled; LCM76/MiSeq169 not substituted.')
    g.add_source_item(SourceItem('source_scope','The1S whole source ordering has30terminalT after a33nttail (canonicalR2 minus itslastT); adding30T afterfullcanonicalR2 would introduceoneextraT.','Wholemodifiedsourcefixtures, no publicationerror asserted.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Fragmentation precedesRT',g.steps['rt_ts'].inputs[0].material,'fragment_product')
    c('OnlyoneSMARTScribeaddition',sum(any(p.material=='smartscribe'for p in step.inputs)for step in g.steps.values()),1)
    c('ActualPCRpool beforeSPRI',g.steps['spri1'].inputs[0].material,'pre_spri_pool')
    c('OtherindexedPCR contributespool',len(g.steps['pool'].inputs),2)
    c('PhysicalPhiXspike',g.steps['spike'].inputs[1].role,R.SPIKE_IN)
    c('NoextraTwholefirstprimer',FIRST,'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC'+'T'*30+'V')
    p=sources_root()/Path(__file__).resolve().parents[1].name/'gr.234807.118_supp_supp_gr.234807.118_Supplemental_File_2.pdf.txt'
    if have(c,Source(p,'https://doi.org/10.1101/gr.234807.118')):
        text=''.join(p.read_text().split())
        for name,seq in SOURCE_ORDERINGS.items():c('Wholeoriginalmodifiedordering '+name,seq in text,True)
    c('Partialnotpublishedready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
