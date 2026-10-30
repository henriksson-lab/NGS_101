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

from nextera import ME,ME_RC
NSR='TCCGATCTGA'+'N'*7
CUSTOM_P5=TRUSEQ_READ1[-23:]+'CT'
OLIGODT=CUSTOM_P5+'T'*30+'VN'
T7_CUSTOM='GATCACTAATACGACTCACTATAGGGTCTCG'
T7_ME=T7_CUSTOM+'TGGGCTCGG'+ME
FINAL_P5=P5+TRUSEQ_READ1[4:]+'CT'
FINAL_P7=P7+'N'*8+'GTCTCGTGGGCTCGGAGATGT'
SOURCE_ORDERINGS={'NSR':NSR,'oligodT':CUSTOM_P5+'T30VN','MErev':'/Phos/'+ME_RC,'T7ME':T7_ME,'customP5':CUSTOM_P5,'T7custom':T7_CUSTOM,'P5':FINAL_P5,'P7':FINAL_P7}
AUTHORITY='Published RAG mainPMC11658833Methods198/214–240, originalTableS2wholeoligos andTableS3guides. Selected3.0single-cell branch, not earlier1/2 orSmart3 comparator.'
def nsr_scene():
    s=Scene();s.strand('nsr',[Segment('source adapter portion','TCCGATCTGA'),Segment('source preceding N','N',placeholder=True,note='Source one N before six-underlined NSR positions; purpose unassigned'),Segment('source NSR six-class window','N'*6,placeholder=True,note='Source underlined non-rRNA priming classes; actual selected pool members unassembled, not a UMI')]);return s

def index_scene():
    s=Scene();s.strand('primer',[Segment('P7',P7),Segment('sample index','N'*8,placeholder=True,feature=feature('rag.i7','sample_index','unknown')),Segment('source right body','GTCTCGTGGGCTCGGAGATGT')]);return s

def loaded_duplex():
    s=Scene();s.strand('top',[Segment('source T7 tail',T7_ME[:-len(ME)]),Segment('source ME',ME)])
    s.anneal('bottom',[Segment('source ME reverse',ME_RC)],to='top',pair=('source ME reverse','source ME'),mod5='P');return s

def protocol_graph():
    b=B('RAG-seq3.0 selected single-cell preparation',AUTHORITY,'Actual isolated selected cell transferred fromPBS/1%BSA; upstream dissociation outside graph',K.BIOLOGICAL);g=b.g
    # Independent loaded stock factory; not a second addition at tagmentation.
    a=g.input('tn5_upper','Source wholeT7-Tn5-MEB',K.OLIGO,state=oligo(T7_ME));r=g.input('tn5_lower','Source5PwholeMErev',K.OLIGO,state=oligo(ME_RC,mod5='P'),chemical_state={'source_ordering':'/Phos/'+ME_RC})
    pair=g.product('tn5_pair','Source annealed19bpMEpair',K.ADAPTER,state=loaded_duplex())
    g.react('tn5_anneal','Sourceequal100uM strands/annealingbuffer94C5min then1Cmin to10C',S.HYBRIDIZATION,((a,R.TEMPLATE),(r,R.STRUCTURAL),b.cond('annealbuffer','Source10mMTris/10mMNaCl')),(pair,))
    stock=g.product('loaded_tn5_stock','Source assembled Tn5 stock, loading formulation/stoichiometry/unbound classes unassembled',K.MOLECULAR)
    g.react('load_tn5','Source transposome loading delegates cited predecessor53; selected original loading unassembled',S.MAGIC,((pair,R.TEMPLATE),b.cat('tn5_enzyme','Actual Tn5 enzyme loading preparation, formulation unassembled')),(stock,),reason='Source loading positive; predecessor original ordinary unread, no SHERRY/Tn5 stock recipe inherited.')
    prior=b.x;b.x=stock;ta=b.aliquot('tn5_input','Source1ulselectedloadedtransposome','Actual1ulloadedtransposome','UntransferredloadedTn5stock, amount/fateunassigned');b.x=prior
    b.op('lysis','Source2ulTriton/RiboLock4U/DNaseI.4Ubuffer20C30min',inputs=(b.cat('dnase','Source.4UDNaseI'),b.cat('lysis_inhibitor','Source4URiboLock'),b.cond('lysis_buffer','Source2ul.5%Triton/1xDNasebuffer')))
    b.op('pk','Source1ul.12UulthermolabilePK25C30min/55C10min',S.DEPROTEINIZATION,(b.cat('pk_enzyme','Source.12UthermolabileProteinaseK'),b.cond('pk_program','Source25C30min then55C10min')))
    nsr=g.input('nsr','Source50pmolNSRclass',K.OLIGO,state=nsr_scene(),chemical_state={'source_ordering':NSR,'underlined_positions_one_based':(12,13,14,15,16,17),'meaning':'priming sequence class, not molecule UMI'})
    dt=g.input('oligodt','Source5fmolP5oligodT',K.OLIGO,state=oligo(OLIGODT),chemical_state={'source_ordering':SOURCE_ORDERINGS['oligodT']})
    b.op('rt_prepare','SourceRNA72C3min/ice withNSR/oligodT/dNTP',S.THERMAL_TREATMENT,((nsr,R.STRUCTURAL),(dt,R.STRUCTURAL),b.stock('rt_nt','Source.8ul25mMdNTP pool, total/perbase convention unassigned'),b.cond('rt_denature','Source72C3min/ice')))
    template=g.product('rt_template','Same preparedRNA/template component',K.MOLECULAR);n=g.product('retained_nsr','Same-tube retainedNSR component',K.OLIGO,state=nsr_scene());d=g.product('retained_dt','Same-tube retainedoligodTcomponent',K.OLIGO,state=oligo(OLIGODT))
    g.react('rt_components','Copresent preparedRNA and primers, not physical fractions or fresh primeradditions',S.GENERIC,((b.x,R.TEMPLATE),),(template,n,d),changes_molecule=False);b.x=template
    b.op('rt','SourceMaxima100U/5ulRTmix25C5/37C15/40C90/70C15min',S.REVERSE_TRANSCRIPTION,((n,R.STRUCTURAL),(d,R.STRUCTURAL),b.cat('maxima','Source100UMaximaHMinus'),b.cat('rt_inhibitor','Source10URiboLock'),b.cond('rt_buffer','Source2ul5xRTbuffer/2ul5Mbetaine/.25ulwater')))
    b.op('tagment','Source10ultagmentmix55C30min',inputs=((ta,R.STRUCTURAL),b.stock('tag_atp','ActualATPcomponent of4ul5xTDsource4.25mM descriptor'),b.cond('td','SameTD50mMTris/25mMMg/50%DMF/PEG4.5ul40%/.5ulwater')),product='Source tagmented RNA/cDNA mixture plus intact/off-target/free stocks; strand transfers and products unassembled')
    b.op('sds_stop','Source1ul.4%SDS55C7min',S.QUENCHING,(b.cond('sds','Source1ul.4%SDS'),))
    b.op('gapfill','SourceBst3.0/Q5gapfill72C15min80C5min',inputs=(b.cat('bst','Source1ul8UulBst3.0'),b.stock('gap_q5','Actual1xQ5mastermix enzyme+nucleotides, whole occupiedpolymerization geometry unassembled')))
    b.pcr('preamp','Sourcefresh3ul2xQ5/4.5pmoleachT7customB+customP5;2low40+2at60+n11singlecell, incrementalextension',oligo(T7_CUSTOM),oligo(CUSTOM_P5))
    # Source guides are a prepared pool; guide-target-only table does not specify complete sgRNA.
    guides=g.input('guide_pool','Actualequal-pooledT7IVTsgRNAstocks forTableS3 27rRNA targets; whole scaffolds/provenance intermediates unassembled',K.OLIGO)
    rn=g.product('dash_rnp','Source Cas9/sgRNAprecomplex, target occupancy and free components unassigned',K.MOLECULAR)
    g.react('dash_assemble','SourceCas9:sgRNA1:2buffer3.137C15min',S.BINDING,((guides,R.TEMPLATE),b.stock('cas9','ActualsourceCas9proteinbindingstock'),b.cond('dash_buffer','Source1xbuffer3.1')),(rn,))
    prior=b.x;b.x=rn;rnp=b.aliquot('dash_input','Source5ulRNPmix selected','Actual5ulRNPmix','UntransferredRNPmix stock');b.x=prior
    b.aliquot('dash_dna','Source50ulPCRproductsselected, otherpreparedPCRvolumeunassigned','Actual50ulPCRproducts','Untransferredpreampstock, amount/fateunassigned')
    b.op('dash','Source5ulRNP+50ulPCR37C2h',inputs=((rnp,R.STRUCTURAL),),product='Source DASH-cleaved and surviving DNA classes, not guaranteed zero rRNA')
    b.op('dash_pk','Source1ulthermolabilePK37C15min/55C10min',S.DEPROTEINIZATION,(b.cat('dash_pk_stock','Source1ulthermolabileProteinaseK'),))
    b.clean('dash_cleanup','Source1xAMPure/elute5ul; RETAINbeads in tube forIVT')
    b.op('ivt','Source5ulIVTmix onto5ulDNA+retainedbeads37C10–16h',S.TRANSCRIPTION,(b.cat('t7','Source1ul50UulT7RNApolymerase'),b.cat('yipp','Source.5ul.1UulYIPP'),b.cat('ivt_inhib','Source.5ulRNaseinhibitor'),b.stock('ivt_nt','Source2ul10mMNTPpool'),b.cond('ivt_buf','Source1ul10xIVTbuffer')),product='SourceIVTRNA plus retainedDNA/beads/freeNTPs/enzymes, products/yieldunassembled')
    b.clean('rna_cleanup','Sourceadd45ul10mMTris to55ultotal then55ulHXP20%PEG/2.5MNaCl; eluteRNA11ul')
    b.op('rt2_prepare','Source11ulRNA+1ul10mMdNTP+1ul10uMP5primer65C5min/ice2min',S.THERMAL_TREATMENT,(b.stock('rt2_nt','Source1ul10mMdNTP'),b.stock('rt2_primer','Source1ul10uMwholeP5primer',K.OLIGO,state=oligo(FINAL_P5)),b.cond('rt2_anneal','Source65C5min/ice2min')))
    tmp=g.product('rt2_template','Same-tube preparedRNAtemplate',K.MOLECULAR);pr=g.product('rt2_retained_primer','Same-tube retainedP5primer',K.OLIGO,state=oligo(FINAL_P5));g.react('rt2_components','CopresentRNA and retained externalprimer, not fractions',S.GENERIC,((b.x,R.TEMPLATE),),(tmp,pr),changes_molecule=False);b.x=tmp
    b.op('rt2','SourceSSIII200U25C5min50C1h70C15min',S.REVERSE_TRANSCRIPTION,((pr,R.STRUCTURAL),b.cat('ssiii','Source1ul200UulSSIII'),b.cat('rt2_inhib','Source1ulRNaseinhibitor'),b.cond('rt2_buffer','Source4ul5xbuffer+1ulwater')))
    b.aliquot('final_template','Selectsource10ulRTcDNA from20ultotal','Actual10ulRTcDNA','Untransferred10ulRTcDNA')
    b.pcr('final_pcr','Source50ulPCR:25ul2xQ5,.5uleach10uMP5/P7,14ulwater;3–5cycles',oligo(FINAL_P5),index_scene())
    b.clean('final_cleanup','Source.8xAMPure/elute20ulTE')
    b.aliquot('qubit','SourceQubitanalyticalallocation amountunassigned',qc=True)
    b.aliquot('loading','SourceNovaSeqloadingallocation amountunassigned')
    b.finish('library','RAG3.0selectedlibrary','SourceNovaSeq6000pairedend; actualsampleindex/completephysicalreadmask/primerlanding unassembled. No claim of exact RNAstrand-transfer endpoints or complete rRNAremoval.')
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();c('MEwholepair exact',revcomp(ME),ME_RC);loaded_duplex()
    c('RetainedRNAprimer not freshlyadded',g.steps['rt2'].inputs[1].material,'rt2_retained_primer')
    c('SinglecellDASHuses actual50ulallocation',g.steps['dash'].inputs[0].material,'dash_dna_selected')
    c('IVTretainsselectedDASHcleanup',g.steps['ivt'].inputs[0].material,'dash_cleanup_kept')
    c('HalfRTstockremains',g.terminals['final_template_remainder'].disposition,Disposition.UNRESOLVED)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'gpbjnl_qzae072_supp_qzae072_supplementary_data_Supplementary_Table_S2.xlsx.txt'
    if have(c,Source(p,'https://doi.org/10.1093/gpbjnl/qzae072')):
        text=''.join(p.read_text().split())
        for name,seq in SOURCE_ORDERINGS.items():c('Wholeoriginalordering '+name,seq in text,True)
    c('Partialnotpublishedready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
