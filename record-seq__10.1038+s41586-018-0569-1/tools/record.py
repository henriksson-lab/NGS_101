"""Published Record-seq/SENECA: array classes and source-local library operations."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
from chemdraw import Scene,Segment,Construct,feature,workflow_from_protocol_graph,workflow_panel
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SequencingPlan,SourceItem
from pcr import PCRGeometry
from illumina import P5,P7,TRUSEQ_READ1,TRUSEQ_READ2,INDEX1_PRIMER
from page import head,info

PAGE_FILE='record.html'
SOURCE='Schmidt et al.2018,10.1038/s41586-018-0569-1; published university-mirrored article Methods1385–1458/ExtendedDataFig2, published SI Notes9–12/PrimerDesignNotes1–2/Tables4–7.'
LIMIT='Private incomplete model. RNA-derived acquisition is supported but DNA-derived spacers are not excluded. Exact RT/integration mechanism, selected array/template ends and full primer-product/read geometry remain unassembled. This population recorder is not a single-cell lineage barcode, a UMI or an assured chronological transcript history.'
DISCREPANCIES='Published main specifies a10µl SENECA ligation but10.38µl input to one firstPCR. It also cites PrimerDesignNote3, while the published supplement lists only Notes1–2. Note1 prints the same P5-arm string under its reverse/P7 schematic. The model retains those ambiguities; it uses actual FS_0911/0963/0964 table entries and Note2 for secondPCR, without substituting a later author-protocol revision.'
I5=feature('record.sample.p5','sample_index','whitelist',whitelist='Published SI Table4 BC1')
I7=feature('record.sample.p7','sample_index','whitelist',whitelist='Published SI Table4 BC2')
CORE=INDEX1_PRIMER
FS0964=Construct([Segment('source4ntoverhang','AAAG'),Segment('pairedcore',CORE)],name='FS_0964')
FS0963=Construct([Segment('universal',TRUSEQ_READ2[:-1])],name='FS_0963')
FS0911=Construct([Segment('source reverse binding arm',TRUSEQ_READ2[:18])],name='FS_0911')
SECOND_F=Construct([Segment('P5',P5),Segment('BC1','N'*8,placeholder=True,feature=I5),Segment('read1arm',TRUSEQ_READ1)],name='Source secondPCR forward index family')
SECOND_R=Construct([Segment('P7',P7),Segment('BC2','N'*8,placeholder=True,feature=I7),Segment('source short read2arm',TRUSEQ_READ2[:-1])],name='Source secondPCR reverse index family')
DR_ARM='CCTAAAAGGAATTGAAAC'


def unknown(label):
    sc=Scene();sc.strand('region',[role(label)]);return sc


def adapter_scene():
    sc=Scene();sc.strand('fs0964',FS0964.segments,label='FS_0964, physical5prime to3prime')
    sc.anneal('fs0963',FS0963.segments,to='fs0964',pair=('universal','pairedcore'),label='FS_0963, antiparallel paired core')
    return sc


@dataclass(frozen=True)
class ArrayClass:
    """Conditional topology of the source figure, not a digest prediction.

    The depicted array cut context and intact copied repeats are assumed for this
    reference comparison. Extra acquired restriction sites, damaged repeats,
    exact templates and empirical amplification are outside the relation.
    """
    spacers:int
    adapter:bool=False
    def __post_init__(self):
        if type(self.spacers)is not int or self.spacers<0:raise ValueError('Explicit nonnegative acquired spacer count required')
    @property
    def full_repeat_sites_after_faqi(self):return self.spacers
    @property
    def eligible_primer_sides(self):return (self.spacers>0,self.adapter)
    def ligate_source_adapter(self):return ArrayClass(self.spacers,True)


@dataclass(frozen=True)
class ArrayPopulation:
    stage:str='recorded'
    def __post_init__(self):
        if self.stage not in {'recorded','ligated'}:raise ValueError('Explicit recorded or ligated array stage required')
    def rows(self):
        sc=Scene()
        for name,count in (('parental',0),('expanded',1)):
            state=ArrayClass(count,self.stage!='recorded')
            parts=[]
            if count:parts.extend([role('intact full direct repeat; bases unassembled'),role('newly acquired spacer; RNA or DNA origin unresolved per molecule')])
            parts.append(role('source last direct repeat' if self.stage=='recorded' else 'source FaqI-truncated last repeat'))
            if state.adapter:parts.append(role('ligated source P7-arm adapter; exact template junction unassembled'))
            sc.strand(name,parts,label=name+' regional DNA strand view; cohort amount unassigned')
            if self.stage!='recorded':sc.note(name,'Source-figure reference primer eligibility: '+('both primer sides present' if all(state.eligible_primer_sides) else 'only adapter/reverse side present')+'; not a yield or perfect selection prediction.')
        sc.strand('other',[role('higher acquisition multiplicities, original DNA/RNA/protein/oligo and reaction remnants; recovery unassigned')],label='Other source material not declared physically absent')
        return sc.rows()


def first_primer_inventory(stagger_length=2):
    if type(stagger_length)is not int or not 2<=stagger_length<=8:raise ValueError('Source stagger family is2–8, not a UMI')
    sc=Scene()
    sc.strand('tail',[Segment('source5prime tail',TRUSEQ_READ1)],label='Components of ONE source forward primer; complete joined model unassembled')
    sc.strand('stagger',[Segment('sequence-diversity stagger','N'*stagger_length,placeholder=True)],label='Diversity, not molecular identity')
    sc.strand('binding',[Segment('source FsArray2 DR binding arm',DR_ARM)],label='Same primer binding component; no separate free reagent')
    return sc


def protocol_graph():
    g=ProtocolGraph('Record-seq — source SENECA population readout',authority=SOURCE);scope=[]
    def inp(key,label,kind=K.CATALYST,state=None,primary=False):return g.input(key,label,kind,state=state,primary=primary)
    def op(key,x,label,kind,reagents=(),state=None,geometry=None,reason=None):
        out=g.product(key+'_out',label+' preparation',K.MOLECULAR,state=state or unknown(label+'; exact molecular geometry unassembled'))
        g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(reagents),(out,),geometry=geometry,**({'reason':reason} if reason else {}));scope.append(key);return out
    def clean(key,x,label):
        stock=inp(key+'_reagents',label,K.SELECTION)
        kept=g.product(key+'_kept',label+' retained fraction',K.MOLECULAR)
        waste=g.product(key+'_removed','Actually nonretained source fractions; purity/recovery unquantified',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(stock,R.SELECTION)),(kept,waste));g.discard(waste,'Source physical fraction transfer,not computational rejection.');scope.append(key);return kept
    cell=inp('cells','Source BL21(DE3)Gold carrying pFS_0235 minimalFsArray2 recorder;construction outside selected graph',K.BIOLOGICAL,primary=True)
    expression=inp('expression','Source FsRT-Cas1/Cas2 expression;TB3ml/100µMIPTG/37°C300rpm until saturation12–14h',K.CONDITION)
    x=op('record',cell,'Source induction and intracellular spacer acquisition;RNA-compatible recording with possible DNA acquisition',S.MAGIC,((expression,R.CONDITION),),ArrayPopulation(),reason='Source RNA-only junction/virus evidence supports RNA acquisition but does not exclude DNA-derived spacers. Exact RT/integration order and per-molecule substrate/product not assembled. No assured temporal/single-cell history.')
    x=clean('miniprep',x,'Harvest2mlculture;standard plasmid miniprep and PicoGreen quantification')
    a=inp('fs0963','Source universal reverse adapter FS_0963',K.OLIGO,FS0963)
    b=inp('fs0964','Source FsArray2 forward adapter FS_0964',K.OLIGO,FS0964)
    anneal=inp('anneal_conditions','95°C5min then20°C at0.12°C/s;sourceNEBuffer2,then1:100TEdilution',K.CONDITION)
    adapter=g.product('adapter','Annealed source adapter with computed donor overhang',K.ADAPTER,state=adapter_scene())
    g.react('anneal_adapter','Anneal source adapter oligos',S.HYBRIDIZATION,((a,R.TEMPLATE),(b,R.STRUCTURAL),(anneal,R.CONDITION)),(adapter,));scope.append('anneal_adapter')
    faq=inp('faq','Source FastDigestFaqI')
    ligase=inp('t7_ligase','SourceT7DNA ligase')
    mix=inp('golden_mix','Source40fmolplasmid/ATP/DTT/SAM;99cycles37°C3min/20°C3min,55°C15min;volume inconsistency retained',K.CONDITION)
    x=op('seneca',x,'Source coupled FaqI digestion/adapter ligation;parental and expanded classes both remain',S.MAGIC,((adapter,R.STRUCTURAL),(faq,R.CATALYST),(ligase,R.CATALYST),(mix,R.CONDITION)),ArrayPopulation('ligated'),reason='Source donor duplex is computed;full plasmid cut sites,ends and covalent junctions not assembled. Preferential PCR eligibility is not perfect physical removal of parental templates.')
    f=inp('first_forward','Source FS_0968–0974 forward family;Table6 binding arm andNote1stagger design as one incomplete primer inventory',K.OLIGO,first_primer_inventory())
    r=inp('first_reverse','Source FS_0911',K.OLIGO,FS0911)
    pol=inp('first_polymerase','SourceNEBNextHighFidelity2xPCRmix;22cycles98/57/72')
    x=op('first_pcr',x,'Source22cycle firstPCR;parental one-primer and expanded two-primer classes not assigned exclusive yields',S.PCR,((f,R.STRUCTURAL),(r,R.STRUCTURAL),(pol,R.CATALYST)),geometry=PCRGeometry.unavailable(Scene(),reason='Source reaction,primer components and relative eligibility known;native complete primer/template geometry and PCR/background fates unassembled.'))
    x=op('first_pool',x,'Source firstPCR pooling;membership prior indexing not assigned beyond this representative preparation',S.MAGIC,reason='Main says pooled before cleanup but scope/membership is not assigned;no invented cross-biosample association or index identity.')
    x=clean('first_cleanup',x,'Source magnetic bead cleanup1:1.6(v/v),recover25µlTE')
    i5=inp('second_forward','Source universal secondPCR BC1 family',K.OLIGO,SECOND_F)
    i7=inp('second_reverse','Source universal secondPCR BC2 family',K.OLIGO,SECOND_R)
    pol2=inp('second_polymerase','Source8cycle98/65/72 secondPCR;5µlfirstPCR input/20µlreaction')
    x=op('index_pcr',x,'Source dual sample-index PCR;no UMI or cell-barcode assignment',S.PCR,((i5,R.STRUCTURAL),(i7,R.STRUCTURAL),(pol2,R.CATALYST)),geometry=PCRGeometry.unavailable(Scene(),reason='Published universal index-primer family known;full joined source amplicon/primer landing unassembled.'))
    x=op('pool_indexed',x,'Source indexed-sample pooling;other sample assignments not expanded',S.MAGIC,reason='Known physical pooling operation;this graph expands one source preparation,not a fabricated multi-sample experiment.')
    x=clean('desalt',x,'Source QIAquickPCRpurification after indexed pooling')
    x=clean('gel',x,'Source2%E-gelEX size selection/QIAquickgel extraction;do not borrow classical211–300bpwindow')
    library=g.product('library','Source indexed SENECA library;exact strands and run geometry unassembled',K.LIBRARY)
    g.react('handoff','Source qPCR quantification and library loading handoff',S.GENERIC,((x,R.TEMPLATE),),(library,));scope.append('handoff')
    g.sequence(library,SequencingPlan(acquisition='unspecified',limitation='SourceMiSeq/NextSeq500,150cyclekits,1.4–1.6pMloading and5–10%PhiX. Kit cycle count is not assigned as a read length or paired/single configuration. Exact index acquisition/primer landing unassembled.'))
    g.declare_reviewed_scope(scope,SOURCE)
    g.add_source_item(SourceItem('unassembled','Exact recording intermediates,plasmid/end/primer geometry,source volume/design-note ambiguities and acquisition','Published main and SI known;no unread later protocol revision substituted and no confirmed absence.'))
    return g


def render_page():
    return head('Record-seq — source-limited model')+'<main class="wrap"><h1>Record-seq</h1><p><a href="01_record.html">Research notes</a></p>'+info(LIMIT)+info(DISCREPANCIES)+workflow_panel(workflow_from_protocol_graph(protocol_graph(),'library'))+'</main>'


def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();parent=ArrayClass(0).ligate_source_adapter();expanded=ArrayClass(1).ligate_source_adapter()
    c('Adapter ligation does not invent an acquired spacer',parent.spacers,0)
    c('Source parental class lacks the intact DR primer side',parent.eligible_primer_sides,(False,True))
    c('Expanded source class retains both primer sides',expanded.eligible_primer_sides,(True,True))
    c('Multiple acquisitions retain multiple intact-repeat sites',ArrayClass(2,True).full_repeat_sites_after_faqi,2)
    try:first_primer_inventory(1);bad=False
    except ValueError:bad=True
    c('Stagger construction rejects outside-source lengths',bad,True)
    g=protocol_graph();c('FirstPCR uses both parental/expanded mixture',g.steps['first_pcr'].inputs[0].material,'seneca_out')
    c('No added library reverse-transcription step',any(s.kind==S.REVERSE_TRANSCRIPTION for s in g.steps.values()),False)
    source=sources_root()/HERE.parent.name/'s41586-018-0569-1_41586_2018_569_MOESM1_ESM.pdf.txt'
    if source.exists():
        text=source.read_text()
        for construct in (FS0963,FS0964,FS0911):c('PublishedTable7 exact '+construct.name,construct.top() in text.replace('\n',''),True)
        c('PublishedTable6 DR binding arm',DR_ARM in text,True)
    else:c.skip('Published oligo originals absent','External sources unavailable')
    c.report()
