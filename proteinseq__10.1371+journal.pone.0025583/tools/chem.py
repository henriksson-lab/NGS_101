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
        g.react(key,label,kind,inputs,(y,));return y
    def clean(key,x,label):
        z=g.input(key+'_selection',label,K.SELECTION)
        y=g.product(key+'_kept','Physically retained source fraction: '+label,K.MOLECULAR)
        w=g.product(key+'_removed','Physically removed fractions: '+label,K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(z,R.SELECTION)),(y,w));g.discard(w,'Source physical workup only; no perfect specificity or purity assertion.');return y
    return op,clean

def end(g,x,label):g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=label))
def incomplete(g):
    g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE)
    g.add_source_item(SourceItem('unassembled','Unreviewed linked originals and full molecular geometry',LIMIT));return g


from protocol_graph import Terminal,Disposition

def save(g,key,x,label):
    take=g.product(key+'_aliquot',label,K.MOLECULAR)
    rest=g.product(key+'_remainder','Remaining preparation; later use unassigned',K.MOLECULAR)
    g.react(key,'Physical aliquot: '+label,S.GENERIC,((x,R.TEMPLATE),),(take,rest))
    g.terminals[rest.id]=Terminal(rest.id,Disposition.STORED,'Remainder retained explicitly; later fate unassigned.')
    return take


def measure(g,key,x,label):
    aliquot=g.product(key+'_aliquot','Separate source analytical aliquot;exact amount unassigned unless stated',K.MOLECULAR)
    keep=g.product(key+'_production','Production preparation retained separately from analytical aliquot',K.MOLECULAR)
    g.react(key+'_split','Separate analytical sample: '+label,S.SPLIT,((x,R.TEMPLATE),),(aliquot,keep))
    measured=g.product(key+'_qc','Analytically assayed aliquot: '+label,K.QC)
    g.react(key,label,S.GENERIC,((aliquot,R.TEMPLATE),),(measured,))
    g.terminals[measured.id]=Terminal(measured.id,Disposition.QC,'Analytical aliquot not the production material pooled for sequencing.')
    return keep

def render_page():
    g=protocol_graph()
    return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_chem.html">Research notes</a></p>'+workflow_panel(workflow_from_protocol_graph(g,ENDPOINT))+'</main>'

from illumina import P5,P7
TITLE='ProteinSeq — published solid-phase PLA and NGS continuation'
SOURCE='Published PMC3183061 Methods194–274; published TableS8(s012 DOCX)'
LIMIT='Private incomplete graph. Source-known NGS operations represented, but ligating probe bodies/splint, exact binding and complete reporter/copy geometry remain unassembled. No ideal affinity selection or ligation yield is predicted. The qPCR dUTP/15-cycle branch is excluded.'
ENDPOINT='pool_product'
INDEX1=P5+'TCTTTCCCTACACGACGCTCTTCCGATC'+'CGTGAT'+'GCGAAACCTGGTCCGGTATC'
REVERSE=P7+'CGATTCGAGAACGTGACTGC'
def index_primer():
    return Construct([Segment('source P5 prefix',P5),Segment('source intervening sequence','TCTTTCCCTACACGACGCTCTTCCGATC'),Segment('source index1','CGTGAT',feature=feature('proteinseq.index1','sample_index','fixed',note='Printed TableS8 member; not an assigned experimental sample.')),Segment('source primer body','GCGAAACCTGGTCCGGTATC')])
def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE);op,clean=operations(g)
    x=g.input('sample','Source diluted protein/plasma preparation',K.BIOLOGICAL,primary=True)
    x=op('capture',x,'Incubate 5µL pooled capture-antibody beads with45µL sample;1.5h RT rotation','Source biotinylated capture-antibody streptavidin beads/PLA buffer;manufacturing predecessor unreviewed',reagent_kind=K.MOLECULAR,reagent_role=R.STRUCTURAL)
    x=clean('capture_wash',x,'Two magnetic PBS/0.05%Tween washes;keep beads and associated classes')
    x=op('probe',x,'50µL paired antibody–DNA probe mix,1.5h RT rotation','Source antibody conjugates carrying5′ or3′ thiol-associated DNA;sequences/conjugation geometry unreviewed',reagent_kind=K.MOLECULAR,reagent_role=R.STRUCTURAL)
    x=clean('probe_wash',x,'Two magnetic washes;keep bead fraction,do not assert perfect binding specificity')
    splint=g.input('splint','100nM source sliding splint;whole sequence unreviewed in deferred TableS7',K.OLIGO)
    y=g.product('ligation_product','Source ligation preparation;associated DNA/protein/beads,unligated probes/splint and reaction remnants retained',K.MOLECULAR)
    g.react('ligation','50µL AmpLigase/NAD/sliding-splint reaction50°C10min;source intended cognate ligation',S.GENERIC,((x,R.TEMPLATE),(splint,R.STRUCTURAL)),(y,));x=y
    x=clean('ligation_wash',x,'One source bead wash after ligation;associated DNA and proteins remain')
    uf=g.input('universal_forward','100nM source universal forward primer;whole body unreviewed',K.OLIGO);ur=g.input('universal_reverse','100nM source universal reverse primer;whole body unreviewed',K.OLIGO)
    y=g.product('universal_pcr_product','New DNA PCR copies plus original bead/protein/DNA classes and reaction components',K.MOLECULAR)
    g.react('universal_pcr','Published NGS universalPCR35cycles;PlatinumTaq/dNTPmix/UNG,dNTPcompositionunassembled,noimportedqPCRdUTPrecipe',S.GENERIC,((x,R.TEMPLATE),(uf,R.STRUCTURAL),(ur,R.STRUCTURAL)),(y,));x=y
    x=save(g,'preamp_split',x,'24µL preamplification preparation for indexing')
    f=g.input('index_forward','Published TableS8 Index1 example',K.OLIGO,state=index_primer());r=g.input('index_reverse','Published TableS8 reverse primer',K.OLIGO,state=Construct([Segment('source P7 prefix',P7),Segment('source reverse body','CGATTCGAGAACGTGACTGC')]))
    y=g.product('indexed','Index PCR preparation with newly copied DNA;original input components not declared removed',K.MOLECULAR)
    g.react('indexing','Add26µL index mix;published stated95/62/72°C program;extra repeat count not invented',S.GENERIC,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL)),(y,))
    x=save(g,'indexed_split',y,'10µL indexed preparation for enrichment')
    ef=g.input('enrich_forward','Published TableS8 enrichment forward primer',K.OLIGO,state=Construct([Segment('source enrichment forward',P7)]));er=g.input('enrich_reverse','Published TableS8 enrichment reverse primer',K.OLIGO,state=Construct([Segment('source enrichment reverse',P5[:25])]))
    y=g.product('enrich_product','Source enrichment PCR new DNA and original template/reagent classes retained',K.MOLECULAR)
    g.react('enrich','Add40µL enrichmentmix;3PCRcycles95/65/72°C',S.GENERIC,((x,R.TEMPLATE),(ef,R.STRUCTURAL),(er,R.STRUCTURAL)),(y,));x=y
    x=clean('gfx',x,'GFX column purification;retained DNA eluate')
    x=clean('egel',x,'2%E-gel selection of199bp product;retain selected DNA fraction')
    x=measure(g,'quantify',x,'Source qPCR quantification of separate analytical aliquot;amount unreviewed')
    x=op('pool',x,'Pool source groups of12indexes(7groups);actual sample/index assignments unassembled')
    end(g,x,'Published GAIIx single-end95bases;6bp sample index. Exact reporter/read-primer junction and site geometry unassembled;no arbitrary canonical run layout.')
    return incomplete(g)
def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Index semantic role distinct from analyte identity',index_primer().segments[2].feature.role,'sample_index')
    c('Indexing uses aliquot,not whole preamp',g.steps['indexing'].inputs[0].material,'preamp_split_aliquot')
    c('Saved preamp remainder retained',g.terminals['preamp_split_remainder'].disposition,Disposition.STORED)
    c('Enrichment follows indexed aliquot',g.steps['enrich'].inputs[0].material,'indexed_split_aliquot')
    c('Probe wash physically precedes ligation',g.steps['ligation'].inputs[0].material,'probe_wash_kept')
    c('Pooled library is retained production preparation,not qPCR tube',g.steps['pool'].inputs[0].material,'quantify_production')
    c('Private incomplete source scope',g.publishable,False)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'journal.pone.0025583_supp_pone.0025583.s012.docx.txt'
    if p.exists():c('Whole printed index and reverse primers',all(v in p.read_text() for v in (INDEX1,REVERSE)),True)
    else:c.skip('Original TableS8','External original absent')
    c.report()
