"""Source-local bounded partial model; no inherited kit recipe or finished library claim."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from dataclasses import dataclass,replace
from chemdraw import Construct,Scene,Segment,feature,revcomp,panel,workflow_from_protocol_graph,workflow_panel
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem
from pcr import PCRGeometry
from page import head,info
PAGE_FILE='chem.html'

def unknown(label):return Construct([Segment(label,'[unassembled]',placeholder=True)])
def cleanup(g,key,x,label,procedure,state=None):
    sel=g.input(key+'_selection',procedure,K.SELECTION)
    kept=g.product(key+'_kept',label,K.MOLECULAR,state=state)
    waste=g.product(key+'_waste','Physically separated '+key+' rejects',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(sel,R.SELECTION)),(kept,waste))
    g.discard(waste,'Source-described physical separation; purity and quantitative recovery are not assumed')
    return kept

def aliquot(g,key,x,label):
    chosen=g.product(key+'_chosen',label,K.MOLECULAR)
    rest=g.product(key+'_remainder','Untransferred '+key+' material',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(chosen,rest))
    g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Source aliquot remainder fate is not specified in the inspected section; no disposal inferred'))
    return chosen

def limit(g,x):
    g.mark_terminal(Terminal(x.id,Disposition.UNRESOLVED,LIMIT))
    g.declare_reviewed_scope(REVIEWED,SOURCE)
    g.add_source_item(SourceItem('unassembled_source','Remaining source library/input chemistry and full strand/read assembly',LIMIT))
    return g

def render_page():
    g=protocol_graph()
    return '\n'.join([head(TITLE),'<main class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="'+NOTE+'">Research notes</a></p>',
        info('Private incomplete source model. '+LIMIT),detail(),
        ''.join(workflow_panel(workflow_from_protocol_graph(g,x)) for x in ENDPOINTS),'</main>'])

def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Unassembled source and endpoint remain unpublished',g.publishable,False)
    c('Private preview has no completion marker','data-protocol-complete=' in render_page(),False)
    p=sources_root()/Path(__file__).resolve().parents[1].name/SOURCE_FILE
    if p.exists():
        source=''.join(p.read_text().split()).upper()
        for seq in TRANSCRIBED:c('Original transcription '+seq,seq.upper() in source,True)
    else:c.skip('Source transcription','External original unavailable')
    source_checks(c,g)
    c.report()
import illumina as il
import nextera as nx
TITLE='Perturb-seq — original expressed barcode proxy'
NOTE='01_perturb_seq.html'
SOURCE='Dixit etal.2016 associated accepted manuscript PMC5181115; vector construction256, single-cell preparation318–379, processing387'
LIMIT='Starts with original10x WTA material. Full historical capture/RT substrate, GBC width/value, full final-strand assembly and paired-read primer landing remain unassembled. The proxy transcript is distinct from sgRNA; sgRNA assignment uses a separately measured dictionary.'
SOURCE_FILE='j.cell.2016.11.038_PMC5181115.xml.txt'
BFP_SITE='TAGCAAACTGGGGCACAAGC'
INDEX='TCGCCTTA'
# Source names this dial-out handle R1; sharing a sequence is not recipe inheritance.
DIAL_OUT=il.P5+INDEX+nx.ADAPTOR_S5+BFP_SITE
TRANSCRIBED=(DIAL_OUT,il.P7)
ENDPOINTS=('gel_kept',)

def dial_primer():return Construct([Segment('P5',il.P5),Segment('sequencing sample index',INDEX,'idx',feature=feature('perturb.sample','sample_index','fixed')),
    Segment('source R1 dial-out handle',nx.ADAPTOR_S5),Segment('BFP-targeting3prime arm',BFP_SITE)])
def identifiers():
    return [Segment('guide barcode proxy','[GBC]',placeholder=True,feature=feature('perturb.guide','guide_barcode','unknown')),
        Segment('cell identity','[CBC]',placeholder=True,feature=feature('perturb.cell','cell_barcode','unknown')),
        Segment('molecular identifier','[UMI]',placeholder=True,feature=feature('perturb.umi','umi','unknown'))]

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    wta=g.input('wta','Source10x whole-transcriptome amplified cDNA, including expressed polyadenylated BFP/GBC transcript products',K.MOLECULAR,primary=True)
    source=aliquot(g,'wta_aliquot',wta,'Approximately5ng WTA used for GBC enrichment; remaining WTA retained separately')
    f=g.input('dial_out','Source indexed dial-out primer, selected printed TCGCCTTA member',K.OLIGO,state=dial_primer())
    r=g.input('p7','Printed P7 Illumina reverse primer',K.OLIGO,state=Construct([Segment('P7 primer',il.P7)]))
    enzyme=g.input('polymerase','Q5 2×mastermix;25cycles with65°C5s annealing',K.CATALYST)
    product=g.product('gbc_product','GBC-enriched products preserving source CBC/UMI linkage; residual WTA/off-products remain in reaction',K.MOLECULAR)
    g.react('dial_pcr','Dial-out PCR enriches barcode-proxy transcript after WTA',S.PCR,((source,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(enzyme,R.CATALYST)),(product,),
        geometry=PCRGeometry.unavailable(unknown('original barcode-proxy cDNA'),reason='Printed primers known; original GBC transcript and historical10x opposite end are not assembled. This is not direct sgRNA capture.'))
    purified=cleanup(g,'gel',product,'Gel-recovered original Perturb-seq GBC library material','2%agarose gel;appropriate band excision')
    return limit(g,purified)

def detail():
    sc=Scene()
    for s in identifiers():sc.strand(s.name,[s],label=s.name+' (role only; not asserted contiguous or assigned width)')
    return panel(sc.rows(),caption='Three distinct identities: proxy guide barcode, cell barcode and UMI')
def source_checks(c,g):
    c('Proxy guide identity is distinct from cell and molecular identity',tuple(x.feature.role for x in identifiers()),('guide_barcode','cell_barcode','umi'))
    c('Dial-out reaction starts from transferred WTA',g.steps['dial_pcr'].inputs[0].material,'wta_aliquot_chosen')
    c('Guide dictionary lookup does not appear as a molecular transformation',tuple(g.steps),('wta_aliquot','dial_pcr','gel'))

REVIEWED=('wta', 'wta_aliquot_chosen', 'wta_aliquot_remainder', 'dial_out', 'p7', 'polymerase', 'gbc_product', 'gel_selection', 'gel_kept', 'gel_waste', 'wta_aliquot', 'dial_pcr', 'gel')
