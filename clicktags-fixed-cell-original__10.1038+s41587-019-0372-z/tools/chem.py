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
import nextera as nx
TITLE='ClickTags — published DNA/probe activity component'
NOTE='01_clicktags.html'
SOURCE='Published SI DOI10.1038/s41587-019-0372-z; MOESM1 SupplementaryFig1c, renderedPDFp2; MOESM4 modified-oligo inventory'
LIMIT='Only published MTZ-DNA/TCO-Cy5 activity assay is modeled. Main published Methods remain unavailable. Actual assay DNA sequence/lane identity, DNA activation workup, fixed-cell tagging, capture, library and run chemistry remain unassembled; BC21 inventory is not asserted to be an assayed lane. Its variable8base sample-tag interpretation is derived from the source oligo family, not a source-assigned assay-lane sequence.'
SOURCE_FILE='s41587-019-0372-z_41587_2019_372_MOESM4_ESM.xlsx.txt'
HANDLE=nx.ADAPTOR_S5[:23]
TAG='AGCAGTTA'
BC21=HANDLE+TAG+'CAGB'+'A'*25
TRANSCRIBED=('/5AmMC6/'+BC21,)
ENDPOINTS=('assay_product',)
@dataclass(frozen=True)
class ProbePath:
    methyltetrazine:bool=True
    cy5_appended:bool=False
    dna_backbone_retained:bool=True
    def __post_init__(self):
        if not self.dna_backbone_retained:raise ValueError('Selected IEDDA coupling does not delete barcode DNA')
        if self.cy5_appended and self.methyltetrazine:raise ValueError('Selected coupled handle is not an unreacted methyltetrazine')
    def couple(self):
        if not self.methyltetrazine or self.cy5_appended:raise ValueError('A fresh complementary reactive handle is required')
        return replace(self,methyltetrazine=False,cy5_appended=True)

def inventory():
    sc=Scene();sc.strand('BC21',[Segment('source handle',HANDLE),Segment('sample tag',TAG,'idx',inferred=True,feature=feature('clicktags.sample','sample_index','fixed')),
        Segment('source CAGB window','CAGB',placeholder=True,note='B source degeneracy retained; no UMI role inferred'),Segment('polyA','A'*25)],mod5='amine-C6',label='BC21 source inventory; not assigned to an assay lane')
    return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    dna=g.input('mtz_dna','Published Fig1c methyltetrazine-activated ClickTag DNA; exact assayed strand unassigned',K.MOLECULAR,primary=True,chemical_state=ProbePath())
    tco=g.input('tco_cy5','TCO-Cy5 as printed in published SI Fig1c',K.MOLECULAR)
    conditions=g.input('conditions','Published Fig1c reaction conditions unassigned; no time/buffer/catalyst invented',K.CONDITION,inferred=True)
    product=g.product('assay_product','Selected Cy5-labeled DNA class plus unreacted dye/DNA fractions; no quantitative yield',K.MOLECULAR,chemical_state={'selected':ProbePath().couple(),'other_species':'unreacted dye/DNA remain until actual PAGE separation'})
    g.react('iedda','Published probe IEDDA coupling; no RNA/base conversion or invented copper catalyst',S.CHEMICAL_CONVERSION,((dna,R.TEMPLATE),(tco,R.STRUCTURAL),(conditions,R.CONDITION)),(product,))
    return limit(g,product)

def detail():return info('Published Fig1a/b yeast assay uses NHS-TCO then MTZ-Cy5; Fig1c instead couples MTZ-DNA to TCO-Cy5. Those assays are not merged into an invented full cell-DNA labeling recipe.')+panel(inventory().rows(),caption='INFERRED — BC21 variable8base sample-tag interpretation; source oligo/modification is verbatim, assay-lane mapping unassigned')
def source_checks(c,g):
    p=ProbePath().couple()
    c('IEDDA labeling retains DNA while consuming selected handle',(p.cy5_appended,p.methyltetrazine,p.dna_backbone_retained),(True,False,True))
    c('BC21 retains constrained B rather than a UMI declaration',next(s.top for s in inventory().strands['BC21'].segs if s.name=='source CAGB window'),'CAGB')
    c('Assay sequence is not silently filled with BC21',g.materials['mtz_dna'].state,None)

REVIEWED=('mtz_dna', 'tco_cy5', 'conditions', 'assay_product', 'iedda')
