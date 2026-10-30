"""Source-version-limited original GESTALT genomic reporter recovery."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature,workflow_from_protocol_graph,workflow_panel,panel
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan
from pcr import PCRGeometry
from page import head,info
PAGE_FILE='chem.html'
TITLE='GESTALT — original genomic UMI-tagged reporter recovery'
SOURCE='Source-version conditional:accepted author manuscript PMC4967023 Results65–87,especially single-extensionUMI paragraph79;finalpublishedmain and detailedSI not retrieved/read'
LIMIT='Private source-limited accepted-manuscript route. Detailed source Methods/primer workbook are gated;enzyme identities,primer ends,cleanup order,full molecular products and sequencing-primer/run geometry remain unassembled. The published article mirror returned HTML rather than PDF. Excluded preprint material does not fill these gaps.'
LINEAGE=feature('gestalt.edited_array','lineage_barcode','unknown',note='Inherited edited recording allele;not a unique present-day cell identifier.')
UMI=feature('gestalt.recovery_umi','umi','random')

@dataclass(frozen=True)
class RecorderContext:
    variant:str='v1 cell-culture array'
    target_sites:int=10
    wildtype_nt:int=257
    def rows(self):
        sc=Scene();sc.strand('array',[Segment('edited source recorder','[source/sample-specific edited array]',placeholder=True,feature=LINEAGE)])
        sc.note('array','Recorder design context:sourcev1 wild-type257bp/10targets downstreamEGFP;edited lengths/alleles not predicted. Array creation/editing occurs before this genomic recovery graph.')
        return sc.rows()

def identifier_inventory():
    sc=Scene();sc.strand('lineage',[Segment('inherited edited array','[edited allele]',placeholder=True,feature=LINEAGE)])
    sc.strand('umi',[Segment('source recovery UMI','N'*10,placeholder=True,feature=UMI)])
    sc.note('lineage','Separate identifier meanings;not an assembled primer or a guarantee of one unique lineage barcode per cell.');return sc

def umi_primer_inventory():
    sc=Scene();sc.strand('umi',[Segment('known primer UMI component','N'*10,placeholder=True,feature=UMI)],label='Component inventory only;full primer/end placement unassembled')
    sc.note('umi','This does not place the UMI relative to unreviewed binding/adaptor regions.');return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    dna=g.input('genomic_dna','Source gDNA from editedv1 cell-culture populations;extraction chemistry unassembled',K.MOLECULAR,primary=True,state=RecorderContext())
    primer=g.input('umi_primer','Source10ntUMI-tagging extension primer;full sequence/end geometry unavailable in this pass',K.OLIGO,state=umi_primer_inventory())
    polymerase=g.input('extension_polymerase','Source polymerase extension;exact enzyme/program not inspected',K.CATALYST)
    tagged=g.product('tagged','Genomic reporter extension reaction with UMI-tagged DNA,original DNA and residual reagents',K.MOLECULAR)
    g.react('umi_extension','Single round of source polymerase extension introduces recoveryUMI;not repeated independent PCR barcoding',S.EXTENSION,((dna,R.TEMPLATE),(primer,R.STRUCTURAL),(polymerase,R.CATALYST)),(tagged,))
    f=g.input('forward','Source reporter-amplification forward primer;detailedSI not read',K.OLIGO)
    r=g.input('reverse','Source reporter-amplification reverse primer;detailedSI not read',K.OLIGO)
    pcrmix=g.input('pcr_mix','Source reporter PCR reagents/cycle count unassembled',K.CATALYST)
    amplified=g.product('amplified','Targeted GESTALT genomic reporter PCR preparation;full duplex unassembled',K.MOLECULAR)
    g.react('reporter_pcr','Source targeted amplification of UMI-tagged reporter',S.PCR,((tagged,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(pcrmix,R.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(identifier_inventory(),reason='Exact source oligos/whole genomic strands and intervening workup order await detailed Methods original. Identifier components do not establish a complete substrate.'))
    g.sequence(amplified,SequencingPlan(acquisition='unspecified',limitation=LIMIT+' Source reports sequencing edited barcodes;no actual read-length/index/primer assignment inferred from later scGESTALT or CARLIN.'))
    g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE)
    g.add_source_item(SourceItem('detailed_methods','Actual defining Methods and primer original',LIMIT));return g

def render_page():
    return '\n'.join([head(TITLE),'<main class="wrap"><h1>'+TITLE+'</h1>',info(LIMIT),'<p><a href="01_gestalt.html">Research notes</a></p>',workflow_panel(workflow_from_protocol_graph(protocol_graph(),'amplified')),panel(RecorderContext().rows(),caption='Recorder context,not a library-construction reaction or an indel prediction'),'</main>'])

def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph();c('Incomplete source model stays private',g.publishable,False)
    c('No finished marker','data-protocol-complete=' in render_page(),False)
    c('Single extension precedes reporterPCR',g.steps['reporter_pcr'].inputs[0].material,'tagged')
    c('Lineage identity separate from recovery UMI',LINEAGE.role!=UMI.role,True)
    c('No later scGESTALT or CARLIN recipe operations',set(g.steps),{'umi_extension','reporter_pcr'})
    p=sources_root()/Path(__file__).resolve().parents[1].name/'science.aaf7907_PMC4967023.xml.txt'
    if p.exists():
        source=' '.join(p.read_text().split())
        c('Selectedoriginal describes10bpUMI singleextension','UMIs; 10 bp' in source and 'single round of polymerase extension' in source,True)
    else:c.skip('Source extension/UMI transcription','External original unavailable')
    c.report()
