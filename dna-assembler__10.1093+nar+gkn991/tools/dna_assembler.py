"""Source-scoped yeast assembly, using prepared neutral example units."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Row,Scene,Segment
from assembly_geometry import CircularOverlapReference
from protocol_graph import ProtocolGraph,Disposition,Terminal,InputRole as R,MaterialKind as K,StepKind as S,MolecularRelation as M
from protocol_model import Drawing,closed_circle
from protocol_views import WorkflowSection,assemble_page
from page import head,info
PAGE_FILE='dna_assembler.html'
TITLE='DNA assembler — illustrative yeast-mediated assembly'
SOURCE='Shao, Zhao and Zhao2009 doi10.1093/nar/gkn991; Methods DNA manipulation/Yeast transformation and Results Design of the DNA assembler method. Prepared-unit/plasmid-assembly branch only; chromosome integration and functional pathways excluded.'
EXAMPLE='Invented neutral prepared units and recipient DNA, not the paper constructs, functional expression cassettes or a replicating plasmid. This models intended reference order, not an experimentally validated assembly.'
OVERLAPS=('ACGTCAGTACGATCGTGACTACGAGATCTACGAGTCTCGAGTAC',
          'TCGATGACTAGCATCGTACGATGCTAGTCAGATCGACTAGCTGA',
          'GATCAGTCGATACGCTAGCATGACTCGATGACATCGTAGCTAGC')
BODIES=('CGATGACCTGATAC','TCAGGATCGATACG','ATCGTAGCTCAGTC')

def references():
    result=[]
    for i,body in enumerate(BODIES):
        left=OVERLAPS[i];right=OVERLAPS[(i+1)%3]
        if min(len(left),len(right))<=40:raise ValueError('example adopts the source greater-than40 overlap constraint')
        con=Construct((Segment('left homology',left,inferred=True),Segment('neutral body',body,inferred=True),Segment('right homology',right,inferred=True)),'prepared toy unit '+str(i+1))
        result.append((('unit_a','unit_b','recipient')[i],con,len(left),len(right)))
    return tuple(result)

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT,representation='illustrative_example',example_description=EXAMPLE)
    refs=references();materials=[]
    for identifier,con,_,_ in refs:
        materials.append(g.input(identifier,'invented prepared '+identifier,K.MOLECULAR,primary=True,state=Scene.duplex(tuple(con)),chemical_state=con,inferred=True,source=EXAMPLE))
    yeast=g.input('yeast','source Saccharomyces cerevisiae assembly/recovery host',K.BIOLOGICAL)
    geom=CircularOverlapReference(refs,'Declared source-style overlap topology on invented prepared units; no enzymatic path is supplied by this reference relation.')
    con=geom.reference('intended toy circular reference')
    final=g.product('final','INFERRED — intended yeast-assembled circular reference',K.MOLECULAR,state=closed_circle('intended toy reference',tuple(con)),chemical_state=con)
    other=g.product('other','unresolved input, alternative products and cellular material',K.MOLECULAR,state=Drawing((Row(chunks=[('Residual substrates, alternative products and cell-associated fates are unexpanded.',None,False)]),)))
    g.react('host_assembly','INFERRED — co-transform prepared units for yeast assembly',S.MAGIC,[(m,R.TEMPLATE) for m in materials]+[(yeast,R.TEMPLATE)],[final,other],geometry=geom,reason='Source in-vivo homologous recombination and recovery are an explicit boundary. Resection, strand exchange, intracellular repair, selection and yield are not inferred from matching overlaps. The toy lacks a replicon.')
    g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'No quantitative fate, complete conversion or purity asserted.'))
    g.constructed(final,'Conditional constructed-DNA reference only; actual construct function, integration and sequencing remain outside scope.')
    g.declare_reviewed_scope(('unit_a','unit_b','recipient','yeast','host_assembly','final'),SOURCE+' Exact units, overlap strings and intended body order are invented. Source OE-PCR preparation is outside this prepared-input scope.')
    return g.require_publishable()

def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_dna_assembler.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1093/nar/gkn991">Shao et al.2009</a>. INFERRED — prepared neutral units, homology strings and intended reference order. This begins after source OE-PCR/purification and recipient linearization; it does not substitute a functional pathway or model chromosome integration.'),
        WorkflowSection('Prepared overlapping units → yeast boundary → intended circular reference'),'</div>'))

def selftest():
    from checks import Check,Source,have
    from chemistry_validation import check_graph,Status
    from storage import sources_root
    c=Check();c('no contradictory declared reference relation',[x for x in check_graph(protocol_graph()) if x.status is Status.CONTRADICTORY],[])
    source=Source(sources_root()/Path(__file__).resolve().parents[1].name/'nar_gkn991_PMC2632897.xml.txt','Published Methods and Design of DNA assembler')
    if have(c,source,label='published design scope'):
        text=source.path.read_text();c('source explicitly uses co-transformation and overlaps greater than40bp','All overlaps are designed to have >40 bp' in text and 'co-transformed into S. cerevisiae' in text,True)
    c.report()
