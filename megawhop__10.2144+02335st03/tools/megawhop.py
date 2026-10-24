"""Abstract-scoped MEGAWHOP draft; product topology deliberately unassigned."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Row
from protocol_model import Drawing
from protocol_graph import ProtocolGraph, Disposition, InputRole as R, MaterialKind as K, StepKind as S
from protocol_views import assemble_page, WorkflowSection
from page import head,info
PAGE_FILE='megawhop.html'
TITLE='MEGAWHOP — abstract-scoped draft'
SOURCE='Miyazaki and Takenouchi2002, doi10.2144/02335st03; defining PubMed abstract PMID12449380 only. Full original methods are deferred.'

def unknown(statement):
    return Drawing((Row(chunks=[(statement,None,False)]),))

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT)
    parent=g.input('parent','source parental plasmid',K.MOLECULAR,primary=True,state=unknown('Parental plasmid; map and exact strands unassigned.'))
    fragment=g.input('megaprimer','source insert-fragment megaprimer',K.MOLECULAR,primary=True,state=unknown('Megaprimer fragment; sequence, strand assignment and binding geometry deferred.'))
    polymerase=g.input('polymerase','whole-plasmid synthesis polymerase; identity deferred',K.CATALYST)
    product=g.product('product','unresolved whole-plasmid synthesis product class',K.MOLECULAR,state=unknown('INFERRED product class only. Topology, nicks, end chemistry and exact sequence are unassigned; no closed circle is drawn.'))
    g.react('synthesis','INFERRED product class — fragment-primed whole-plasmid synthesis',S.EXTENSION,[(parent,R.TEMPLATE),(fragment,R.STRUCTURAL),(polymerase,R.CATALYST)],[product])
    g.constructed(product,'Abstract-scoped DNA synthesis endpoint only. Full original reaction, parental-template handling and recovered clone are not claimed; catalogue status remains notes.')
    g.declare_reviewed_scope(('parent','megaprimer','polymerase','synthesis','product'),SOURCE)
    return g.require_publishable()

def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_megawhop.html">Research notes</a></p>',
        info('<a href="https://pubmed.ncbi.nlm.nih.gov/12449380/">Defining abstract, Miyazaki and Takenouchi (2002)</a>.'),
        info('DRAFT — INFERRED generic synthesis-product class. Complete original methods, exact sequences, product topology and downstream processing remain unassigned.'),
        WorkflowSection('Fragment-primed synthesis: unresolved endpoint'),'</div>'))

def selftest():
    from checks import Check
    from chemistry_validation import Status,check_graph
    check=Check()
    check('no contradictory molecular claims',[c for c in check_graph(protocol_graph()) if c.status is Status.CONTRADICTORY],[])
    check.report()
