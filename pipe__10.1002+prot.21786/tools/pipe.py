"""PIPE: a population-level assembly boundary, with no invented end coordinates."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Row
from protocol_graph import ProtocolGraph, Disposition, InputRole as R, MaterialKind as K, StepKind as S
from protocol_model import Drawing, role, closed_circle, unsealed_circular_duplex
from protocol_views import WorkflowSection, assemble_page
from page import head, info

PAGE_FILE='pipe.html'
TITLE='PIPE — incomplete-extension fragment populations'
URL='https://onlinelibrary.wiley.com/doi/10.1002/prot.21786'
SOURCE='Klock et al., Proteins71:982–994, DOI10.1002/prot.21786; Methods PIPE overview and PIPE cloning; Results Application of PIPE; Figure1 caption. Browser-readable defining original.'


def parts():
    return (role('insert',inferred=True),role('backbone',inferred=True))


def population(label):
    return Drawing((Row(chunks=[(label+': heterogeneous partial duplex population; individual strand ends and base coordinates unassigned.',None,False)]),))


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT)
    vector=g.input('vector_population','vector PCR population',K.MOLECULAR,primary=True,state=population('vector'))
    insert=g.input('insert_population','insert PCR population',K.MOLECULAR,primary=True,state=population('insert'))
    hybrid=g.product('hybrid','INFERRED — intended unsealed vector/insert hybrid class',K.MOLECULAR,
        state=unsealed_circular_duplex(parts(),statement='Candidate circular assembly only; variable gaps and nicks remain unassigned. No sealed strand is implied.'))
    g.react('anneal','Mix and anneal complementary fragment populations',S.HYBRIDIZATION,
        [(vector,R.TEMPLATE),(insert,R.STRUCTURAL)],[hybrid])
    host=g.input('host','recipient host',K.BIOLOGICAL)
    repaired=g.product('repaired','INFERRED — intended repaired construct class',K.MOLECULAR,
        state=closed_circle('intended host-repaired construct',parts()))
    g.react('recover','Host recovery',S.MAGIC,[(hybrid,R.TEMPLATE),(host,R.TEMPLATE)],[repaired],
        reason='Host repair is a source-qualified boundary; intracellular reactions and clone identity are not resolved.')
    g.constructed(repaired,'Intended constructed-DNA class; no genotype or clone purity guarantee.')
    g.declare_reviewed_scope(('vector_population','insert_population','anneal','host','recover','repaired'),SOURCE)
    return g.require_publishable()


def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap">','<h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_pipe.html">Research notes</a></p>',
        info('<a href="'+URL+'">Klock et al., published original</a>.'),
        info('INFERRED: generic roles and candidate topology. PCR products have primer-defined 5′ sequences and variable 3′ extensions; no fixed cohesive-end geometry is assigned. Unpurified populations are mixed and transformed; gaps and nicks are repaired in the host.'),
        WorkflowSection('Population assembly boundary'),
        info('The diagnostic nuclease/ligation experiment is excluded from this cloning branch. Exact oligos, target maps and sequencing readouts are outside this diagram’s scope.'),'</div>'))


def selftest():
    from checks import Check,Source,have
    from chemistry_validation import Status,check_graph
    from storage import sources_root
    check=Check()
    check('no contradictory declared molecular relation',[c for c in check_graph(protocol_graph()) if c.status is Status.CONTRADICTORY],[])
    source=Source(sources_root()/Path(__file__).resolve().parents[1].name/'prot.21786.txt','published Wiley original; browser read, ordinary fetch403')
    if have(check,source,label='local original source availability'):
        check('original DOI identity','10.1002/prot.21786' in source.path.read_text(),True)
    check.report()
