"""CPEC: source-scoped generic single-insert topology, not a target design."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene, complement_segments
from protocol_graph import Disposition, ProtocolGraph, InputRole as R, MaterialKind as K, StepKind as S
from protocol_model import role, closed_circle, unsealed_circular_duplex
from protocol_views import WorkflowSection, assemble_page
from page import head, info

PAGE_FILE='cpec.html'
TITLE='CPEC — circular polymerase extension cloning'
SOURCE='Quan and Tian2009, doi10.1371/journal.pone.0006441; published Results Single gene cloning using CPEC and Figure1A. Generic single-insert mechanism only.'

def duplex(label):
    parts=(role(label+'; exact overlap and internal bases unassigned',inferred=True),)
    bottom=complement_segments(parts)
    sc=Scene(); sc.strand('top',parts,label=label)
    sc.anneal('bottom',bottom,to='top',pair=(bottom[0].name,parts[0].name),unpaired=tuple(s.name for s in bottom))
    return sc

def parts():
    return (role('backbone; unassigned map',inferred=True),role('insert; unassigned sequence',inferred=True))

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT)
    vector=g.input('vector','linear duplex backbone with source-required terminal overlaps',K.MOLECULAR,primary=True,state=duplex('linear backbone'))
    insert=g.input('insert','linear duplex insert with matching terminal overlaps',K.MOLECULAR,primary=True,state=duplex('linear insert'))
    polymerase=g.input('polymerase','source high-fidelity DNA polymerase; identity outside this scope',K.CATALYST)
    nicked=g.product('nicked','INFERRED unsealed circular duplex; one nick per strand',K.MOLECULAR,state=unsealed_circular_duplex(parts(),statement='Figure1A defines mutual fragment-templated extension and one nick per strand. Exact overlap bases, strand ends and nick coordinates are not assigned.'))
    g.react('extension','INFERRED topology — overlapping fragments prime extension on each other',S.EXTENSION,[(vector,R.TEMPLATE),(insert,R.TEMPLATE),(polymerase,R.CATALYST)],[nicked])
    host=g.input('host','source transformation host; intracellular processes outside scope',K.BIOLOGICAL)
    recovered=g.product('construct','INFERRED intended recovered insert-containing plasmid',K.MOLECULAR,state=closed_circle('intended recovered reference strand',parts()))
    g.react('recovery','INFERRED intended construct — host recovery after transformation',S.MAGIC,[(nicked,R.TEMPLATE),(host,R.TEMPLATE)],[recovered],reason='Source reports recovered sequence-verified clones. Intracellular repair and clone selection are an explicit unexpanded boundary; extension alone is not covalent closure.')
    g.constructed(recovered,'Intended generic construct class only; clone identity, efficiency and alternatives are not predicted.')
    g.declare_reviewed_scope(('vector','insert','polymerase','extension','host','recovery','construct'),SOURCE)
    return g.require_publishable()

def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_cpec.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1371/journal.pone.0006441">Quan and Tian (2009)</a>, Results and Figure1A.'),
        info('INFERRED: generic backbone/insert roles, no selected target, complete map, overlap bases or nick coordinates. The unsealed extension intermediate is distinct from the intended host-recovered construct.'),
        WorkflowSection('Linear fragments → nicked circular duplex → host recovery'),
        info('Source fragment-templated extension uses no added PCR primers. This scope begins with prepared linear fragments and excludes their preparation, multi-fragment/library variants and sequencing.'),'</div>'))

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    from chemistry_validation import Status,check_graph
    check=Check()
    check('no contradictory molecular claims',[c for c in check_graph(protocol_graph()) if c.status is Status.CONTRADICTORY],[])
    src=Source(sources_root()/Path(__file__).resolve().parents[1].name/'journal.pone.0006441_PMC2713398.xml.txt','Published article Results and Figure1A')
    if have(check,src,label='published CPEC mechanism'):
        text=' '.join(src.path.read_text().split())
        check('original specifies nicked duplex','leaving only one nick in each strand' in text,True)
        check('original excludes added PCR primers','without adding any PCR primers' in text,True)
    check.report()
