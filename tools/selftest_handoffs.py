"""Causal contract checks for versioned protocol material interfaces."""
from pathlib import Path
import sys
from dataclasses import replace
sys.path[:0]=[str(Path(__file__).resolve().parents[1]/'lib')]
from checks import Check
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Disposition,Terminal
from protocol_handoff import export_material,import_material,declare_input_origin
from chemdraw import Construct,Segment,feature
from protocol_views import handoff_section


def producer(*, drawable=True):
    g=ProtocolGraph('Synthetic source',authority='Synthetic fixture, not laboratory chemistry',required_output=Disposition.STORED)
    x=g.input('sample','input',K.MOLECULAR,state=Construct([Segment('body','ACGT')]),primary=True)
    y=g.product('cdna','unfragmented DNA',K.MOLECULAR,state=Construct([Segment('body','ACGT'),Segment('umi','NNN',placeholder=True,feature=feature('fixture.umi','umi','random'))])if drawable else None,chemical_state={'stage':'full_length','fragmented':False},inferred=True,source='Explicit synthetic inferred example')
    enzyme=g.input('enzyme','synthetic fixture enzyme',K.CATALYST)
    g.react('prepare','Synthetic preparation',S.EXTENSION,((x,R.TEMPLATE),(enzyme,R.CATALYST)),(y,))
    g.mark_terminal(Terminal(y.id,Disposition.STORED,'Intermediate fixture'))
    g.declare_reviewed_scope(('sample','enzyme','prepare','cdna'),'Synthetic fixture')
    e=export_material(g,'full_length',y,protocol='fixture-source',version='Rev A',page='source.html',source='Synthetic fixture',limitation='Synthetic inference retained; not a chemistry validation claim')
    return g,e


def consumer(e):
    g=ProtocolGraph('Synthetic consumer',authority='Synthetic fixture',required_output=Disposition.STORED)
    x=import_material(g,'input',e,expected_protocol='fixture-source',expected_version='Rev A',source='Synthetic affirmative handoff',primary=True,requirements={'fragmented':False})
    y=g.product('output','retained input',K.MOLECULAR,state=x.state,chemical_state=x.chemical_state,inferred=x.inferred)
    g.react('retain','Synthetic retention',S.GENERIC,((x,R.TEMPLATE),),(y,),changes_molecule=False)
    g.mark_terminal(Terminal(y.id,Disposition.STORED,'Intermediate fixture'))
    g.declare_reviewed_scope(('input','retain','output'),'Synthetic fixture')
    return g


def run(c):
    p,e=producer();g=consumer(e);g.require_publishable()
    c('handoff preserves exact declared bases',g.materials['input'].state.top(),e.material.state.top())
    c('handoff preserves identifier semantics',g.materials['input'].state.get('umi').feature,e.material.state.get('umi').feature)
    c('handoff cannot promote inferred state',g.materials['input'].inferred,True)
    c('handoff copies mutable declarations',g.materials['input'].state is e.material.state,False)
    c('page renders exact upstream revision and anchor','Rev A'in handoff_section(g)and'export-full_length'in handoff_section(g),True)
    chemical_only=consumer(producer(drawable=False)[1]);chemical_only.require_publishable()
    c('chemical-only handoff renders without inventing sequence geometry','fixture-source'in handoff_section(chemical_only),True)
    c.raises('wrong source version rejected',lambda:import_material(ProtocolGraph('bad'), 'x',e,expected_protocol='fixture-source',expected_version='Rev B',source='test'))
    c.raises('wrong source protocol rejected',lambda:import_material(ProtocolGraph('bad'), 'x',e,expected_protocol='another-source',expected_version='Rev A',source='test'))
    c.raises('fragmented/wrong-stage contracts rejected',lambda:import_material(ProtocolGraph('bad'),'x',e,expected_protocol='fixture-source',expected_version='Rev A',source='test',requirements={'fragmented':True}))
    c.raises('a missing downstream primer site rejects the handoff',lambda:import_material(ProtocolGraph('bad'),'x',e,expected_protocol='fixture-source',expected_version='Rev A',source='test',required_sites={'missing site':'GAGAGAGAGAGA'}))
    c.raises('a variable UMI cannot prove a concrete primer site',lambda:import_material(ProtocolGraph('bad'),'x',e,expected_protocol='fixture-source',expected_version='Rev A',source='test',required_sites={'unknown site':'AAA'}))
    c.raises('unprocessed inputs cannot masquerade as exports',lambda:export_material(p,'bad','sample',protocol='fixture-source',version='Rev A',page='source.html',source='test',limitation='fixture'))
    c.raises('undeclared export rejected',lambda:export_material(p,'missing','absent',protocol='fixture-source',version='Rev A',page='source.html',source='test',limitation='fixture'))
    g.materials['input'].state.segments[0]=Segment('body','AAAA')
    c.raises('consumer mutation is detected',g.require_publishable)
    p2,e2=producer();g2=consumer(e2);p2.materials['cdna'].state.segments[0]=Segment('body','AAAA')
    c.raises('producer mutation is detected',g2.require_publishable)
    p3,e3=producer();g3=consumer(e3);del p3.exports['full_length']
    c.raises('retired export is detected',g3.require_publishable)
    p4,e4=producer();g4=consumer(e4)
    declare_input_origin(p4,'sample',protocol='named-source',stage='source-named input',page='named.html',source='Synthetic located claim',limitation='No exact version binding')
    c('source-named reference does not import a state',len(p4.imports),0)
    c('source-named reference remains qualified','no imported molecular state'in handoff_section(p4),True)
    c.raises('source-named reference cannot target a produced output',lambda:declare_input_origin(p4,'cdna',protocol='source',stage='stage',page='page.html',source='test',limitation='test'))
    external=consumer(producer()[1]); external.imports.clear()
    origin=declare_input_origin(external,'input',title='Unmodelled precursor',url='https://doi.org/10.1/fixture',stage='Prepared DNA before fragmentation',source='Synthetic located claim',limitation='No producer model; no state inheritance')
    c('external precursor retains original URL and stage',origin.href,'https://doi.org/10.1/fixture')
    c('external precursor renders its method name','Unmodelled precursor'in handoff_section(external),True)
    c('external precursor does not invent a producer model',origin.protocol,'')
    c.raises('mixed external and local provenance rejected',lambda:declare_input_origin(consumer(producer()[1]),'input',protocol='source',page='page.html',title='source',url='https://doi.org/10.1/fixture',stage='stage',source='source',limitation='limit'))
    raw_external=ProtocolGraph('URL fixture');raw_external.input('x','x',K.MOLECULAR)
    c.raises('unsafe upstream URL rejected',lambda:declare_input_origin(raw_external,'x',title='source',url='javascript:alert(1)',stage='stage',source='source',limitation='limit'))
    # An actual dependency cycle must be caught without unbounded recursion.
    g4.exports['other']=replace(e4,graph=g4,reference=replace(e4.reference,name='other',material='output',producer='retain'))
    p4.imports['sample']=replace(g4.imports['input'],export=g4.exports['other'],expected_protocol=e4.reference.protocol)
    c('cross-protocol dependency cycle is reported',any('cycle'in x for x in p4.publication_problems()),True)
    raw=ProtocolGraph('raw',authority='test',required_output=Disposition.STORED);x=raw.input('x','raw input',K.MOLECULAR,primary=True,state=Construct([Segment('x','ACGT')]))
    raw.mark_terminal(Terminal(x.id,Disposition.STORED,'not processed'));raw.declare_reviewed_scope(('x',),'test')
    c.raises('stored boundary requires an actual producer',raw.require_publishable)
    ngs=consumer(producer()[1]);ngs.required_output=Disposition.SEQUENCING
    c.raises('stored intermediate cannot satisfy an NGS endpoint',ngs.require_publishable)

if __name__=='__main__':
    c=Check();run(c);c.report()
