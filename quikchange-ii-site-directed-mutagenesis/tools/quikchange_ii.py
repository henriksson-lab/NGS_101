"""Agilent QuikChangeII RevisionE1: complementary primer extension and host nick repair."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Segment, revcomp
from protocol_graph import Disposition, Terminal, ProtocolGraph, InputRole as R, MaterialKind as K, StepKind as S
from protocol_model import closed_circle, role, unsealed_circular_duplex
from protocol_views import OligoEntry, OligoSection, WorkflowSection, assemble_page
from page import head,info

PAGE_FILE='quikchange_ii.html'
TITLE='QuikChange II — complementary primers and nick repair'
URL='https://www.agilent.com/cs/library/usermanuals/public/200523.pdf'
SOURCE='Agilent QuikChangeII200523/200524 instruction manual RevisionE1; printed materials p1, Introduction pp2–3/Fig1, control p4 and primer design p5; browser-readable original, local PDF fetch blocked'
FORWARD='CCATGATTACGCCAAGCGCGCAATTAACCCTCAC'
REVERSE='GTGAGGGTTAATTGCGCGCTTGGCGTAATCATGG'


def primer(name,sequence):
    if revcomp(FORWARD)!=REVERSE:
        raise ValueError('the independently printed control primers must be complementary')
    return Construct((Segment(name,sequence),),name)


def candidate_parts():
    return (Segment('source control-primer-written region',FORWARD,inferred=True),
        role('remaining control-plasmid interval; map and length unassigned',inferred=True))


def nicked():
    return unsealed_circular_duplex(candidate_parts(),statement='RevisionE1 states a circular duplex with staggered nicks after PfuUltra strand synthesis. Complete strands, nick positions and parental/new-strand mixtures remain unassigned; this is not a covalently closed strand.')


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT)
    parent=g.input('parent','source methylated parental control plasmid',K.MOLECULAR,primary=True,
        state=closed_circle('source parental plasmid',(role('parental control plasmid; full sequence unassigned'),)))
    f=g.input('forward','source control oligo1',K.OLIGO,state=primer('control oligo1',FORWARD))
    r=g.input('reverse','source control oligo2',K.OLIGO,state=primer('control oligo2',REVERSE))
    pol=g.input('pfu','PfuUltra high-fidelity polymerase',K.CATALYST)
    product=g.product('nicked','INFERRED intended staggered-nick circular-duplex class',K.MOLECULAR,state=nicked())
    carry=g.product('carryover','parental-template carryover in the same reaction',K.MOLECULAR,state=parent.state)
    g.react('synthesis','Complementary-primer whole-plasmid strand synthesis',S.EXTENSION,
        [(parent,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(pol,R.CATALYST)],[product,carry])
    dpni=g.input('dpni','DpnI',K.CATALYST)
    ready=g.product('ready','INFERRED nicked product in DpnI-treated mixture',K.MOLECULAR,state=nicked())
    remnants=g.product('remnants','digested parental-template class',K.MOLECULAR)
    g.react('digest_parent','Digest parental methylated/hemimethylated template',S.DIGESTION,
        [(product,R.TEMPLATE),(carry,R.TEMPLATE),(dpni,R.CATALYST)],[ready,remnants])
    g.mark_terminal(Terminal(remnants.id,Disposition.UNRESOLVED,'Digested-parent fragments remain in the reaction mixture. No separate physical cleanup, exact genotype or zero background is asserted.'))
    cells=g.input('cells','source competent-cell host for nick repair',K.BIOLOGICAL)
    repaired=g.product('candidate','INFERRED intended host-repaired construct class',K.MOLECULAR,
        state=closed_circle('intended repaired control construct',candidate_parts()))
    g.react('host_repair','Recover an intended nick-repaired construct in the host',S.MAGIC,
        [(ready,R.TEMPLATE),(cells,R.TEMPLATE)],[repaired],
        reason='The vendor attributes nick repair to cells after transformation. Intracellular repair, clone selection and genotype are not expanded; no in-vitro ligase closure is substituted.')
    g.constructed(repaired,'Source-expected repaired construct class with an unknown plasmid interval; not a verified clone or NGS library.')
    g.declare_reviewed_scope(('parent','forward','reverse','pfu','synthesis','dpni','digest_parent','cells','host_repair','candidate'),SOURCE)
    return g.require_publishable()


def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap">','<h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_quikchange_ii.html">Research notes</a></p>',
        info('<a href="'+URL+'">Agilent QuikChangeII manual, RevisionE1</a>.'),
        info('INFERRED: control-primer region with unassigned plasmid body and nick positions. The pre-host intermediate is an unsealed circular duplex; only the intended host-repaired endpoint is shown as a closed reference strand. Complementary overlapping primers are kept distinct from back-to-back inverse-PCR geometry.'),
        OligoSection('Published control oligos',(OligoEntry('control oligo1',tuple(primer('oligo1',FORWARD))),OligoEntry('control oligo2',tuple(primer('oligo2',REVERSE))))),
        WorkflowSection('Strand synthesis, parental digestion and host repair'),
        info('No sample index or sequencing read is assigned. Source oligos and complementarity do not establish genotype, reaction yield or quantitative template removal.'),'</div>'))


def selftest():
    from checks import Check,Source,have
    from chemistry_validation import Status,check_graph
    from storage import sources_root
    check=Check()
    check('no contradictory source-scoped molecular claims',[c for c in check_graph(protocol_graph()) if c.status is Status.CONTRADICTORY],[])
    source=Source(sources_root()/Path(__file__).resolve().parents[1].name/'200523.pdf.txt','Agilent original official200523 manual RevisionE1; ordinary local fetch currently403')
    if have(check,source,label='original vendor primer transcription'):
        text=source.path.read_text();compact=''.join(text.split())
        check('reviewed original revision retained','Revision E1' in text,True)
        check('control oligo1 in original',FORWARD in compact,True)
        check('control oligo2 in original',REVERSE in compact,True)
    check.report()
