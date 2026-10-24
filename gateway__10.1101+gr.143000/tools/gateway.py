from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Row, Scene, Segment, revcomp
from cohesive_assembly import CohesiveFragment, CircularCohesiveAssembly
from assembly_geometry import StrandRetention
from protocol_graph import ProtocolGraph, Disposition, Terminal, InputRole as R, MaterialKind as K, StepKind as S, MolecularRelation as M
from protocol_model import Drawing, closed_circle
from protocol_views import WorkflowSection, assemble_page
from page import head,info

REPRESENTATION='illustrative_example'
EXAMPLE='Invented short, nonfunctional DNA parts demonstrate the declared molecular relationships. They were not used in the cited experiments and are not a working vector or validated laboratory design.'

def base(name,sequence):
    return Segment(name,sequence,inferred=True)

def duplex(top,bottom=None):
    topseq=top.top(); bottomseq=revcomp(topseq) if bottom is None else bottom
    sc=Scene();sc.strand('top',(base('top reference',topseq),))
    sc.anneal('bottom',(base('bottom reference',bottomseq),),to='top',pair=('bottom reference','top reference'))
    return sc

from protocol_model import role
PAGE_FILE='gateway.html'
TITLE='Gateway — illustrative BP/LR site-identity example'
SOURCE='Hartley et al.2000 doi10.1101/gr.143000 (defining original retrieval blocked); mechanism scope follows Reece-Hoyes and Walhout2018 Gateway Recombinational Cloning doi10.1101/pdb.top094912, Basic principles and Figures1–4 text. This later explanatory source is explicitly distinct from the original experiments.'
INSERT='CGATGACCTG'
DONOR='TCAGGATCGA'
DESTINATION='ATCGTACGCA'
DROPOUT='CATGCAGTAC'

def site(name):
    return role(name+'; exact recognition bases unassigned',inferred=True)

def cassette(first,body,last,backbone=None):
    parts=(site(first),base('invented central body',body),site(last))
    if backbone is not None:parts+= (base('invented backbone body',backbone),)
    return Construct(parts,'symbolic site cassette')

def circle(con,label):
    return closed_circle(label,tuple(con))

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT,representation=REPRESENTATION,example_description=EXAMPLE+' Actual att-site sequences are unassigned; their symbolic identities are not invented functional recognition sites.')
    insert_con=cassette('attB1',INSERT,'attB2')
    donor_con=cassette('attP1',DROPOUT,'attP2',DONOR)
    destination_con=cassette('attR1',DROPOUT,'attR2',DESTINATION)
    insert=g.input('insert','toy linear attB-flanked insert',K.MOLECULAR,primary=True,state=Scene.duplex(tuple(insert_con),unpaired=tuple(s.name+"'" for s in insert_con if s.is_role_token())),inferred=True,source=EXAMPLE)
    donor=g.input('donor','toy circular attP donor',K.MOLECULAR,primary=True,state=circle(donor_con,'toy donor circle'),inferred=True,source=EXAMPLE)
    bp=g.input('bp','source BP recombination enzyme mixture',K.CATALYST)
    entry_con=cassette('attL1',INSERT,'attL2',DONOR)
    entry=g.product('entry','INFERRED — intended attL entry reference',K.MOLECULAR,state=circle(entry_con,'symbolic entry reference'))
    bp_by=g.product('bp_byproduct','INFERRED — attR-containing byproduct class; exact strands unassigned',K.MOLECULAR,state=Drawing((Row(chunks=[("Byproduct exact strands and topology unassigned; retained in the conceptual reaction mixture.",None,False)]),)))
    g.react('bp_reaction','INFERRED — BP site-specific recombination: B/P → L/R',S.SITE_SPECIFIC_RECOMBINATION,[(insert,R.TEMPLATE),(donor,R.TEMPLATE),(bp,R.CATALYST)],[entry,bp_by],reason='Only source-declared site identities and intended cassette transfer are represented. Exact att sequences, crossover coordinates, byproduct strands and alternative products are unresolved; no sequence-resolved recombination validation is claimed.')
    g.mark_terminal(Terminal(bp_by.id,Disposition.UNRESOLVED,'attR-containing byproduct class remains unexpanded; no purification or complete conversion is claimed.'))
    destination=g.input('destination','toy circular attR destination',K.MOLECULAR,primary=True,state=circle(destination_con,'toy destination circle'),inferred=True,source=EXAMPLE)
    lr=g.input('lr','source LR recombination enzyme mixture',K.CATALYST)
    final_con=cassette('attB1',INSERT,'attB2',DESTINATION)
    final=g.product('final','INFERRED — intended attB destination reference',K.MOLECULAR,state=circle(final_con,'symbolic final reference'))
    lr_by=g.product('lr_byproduct','INFERRED — attP-containing byproduct class; exact strands unassigned',K.MOLECULAR,state=Drawing((Row(chunks=[("Byproduct exact strands and topology unassigned; retained in the conceptual reaction mixture.",None,False)]),)))
    g.react('lr_reaction','INFERRED — LR site-specific recombination: L/R → B/P',S.SITE_SPECIFIC_RECOMBINATION,[(entry,R.TEMPLATE),(destination,R.TEMPLATE),(lr,R.CATALYST)],[final,lr_by],reason='Symbolic site-identity bookkeeping only; exact strand exchange and biochemical intermediates require actual recognition sites and source maps.')
    g.mark_terminal(Terminal(lr_by.id,Disposition.UNRESOLVED,'attP-containing byproduct class is unexpanded. Selection and residual parental molecules are outside this intended-product example.'))
    g.constructed(final,'Symbolic illustrative destination construct, not a complete functional vector, verified experimental clone or NGS library.')
    g.declare_reviewed_scope(('insert','donor','bp','bp_reaction','destination','lr','lr_reaction','final'),SOURCE+' DNA bodies are invented, exact recognition bases unassigned.')
    return g.require_publishable()

def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_gateway.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1101/gr.143000">Defining method (2000)</a>; <a href="https://doi.org/10.1101/pdb.top094912">explicit explanatory source (2018)</a>. INFERRED — invented DNA bodies and symbolic att-site roles. Exact recombination bases and strand exchange are unassigned, rather than invented.'),
        WorkflowSection('BP entry formation → LR destination transfer'),'</div>'))

def selftest():
    from checks import Check
    from chemistry_validation import check_graph,Status
    c=Check();c('no contradictory declared molecular relations',[x for x in check_graph(protocol_graph()) if x.status is Status.CONTRADICTORY],[]);c.report()
