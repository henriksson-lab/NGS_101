"""Original two-step EMP cloning; no inherited RF/one-step recipe."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, complement_segments
from circular import CircularizationGeometry
from pcr import PCRGeometry
from protocol_graph import ProtocolGraph, Disposition, InputRole as R, MaterialKind as K, StepKind as S
from protocol_model import role, closed_circle
from protocol_views import WorkflowSection, assemble_page
from page import head, info

PAGE_FILE='emp_cloning.html'
TITLE='EMP cloning — megaprimer-assisted inverse PCR'
SOURCE='Ulrich et al.2012, doi:10.1371/journal.pone.0053360; original Figure1A, Methods Exponential Megapriming PCR, synthesis, insertion, in-vitro ligation/transformation; published TableS1'


def part(name):
    return role(name,inferred=True)


def insert_parts():
    return (part('insert 5prime / F1 site'),part('insert body'),part('insert 3prime / R1 core'))


def recipient_parts():
    # Cut the circular display immediately downstream of the intended insertion.
    return (part('B downstream landing site'),part('retained backbone'),part('upstream R2 site'))


def product_parts():
    return (*insert_parts(),*recipient_parts())


def duplex(parts,label,phosphorylated=False):
    sc=Scene()
    bottom=complement_segments(parts)
    mod='P' if phosphorylated else ''
    sc.strand('top',parts,label=label,mod5=mod)
    sc.anneal('bottom',bottom,to='top',pair=(bottom[0].name,parts[-1].name),
              unpaired=tuple(s.name for s in bottom),mod5=mod)
    return sc


def landing_view():
    sc=Scene()
    receiver=complement_segments(recipient_parts())
    mega=(*insert_parts(),part('B downstream landing site'))
    sc.strand('recipient strand',receiver,label='recipient strand at the insertion boundary')
    sc.anneal('megaprimer active strand',mega,to='recipient strand',
              pair=('B downstream landing site',"B downstream landing site'"),
              unpaired=tuple(s.name for s in mega))
    sc.arrow('megaprimer active strand','extension')
    return sc


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT)
    source=g.input('insert_source','insert template; neutral unassigned payload',K.MOLECULAR,
                   state=duplex(insert_parts(),'insert template'),primary=True)
    recipient=g.input('recipient','circular recipient plasmid',K.MOLECULAR,
                      state=closed_circle('recipient',recipient_parts()),primary=True)
    f1=g.input('f1','F1: insert forward primer, reused in step2',K.OLIGO,
               state=Construct((part('insert 5prime / F1 site'),),'F1'))
    r1=g.input('r1','R1: 5prime recipient homology + insert reverse core',K.OLIGO,
               state=Construct((part('reverse-complement B tail'),part('reverse-complement insert 3prime')),'R1'))
    r2=g.input('r2','R2: reverse primer upstream of insertion',K.OLIGO,
               state=Construct((part('reverse-complement upstream R2 site'),),'R2'))
    pol=g.input('phusion','Phusion DNA polymerase',K.CATALYST)
    pnk=g.input('pnk','T4 polynucleotide kinase',K.CATALYST)
    ligase=g.input('ligase','T4 DNA ligase',K.CATALYST)
    dpni=g.input('dpni','DpnI',K.CATALYST)
    host=g.input('host','E.coli colony recovery and screening',K.SELECTION)
    mega_parts=(*insert_parts(),part('B downstream landing site'))
    mega=g.product('megaprimer','INFERRED — duplex insert plus recipient homology B',K.MOLECULAR,
                   state=duplex(mega_parts,'first-PCR megaprimer product'))
    g.react('first_pcr','Amplify insert and add the recipient homology tail',S.PCR,
            [(source,R.TEMPLATE),(f1,R.STRUCTURAL),(r1,R.STRUCTURAL),(pol,R.CATALYST)],[mega],
            geometry=PCRGeometry.unavailable(Construct(insert_parts(),'insert template'),
              reason='Figure1A role architecture only. A particular insert, primer sequences and exact binding coordinates are not assigned.'))
    purified=g.product('purified_megaprimer','INFERRED — purified megaprimer',K.MOLECULAR,state=mega.state)
    removed=g.product('first_cleanup_waste','material removed during PCR-product purification',K.WASTE)
    g.react('first_cleanup','Purify the first-PCR product',S.CLEANUP,[(mega,R.TEMPLATE)],[purified,removed],changes_molecule=False)
    g.discard(removed,'Purification discard; exact composition and recovery not assigned.')
    linear=g.product('linear','INFERRED — linear insert–backbone product',K.MOLECULAR,
                     state=duplex(product_parts(),'linear second-PCR product'))
    g.react('second_pcr','Megaprimer initiates; F1 and R2 amplify the linear product',S.PCR,
            [(recipient,R.TEMPLATE),(purified,R.STRUCTURAL),(f1,R.STRUCTURAL),(r2,R.STRUCTURAL),(pol,R.CATALYST)],[linear],
            geometry=PCRGeometry.unavailable(landing_view(),
              reason='The active megaprimer strand binds B downstream of insertion; R2 binds upstream on the circle. Later F1/R2 cycles use the linear product. Exact target maps are not assigned; no invented inward primer pair is imposed on the starting circle.'))
    clean=g.product('purified_linear','INFERRED — purified linear second-PCR product',K.MOLECULAR,state=linear.state)
    waste2=g.product('second_cleanup_waste','material removed in second product purification',K.WASTE)
    g.react('second_cleanup','Purify second-PCR DNA',S.CLEANUP,[(linear,R.TEMPLATE)],[clean,waste2],changes_molecule=False)
    g.discard(waste2,'Purification discard; parental-template clearance is not assumed.')
    phosphorylated=g.product('phosphorylated','INFERRED — 5prime-phosphorylated linear product',K.MOLECULAR,
                             state=duplex(product_parts(),'PNK-treated product',True))
    g.react('phosphorylate','Phosphorylate product ends',S.PHOSPHORYLATION,
            [(clean,R.TEMPLATE),(pnk,R.CATALYST)],[phosphorylated])
    circular=g.product('circular','INFERRED — expected circular insert–backbone candidate',K.MOLECULAR,
                       state=closed_circle('ligated candidate',product_parts(),closure_label='ligation at insert-start / upstream-backbone boundary'))
    g.react('ligate','Circularize phosphorylated product with T4 DNA ligase',S.CIRCULARIZATION,
            [(phosphorylated,R.TEMPLATE),(ligase,R.CATALYST)],[circular],
            geometry=CircularizationGeometry.unavailable('The source establishes in-vitro phosphorylation and ligation, but this role-only instance does not assign exact end bases, terminal structures or junction coordinates.'))
    digested=g.product('after_dpni','INFERRED — candidate after parental-template digestion',K.MOLECULAR,state=circular.state)
    parental=g.product('parental_fragments','parental-template digest fragments',K.MOLECULAR,
                       state=Construct((part('parental digest fragments'),),'parental fragments'))
    g.react('digest_parent','DpnI digest remaining parental plasmid after ligation',S.DIGESTION,
            [(circular,R.TEMPLATE),(dpni,R.CATALYST)],[digested,parental])
    candidate=g.product('candidate','INFERRED — recovered construct candidate',K.MOLECULAR,state=circular.state)
    rejected=g.product('unselected','material not carried into recovered colonies',K.WASTE)
    g.react('recover','Recover and screen candidate colonies',S.MAGIC,
            [(digested,R.TEMPLATE),(parental,R.TEMPLATE),(host,R.SELECTION)],[candidate,rejected],
            reason='Colony recovery and screening are not a molecular proof of genotype or insert function.')
    g.discard(rejected,'Material not carried into recovered candidates; no yield or detailed molecular fate assigned.')
    g.constructed(candidate,'Expected two-step EMP product. No complete plasmid map, particular experimental clone, insert function or NGS run is assigned.')
    g.declare_reviewed_scope(('insert_source','recipient','f1','r1','r2','phusion','pnk','ligase','dpni','host','first_pcr','first_cleanup','second_pcr','second_cleanup','phosphorylate','ligate','digest_parent','recover','candidate'),SOURCE)
    return g.require_publishable()


def render_page():
    return assemble_page(protocol_graph(),(
        head(TITLE),'<div class="wrap">',f'<h1>{TITLE}</h1>',
        '<p class="research-notes"><a href="01_emp_cloning.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1371/journal.pone.0053360">Ulrich et al., PLOS ONE (2012)</a>; published Figure1A, Methods and TableS1.'),
        info('Original two-step EMP branch. The insertion topology is source-defined; neutral payload, primer sites and backbone are inferred role placeholders. No target sequence or experimental clone is reconstructed. The first product is duplex DNA; its active strand exposes the recipient homology after denaturation, rather than carrying an asserted persistent single-stranded overhang.'),
        info('Published discrepancy: Methods, second EMP PCR (p3), specifies 25 cycles; TableS1 second-EMP cycle scheme specifies 35. Its adjacent lower-left heading reads “2nd RF PCR”. The model follows the Methods reaction identity and records 25 cycles in the note; no molecular claim depends on the disputed cycle count. The intended correction is unconfirmed.'),
        WorkflowSection('Two-step EMP assembly',paths=(('purified_megaprimer',()),('candidate',()))),
        '<h2>Output</h2>',info('Constructed DNA. One-step EMP, multi-insert EMP and restriction-free cloning comparisons are separate branches outside this schematic. No NGS index, library or run is assigned.'),'</div>'))


def selftest():
    from checks import Check,Source,have
    from chemistry_validation import check_graph,Status
    from storage import sources_root
    c=Check()
    c('no contradictory scoped molecular claims',[x for x in check_graph(protocol_graph()) if x.status is Status.CONTRADICTORY],[])
    s=Source(sources_root()/Path(__file__).resolve().parents[1].name/'journal.pone.0053360_PMC3534072.pdf.txt','Published Methods in-vitro ligation/transformation')
    if have(c,s,label='published EMP reagent order'):
        t=s.path.read_text()
        order=[t.find(x,t.find('In vitro ligation and transformation of EMP')) for x in ('5 U T4 PNK','200 U of T4 DNA ligase','10 U DpnI')]
        c('source orders kinase before ligase before DpnI',all(x>=0 for x in order) and order==sorted(order))
    c.report()
