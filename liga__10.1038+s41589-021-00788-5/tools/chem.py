from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, feature, workflow_from_protocol_graph, workflow_panel
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan
from page import head,info
PAGE_FILE='chem.html'

def operations(g):
    def op(key,x,label,reagent=None,state=None,kind=S.GENERIC,reagent_kind=K.CONDITION,reagent_role=R.CONDITION):
        inputs=[(x,R.TEMPLATE)]
        if reagent:
            z=g.input(key+'_reagent',reagent,reagent_kind);inputs.append((z,reagent_role))
        y=g.product(key+'_product',label+'; preparation includes unremoved reaction components',K.MOLECULAR,state=state)
        g.react(key,label,kind,inputs,(y,));return y
    def clean(key,x,label):
        z=g.input(key+'_selection',label,K.SELECTION)
        y=g.product(key+'_kept','Physically retained source fraction: '+label,K.MOLECULAR)
        w=g.product(key+'_removed','Physically removed fractions: '+label,K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(z,R.SELECTION)),(y,w));g.discard(w,'Source physical workup only; no perfect specificity or purity assertion.');return y
    return op,clean

def end(g,x,label):g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=label))
def incomplete(g):
    g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE)
    g.add_source_item(SourceItem('unassembled','Unreviewed linked originals and full molecular geometry',LIMIT));return g


from protocol_graph import Terminal,Disposition

def save(g,key,x,label):
    take=g.product(key+'_aliquot',label,K.MOLECULAR)
    rest=g.product(key+'_remainder','Remaining preparation; later use unassigned',K.MOLECULAR)
    g.react(key,'Physical aliquot: '+label,S.GENERIC,((x,R.TEMPLATE),),(take,rest))
    g.terminals[rest.id]=Terminal(rest.id,Disposition.STORED,'Remainder retained explicitly; later fate unassigned.')
    return take


def measure(g,key,x,label):
    aliquot=g.product(key+'_aliquot','Separate source analytical aliquot;exact amount unassigned unless stated',K.MOLECULAR)
    keep=g.product(key+'_production','Production preparation retained separately from analytical aliquot',K.MOLECULAR)
    g.react(key+'_split','Separate analytical sample: '+label,S.SPLIT,((x,R.TEMPLATE),),(aliquot,keep))
    measured=g.product(key+'_qc','Analytically assayed aliquot: '+label,K.QC)
    g.react(key,label,S.GENERIC,((aliquot,R.TEMPLATE),),(measured,))
    g.terminals[measured.id]=Terminal(measured.id,Disposition.QC,'Analytical aliquot not the production material pooled for sequencing.')
    return keep

def render_page():
    g=protocol_graph()
    return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_chem.html">Research notes</a></p>'+workflow_panel(workflow_from_protocol_graph(g,ENDPOINT))+'</main>'

TITLE='LiGA — version-qualified glycan-phage selection and NGS interface'
SOURCE='Published SI MOESM53 Fig1–3;associated accepted manuscript PMC8380037 Methods1.8–1.15'
LIMIT='Private incomplete source graph. Operational route below is explicitly associated accepted-manuscript production/bead-capture scope,not verified whole-publication equivalence. Published SI Fig1 viability has intermediate Zeba purification and8h glycosylation;that experiment is not merged into the production route. Full phage genome/chemical sites/primer geometry unassembled;no quantitative coupling or binding prediction.'
ENDPOINT='library_product'
def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE);op,clean=operations(g)
    x=g.input('phage','Source preconstructed silently DNA-barcoded M13 phage clones;genome sequence and manufacture unassembled',K.BIOLOGICAL,primary=True)
    x=op('acylation',x,'Associated MS production1.8:DBCO-NHS0.1–1.5mM,45min incubation;phage modification distribution not computed','DBCO-NHS in DMF/PBS;hydrolyzed and unreacted reagent may remain')
    x=op('glycosylation',x,'Associated MS production1.8:add glycan-azide directly to2–3mM;overnight4°C','Glycan-azide;source SPAAC product ensemble not guaranteed uniform')
    x=clean('zeba',x,'Source production Zeba Spin purification;retain glycan-phage fraction')
    x=measure(g,'maldi',x,'Source MALDI-TOF analysis;separate2µL analytical sample exposed to matrix/organic solvent')
    x=op('mix',x,'Mix desired glycan-phage clones,108PFU each;physical mixture distinct from external barcode→glycan dictionary')
    x=op('bind',x,'Selected MS1.11.3 CD22-Fc/ProteinG magnetic-bead binding;other plate/cell/lectin recipes not merged','Source prepared CD22-Fc ProteinG beads;manufacturer preparation not reviewed')
    x=clean('wash',x,'KingFisher2HEPES+1%Tween washes then1HEPES wash;retain bead-associated preparation')
    x=op('water',x,'Transfer beads to30µL nuclease-free water')
    t=g.product('titer_aliquot','2µL suspension for titering;not the PCR aliquot',K.MOLECULAR);rest=g.product('heat_input','Remaining suspension after titer aliquot',K.MOLECULAR)
    g.react('titer_split','Remove2µL titer aliquot from30µL',S.GENERIC,((x,R.TEMPLATE),),(t,rest));g.terminals[t.id]=Terminal(t.id,Disposition.STORED,'Source titer aliquot;not NGS library.')
    x=op('heat',rest,'Remaining preparation90°C15min;DNA/particle remnants and beads not magically removed')
    y=g.product('heat_supernatant','Source25µL supernatant PCR input',K.MOLECULAR);rem=g.product('heat_residue','Residual liquid/beads/debris;not source PCR fraction',K.MOLECULAR)
    g.react('supernatant','Source physically collect25µL supernatant',S.CLEANUP,((x,R.TEMPLATE),),(y,rem));g.terminals[rem.id]=Terminal(rem.id,Disposition.UNRESOLVED,'Residual fraction not automatically destroyed or sequenced.')
    x=op('pcr',y,'Associated MS1.15:50µL Phusion PCR with indexed forward/reverse primers;new dsDNA copies distinct from input phage material','Source barcoded primers unassembled;Phusion/dNTP/MgCl2 buffer mixture;thermal details deferred')
    x=measure(g,'agarose',x,'Source2%agarose analytical quantification;not the preparative eGel step')
    x=op('pool',x,'Pool different indexed PCR products,10ng each;poolmembership unassigned')
    x=clean('egel',x,'Source eGel purification;retain selected DNA')
    x=measure(g,'qubit',x,'Source Qubit analytical quantification')
    x=op('library',x,'Sequencing-ready source-selected pooled DNA;full end/read geometry unassembled')
    end(g,x,'Associated MS1.15 NextSeq500/550 HighOutputv2.5 paired-end2×75;full indexed primer geometry/source sequence windows not assembled. Silent genome identity is not a per-molecule UMI.')
    return incomplete(g)
def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Production glycan addition directly follows acylation',g.steps['glycosylation'].inputs[0].material,'acylation_product')
    c('Binding library mixture uses production remainder,not MALDI-treated sample',g.steps['mix'].inputs[0].material,'maldi_production')
    c('MALDI aliquot remains independent QC',g.terminals['maldi_qc'].disposition,Disposition.QC)
    c('Production cleanup after glycosylation',g.steps['zeba'].inputs[0].material,'glycosylation_product')
    c('Titer aliquot preserved separately',g.terminals['titer_aliquot'].disposition,Disposition.STORED)
    c('Heating uses remainder not saved titer',g.steps['heat'].inputs[0].material,'heat_input')
    c('PCR uses physically collected supernatant',g.steps['pcr'].inputs[0].material,'heat_supernatant')
    c('Residual beads not converted into PCR fraction',g.terminals['heat_residue'].disposition,Disposition.UNRESOLVED)
    c('Source-limited version scope stays private',g.publishable,False)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'s41589-021-00788-5_PMC8380037.xml.txt'
    if p.exists():c('Selected original explicitly supplies collection and run',all(s in p.read_text() for s in ('25 μL of the supernatant','2×75 Cycles')),True)
    else:c.skip('Associated accepted manuscript','External original absent')
    c.report()
