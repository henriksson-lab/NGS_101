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


from dataclasses import dataclass,replace
from protocol_graph import Terminal,Disposition

def join(g,key,inputs,label):
    y=g.product(key+'_product',label+';unremoved companion species remain',K.MOLECULAR)
    g.react(key,label,S.GENERIC,tuple((x,R.TEMPLATE if i==0 else R.STRUCTURAL) for i,x in enumerate(inputs)),(y,));return y

def render_page():
    return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_chem.html">Research notes</a></p>'+workflow_panel(workflow_from_protocol_graph(protocol_graph(),ENDPOINT))+'</main>'

TITLE='NGGM — associated-manuscript direct glycan–DNA preparation'
SOURCE='Associated accepted manuscript PMC13016201 Methods103–126;publishedSI unavailable/unreviewed'
LIMIT='Private incomplete source graph. Selected associated author-manuscript operations only;whole final-publication equivalence unverified. Actual sequences/index/adapter geometry deferred to sourceSI. DNA reporters measure selected glycan ligand classes,not cell identity. No perfect conjugation,binding or fraction recovery predicted.'
ENDPOINT='library_kept'
@dataclass(frozen=True)
class CodeGlycan:
    code_id:str
    stage:str='azide-DNA'
    def __post_init__(self):
        if not self.code_id or self.stage not in {'azide-DNA','glycan-DBCO','glycan-DNA','quenched'}:raise ValueError('Explicit source code identity/stage needed')
    def click(self):
        if self.stage!='azide-DNA':raise ValueError('Azide DNA needed')
        return replace(self,stage='glycan-DNA')
    def quench(self):
        if self.stage!='glycan-DNA':raise ValueError('Selected conjugated class needed')
        return replace(self,stage='quenched')

def branch_graph(branch='lectin'):
    if branch not in {'lectin','antibody','bacteria'}:raise ValueError('Source branch required')
    g=ProtocolGraph(TITLE+' ('+branch+')',authority=SOURCE);op,clean=operations(g)
    gly=g.input('glycan','Source defined AEAB-glycan mixture;source structures/code dictionary not instantiated',K.MOLECULAR,primary=True)
    gly=op('dbco',gly,'Defined library:AEABglycan1eq+DBCO-PEG4-NHS2eq,phosphatepH8.5,RT2h','DBCO-PEG4-NHS;coupling ensemble/yield uncomputed')
    gly=clean('hplc',gly,'Source reversed-phase HPLC glycan-PEG4-DBCO purification/quantification')
    a=g.input('azide_oligo','Source5′azideDNA,18nt primer/20nt code/18nt primer;whole sequence unreviewed',K.OLIGO)
    b=g.input('complement','Source unmodified complementary oligo;whole sequence unreviewed',K.OLIGO)
    dna=join(g,'anneal',(a,b),'Equal molar strands95°C5min,cool0.1°C/s to10°C;source azide dsDNA')
    x=join(g,'click',(gly,dna),'Defined library:1µL0.1mMdsDNA+1µL0.5mMglycan-PEG4-DBCO,RT3h')
    x=op('quench',x,'Add2µL0.1Mazidoethanol;quench remaining click-reactive reagent,store−20°C','Azidoethanol;unreacted compounds/reaction products remain')
    if branch=='bacteria':
        cells=g.input('bacteria','Source400µL cultureOD~1',K.BIOLOGICAL)
        cells=clean('cell_pellet',cells,'Source3000g1min;retain bacteria pellet')
        x=join(g,'bind',(x,cells),'Resuspend/incubate bacteria+glycanDNA25fmol each in100µLTSM/BSA/Triton/salmonDNA1h')
        x=clean('binding_wash',x,'Centrifuge/three100µLTSM washes;keep cell pellet')
        x=op('resuspend',x,'Resuspend washed cells in100µLPBS;no boil or clarified-supernatant step inherited')
        fraction='Source1µL cell suspension forPCR'
    else:
        protein=g.input('binding_protein','Source biotin lectin1µg' if branch=='lectin' else 'Source blood-type antibody1µL',K.MOLECULAR)
        x=join(g,'bind',(x,protein),'Incubate glycanDNA1.25fmol each+bindingprotein in10µLTSM/BSA/Triton/salmonDNA2.5h')
        if branch=='antibody':
            beads=g.input('beads','Source SA magneticbeads preincubated with biotin antimouseIgM2h;not direct antibody biotin inference',K.MOLECULAR)
        else:beads=g.input('beads','Source1µL1mg/mL SA magneticbeads',K.MOLECULAR)
        x=join(g,'capture',(x,beads),'Source30min bead capture after solution-phase ligand binding')
        x=clean('binding_wash',x,'Magnetic separation/three10µLTSM washes;keep beads')
        x=op('boil',x,'Resuspend beads10µLwater;95°C10min;beads/reaction material retained until supernatant sampling')
        fraction='Source1µL supernatant forPCR;residual beads/liquid independent'
    chosen=g.product('pcr_aliquot',fraction,K.MOLECULAR);rest=g.product('remainder','Remaining source suspension or beads/liquid;laterfate unassigned',K.MOLECULAR)
    g.react('allocate',fraction,S.SPLIT,((x,R.TEMPLATE),),(chosen,rest));g.terminals[rest.id]=Terminal(rest.id,Disposition.STORED,'Source residual retained,not automatically destroyed or selected.')
    f=g.input('index_forward','Source firstPCR indexedforward primer;whole sequence unreviewed',K.OLIGO);r=g.input('index_reverse','Source firstPCR indexedreverse primer;whole sequence unreviewed',K.OLIGO)
    x=join(g,'pcr1',(chosen,f,r),'Source20–30cyclePCR;two sample-index sequences,newDNAcopies distinct from conjugate template')
    x=op('pool',x,'Mix comparable amounts of different first-roundPCR products;membership unassigned')
    x=clean('purify1',x,'PureLinkQuickPCR purification')
    analytical=g.product('quantification_aliquot','Source analytical quantification aliquot;amount unassigned',K.MOLECULAR);production=g.product('quantified_production','ProductionDNA retained separately from analytical aliquot',K.MOLECULAR)
    g.react('quantify_split','Separate source quantification aliquot;amount unreviewed',S.SPLIT,((x,R.TEMPLATE),),(analytical,production))
    qc=g.product('quantification_qc','Analytically assayed DNA fraction;method unreviewed',K.QC);g.react('quantify','Source analytical quantification',S.GENERIC,((analytical,R.TEMPLATE),),(qc,));g.terminals[qc.id]=Terminal(qc.id,Disposition.QC,'Assayed fraction not automatically amplified as production library.');x=production
    x=op('pcr2',x,'Source printed10µg mixture,10cyclePCR addingNGSadapters;literal quantity retained','Source adapter primers/polymerase mixture;exact oligo/end geometry unassembled')
    x=clean('library',x,'Source secondPCR product purification beforeAdmeraHealth submission')
    end(g,x,'SourceNGS service submission explicit;actualinstrument/cycles/primers/readlayout unreviewed. Barcode counts computational,not physical binding amount guarantee.')
    return incomplete(g)
def protocol_graph():return branch_graph()
def selftest():
    from checks import Check
    c=Check();a=branch_graph();b=branch_graph('bacteria');ab=branch_graph('antibody')
    c('Lectin capture follows solution binding',a.steps['capture'].inputs[0].material,'bind_product')
    c('BacterialPCR selects actual cell suspension','suspension' in b.materials['pcr_aliquot'].label,True)
    c('Bacterial branch does not inherit boiling','boil' in b.steps,False)
    c('Antibody capture has source secondary-antibody beads','antimouseIgM' in ab.materials['beads'].label,True)
    c('Adapter PCR uses retained productionDNA',a.steps['pcr2'].inputs[0].material,'quantified_production')
    c('Saved source remainder remains separate',a.terminals['remainder'].disposition,Disposition.STORED)
    n=CodeGlycan('selected_source_code');c('Conjugation retains code identity',n.click().quench().code_id,n.code_id)
    try:n.quench();bad=False
    except ValueError:bad=True
    c('Cannot claim quenched conjugate before declared coupling',bad,True)
    c('Incomplete final source geometry prevents publication',a.publishable,False);c.report()
