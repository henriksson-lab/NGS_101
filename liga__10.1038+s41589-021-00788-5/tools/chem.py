from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, feature, workflow_from_protocol_graph, workflow_panel
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan
from page import head,info
from pcr import PCRGeometry
PAGE_FILE='chem.html'

def operations(g):
    def op(key,x,label,reagent=None,state=None,kind=S.MAGIC,reagent_kind=K.MOLECULAR,reagent_role=R.STRUCTURAL):
        inputs=[(x,R.TEMPLATE)]
        if reagent:
            z=g.input(key+'_reagent',reagent,reagent_kind);inputs.append((z,reagent_role))
        y=g.product(key+'_product',label+'; preparation includes unremoved reaction components',K.MOLECULAR,state=state)
        kw={'reason':'Actualsourceoperation; occupiedchemicalstates/internalgeometryunassembled, nofamilyrecipe orquantitativeyield claim.'} if kind==S.MAGIC else {}
        if kind==S.GENERIC:kw['changes_molecule']=False
        if kind==S.CHEMICAL_CONVERSION:
            cond=g.input(key+'_conditions',label+'; actualsourceconditions, notextra reagentaddition',K.CONDITION);inputs.append((cond,R.CONDITION))
        g.react(key,label,kind,inputs,(y,),**kw);return y
    def clean(key,x,label):
        z=g.input(key+'_selection',label,K.SELECTION)
        y=g.product(key+'_kept','Physically retained source fraction: '+label,K.MOLECULAR)
        w=g.product(key+'_removed','Physically removed fractions: '+label,K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(z,R.SELECTION)),(y,w));g.discard(w,'Source physical workup only; no perfect specificity or purity assertion.');return y
    return op,clean

def end(g,x,label):g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=label))
def incomplete(g):
    g.declare_reviewed_scope(('phage','acylation_reagent','acylation_conditions','acylation','glycosylation_reagent','glycosylation_conditions','glycosylation','zeba_selection','zeba','maldi_split','maldi','mix_allocate','mix_other','mix','bind_reagent','bind','wash_selection','wash','water_reagent','water','titer_split','heat_reagent','heat','supernatant','pcr_supplies','pcr_components','pcr_forward','pcr_reverse','pcr_enzyme','pcr_nt','pcr_conditions','pcr','agarose_split','agarose','pool_allocate','pool_other','pool','egel_selection','egel','qubit_split','qubit','loading_allocate','library'),SOURCE+'; fixedactualsourceproduction/capture/NGSinterface scope, not fullgenome/primer/linked-original review.')
    g.add_source_item(SourceItem('unassembled','Unreviewed linked originals and full molecular geometry',LIMIT));return g


from protocol_graph import Terminal,Disposition

def save(g,key,x,label):
    take=g.product(key+'_aliquot',label,K.MOLECULAR)
    rest=g.product(key+'_remainder','Remaining preparation; later use unassigned',K.MOLECULAR)
    g.react(key,'Physical aliquot: '+label,S.SPLIT,((x,R.TEMPLATE),),(take,rest))
    g.terminals[rest.id]=Terminal(rest.id,Disposition.UNRESOLVED,'Remainder retained explicitly; later fate unassigned.')
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
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance');op,clean=operations(g)
    def selected(k,x,label):
        take=g.product(k+'_selected',label,K.MOLECULAR);rest=g.product(k+'_remainder','Untransferredactualstock; amount/fateunassigned',K.MOLECULAR)
        g.react(k+'_allocate',label,S.SPLIT,((x,R.TEMPLATE),),(take,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Untransferredsource stock.'));return take
    def pool(k,x,label,ownlabel,otherlabel):
        own=selected(k,x,ownlabel);other=g.input(k+'_other',otherlabel,K.MOLECULAR);out=g.product(k+'_product',label,K.MOLECULAR)
        g.react(k,label,S.POOL,((own,R.TEMPLATE),(other,R.TEMPLATE)),(out,));return out
    x=g.input('phage','Source preconstructed silently DNA-barcoded M13 phage clones;genome sequence and manufacture unassembled',K.BIOLOGICAL,primary=True)
    x=op('acylation',x,'Associated MS production1.8:DBCO-NHS0.1–1.5mM,45min incubation;phage modification distribution not computed','ActualDBCO-NHS inDMF/PBS;hydrolyzed/unreactedreagent mayremain',kind=S.CHEMICAL_CONVERSION)
    x=op('glycosylation',x,'Associated MS production1.8:add glycan-azide directly to2–3mM;overnight4°C','Actualglycan-azide;sourceSPAACproductensemble notguaranteeduniform',kind=S.CHEMICAL_CONVERSION)
    x=clean('zeba',x,'Source production Zeba Spin purification;retain glycan-phage fraction')
    x=measure(g,'maldi',x,'Source MALDI-TOF analysis;separate2µL analytical sample exposed to matrix/organic solvent')
    x=pool('mix',x,'Actualsource desired glycan-phage clonepool; dictionaryassociation notphysicaltransformation','Actual10^8PFU ownclonecontribution; untransferredsourceclone stockretained','Actual10^8PFUeachotherdesired glycan-phageclonecontributions; membershipunexpanded, notduplicatesofownclone')
    x=op('bind',x,'Selected MS1.11.3 CD22-Fc/ProteinG magnetic-bead binding;other plate/cell/lectin recipes not merged','Actual preparedCD22-FcProteinGbeads; manufacturerpreparationunreviewed',reagent_kind=K.SELECTION,reagent_role=R.SELECTION)
    x=clean('wash',x,'KingFisher2HEPES+1%Tween washes then1HEPES wash;retain bead-associated preparation')
    x=op('water',x,'Transfer beads to30µL nuclease-free water','Actualsourcewaterresuspension',kind=S.GENERIC,reagent_kind=K.CONDITION,reagent_role=R.CONDITION)
    t=g.product('titer_aliquot','2µL suspension for titering;not the PCR aliquot',K.MOLECULAR);rest=g.product('heat_input','Remaining suspension after titer aliquot',K.MOLECULAR)
    g.react('titer_split','Remove2µL titer aliquot from30µL',S.SPLIT,((x,R.TEMPLATE),),(t,rest));g.terminals[t.id]=Terminal(t.id,Disposition.UNRESOLVED,'Source titer aliquot retainedforanalysis; no inventedstorage/NGS use.')
    x=op('heat',rest,'Remaining preparation90°C15min;DNA/particle remnants and beads not magically removed','Actualsource90C15min',kind=S.THERMAL_TREATMENT,reagent_kind=K.CONDITION,reagent_role=R.CONDITION)
    y=g.product('heat_supernatant','Source25µL supernatant PCR input',K.MOLECULAR);rem=g.product('heat_residue','Residual liquid/beads/debris;not source PCR fraction',K.MOLECULAR)
    g.react('supernatant','Source physically collect25µL supernatant',S.CLEANUP,((x,R.TEMPLATE),),(y,rem));g.terminals[rem.id]=Terminal(rem.id,Disposition.UNRESOLVED,'Residual fraction not automatically destroyed or sequenced.')
    supplies=g.input('pcr_supplies','Actual source50ulPCR supplies; individualstockpackaging/thermalprogramunassembled',K.MOLECULAR)
    f=g.product('pcr_forward','Actual1uMforwardindexedprimerclass; wholebodyunassembled',K.OLIGO);r=g.product('pcr_reverse','Actual1uMreverseindexedprimerclass; wholebodyunassembled',K.OLIGO)
    cat=g.product('pcr_enzyme','Actualsource1UPhusion polymerase component',K.CATALYST);nt=g.product('pcr_nt','Actualsource50uMeachdNTPsubstratecomponent',K.MOLECULAR);buf=g.product('pcr_conditions','Actual1xPhusionbuffer/500uMMgCl2; exactsourcecyclingdeferred',K.CONDITION)
    g.react('pcr_components','ActualPCRsupplyfunctions, notpremixclaimorextra freshstockdoses',S.GENERIC,((supplies,R.TEMPLATE),),(f,r,cat,nt,buf),changes_molecule=False)
    amp=g.product('pcr_product','ActualPCRDNAcopies plusoriginalphagegenome/remainingreagents, notcopiedglycan/protein',K.MOLECULAR)
    g.react('pcr','AMMS1.15 actual50ulindexedPhusionPCR',S.PCR,((y,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST),(nt,R.STRUCTURAL),(buf,R.CONDITION)),(amp,),geometry=PCRGeometry.unavailable(Scene(),reason='ActualsourcePCRproducer; exactsourceprimerbodies/fullphagewindow/thermalprogramunassembled, no linkagefromprotein/glycantosequence invented.'));x=amp
    x=measure(g,'agarose',x,'Source2%agarose analytical quantification;not the preparative eGel step')
    x=pool('pool',x,'ActualdifferentindexedPCRproductpool','Actual10ngown indexedPCRcontribution','Actual10ngeachotherindexedPCRproductcontributions; membership/sourcepreparationsunexpanded')
    x=clean('egel',x,'Source eGel purification;retain selected DNA')
    x=measure(g,'qubit',x,'Source Qubit analytical quantification')
    x=selected('loading',x,'ActualselectedpooledPCRDNAloadingcontribution; volumeunassigned')
    x=op('library',x,'Sequencing-ready source-selected pooled DNA;full end/read geometry unassembled',kind=S.GENERIC)
    end(g,x,'Associated MS1.15 NextSeq500/550 HighOutputv2.5 paired-end2×75;full indexed primer geometry/source sequence windows not assembled. Silent genome identity is not a per-molecule UMI.')
    return incomplete(g)
def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Production glycan addition directly follows acylation',g.steps['glycosylation'].inputs[0].material,'acylation_product')
    c('Binding library mixture uses production remainder,not MALDI-treated sample',g.steps['mix_allocate'].inputs[0].material,'maldi_production')
    c('MALDI aliquot remains independent QC',g.terminals['maldi_qc'].disposition,Disposition.QC)
    c('Production cleanup after glycosylation',g.steps['zeba'].inputs[0].material,'glycosylation_product')
    c('Titer aliquot preserved separately',g.terminals['titer_aliquot'].disposition,Disposition.UNRESOLVED)
    c('Heating uses remainder not saved titer',g.steps['heat'].inputs[0].material,'heat_input')
    c('PCR uses physically collected supernatant',g.steps['pcr'].inputs[0].material,'heat_supernatant')
    c('Residual beads not converted into PCR fraction',g.terminals['heat_residue'].disposition,Disposition.UNRESOLVED)
    c('PCR is actualDNA-producingstep',g.steps['pcr'].kind,S.PCR)
    c('Ownindexedpoolcontribution selectedfromproductionstock',g.steps['pool_allocate'].inputs[0].material,'agarose_production')
    c('Actualindexedpoolhasdistinctothercontributions',len(g.steps['pool'].inputs),2)
    c('PCRnucleotides come fromsamesourceinventory',g.steps['pcr'].inputs[4].material,'pcr_nt')
    c('Loadingretainsuntransferredlibrary',g.terminals['loading_remainder'].disposition,Disposition.UNRESOLVED)
    c('Source-limited version scope stays private',g.publishable,False)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'s41589-021-00788-5_PMC8380037.xml.txt'
    if p.exists():c('Selected original explicitly supplies collection and run',all(s in p.read_text() for s in ('25 μL of the supernatant','2×75 Cycles')),True)
    else:c.skip('Associated accepted manuscript','External original absent')
    c.report()
