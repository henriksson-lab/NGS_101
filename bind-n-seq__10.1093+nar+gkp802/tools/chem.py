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
        g.react(key,label,kind,inputs,(y,),changes_molecule=False);return y
    def clean(key,x,label):
        z=g.input(key+'_selection',label,K.SELECTION)
        y=g.product(key+'_kept','Physically retained source fraction: '+label,K.MOLECULAR)
        w=g.product(key+'_removed','Physically removed fractions: '+label,K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(z,R.SELECTION)),(y,w));g.discard(w,'Source physical workup only; no perfect specificity or purity assertion.');return y
    return op,clean

def end(g,x,label):g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=label))
def incomplete(g,branch):
    g.declare_reviewed_scope(COMMON_SCOPE+(RESIN_SCOPE if branch=='resin' else GEL_SCOPE),SOURCE)
    g.add_source_item(SourceItem('unassembled','Unreviewed linked originals and full molecular geometry',LIMIT));return g


from dataclasses import dataclass,replace
from protocol_graph import Terminal,Disposition

def join(g,key,inputs,label,kind=S.BINDING):
    y=g.product(key+'_product',label+';unremoved companion species remain',K.MOLECULAR)
    g.react(key,label,kind,tuple((x,R.TEMPLATE if i==0 else R.STRUCTURAL) for i,x in enumerate(inputs)),(y,),reason='Source-positive primer extension with named TaqPro Complete stock; complete primer/template and stock chemistry unassembled.' if kind is S.MAGIC else '');return y

def split(g,key,x,labels):
    out=[g.product(key+'_'+str(i+1),label,K.MOLECULAR) for i,label in enumerate(labels)]
    g.react(key,'Physical source allocation;amounts explicit only where inspected',S.SPLIT,((x,R.TEMPLATE),),tuple(out));return out

def render_page():
    return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_chem.html">Research notes</a></p>'+workflow_panel(workflow_from_protocol_graph(protocol_graph(),ENDPOINT))+'</main>'

TITLE='Bind-n-Seq — source DNA selection and sequencing'
SOURCE='Published PMC2794170 Methods25–59/Results63–77,Fig1originalPDFp3;published SI TableS1A/B'
LIMIT='Private incomplete graph;whole93mer andsourcePr4-f/Pr3-r primers/read sites unassembled. Selected resin/gel workups kept distinct,no predicted affinity/yield/ideal purity. Real-timePCR of DNA is quantification,not RNA reverse transcription. Optional second-round selection not inserted into baseline route.'
ENDPOINT='pool_product'
@dataclass(frozen=True)
class SourceTarget:
    barcode:str='AAC'
    def __post_init__(self):
        if len(self.barcode)!=3 or set(self.barcode)-set('ACGT'):raise ValueError('Source3ntbarcode required')
    def encoded_region(self):
        return Construct([Segment('source5primeAAleader','AA'),Segment('source condition barcode',self.barcode,feature=feature('bindnseq.condition.'+self.barcode.lower(),'sample_index','fixed',note='Reaction condition identity,not molecule or binding motif.')),Segment('source21randompositions','N'*21)])

def branch_graph(branch='resin'):
    if branch not in {'resin','gel'}:raise ValueError('Published resin or gel enrichment branch required')
    g=ProtocolGraph(TITLE+' ('+branch+')',authority=SOURCE);op,clean=operations(g)
    dna=g.input('ss_template','Source synthetic93mer with2ntAAleader/3ntconditionbarcode/N21/Illumina flanks;fullsequenceunassembled',K.OLIGO,primary=True)
    rp=g.input('extension_primer','Source reverse primer88µM;wholebodyunreviewed',K.OLIGO)
    supply=g.input('extension_supplies','Actual source1xTaqProComplete contributions including stated2mMMg2+; stockinternals/primer-end geometry unassembled',K.MOLECULAR)
    x=join(g,'second_strand',(dna,rp,supply),'SourceTaqProCompletestrandextension95°C2min63°C1min72°C4min;dsDNAcopies+originalssDNA/primer/reagents retained',kind=S.MAGIC)
    protein=g.input('binding_protein','Source purifiedMBP-zincfingerDBD;actualconcentration/saltcondition unassigned',K.MOLECULAR)
    x=join(g,'binding',(x,protein),'Source50µL protein/DNA equilibrium incubation23°C2h;herringDNA/Zn/DTT/BSA/KClconditions')
    if branch=='resin':
        resin=g.input('amylose','Source50µLamylose resin;two water/two washbufferprewashes',K.MOLECULAR)
        x=join(g,'resin_capture',(x,resin),'Source30min23°CMBP-amylosecapture,mixeach10min;unboundDNA/protein coexist')
        x=clean('wash',x,'Three10min sourceparameter-specific bufferwashes;keep resin-associated preparation')
        x=op('maltose',x,'Source50µL10mMmaltose10min elution;releasedMBPprotein/DNAclasses+resinremainuntilpellet','10mMmaltose/10mMTrispH8.5',reagent_kind=K.MOLECULAR,reagent_role=R.STRUCTURAL)
        y=g.product('eluate','Source collected supernatant;protein/DNA association may persist,carryover notperfect',K.MOLECULAR);r=g.product('resin_residue','Source pelletedresinandretainedmaterial',K.WASTE)
        g.react('pellet_collect','Source resinpellet/supernatanttransfer,thenrepeatpellet/transfer',S.CLEANUP,((x,R.TEMPLATE),),(y,r));g.discard(r,'Actual source resinremoval only;not quantitative separation.');x=y
    else:
        x=op('emsa',x,'Source5%TBEpolyacrylamidegel120V30min;FAMcontrol indicatescomplexmobility','FAM-labelledcontrolbindingreaction is separate marker,notlabeloneverylibrarymolecule')
        x=clean('band',x,'Exciseprotein–DNAcomplexband;unbound/othergelregions not baseline sequencing fraction')
        x=op('elute_gel',x,'Sourceovernight23°CgelDNAelution400µL ammoniumacetate/Mgacetate/0.1%SDS;protein/remnants notmagicallyremoved')
        x=clean('gel_liquid',x,'Transfer source liquid awayfromgel')
        x=clean('ethanol',x,'Source2volumesethanol−80°C20min,pellet/wash500µLethanol/dry/resuspend50µL')
    analytical,production=split(g,'qc_allocate',x,('Source1µLanalyticalreal-timePCRinput','Retained production preparation,analyticaltube notsequencinglibrary'))
    q=g.product('qpcr_result','SourcequantitativePCRdiagnostic;Pr4-f/Pr3-r/SYBR,notRT',K.QC);qa=g.input('qpcr_supplies','Actual source20ulanalyticalSYBRGreenPCRMasterKit/0.5uMPr4-f/Pr3-r contributions; wholeprimerorders and reagentformulation unassembled',K.MOLECULAR);g.react('qpcr','Sourceanalytical40-cyclethermalprogram+meltingcurve+standardcurve',S.MAGIC,((analytical,R.TEMPLATE),(qa,R.STRUCTURAL)),(q,),reason='Source actual real-time DNA PCR analytical assay; detailed stock/primer/input-product geometry unassembled. This QC tube is not a production library.');g.terminals[q.id]=Terminal(q.id,Disposition.QC,'AnalyticalPCR tube retainedasQC,notpooledassourcelibrary.')
    f=g.input('pcr_forward','SourcePr4-f wholebodyunreviewed',K.OLIGO);r=g.input('pcr_reverse','SourcePr3-r wholebodyunreviewed',K.OLIGO)
    normalized=g.product('normalized','Source normalized/amplified outputDNA forpooling;newcopies distinct fromoriginalprotein-associatedDNA,fullworkupunassembled',K.MOLECULAR)
    g.react('normalization','Sourceproductionnormalization/amplificationboundaryguidedbyquantitativePCR',S.MAGIC,((production,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL)),(normalized,),reason='Main specifies analyticalPCR andsource goal ofequalconcentration sequencingoutputs;exactproductioncycles/workupnotfullyassembled,noanalyticaltubeconvertedintolibrary.')
    x=op('pool',normalized,'SourceapproximatelyequalDNAconcentrationsfromdifferentbarcodedexperiments pooled')
    end(g,x,'PublishedIlluminaGenomeAnalyzerGAI(run1)/GAII(run2),version1reagents;actualrun-specific readcycles/primergeometryunassembled. N21 isbindingtarget,notUMI;conditionbarcode3nt.')
    return incomplete(g,branch)

# Independently curated source inventory: main Methods25-59, with actual opaque
# production/primer/linked-original limits still behind the unassembled item.
COMMON_SCOPE=tuple('ss_template extension_primer extension_supplies second_strand binding_protein binding qc_allocate qpcr_supplies qpcr pcr_forward pcr_reverse normalization pool pool_product'.split())
RESIN_SCOPE=tuple('amylose resin_capture wash_selection wash maltose_reagent maltose pellet_collect'.split())
GEL_SCOPE=tuple('emsa_reagent emsa band_selection band elute_gel gel_liquid_selection gel_liquid ethanol_selection ethanol'.split())

def protocol_graph():return branch_graph()
def selftest():
    from checks import Check
    c=Check();g=protocol_graph();gel=branch_graph('gel')
    c('DNA binding is not unspecified generic chemistry',g.steps['binding'].kind,S.BINDING)
    c('TaqPro source stock participates in primer extension',g.steps['second_strand'].inputs[-1].material,'extension_supplies')
    c('Physical elution/pooling preserve molecular identity',(g.steps['maltose'].changes_molecule,g.steps['pool'].changes_molecule),(False,False))
    c('Analytical PCR receives actual assay supplies',g.steps['qpcr'].inputs[-1].material,'qpcr_supplies')
    c('Source reaction barcode is structuredsampleidentity',SourceTarget().encoded_region().segments[1].feature.role,'sample_index')
    try:SourceTarget('AAAA');bad=False
    except ValueError:bad=True
    c('Sourceconditionbarcode rejectswrongwidth',bad,True)
    c('Maltoseelution precedes actual resin separation',g.steps['pellet_collect'].inputs[0].material,'maltose_product')
    c('Production excludes qPCR tube',g.steps['normalization'].inputs[0].material,'qc_allocate_2')
    c('Quantitative PCR retainedasQC',g.terminals['qpcr_result'].disposition,Disposition.QC)
    c('Gel fraction eluted before ethanol cleanup',gel.steps['ethanol'].inputs[0].material,'gel_liquid_kept')
    c('Baseline doesnotrepeatselection','second_selection' in g.steps,False)
    c('Sourcegeometry limits preventpublication',g.publishable,False);c.report()
