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
    y=g.product(key+'_product',label+'; unremoved companion species remain',K.MOLECULAR)
    g.react(key,label,S.GENERIC,tuple((x,R.TEMPLATE if i==0 else R.STRUCTURAL) for i,x in enumerate(inputs)),(y,));return y

def render_page():
    return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_chem.html">Research notes</a></p>'+workflow_panel(workflow_from_protocol_graph(protocol_graph(),ENDPOINT))+'</main>'

TITLE='SEQRS — associated-manuscript iterative RNA selection'
SOURCE='Associated accepted manuscript PMC3375920 Approach58–60,Methods132–142,Fig1caption384–390; actual chemistry supplement gated'
LIMIT='Private incomplete source graph. Associated manuscript operations explicitly version-qualified; full published equivalence and TableS1 polymers/read primer remain unreviewed. Random20nt RNA is binding analyte, six-base sample index identifies selected libraries. No predicted affinity, purity or successful retention of every RNA.'
ENDPOINT='pool_product'
@dataclass(frozen=True)
class RoundSample:
    completed:int=0
    polymer:str='RNA'
    def __post_init__(self):
        if type(self.completed)is not int or self.completed<0 or self.polymer not in {'RNA','DNA'}:raise ValueError('Explicit source stage and polymer required')
    def rt_copy(self):
        if self.polymer!='RNA':raise ValueError('RT copies RNA input')
        return replace(self,polymer='DNA')
    def retranscribe(self):
        if self.polymer!='DNA':raise ValueError('T7 transcription needs DNA template')
        return replace(self,completed=self.completed+1,polymer='RNA')

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE);op,clean=operations(g)
    x=g.input('dna_template','Source dsDNA/T7 template encoding N20 between source constant regions; full promoter/strand geometry unreviewed',K.OLIGO,primary=True)
    libraries=[]
    for cycle in range(1,6):
        k='round'+str(cycle)
        x=op(k+'_ivt',x,'Source T7 transcription; initial input1µg dsDNA, later template allocation/workup unassembled','AmpliScribe T7-Flash; new RNA plus original DNA/reagents retained')
        if cycle==1:
            x=op(k+'_dnase',x,'Source initial-library RNase-free DNase treatment of residual template DNA; digestion is not physical removal','Source RNase-free DNase',reagent_kind=K.CATALYST,reagent_role=R.CATALYST)
            x=clean(k+'_rna_purify',x,'Source initial-library GeneJET RNA purification; physically retain RNA preparation')
        else:
            prepared=g.product(k+'_rna_prepared','Source regenerated RNA preparation for next binding; DNA/template/enzyme/remnant inventory unassembled',K.MOLECULAR)
            g.react(k+'_regeneration_workup','Unreviewed later-cycle preparation/workup after T7 regeneration',S.MAGIC,((x,R.TEMPLATE),),(prepared,),reason='Approach/Fig1 establish iterative transcription and selection, but Methods specifies DNase/GeneJET for initial input only. No identical later-round cleanup or perfect template removal inferred.');x=prepared
        protein=g.input(k+'_protein','Source prepared GST-fusion RBP on magnetic GST-agarose; selected source protein identity/occupancy unassigned',K.MOLECULAR)
        x=join(g,k+'_bind',(x,protein),'Source150ng RNA+50–100nmol fusion protein in100µL buffer/tRNA competitor/inhibitor;30min ambient')
        x=clean(k+'_wash',x,'Source magnetic block/aspirate reaction; four200µL ice-cold washes, retain resin-associated population')
        primer=g.input(k+'_rt_primer','Source TableS1 RT primer; whole polymer unreviewed',K.OLIGO)
        x=join(g,k+'_heat',(x,primer),'Resuspend resin in1mMTrispH8+10pmolRTprimer;65°C10min then ice; resin/protein/RNA not automatically separated')
        aliquot=g.product(k+'_rt_input','Source5µL of heated/cooled sample for10µL ImProm-II RT; possible resin companions unassembled',K.MOLECULAR)
        rest=g.product(k+'_elution_remainder','Remaining source sample/resin after RT allocation; later fate unassigned',K.MOLECULAR)
        g.react(k+'_rt_allocate','Physical5µL source sample transfer to RT, not invented clear-supernatant collection',S.SPLIT,((x,R.TEMPLATE),),(aliquot,rest));g.terminals[rest.id]=Terminal(rest.id,Disposition.UNRESOLVED,'Actual remaining preparation not forced into RT or discarded.')
        x=op(k+'_rt',aliquot,'Source ImProm-II makes new ssDNA copies; original RNA/primer/protein companions remain until later workup','ImProm-II RT system; exact whole product geometry unassembled')
        f=g.input(k+'_pcr_forward','Source PCR primer regenerating T7 promoter; full TableS1 body unreviewed',K.OLIGO)
        r=g.input(k+'_pcr_reverse','Source PCR reverse primer; full body unreviewed',K.OLIGO)
        x=join(g,k+'_pcr',(x,f,r),'Source PCR generates dsDNA templates for next selection cycle; copies distinct from original RNA/protein')
        sample=g.product(k+'_library_input','Source dsDNA aliquot for this round library; amount unreviewed',K.MOLECULAR)
        onward=g.product(k+'_continuation','Remaining source dsDNA preparation; subsequent template transfer/workup unassembled',K.MOLECULAR)
        g.react(k+'_sample','Source round sampling for sequencing versus continued selection',S.SPLIT,((x,R.TEMPLATE),),(sample,onward))
        af=g.input(k+'_adapter_forward','Source sequencing-adapter/index primer; exact6nt member/body unreviewed',K.OLIGO)
        ar=g.input(k+'_adapter_reverse','Source sequencing-adapter primer; exact body unreviewed',K.OLIGO)
        lib=join(g,k+'_library_pcr',(sample,af,ar),'Source sampled-round PCR adds Illumina adapters and unique6bp condition/round index, not a UMI')
        diagnostic=g.product(k+'_qc','Separate source electrophoresis purity-assay material; actual aliquot amount unreviewed',K.QC)
        kept=g.product(k+'_production','Source production library retained separately from analytical assays',K.MOLECULAR)
        g.react(k+'_qc_allocate','Separate analytical electrophoresis purity assessment',S.SPLIT,((lib,R.TEMPLATE),),(kept,diagnostic));g.terminals[diagnostic.id]=Terminal(diagnostic.id,Disposition.QC,'Analytical material not automatically pooled as sequencing substrate.')
        lib=clean(k+'_library_clean',kept,'Source individual WizardSV Gel/PCRClean-Up purification')
        quantified=g.product(k+'_quantified_production','Source purified library retained separately from PicoGreen aliquot',K.MOLECULAR);quantity=g.product(k+'_picogreen','Separate source PicoGreen fluorescence quantification aliquot',K.QC)
        g.react(k+'_quantify','Source concentration assessment after individual purification',S.SPLIT,((lib,R.TEMPLATE),),(quantified,quantity));g.terminals[quantity.id]=Terminal(quantity.id,Disposition.QC,'PicoGreen assay material not automatically pooled for sequencing.')
        libraries.append(quantified);x=onward
    g.terminals[x.id]=Terminal(x.id,Disposition.UNRESOLVED,'Remaining DNA after selected typical five-cycle workflow; further selection not assumed.')
    allocated=g.product('pool_inputs','Source barcoded production-library fractions assigned to this batch; exact round/sample membership and allocation unassembled',K.MOLECULAR)
    allocated_other=g.product('pool_inputs_other','Other allocated barcoded source-library fractions; count and round/sample identities unassigned',K.MOLECULAR)
    unallocated=g.product('unallocated_libraries','Other prepared round libraries/remainders outside this source pool; exact identities/amounts unassigned',K.MOLECULAR)
    g.react('batch_allocation','Unreviewed source round/sample-to-sequencing-batch allocation',S.MAGIC,tuple((x,R.TEMPLATE) for x in libraries),(allocated,allocated_other,unallocated),reason='Source samples multiple rounds and pools barcoded DNA, but selected main does not assign exact lane/batch membership. Do not force every prepared round into one pool.');g.terminals[unallocated.id]=Terminal(unallocated.id,Disposition.UNRESOLVED,'Unallocated library fractions retained independently; no forced discard or acquisition.')
    pooled=g.product('pool_product','Approximately equal allocated barcoded DNA libraries; actual round/sample membership unassigned',K.MOLECULAR)
    g.react('pool','Source concentration-guided multiplex pooling after unresolved batch allocation',S.POOL,((allocated,R.TEMPLATE),(allocated_other,R.TEMPLATE)),(pooled,))
    end(g,pooled,'Associated source3pmol pool onIlluminaHiSeq2000 with custom primer; exact read cycles, index members/primer landing unreviewed. Random binding RNA sequences are not molecular identifiers.')
    return incomplete(g)
def selftest():
    from checks import Check
    c=Check();g=protocol_graph();r=RoundSample().rt_copy().retranscribe()
    c('Iteration creates RNA from a new DNA intermediate',(r.polymer,r.completed),('RNA',1))
    try:RoundSample().retranscribe();bad=False
    except ValueError:bad=True
    c('Cannot treat native RNA as a T7 DNA template',bad,True)
    c('DNase remnants remain until physical purification',g.steps['round1_rna_purify'].inputs[0].material,'round1_dnase_product')
    c('RT receives actual transferred sample',g.steps['round1_rt'].inputs[0].material,'round1_rt_input')
    c('Later-round workup does not inherit initial DNase/GeneJET','round2_dnase' in g.steps or 'round2_rna_purify' in g.steps,False)
    c('Regenerated RNA reaches binding through explicit unknown workup',g.steps['round2_bind'].inputs[0].material,'round2_rna_prepared')
    c('Remaining elution sample has unresolved later fate',g.terminals['round1_elution_remainder'].disposition,Disposition.UNRESOLVED)
    c('Next round starts with continued PCR preparation',g.steps['round2_ivt'].inputs[0].material,'round1_continuation')
    c('Each pooled branch uses retained production sample',g.steps['round1_library_clean'].inputs[0].material,'round1_production')
    c('Quantification follows source individual cleanup',g.steps['round1_quantify'].inputs[0].material,'round1_library_clean_kept')
    c('Batch allocation excludes PicoGreen assay material',g.steps['batch_allocation'].inputs[0].material,'round1_quantified_production')
    c('Pool receives allocated fractions,not forced every round',g.steps['pool'].inputs[0].material,'pool_inputs')
    c('Final continuation not silently discarded',g.terminals['round5_continuation'].disposition,Disposition.UNRESOLVED)
    c('Unreviewed whole geometry prevents promotion',g.publishable,False);c.report()
