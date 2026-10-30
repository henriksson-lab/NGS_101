"""SLIC-CAGE carrier removal after a source-qualified nAnT-iCAGE boundary."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import illumina as il
from chemdraw import Construct,Segment,Scene
from protocol_graph import (ProtocolGraph,MaterialKind as K,InputRole as R,
    StepKind as S,SequencingPlan,SourceItem,Terminal,Disposition)
from pcr import PCRGeometry

SOURCE='PublishedPMC6280763 Methods169–201; publishedSupplementalMaterial TablesS1/2/5/6 andFigS20. TableS6 originalp8 rendered. ExactMurata2014precedessor internally unreviewed,not imported.'
FWD=il.P5[:20]
REV=il.P7[:21]
CARRIER_F='GCGGCAGCGTTCGCTATAAC'
SITES=('TAGGGATAACAGGGTAAT','TAACTATAACGGTCCTAAGGTAGCGAA')

@dataclass(frozen=True)
class CarrierInventory:
    origin:str='Manufactured carrier RNA, distinct from biological RNA of interest'
    cap_mix:str='Main nominal90percentuncapped/10percentcapped; TableS5approximately44ug+5ug,not exact molecule fractions'
    template_range_bp:tuple=(386,1034)
    motifs:tuple=SITES
    random_ends:str='Source PCR templates have5random5prime and6random3prime nucleotides;identifier function unassigned,not UMI by default'
    geometry:str='Full pJ241/template/IVT/cap/RT/library molecules unassembled; PCR-template bp are not asserted RNA transcript lengths'

@dataclass(frozen=True)
class SelectedCuts:
    """Explicit selected strand cuts; no enzyme offsets, recognition or yield defaults."""
    sequence:str
    positions:tuple
    inferred:bool=True
    def __post_init__(self):
        if not self.sequence or set(self.sequence)-set('ACGT') or tuple(sorted(set(self.positions)))!=self.positions or any(type(i)is not int or not 0<i<len(self.sequence) for i in self.positions):
            raise ValueError('ExplicitDNA and sorted unique internal cut coordinates required')
    @property
    def pieces(self):
        points=(0,)+self.positions+(len(self.sequence),)
        return tuple(self.sequence[a:b] for a,b in zip(points,points[1:]))

def oligo(name,seq):return Construct([Segment(name,seq)])
def unknown_scene(label):
    s=Scene();s.strand('population',[Segment('unassembled','['+label+']',placeholder=True)]);return s
def opaque(g,key,x,label):
    y=g.product(key+'_product',label+'; molecular inventory and companions remain unassembled',K.MOLECULAR)
    g.react(key,label,S.MAGIC,((x,R.TEMPLATE),),(y,),reason='Source-known continuation; complete source component chemistry/conditions/geometry not assembled.');return y
def cleanup(g,key,x,label):
    z=g.input(key+'_selection',label,K.SELECTION)
    kept=g.product(key+'_kept','Actual source-retained DNA preparation; purity/yield not predicted',K.MOLECULAR)
    waste=g.product(key+'_removed','Source physically separated fractions; carrier/target purity not assumed',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(z,R.SELECTION)),(kept,waste));g.discard(waste,'Actual source workup only; no perfect selectivity claim.');return kept
def digest(g,key,x):
    enzymes=g.input(key+'_enzymes','Source I-SceI5U/I-CeuI5U/CutSmart; genomic off-target sites and carrier copy distribution unassigned',K.CATALYST)
    y=g.product(key+'_mixture','Selected cut carrier/other site-bearing DNA fragments plus uncut/partly cut DNA, enzymes and buffers; no chemical deletion or perfect depletion',K.MOLECULAR)
    g.react(key,'Source carrier digestion37C3h',S.DIGESTION,((x,R.TEMPLATE),(enzymes,R.CATALYST)),(y,))
    y=opaque(g,key+'_heat',y,'Source65C20min enzyme heat-inactivation; fragments remain until beads')
    return y
def diagnostic(g,key,x,label):
    keep=g.product(key+'_production','Source retained production preparation independently from analytical aliquots',K.MOLECULAR)
    assay=g.product(key+'_assay',label,K.QC)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(keep,assay));g.mark_terminal(Terminal(assay.id,Disposition.QC,'Analytical material does not become the preparative DNA.'));return keep

def graph_for_cleanup(short_fragments=None):
    if short_fragments is not None and type(short_fragments)is not bool:raise ValueError('Explicit QC outcome or unassigned outcome required')
    g=ProtocolGraph('SLIC-CAGE — source10ngrecipe and conditional carrier depletion',authority=SOURCE,representation='published_instance')
    target=g.input('target_rna','Source illustrative10ng RNAofinterest; exact study specimen and transcript/cap composition unassigned',K.MOLECULAR,primary=True)
    carrier=g.input('carrier_rna','Source4990ng carrier mix for total5ug RNA; upstream PCR/T7/partialcapping/diagnostic preparation unassembled',K.MOLECULAR,chemical_state=CarrierInventory())
    x=g.product('mixed_rna','Physical target+manufactured carrier RNA, not a covalent assembly',K.MOLECULAR)
    g.react('mix','Source10ng+4990ng illustrative input recipe',S.POOL,((target,R.TEMPLATE),(carrier,R.TEMPLATE)),(x,))
    x=opaque(g,'cage',x,'Source Murata2014nAnT-iCAGE RT/cap-trapping/library steps; original not yet read, no recipe constructor imported. Child samples stay separate; no complete centrifugal drying')
    x=opaque(g,'cage_concentrate',x,'Source final nAnTAMPure product concentrated/adjusted15ul')
    x=diagnostic(g,'initial_bioanalyzer',x,'Source1ulHS DNA Bioanalyzer from15ul; nominal14ul remains')
    x=digest(g,'digest1',x)
    x=cleanup(g,'clean1',x,'Source1.8xAMPure;42ulwater elution thenconcentrate20ul')
    # Separate diagnostic reactions; never place all three primers in one tube.
    prod=g.product('pcr_input','Source18ul production DNA allocated to100ul amplification',K.MOLECULAR)
    sample=g.product('qpcr_input','Source analytical DNA allocation;1ul per qPCR, whole-library and carrier-selective assays; technical replicates/scaling unassigned',K.QC)
    g.react('qpcr_allocate','Source20ul preparation/qPCR/18ul PCR allocation; no invented replicate count',S.SPLIT,((x,R.TEMPLATE),),(prod,sample))
    assay_inputs=[g.product(k+'_input','Source separate1ul '+k+'qPCR input; technical replicate counts/scaling unassigned',K.QC) for k in ('whole','carrier')]
    g.react('diagnostic_allocate','Allocate separate whole-library and carrier-selective qPCR aliquots',S.SPLIT,((sample,R.TEMPLATE),),tuple(assay_inputs))
    for key,substrate,forward in zip(('whole','carrier'),assay_inputs,(FWD,CARRIER_F)):
        qf=g.input(key+'_f','PublishedTableS6 separate diagnostic forward primer',K.OLIGO,state=oligo('source forward',forward))
        qr=g.input(key+'_r','PublishedTableS6 separate diagnostic adapter_r1 aliquot',K.OLIGO,state=oligo('source reverse',REV))
        qe=g.input(key+'_kapa','Source KAPA SYBR FAST qPCR kit;TableS7 conditions',K.CATALYST)
        assay=g.product(key+'_assay','Source '+key+'qPCR products; actualCt not simulated',K.QC)
        g.react(key+'_qpcr','Separate source '+key+'qPCR assay; source1ul input/10ul reaction',S.PCR,((substrate,R.TEMPLATE),(qf,R.STRUCTURAL),(qr,R.STRUCTURAL),(qe,R.CATALYST)),(assay,),geometry=PCRGeometry.unavailable(unknown_scene('Diagnostic full-template geometry unassembled'),reason='Source diagnostic primer bodies known; full source library binding geometry unassembled.'))
        g.mark_terminal(Terminal(assay.id,Disposition.QC,'Diagnostic products not mixed into the sequencing preparation.'))
    f=g.input('pcr_f','Publishedadapter_f1 production primer, separate stock aliquot',K.OLIGO,state=oligo('adapter_f1',FWD));r=g.input('pcr_r','Publishedadapter_r1 production primer, separate stock aliquot',K.OLIGO,state=oligo('adapter_r1',REV))
    enz=g.input('kapa','Source KAPA HiFi HSReadyMix; cycle count selected from diagnosticCt, no fixed universalcycles',K.CATALYST)
    x=g.product('pcr_product','New DNA copies plus original cut/uncut templates and primers/reagents; no perfectly carrier-free population assumed',K.MOLECULAR)
    g.react('pcr','Source100ulPCR with18ulDNA; qPCR-selectedcycles',S.PCR,((prod,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(enz,R.CATALYST)),(x,),geometry=PCRGeometry.unavailable(unknown_scene('Full nAnT/copy/primer-landing geometry unassembled'),reason='Source primer bodies known; inherited library/full product geometry not assembled.'))
    x=cleanup(g,'pcr_clean',x,'Source1.8xAMPure;42ulelution thenconcentrate14ul')
    x=digest(g,'digest2',x)
    x=cleanup(g,'clean2',x,'Source stringent1xAMPure;42ulelution thenconcentrate12ul')
    x=diagnostic(g,'second_bioanalyzer',x,'SourceHS DNA Bioanalyzer aliquot determines whether shortfragments<300bp remain; aliquotamount unassigned')
    if short_fragments is True:
        x=cleanup(g,'optional_size',x,'SourcepositiveQC branch:adjust30ul/add0.8xAMPure;42ulelution thenconcentrate12ul')
        x=diagnostic(g,'repeat_bioanalyzer',x,'Source repeatBioanalyzer analytical aliquot after optional size selection')
    elif short_fragments is None:
        x=opaque(g,'optional_size',x,'Source conditional0.8xAMPure/repeatBioanalyzer if<300bp DNA present; actual QC outcome and intervening fractions unassigned')
    x=diagnostic(g,'final_qc',x,'Source final leftover-carrier qPCR and FigS20PicoGreen concentration measurements; actual aliquot counts/amounts unassigned')
    others=g.input('other_libraries','Other source independently prepared/QCchecked libraries; exact batchmembership/amounts unassigned',K.MOLECULAR)
    pool=g.product('sequencing_pool','Source sequencing library pool, no covalent joining of independent samples',K.LIBRARY)
    g.react('pool','Source child libraries pooled only before sequencing',S.POOL,((x,R.TEMPLATE),(others,R.TEMPLATE)),(pool,))
    g.sequence(pool,SequencingPlan(limitation='SourceMiSeqforS.cerevisiae/HiSeq2500forM.musculus;SE50. Exact source linker/index/read-primer geometry and actual batch QC outcomes unassembled. Carrier qPCR/read filtering is computational observation, not physical DNA removal.'))
    g.add_source_item(SourceItem('pending','Direct Murata2014nAnT-iCAGE original, full carrier/plasmid/cap/RT/library geometry','Publishedchild explicitly inherits named2014component with local pooling/drying exceptions; no whole recipe imported. Additional source records/code/results remain uninspected, not declared irrelevant or chemically absent.'))
    return g
def protocol_graph():return graph_for_cleanup()

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();cut=SelectedCuts('ACGTACGTACGT',(2,9))
    c('Selected digestion conserves all source bases',''.join(cut.pieces),cut.sequence)
    c('Reference not presented as source genomic prediction',cut.inferred,True)
    c.raises('Duplicate cut coordinates rejected',lambda:SelectedCuts('ACGT',(2,2)))
    c.raises('Terminal cuts cannot create nominal empty fragments',lambda:SelectedCuts('ACGT',(0,)))
    g=protocol_graph()
    c('Biological RNA distinct from independent carrier stock',{i.material for i in g.steps['mix'].inputs},{'target_rna','carrier_rna'})
    c('Production PCR excludes analytical qPCR DNA',g.steps['pcr'].inputs[0].material,'pcr_input')
    c('Whole and carrier qPCR have separate physical templates',(g.steps['whole_qpcr'].inputs[0].material,g.steps['carrier_qpcr'].inputs[0].material),('whole_input','carrier_input'))
    c('Carrier fragments remain through first heat treatment',g.steps['digest1_heat'].inputs[0].material,'digest1_mixture')
    c('Second digestion follows actual PCR cleanup',g.steps['digest2'].inputs[0].material,'pcr_clean_kept')
    c('Unknown QC outcome remains unassigned',g.steps['optional_size'].kind,S.MAGIC)
    c('Positive QC branch has actual size-selection fractions',graph_for_cleanup(True).steps['optional_size'].kind,S.CLEANUP)
    c('Negative QC branch does not invent size selection','optional_size' in graph_for_cleanup(False).steps,False)
    c('Unassembled predecessor geometry prevents promotion',bool(g.publication_problems()),True)
    src=Source(sources_root()/'slic-cage__10.1101+gr.235937.118'/'gr.235937.118_supp_supp_gr.235937.118_Supplemental_Material.pdf.txt','PublishedTableS6p8original rendered;sameprimersusedforproductionPCR')
    if have(c,src,label='Published primer inventory'):
        t=src.path.read_text();c('Whole three source primer bodies transcribed',all(x in t for x in (FWD,REV,CARRIER_F)),True)
        c('Both source rare-cutter motifs retained',all(x in t for x in SITES),True)
    c.report()
if __name__=='__main__':selftest()
