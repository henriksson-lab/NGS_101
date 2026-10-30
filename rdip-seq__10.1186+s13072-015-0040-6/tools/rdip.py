"""Published RDIP source scope; native hybrid and directional-copy limits stay local."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
from chemdraw import Scene, Segment, feature
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SourceClaim, SourceDiscrepancy, SequencingPlan, Terminal, Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from chemdraw import workflow_from_protocol_graph, workflow_panel
from page import head, info

PAGE_FILE='rdip.html'
SOURCE='Nadel et al.2015,10.1186/s13072-015-0040-6; published main Methods202–216; SI1pp19–24; SI2TableS4'
SAMPLE=feature('rdip.sample','sample_index','whitelist',whitelist='Published SI2 TableS4',note='Published TableS4 prints six fixed bases plus a parenthesized seventh. Physical placement/acquired cycle count unassigned.')
INDEX_ROWS=(('IMR-90 Input',5,'ACAGTG(A)'),('HEK 293T Input',6,'GCCAAT(A)'),('HEK 293T RDIP',12,'CTTGTA(A)'),('IMR-90 RDIP',19,'GTGAAA(C)'))
LIMIT='Private incomplete model. Original hybrid strands, RNaseH-primed replacement, exact historical adapter/primer junctions and acquired index geometry are unassembled. Published reagent/table inconsistencies remain explicit; no ordinary RNA-seq RT or silently corrected recipe is substituted.'


def unknown(label):
    sc=Scene();sc.strand('population',[Segment(label,'[source population; positions and ends unassigned]',placeholder=True)]);return sc


@dataclass(frozen=True)
class UracilFates:
    """Qualitative selected processing class, not a yield or all-molecule prediction."""
    processed: bool=False
    heated: bool=False
    original_dna: bool=True
    rna_remnants: bool=True

    def rows(self):
        sc=Scene()
        if self.original_dna:sc.strand('original',[Segment('original genomic DNA','[original DNA; covalent continuity unassembled]',placeholder=True)])
        label='dUTP-marked copied DNA; positions/extent unassembled' if not self.processed else 'copied DNA with processed AP sites and possible residual uracil; no whole-strand disappearance'
        if self.heated:label+='; heat-derived strand/fragment distribution unassigned'
        sc.strand('copy',[Segment(label,'[copied DNA/backbone class]',placeholder=True)])
        if self.rna_remnants:sc.strand('rna',[Segment('source RNA/cleavage remnants; recovery unknown','[RNA/remnants]',placeholder=True)])
        if self.processed:sc.strand('released',[Segment('uracil nucleobases released by UNG; still in reaction','[free uracil]',placeholder=True)])
        return sc.rows()


def selected_ung_event(state):
    if state.processed:raise ValueError('Selected UNG event already recorded')
    return replace(state,processed=True)


def heat_treatment(state):
    return replace(state,heated=True)


def protocol_graph():
    g=ProtocolGraph('RDIP-seq — published directional hybrid library',authority=SOURCE);scope=[]
    def inp(key,label,kind=K.CATALYST,state=None,primary=False):return g.input(key,label,kind,state=state,primary=primary)
    def op(key,x,label,kind,reagents=(),state=None,geometry=None,reason=None):
        y=g.product(key+'_out',label+' preparation',K.MOLECULAR,state=state or unknown(label))
        g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(reagents),(y,),geometry=geometry,**({'reason':reason} if reason else {}));scope.append(key);return y
    def clean(key,x,label):
        reagent=inp(key+'_reagents','Source extraction/cleanup reagents',K.SELECTION)
        keep=g.product(key+'_kept',label+' retained preparation',K.MOLECULAR,state=unknown('Retained nucleic-acid preparation; class purity not assumed'))
        waste=g.product(key+'_waste','Actually nonretained extraction/wash fractions; composition unresolved',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(reagent,R.SELECTION)),(keep,waste));g.discard(waste,'Actual source phase transfer, wash or column separation.');scope.append(key);return keep
    cell=inp('cells','Source IMR-90 cells; fixation/FISH/EMSA and bisulfite validation excluded',K.BIOLOGICAL,primary=True)
    x=op('wcna',cell,'Source overnight37°C proteinaseK/SDS lysis, salt-out/centrifugation, ethanol hook recovery/wash and phenol-chloroform WCNA extraction',S.MAGIC,reason='Operational source sequence known; exact native hybrid recovery and lysis proteinaseK quantity ambiguity not resolved.')
    shear=inp('covaris','Covaris source five10s on/1min off cycles; main400–600bp versus SI approximately400bp',K.CONDITION)
    x=op('sonicate',x,'Sonicate whole-cell nucleic acids before capture',S.FRAGMENTATION,((shear,R.CONDITION),))
    rnase=inp('rnase_i','RNaseI,37°C1.5h')
    x=op('rnase_i_digest',x,'Source RNaseI pretreatment; cleavage products remain in mixture',S.DIGESTION,((rnase,R.CATALYST),))
    x=clean('pre_ip_cleanup',x,'Phenol-chloroform/ethanol;resuspend EB,not TE')
    control=g.product('input_control','Reserved sonicated/RNaseI-treated input;not RNaseI-untreated cells',K.MOLECULAR)
    ip=g.product('ip_aliquot','Source3µg WCNA for IP',K.MOLECULAR)
    g.react('allocate_input','Reserve input control separately from IP aliquot',S.SPLIT,((x,R.TEMPLATE),),(ip,control));scope.append('allocate_input')
    g.mark_terminal(Terminal(control.id,Disposition.UNRESOLVED,'Source makes matched directional input library by the same downstream program;this graph expands IP library only. Not discarded or silently sequenced without preparation.'))
    antibody=inp('s96','SourceS9.6 antibody;overnight4°C')
    x=op('bind_antibody',ip,'Incubate nucleic acids with S9.6 antibody',S.GENERIC,((antibody,R.CATALYST),))
    beads=inp('beads','M280 sheep anti-mouse Dynabeads;2h4°C',K.SELECTION)
    bound=g.product('bound','Magnetic retained antibody/bead-associated hybrid preparation',K.MOLECULAR,state=unknown('Bound source hybrid-containing mixture'))
    waste=g.product('ip_wash','Magnet-removed supernatant and three washes',K.WASTE)
    g.react('capture','Capture immune complexes then magnet-remove supernatant/wash3times',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(bound,waste));g.discard(waste,'Actual source magnetic transfers;not perfect specificity/purity.');scope.append('capture')
    pk=inp('proteinase_k','100µg proteinaseK37°C2h,then60µg37°C1h')
    x=op('proteolysis',bound,'Digest antibody/protein;beads and remnants are not automatically removed',S.DIGESTION,((pk,R.CATALYST),))
    x=clean('ip_release',x,'Discard beads;phenol-chloroform supernatant;Ethachinmate/ethanol recovery into40µlEB')
    synthesis=inp('second_mix','SI E.coli DNA polymerase+ligase,RNaseH,dUTP,second-strand buffer;16°C2h')
    x=op('second_strand',x,'RNaseH-dependent priming/replacement and dUTP-marked DNA copying from native hybrid;not first-strand RT',S.MAGIC,((synthesis,R.CATALYST),),state=UracilFates(),reason='Main operation/strand attribution and SI reagent order known;primer fragments,replacement extent,nucleotide-mixture formulation and end geometry unassembled.')
    polish=inp('polish_mix','T4 DNApolymerase16°C10min,75°C20min,37°C15min')
    x=op('polish',x,'Source T4DNApolymerase processing and thermal program',S.MAGIC,((polish,R.CATALYST),),reason='Exact end processing across original/copied/RNA-remnant classes unassembled.')
    x=clean('second_cleanup',x,'ZymoDNA Clean&Concentrator cleanup')
    repair=inp('repair_mix','T4DNApolymerase,T4PNK,dNTP,T4ligasebuffer+ATP;1h roomtemperature')
    x=op('end_repair',x,'Source end repair/phosphorylation',S.MAGIC,((repair,R.CATALYST),),reason='Source operational end repair known;exact substrate ends and sealed/phosphorylated products unassembled.')
    x=clean('repair_cleanup',x,'Zymo cleanup after end repair')
    tail=inp('tail_mix','Source AddA step:exo-minusKlenow/NEB2/printed1mMATP;37°C30min')
    x=op('a_step',x,'Source AddA-overhang operation;printed ATP reagent retained unresolved',S.MAGIC,((tail,R.CATALYST),),reason='SI prints ATP,not dATP;no silent nucleotide correction or established A-tailed DNA product.')
    x=clean('tail_cleanup',x,'Zymo cleanup after source AddA step')
    index_state=Scene();index_state.strand('index_component',[Segment('TableS4 sample index class;not a joined whole adapter','[printed index window; full placement/width unassigned]',placeholder=True,feature=SAMPLE)])
    adapters=inp('truseq_adapters','Historical TruSeq indexed adapters;complete polymers/ends unassembled',K.ADAPTER,state=index_state)
    ligase=inp('quick_ligase','Quickligase/2xligationbuffer;roomtemperature15min')
    geom=AdapterLigationGeometry('unavailable','unavailable',(unknown('Source substrate/adapter junctions unassembled'),),'Source ligation reported;exact ends/historical polymers and ATP-step outcome unassembled,not claimed unreported.')
    x=op('adapter_ligation',x,'Source indexed adapter ligation',S.LIGATION,((adapters,R.STRUCTURAL),(ligase,R.CATALYST)),geometry=geom)
    x=clean('ligation_cleanup',x,'Zymo cleanup after adapter ligation')
    f=inp('forward','Source10µMforwardPCRprimer;full polymer unreviewed',K.OLIGO)
    r=inp('reverse','Source10µMreversePCRprimer;full polymer unreviewed',K.OLIGO)
    ung=inp('ung_enzyme','UNG37°C15min')
    ungstate=selected_ung_event(UracilFates())
    x=op('ung',x,'UNG uracil-base excision;AP/backbone/RNA/remnant classes and free bases remain',S.DIGESTION,((ung,R.CATALYST),(f,R.STRUCTURAL),(r,R.STRUCTURAL)),state=ungstate)
    heat=inp('heat','98°C10min',K.CONDITION)
    x=op('ung_heat',x,'Source post-UNG thermal treatment;not physical removal',S.THERMAL_TREATMENT,((heat,R.CONDITION),),state=heat_treatment(ungstate))
    pol=inp('phusion','SourcePhusionHF:98°C2min,18cycles98/60/72,72°C10min')
    x=op('pcr',x,'Source18cycle library PCR;strand survival/recovery not predicted',S.PCR,((f,R.STRUCTURAL),(r,R.STRUCTURAL),(pol,R.CATALYST)),geometry=PCRGeometry.unavailable(Scene(),reason='Historical primer/product geometry and AP/heat-product amplification limits unassembled. No direct whole-strand disposal attributed to UNG.'))
    x=clean('final_cleanup',x,'Immediate Zymo final-library cleanup;Qubit/BioanalyzerQC')
    library=g.product('library','Source indexed IP library;conditional adapter-dimer gel cleanup not inserted as universal step',K.LIBRARY,state=unknown('Source sequenced library;full strands unassembled'))
    g.react('library_handoff','Reported library acquisition handoff',S.GENERIC,((x,R.TEMPLATE),),(library,));scope.append('library_handoff')
    g.sequence(library,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':150,'Read 2':150},limitation='Source HiSeq2500PE150;actual primer landing,index cycles and strand geometry remain unassembled. TableS4 fixed sixbases plus parentheticalseventh retained verbatim,not silently cropped or reverse-complemented.'))
    g.annotate(SourceDiscrepancy('a_step','Printed adenylation reagent', (SourceClaim('SI1 AddAoverhang step,pp22–23','Heading says AddAoverhangs;reagent prints1mMATP5µl.'),SourceClaim('Main directional RDIP Methods206','Reports adenosine tails before adapter ligation.')), 'Keep operational handoff unresolved;do not change ATP to dATP.', 'Published reagent does not specify dATP despite the stated DNA-tail outcome.', 'Actual nucleotide formulation and exact product ends remain unassembled.'))
    g.declare_reviewed_scope(scope,SOURCE)
    g.add_source_item(SourceItem('unreviewed_geometry','Historical buffers,adapter/primer/UNG-product geometry and cited precursor method','Main reference49 and source kit identities not followed under cap;no confirmed absence.'))
    return g


def render_page():
    return head('RDIP-seq — source-limited model')+'<main class="wrap"><h1>RDIP-seq</h1><p><a href="01_rdip.html">Research notes</a></p>'+info(LIMIT)+info('Published SI contradictions: sonication200ng/ml versus RNaseI100µl20µg;lysis0.5mg proteinaseK versus2.5µl20mg/ml;AddA step printsATP rather than dATP. Main calls TableS4 adapter sequences,while the workbook provides index windows only. Main describes a dUTP-containing RNA strand;the model distinguishes newly copied DNA from original RNA and does not attribute RNA digestion to UNG. No quantities/nucleotides/whole adapters are silently repaired.')+workflow_panel(workflow_from_protocol_graph(protocol_graph(),'library'))+'</main>'


def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();initial=UracilFates();processed=selected_ung_event(initial);heated=heat_treatment(processed)
    c('Selected UNG processing preserves originalDNA and RNA-remnant classes',(processed.original_dna,processed.rna_remnants),(True,True))
    c('Heating records treatment without purifying or reverting AP processing',(heated.processed,heated.original_dna,heated.rna_remnants),(True,True,True))
    try:selected_ung_event(processed);bad=False
    except ValueError:bad=True
    c('Selected processing event cannot be silently duplicated',bad,True)
    g=protocol_graph();c('Physical postligation cleanup precedes UNG',g.steps['ung'].inputs[0].material,'ligation_cleanup_kept')
    c('No automatic disposal at UNG output',g.steps['ung_heat'].inputs[0].material,'ung_out')
    c('Native hybrid copying is not imported ordinary first-strandRT',any(s.kind==S.REVERSE_TRANSCRIPTION for s in g.steps.values()),False)
    folder=sources_root()/'rdip-seq__10.1186+s13072-015-0040-6'
    tables=folder/'s13072-015-0040-6_supp_13072_2015_40_MOESM2_ESM.xlsx.txt'
    if tables.exists():
        block=tables.read_text().split('## sheet Table S4\n',1)[1].split('## sheet Table S5',1)[0]
        for name,num,raw in INDEX_ROWS:c('SourceTableS4 exact printed member '+name,(name+'\t'+str(num)+'\t'+raw) in block,True)
    else:c.skip('Original published TableS4 absent')
    c.report()
