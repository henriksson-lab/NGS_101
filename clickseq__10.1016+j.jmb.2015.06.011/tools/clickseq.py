"""Source-local azido-cDNA and triazole-linked template classes for ClickSeq."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
from chemdraw import Scene,Segment,workflow_from_protocol_graph,workflow_panel
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan
from pcr import PCRGeometry
from page import head,info

PAGE_FILE='clickseq.html'
SOURCE='Routh et al.2015 DOI10.1016/j.jmb.2015.06.011;associated author manuscript PMC4523409 main30–34/38–58 andFig1–2. ActualNIHMS705616-supplement.docx/TableS1 is gated;finalVOR applicability unverified.'
LIMIT='Private incomplete source model. Stochastic azido termination shortens newly copied cDNA,not the original RNA. Successful click joining makes a nonphosphodiester triazole template;PCR copying does not copy that chemical bond onto every daughter. Positions,yields,full modified oligos,enzymes/conditions and source-specific library/read geometry remain unassembled. INFERRED — PCR nucleotide supply is a functional requirement, not a disclosed recipe or dose. The selected source100ngFHV RNA library is distinct from the independentCrPV library and targetedPrimer5 validation;hostmRNA mapping comes from the FHV dataset.'


def validate_sequence(sequence):
    if sequence is not None and (not sequence or set(sequence)-set('ACGT')):raise ValueError('Explicit finite DNA or an unassembled sequence is required')


@dataclass(frozen=True)
class AzidoCDNA:
    sequence:str|None=None
    reactive_azide:bool=True
    inferred:bool=True
    def __post_init__(self):validate_sequence(self.sequence)
    def segment(self):
        return role('copied cDNA including its terminal azido-dideoxynucleotide;positions/length unassembled') if self.sequence is None else Segment('selected copied cDNA',self.sequence,inferred=self.inferred)


@dataclass(frozen=True)
class HexynylAdapter:
    sequence:str|None=None
    reactive_alkyne:bool=True
    inferred:bool=True
    def __post_init__(self):validate_sequence(self.sequence)
    def segment(self):
        return role('source Oligo2 adapter DNA body;whole polymer unassembled') if self.sequence is None else Segment('selected source-independent adapter DNA',self.sequence,inferred=self.inferred)


@dataclass(frozen=True)
class ClickTemplate:
    cdna:AzidoCDNA
    adapter:HexynylAdapter
    def __post_init__(self):
        if not self.cdna.reactive_azide or not self.adapter.reactive_alkyne:raise ValueError('Selected click relation needs both supplied reactive termini')
    @property
    def nucleotide_count(self):
        return None if self.cdna.sequence is None or self.adapter.sequence is None else len(self.cdna.sequence)+len(self.adapter.sequence)
    @property
    def reactive_termini(self):return (False,False)
    def parts(self):
        return [self.cdna.segment(),role('nonphosphodiester triazole bridge;not a DNA base'),self.adapter.segment()]
    def rows(self):
        sc=Scene();sc.strand('joined',self.parts(),label='One selected successful covalent click-template class,5prime to3prime;not all source molecules')
        sc.note('joined','The terminal azido nucleotide is retained. The consumed azide/alkyne ends are represented by the triazole junction;no ordinary DNA ligase bond is substituted.')
        return sc.rows()


@dataclass(frozen=True)
class ClickMixture:
    stage:str
    def __post_init__(self):
        if self.stage not in {'rt','click','pcr'}:raise ValueError('Explicit RT,click orPCR source stage required')
    def rows(self):
        sc=Scene()
        if self.stage in {'click','pcr'}:
            sc.strand('clicked',ClickTemplate(AzidoCDNA(),HexynylAdapter()).parts(),label='Triazole-linked template class;whole geometry unassembled')
        sc.strand('unclicked',[AzidoCDNA().segment()],label='Unclicked azido-cDNA class;amount/recovery unassigned',mod3='3′azido block,not3′OH')
        sc.strand('others',[role('possible nonterminated cDNA,free primers,reagents and originalRNA/carryover;physical recovery unassigned')],label='Other material is not automatically disposed at copying/joining')
        if self.stage in {'click','pcr'}:sc.strand('free_adapter',[HexynylAdapter().segment()],label='Possible excess source click adapter/reagents',mod5='5′hexynyl')
        if self.stage=='pcr':sc.strand('copies',[role('reported PCR-derived DNA;sequence/error/length/yield unassembled;native template bond not copied as a chemical triazole')],label='Newly synthesized DNA class,not chemical modification transfer')
        return sc.rows()


def primer_inventory():
    sc=Scene();sc.strand('primer',[role('source constant Illumina3prime-adapter portion;fullTableS1 unread'),Segment('random RNA-priming hexamer','N'*6)],label='One source semirandom RT primer;actual degenerate N6 synthesis, not an unassigned index or UMI');return sc


def protocol_graph():
    g=ProtocolGraph('ClickSeq — source FHV-RNA chemical joining route',authority=SOURCE)
    def inp(key,label,kind=K.CATALYST,state=None,primary=False):return g.input(key,label,kind,state=state,primary=primary)
    def op(key,x,label,kind,reagents=(),state=None,geometry=None,reason=None):
        out=g.product(key+'_out',label+' preparation',K.MOLECULAR,state=state)
        g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(reagents),(out,),geometry=geometry,**({'reason':reason} if reason else {}));return out
    def clean(key,x,label):
        reagent=inp(key+'_reagent',label,K.SELECTION)
        kept=g.product(key+'_kept','Source retained preparation;exact class recovery unresolved',K.MOLECULAR,state=x.state)
        waste=g.product(key+'_removed','Source physically removed fractions;composition/purity not assigned',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(reagent,R.SELECTION)),(kept,waste));g.discard(waste,'Actual source physical workup;not a guarantee of complete class removal.');return kept
    rna=inp('rna','Source100ng purified FHV-RNA preparation;hostRNA not declared absent. CrPV is a separate source library,not pooled here;upstream extraction outside graph',K.MOLECULAR,primary=True)
    primer=inp('rt_primer','Source Oligo1 semirandom RT primer;complete modified polymer unassembled',K.OLIGO,primer_inventory())
    ntps=inp('az_ntps','Source FHV/CrPV library1:20AzNTP:dNTP ratio;3primeazido2prime,3primedideoxy analogs. Absolute mix and per-site choices unassembled',K.MOLECULAR)
    rt=inp('reverse_transcriptase','Source reverse transcription;exact enzyme/buffer/thermal program gated in supplement')
    x=op('rt',rna,'Source random priming and stochastic azido-nucleotide termination;no original-RNA fragmentation step',S.REVERSE_TRANSCRIPTION,((primer,R.STRUCTURAL),(ntps,R.STRUCTURAL),(rt,R.CATALYST)),ClickMixture('rt'))
    x=clean('rt_cleanup',x,'Source DNA cleaning/desalting columns remove unincorporated AzNTPs')
    sc=Scene();sc.strand('adapter',[HexynylAdapter().segment()],label='Source Oligo2 modified DNA body',mod5='5′hexynyl')
    adapter=inp('click_adapter','Source excess5primehexynyl modified Illumina5primeadapter,Oligo2',K.ADAPTER,sc)
    copper=inp('copper','Source active Cu(I)-TBTA click catalyst;stock formulation and reaction geometry unresolved')
    reductant=inp('ascorbate','Actual source ascorbate molecular reducing reagent;dose/stock formulation unreviewed',K.MOLECULAR)
    x=op('click_join',x,'Source Cu-catalyzed azide–alkyne joining produces triazole-linked DNA;unclicked classes remain',S.MAGIC,((adapter,R.STRUCTURAL),(copper,R.CATALYST),(reductant,R.STRUCTURAL)),ClickMixture('click'),reason='Selected reactive termini and nonphosphodiester product class known;full polymer/junction geometry and conversion fractions unassembled. Not a conventional phosphodiester ligation or perfect capture.')
    f=inp('pcr_f','Source PCRPrimer3;fullTableS1 polymer/member unreviewed',K.OLIGO)
    r=inp('pcr_r','Source PCRPrimer4;fullTableS1 polymer/member unreviewed',K.OLIGO)
    pol=inp('pcr_polymerase','Source common-polymerase read-through sufficient for PCR;precise library enzyme/cycles unreviewed. No efficiency prediction.')
    nts=g.input('pcr_nucleotides','INFERRED — nucleotide supply required for reported PCR copying; actual source formulation/dose unreviewed',K.MOLECULAR,inferred=True)
    x=op('pcr',x,'INFERRED — functional nucleotide inventory for source directPCR; new DNA coexists with original modified templates/remnants',S.PCR,((f,R.STRUCTURAL),(r,R.STRUCTURAL),(pol,R.CATALYST),(nts,R.STRUCTURAL)),ClickMixture('pcr'),PCRGeometry.unavailable(Scene(),reason='SourcePCR reported but full primer/product geometry and polymerase behavior at triazole junction unassembled;no guaranteed read-through or covalent-bond transfer.'))
    x=clean('bead_size',x,'Source AMPure bead fragment-size selection after PCR')
    x=clean('gel_size',x,'Source final300–500nt library selection and gel extraction;precise workup details unreviewed')
    library=g.product('library','Reported source ClickSeq library;exact full product/primer/strand geometry unassembled',K.LIBRARY,state=x.state)
    g.react('library_handoff','Source sequencing-library handoff; no new molecular reaction',S.GENERIC,((x,R.TEMPLATE),),(library,),changes_molecule=False)
    g.sequence(library,SequencingPlan(acquired_roles=('Read 1',),read_lengths={'Read 1':100},limitation='Associated main40 reports independent100ngFHV andCrPV libraries,1:20AzNTP:dNTP,300–500ntselection,HiSeq single-end100nt. This graph selectsFHV context;hostmRNA mapping is from that dataset,not a separate mRNA preparation. Full index/primer/end geometry and finalVOR applicability remain unassembled;targetedPrimer5 validation settings not substituted.'))
    # Curator-written main-supported inventory, not accumulated by helper calls.
    g.declare_reviewed_scope(('rna','rt_primer','az_ntps','reverse_transcriptase','rt','rt_cleanup_reagent','rt_cleanup','click_adapter','copper','ascorbate','click_join','pcr_f','pcr_r','pcr_polymerase','pcr','bead_size_reagent','bead_size','gel_size_reagent','gel_size','library_handoff'),SOURCE+' Main30–40 operational scope only; exact supplement recipe/end/library scope remains unreviewed.')
    g.declare_reviewed_scope(('pcr_nucleotides',),SOURCE+' Main34 reports PCR; nucleotide requirement is explicitly inferred, not a published formulation.')
    g.add_source_item(SourceItem('supplement','Actual chemistry supplement/TableS1','Main JATS points NIHMS705616-supplement.docx;standard declared retrieval hit PMCgate. Wholeoligos,preciseenzymes/ratios/rungeometry and finalVORapplicability unreviewed;no absence claim.'))
    return g


def render_page():
    return head('ClickSeq — source-limited model')+'<main class="wrap"><h1>ClickSeq</h1><p><a href="01_clickseq.html">Research notes</a></p>'+info(LIMIT)+workflow_panel(workflow_from_protocol_graph(protocol_graph(),'library'))+'</main>'


def selftest():
    from checks import Check
    c=Check();left=AzidoCDNA('ACGT');right=HexynylAdapter('TGCA');joined=ClickTemplate(left,right)
    c('Click joining preserves supplied DNA nucleotide count',joined.nucleotide_count,len(left.sequence)+len(right.sequence))
    c('Terminal azido nucleotide is not deleted',joined.cdna.sequence[-1],left.sequence[-1])
    c('Successful selected junction consumes reactive termini',joined.reactive_termini,(False,False))
    c('Unknown body length remains unknown',ClickTemplate(AzidoCDNA(),right).nucleotide_count,None)
    try:ClickTemplate(AzidoCDNA('ACGT',False),right);bad=False
    except ValueError:bad=True
    c('No click junction inferred without an azido reactant',bad,True)
    try:ClickTemplate(left,HexynylAdapter('TGCA',False));bad=False
    except ValueError:bad=True
    c('No click junction inferred without an alkyne reactant',bad,True)
    g=protocol_graph();c('DirectPCR uses click reaction rather than invented intermediate cleanup',g.steps['pcr'].inputs[0].material,'click_join_out')
    c('Degenerate priming N6 is not an unassigned placeholder',primer_inventory().strands['primer'].segs[-1].placeholder,False)
    c('Ascorbate is a molecular reaction contributor',g.steps['click_join'].inputs[-1].role,R.STRUCTURAL)
    c('Selected library preserves PCR/template population through actual cleanup',g.materials['library'].state.stage,'pcr')
    c('No original-RNA fragmentation operation inserted',any(s.kind==S.FRAGMENTATION for s in g.steps.values()),False)
    c.report()
