"""Source-qualified irCLIP adapter preparation and single-strand ligation."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Construct, Segment, Scene, feature, revcomp, panel
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import ProtocolGraph, MaterialKind, InputRole, StepKind, Terminal, Disposition, SourceDiscrepancy, SourceClaim, DetailView, Highlight
from illumina import PE_PCR_PRIMER_2
from circular import CircularizationGeometry
from pcr import PCRGeometry
from protocol_views import WorkflowSection,assemble_page
from page import head,info
from protocol_graph import SequencingPlan,AnnotationSlot
from spacer_backbone import DNA, VariableDNA, C18, Backbone, StopHypothesis, bind_suffix, extend_to_linker
SOURCE=('Zarnegar et al.2016 doi:10.1038/nmeth.3840, accepted manuscript PMC5477425 '
        'Online Methods JATS55,63,67,71–77; final-VOR procedural applicability unverified')
ARM='AGATCGGAAGAGCGGTTCAG'
ADAPTER_WRITTEN='5′-OH-AGATCGGAAGAGCGGTTCAGAAAAAAAAAAAA/iAzideN/AAAAAAAAAAAA/3Bio/-3′'
RT_WRITTEN='/5phos/WWWNNNATCACGNNNNNTACCCTTCGCTTCACACACAAG/iSp18/GGATCC/iSp18/TACTGAACCGC'
P3SHORT='CTGAACCGCTCTTCCGATCT'
P3TALL='GCATTCCTGCTGAACCGCTCTTCCGATCT'
P6TALL='TTTCCCCTTGTGTGTGAAGCGAAGGGTA'
P3SOLEXA=PE_PCR_PRIMER_2
P6SOLEXA='AATGATACGGCGACCACCGAGATCTACACTCTTTCCCCTTGTGTGTGAAGCGAAGGGTA'
P6SEQUENCING='CACTCTTTCCCCTTGTGTGTGAAGCGAAGGGTA'
BOUNDARY=('The spacer-bearing RT primer and cDNA circle are retained as typed ordered '
          'backbones. The first extension is conditional on an explicitly inferred C18 '
          'polymerase-stop hypothesis. Actual source two-PCR finishing is declared at a coarse product-class level; full occupied products, acquired read geometry and final-version applicability remain unresolved.')


def rt_primer_backbone():
    """Manuscript barcode1 oligo, with actual non-nucleotide linker nodes."""
    def umi(part):
        return feature('irclip_degenerate_identifier_'+part.lower(),'umi','degenerate',
                       group='irclip_degenerate_identifier',part=part,
                       note='Composite WWWNNN…NNNNN; molecular identifier role inferred from degenerate synthesis and duplicate handling, not printed as UMI.')
    return Backbone((
        DNA('degenerate W3','WWW',umi('W3')),
        DNA('random N3','NNN',umi('N3')),
        DNA('inline barcode1','ATCACG',feature('irclip_sample','sample_index','fixed')),
        DNA('random N5','NNNNN',umi('N5')),
        DNA('P6 handle','TACCCTTCGCTTCACACACAAG'),
        C18('first C18'), DNA('inter-linker BamHI motif','GGATCC'), C18('second C18'),
        DNA('RT tail','TACTGAACCGC'),
    ))


def linear_cdna_backbone():
    """Conditional exact adapter-arm copying, with unassigned captured RNA sequence.

    The longest exact terminal suffix of the printed RT tail pairs the specified
    adapter arm. No binding-affinity/yield claim follows from that sequence match.
    """
    primer=rt_primer_backbone()
    adapter_arm=Backbone((DNA('adapter arm',ARM),),five_prime='OH')
    landing=bind_suffix(adapter_arm,primer.nodes[-1].sequence,block_name='adapter arm')
    upstream_copy=revcomp(ARM[:landing.start])
    return Backbone((*primer.nodes,DNA('copied adapter arm',upstream_copy),
                     VariableDNA('captured RNA cDNA','captured RNA',True)),inferred=True)


def conditional_first_extension(circle=None):
    circle=linear_cdna_backbone().close() if circle is None else circle
    landing=bind_suffix(circle,P6TALL,block_name='P6 handle')
    return extend_to_linker(circle,landing,hypothesis=StopHypothesis(
        frozenset({'/iSp18/'}),
        'INFERRED mechanism: Phusion extension stops before a C18 linker. The inspected '
        'accepted manuscript says the two spacers enable linear PCR products from a '
        'circular cDNA; it does not specify these exact traversal/stop coordinates.'))


@dataclass(frozen=True)
class SpacerMixture:
    backbone: Backbone
    other_material: tuple[str,...]


@dataclass(frozen=True)
class Adapter:
    five_prime: str='OH'
    internal: str='azide'
    three_prime: str='biotin'
    def __post_init__(self):
        if self.five_prime not in {'OH','phosphate','App'} or self.internal not in {'azide','IR800CW conjugate'}:
            raise ValueError('unsupported source adapter state')
        if self.three_prime != 'biotin':
            raise ValueError('source 3-prime biotin must be retained')
    def phosphorylate(self):
        if self.five_prime!='OH': raise ValueError('OptiKinase substrate must have 5-prime OH')
        return replace(self,five_prime='phosphate')
    def adenylate(self):
        if self.five_prime!='phosphate': raise ValueError('adenylation requires phosphate')
        return replace(self,five_prime='App')
    def click(self):
        if self.internal!='azide': raise ValueError('azide conjugation requires unreacted azide')
        return replace(self,internal='IR800CW conjugate')
    def scene(self):
        sc=Scene()
        sc.strand('adapter', [Segment('L3 DNA arm',ARM),Segment('A12 upstream','A'*12)],
                  label='irCLIP adapter; first nucleotide tract',mod5=self.five_prime)
        sc.strand('adapter_after_site',[Segment('A12 downstream','A'*12)],
                  label='same covalent adapter; tract after internal '+self.internal,mod3=self.three_prime)
        sc.note('adapter','ONE covalent source oligo: tracts linked through /iAzideN/ chemistry ('+self.internal+'). No direct nucleotide adjacency or extra base is invented at that modification. Source calls the site nucleotide position32; atom map/sugar convention unassembled.')
        sc.ligation_end('adapter','L3 DNA arm',kind='single-strand',polarity='5-prime',owner='adapter')
        return sc

def protocol_graph():
    g=ProtocolGraph('irCLIP — adapter activation and RNA ligation',authority=SOURCE,
                    representation='published_instance')
    K,R,S=MaterialKind,InputRole,StepKind
    def select(i,x,label,kind=None,disposition=Disposition.UNRESOLVED):
        y=g.product(i+'_selected',label,kind or x.kind,state=x.state,chemical_state=x.chemical_state)
        rest=g.product(i+'_remaining','Other source stock fraction; quantity/future use unassembled',x.kind,state=x.state,chemical_state=x.chemical_state)
        g.react(i,label,S.SPLIT,((x,R.TEMPLATE),),(y,rest));g.mark_terminal(Terminal(rest.id,disposition,'Untransferred actual stock; not automatically discarded.'))
        return y
    def clean(i,x,label,carrier=None,compound=False):
        c=g.input(i+'_conditions',label,K.CONDITION)
        y=g.product(i+'_retained','Retained '+label,x.kind,state=x.state,chemical_state=x.chemical_state)
        w=g.product(i+'_removed','Physical separation reject; composition unassigned',K.WASTE)
        ports=[(x,R.TEMPLATE),(c,R.CONDITION)]
        if carrier:ports.append((carrier,R.STRUCTURAL))
        g.react(i,label,S.MAGIC if compound else S.CLEANUP,ports,(y,w),reason=('Source array cleanup/eluate pooling; individual recovery/yield unassembled.' if compound else ''))
        g.discard(w,'Source physical fractionation only; no quantitative recovery guarantee.')
        return y
    def qc(i,x,label):
        q=g.product(i+'_qc',label,K.QC,state=x.state)
        rest=g.product(i+'_remaining','Remaining source stock after measured aliquot',x.kind,state=x.state,chemical_state=x.chemical_state)
        g.react(i,label,S.SPLIT,((x,R.TEMPLATE),),(q,rest));g.mark_terminal(Terminal(q.id,Disposition.QC,'Actual source measurement aliquot; not a sequencing contribution.'))
        return rest
    state=Adapter()
    stock=g.input('adapter_oh','Source20nmol substrate '+ADAPTER_WRITTEN,MaterialKind.ADAPTER,state=tuple(state.scene().rows()),chemical_state=state)
    prior=stock
    for ident,label,enzyme,op,kind in (
        ('phosphorylate','OptiKinase + ATP (37 °C, 1 h)','OptiKinase',Adapter.phosphorylate,StepKind.PHOSPHORYLATION),
        ('adenylate','Mth RNA ligase + ATP (65 °C, 2 h)','Mth RNA ligase',Adapter.adenylate,StepKind.CHEMICAL_CONVERSION),
        ('click','IR800CW–DBCO conjugation (Methods:37 °C, 2 h; FigS1a:1 h)','Source20µl10mMIR800CW–DBCO',Adapter.click,StepKind.CHEMICAL_CONVERSION)):
        reagent_kind=MaterialKind.MOLECULAR if ident=='click' else MaterialKind.CATALYST
        reagent=g.input(ident+'_reagent',enzyme,reagent_kind)
        role=InputRole.STRUCTURAL if ident=='click' else InputRole.CATALYST
        state=op(state)
        product=g.product(ident+'_adapter',label,MaterialKind.ADAPTER,state=tuple(state.scene().rows()),chemical_state=state)
        ports=[(prior,R.TEMPLATE),(reagent,role)]
        if ident=='click':
            click_condition=g.input('click_conditions','SourcePBS,37C2h; AMscope, publisherFigS1a1hr comparison retained',K.CONDITION)
            ports.append((click_condition,R.CONDITION))
        if ident!='click':
            atp=g.input(ident+'_atp',('Source1mMATP in0.2ml phosphorylation' if ident=='phosphorylate' else 'Source75µl1mMATP added for adenylation'),K.MOLECULAR)
            ports.append((atp,R.STRUCTURAL))
        g.react(ident,label,kind,ports,(product,))
        prior=product
        if ident=='phosphorylate':prior=clean('phospho_precipitation',prior,'Entire source phosphorylated stock acetate/ethanol recovery;0.5mlwater')
        elif ident=='adenylate':
            prior=select('preadapter_half',prior,'Source half0.375ml preadenylated reaction; other half stock retained')
            prior=clean('preadapter_precipitation',prior,'Source half-reaction acetate/ethanol recovery;0.18mlPBS')
        else:prior=clean('adapter_columns',prior,'Source200µlclick+4.8mlPNI→20columns250µleach→50µleluates pooled',compound=True)
    prior=select('adapter_ligation_aliquot',prior,'Source1pmol prepared adapter selected for ligation; parent working stock retained')
    rna=g.input('rna_3oh','RNase-A/1-treated and PNK-dephosphorylated crosslinked RNA on beads',
                MaterialKind.MOLECULAR,primary=True,
                state=Construct([Segment('crosslinked RNA','[variable RNA]',placeholder=True)]))
    ligase=g.input('t4_rnl1','T4 RNA ligase 1',MaterialKind.CATALYST)
    joined=g.product('ligated_rna','RNA–adapter junction with internal IR800CW and 3-prime biotin retained',
                     MaterialKind.MOLECULAR)
    g.react('ligate','Ligate preadenylated DNA adapter to RNA (16 °C overnight, PEG400)',
            StepKind.LIGATION,((rna,InputRole.TEMPLATE),(prior,InputRole.STRUCTURAL),
                               (ligase,InputRole.CATALYST)),(joined,),
            geometry=AdapterLigationGeometry.known('single-strand',state.scene(),
                         statement='Online Methods55 and63: activated adapter 5-prime end ligates RNA 3-prime OH.'))
    pk=g.input('proteinase_k','proteinase K / SDS recovery',MaterialKind.CATALYST)
    recovered=g.product('recovered_rna','gel/membrane-recovered peptide-bearing RNA–adapter',MaterialKind.MOLECULAR)
    g.react('recover_rna','SDS-PAGE/nitrocellulose and proteinase K recovery/extraction',
            StepKind.DEPROTEINIZATION,((joined,InputRole.TEMPLATE),(pk,InputRole.CATALYST)),(recovered,))
    recovered=qc('proteinase_qc',recovered,'Source1µlproteinase reaction dot-blot check')
    recovered=clean('rna_extraction',recovered,'Phenol/phase-lock plus two actual chloroform extractions; aqueous phase retained')
    lpa=g.input('linear_polyacrylamide','Source3µllinearpolyacrylamide carrier',K.MOLECULAR)
    recovered=clean('rna_precipitation',recovered,'NaCl/ethanol precipitation and75percent wash;15µlresuspension',lpa)
    recovered=qc('rna_recovery_qc',recovered,'Source1µlof15µlRNA dot blot')
    recovered=select('rna_rt_aliquot',recovered,'Source13.75µl RNA for RT; residual actual stock unexpanded')
    primer=rt_primer_backbone()
    primer_input=g.input('rt_primer','manuscript barcode1 RT primer: phosphorylated DNA with two C18 links',
                         MaterialKind.OLIGO,state=tuple(primer.rows()),chemical_state=primer)
    annealed=g.product('annealed_rna','Same-tube RNA template component after70C10min/4C',K.MOLECULAR)
    retained=g.product('retained_rt_primer','Retained same-tube RTprimer component, not new dose',K.OLIGO,state=primer_input.state,chemical_state=primer)
    g.react('anneal_rt','Source RNA and primer annealed before RT enzyme',S.HYBRIDIZATION,((recovered,R.TEMPLATE),(primer_input,R.STRUCTURAL)),(annealed,retained))
    tgirt=g.input('tgirt','Source0.2µlTGIRT-III plus sourcebuffer;roomtemperature30min before dNTP addition',MaterialKind.CATALYST)
    nts=g.input('rt_nucleotides','Source1µl25mMdNTP added after enzyme preincubation;60C2h',K.MOLECULAR)
    preincubated=g.product('preincubated_rna','Same-tube template component after sourceRTenzyme/buffer incubation',K.MOLECULAR)
    primer_carried=g.product('preincubated_primer','Retained primer component after enzyme incubation',K.OLIGO,state=retained.state,chemical_state=retained.chemical_state)
    enzyme_carried=g.product('retained_tgirt','Retained same-tube TGIRT component, not fresh added enzyme',K.CATALYST)
    g.react('rt_preincubation','Source adds TGIRT/buffer; roomtemperature30min before dNTP addition',S.HYBRIDIZATION,((annealed,R.TEMPLATE),(retained,R.STRUCTURAL),(tgirt,R.CATALYST)),(preincubated,primer_carried,enzyme_carried))
    linear=linear_cdna_backbone()
    hybrid=g.product('rt_hybrid','INFERRED adapter-arm landing; RNA–cDNA hybrid with spacer-bearing cDNA',
                     MaterialKind.MOLECULAR,state=tuple(linear.rows()),
                     chemical_state=SpacerMixture(linear,('RNA template with biotin adapter retained',)),inferred=True)
    g.react('reverse_transcribe','TGIRT-III RT (60 °C, 2 h); exact landing/copy inferred from printed arm',
            StepKind.REVERSE_TRANSCRIPTION,((preincubated,InputRole.TEMPLATE),(primer_carried,InputRole.STRUCTURAL),(nts,InputRole.STRUCTURAL),
                                           (enzyme_carried,InputRole.CATALYST)),(hybrid,))
    beads=g.input('streptavidin','MyOne Streptavidin C1 Dynabeads',MaterialKind.SELECTION)
    captured=g.product('captured_hybrid','biotin-mediated captured RNA–cDNA hybrid',MaterialKind.MOLECULAR,
                       chemical_state=hybrid.chemical_state,state=hybrid.state,inferred=True)
    flowthrough=g.product('capture_flowthrough','unretained material and bead washes',MaterialKind.WASTE)
    g.react('capture_hybrid','Streptavidin capture and washes',StepKind.AFFINITY_CAPTURE,
            ((hybrid,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(captured,flowthrough))
    g.discard(flowthrough,'Physical capture flow-through and washes; no exact composition assigned')
    p3=g.input('p3short','P3short cDNA elution oligo',MaterialKind.OLIGO,
               state=Construct([Segment('P3short',P3SHORT)]))
    heat=g.input('heat_manganese','P3short/MnCl2 heating and cooling',MaterialKind.CONDITION)
    eluate=g.product('cdna_eluate','INFERRED cDNA backbone in RNA-degradation/elution mixture',MaterialKind.MOLECULAR,
                     chemical_state=SpacerMixture(linear,('RNA degradation products','P3short; occupancy unassigned')),
                     state=tuple(linear.rows()),inferred=True)
    g.react('release_cdna','Source RNA degradation/cDNA elution; do not assign debris to waste',
            StepKind.CHEMICAL_CONVERSION,((captured,InputRole.TEMPLATE),(p3,InputRole.STRUCTURAL),
                                        (heat,InputRole.CONDITION)),(eluate,))
    circligase=g.input('circligase','CircLigase II',MaterialKind.CATALYST)
    circle=linear.close()
    closed=g.product('cdna_circle','INFERRED closed cDNA backbone; both C18 links and GGATCC retained',
                     MaterialKind.MOLECULAR,state=tuple(circle.rows()),
                     chemical_state=SpacerMixture(circle,eluate.chemical_state.other_material),inferred=True)
    g.react('circularize','CircLigase II (60 °C, 1.5 h); preserve all covalent backbone nodes',
            StepKind.CIRCULARIZATION,((eluate,InputRole.TEMPLATE),(circligase,InputRole.CATALYST)),(closed,),
            geometry=CircularizationGeometry.known('single-strand intramolecular ligation',
                statement='The manuscript prints a phosphorylated RT primer and cDNA circularization. C18 links remain in the substrate.'))
    clean_circle=g.product('recovered_circle','INFERRED circle after AMPure/isopropanol recovery; residual mixture unquantified',
                           MaterialKind.MOLECULAR,state=closed.state,
                           chemical_state=SpacerMixture(circle,('residual material partition unquantified',)),inferred=True)
    cleanup_waste=g.product('circle_cleanup_waste','AMPure/isopropanol recovery supernatant/washes',MaterialKind.WASTE)
    g.react('recover_circle','Physical cDNA-circle recovery before PCR1',StepKind.CLEANUP,
            ((closed,InputRole.TEMPLATE),),(clean_circle,cleanup_waste))
    g.discard(cleanup_waste,'Physical cleanup reject; RNA degradation products were molecular material until separation')
    # Keep the mechanistic hypothesis as a nonphysical analytical reference.
    # It neither allocates another circle aliquot nor performs a second PCR.
    extension=conditional_first_extension(circle)
    first=g.product('conditional_first_strand','INFERRED analytical first-copy reference; not an observed separated material',K.MOLECULAR,state=tuple(extension.product.rows()),chemical_state=extension,inferred=True,source=SOURCE+'; conditional C18-stop hypothesis only')
    g.react('conditional_extend','INFERRED analytical reference only: P6 copying under explicit C18-stop hypothesis',S.GENERIC,((clean_circle,R.TEMPLATE),),(first,),changes_molecule=False,reason='Nonphysical reference calculation, not an additional PCR/aliquot or evidence of a unique observed product.')
    g.mark_terminal(Terminal(first.id,Disposition.UNRESOLVED,BOUNDARY))
    stock1=g.input('pcr1_primer_stock','Source0.5µl30µMP3/P6 PCR1 oligonucleotide mix',K.MOLECULAR)
    p3tall=g.product('p3tall','Copresent source P3tall primer',K.OLIGO,state=Construct([Segment('P3tall',P3TALL)]))
    p6=g.product('p6tall','Copresent source P6tall primer',K.OLIGO,state=Construct([Segment('P6tall',P6TALL)]))
    g.react('pcr1_primer_components','Components of one actual source PCR1 primer mix',S.MAGIC,((stock1,R.TEMPLATE),),(p3tall,p6),changes_molecule=False,reason='Source explicit primer mix; not two fresh free doses.')
    def mastermix(i,label):
        stock=g.input(i+'_stock',label,K.MOLECULAR)
        enz=g.product(i+'_enzyme','Copresent PhusionHF polymerase component',K.CATALYST)
        nt=g.product(i+'_nucleotides','Copresent nucleotide component',K.MOLECULAR)
        buf=g.product(i+'_buffer','Copresent prepared-mix buffer',K.CONDITION)
        g.react(i+'_components','Inventory of one actual source prepared master mix',S.MAGIC,((stock,R.TEMPLATE),),(enz,nt,buf),changes_molecule=False,reason='Source names actual2xPhusionHFmastermix; functional stock components only, not fresh separate doses.')
        return enz,nt,buf
    enzyme1,nt1,buf1=mastermix('pcr1_mix','Source15µl2xPhusionHF-PCRMasterMix')
    dye=g.input('sybr_green','Source0.5µl15xSYBRGreenI actual molecular qPCR dye',K.MOLECULAR)
    product1=g.product('pcr1_product','Actual source linear PCR1 product class; full occupied products/C18fates unassembled',K.MOLECULAR,state=Construct([Segment('source PCR1 product','[unassembled occupied product]',placeholder=True,inferred=True)]))
    template_scene=Scene();template_scene.strand('template',[Segment('source circle','[spacer-bearing source circle]',placeholder=True,inferred=True)],label='full occupied PCR template unassembled')
    g.react('pcr1','Source real-time PCR1: programmed15cycles but stop at exponential-growth knee',S.PCR,((clean_circle,R.TEMPLATE),(p3tall,R.STRUCTURAL),(p6,R.STRUCTURAL),(nt1,R.STRUCTURAL),(enzyme1,R.CATALYST),(buf1,R.CONDITION),(dye,R.STRUCTURAL)),(product1,),geometry=PCRGeometry.unavailable(template_scene,reason='Source actual two-primer amplification retained; full cycle-dependent occupied C18 products and binding geometry unassembled. Analytical first-extension hypothesis is not substituted for whole PCR.'))
    selected1=clean('pcr1_size1',product1,'Actual first size selection:4.5µlisopropanol+54µlAMPure;10µlwater')
    selected2=clean('pcr1_size2',selected1,'Actual second size selection:1.5µlisopropanol+18µlAMPure; elution step follows')
    pcr2stock=g.input('pcr2_elution_primers','Source10µl500nMP3solexa/P6solexa oligo mix used as size-selection eluent',K.MOLECULAR)
    eluate2=g.product('pcr2_template','Source selected PCR1 template component in primer-containing eluent',K.MOLECULAR,state=product1.state)
    p3s=g.product('p3solexa','Source retained P3solexa component from actual eluent',K.OLIGO,state=Construct([Segment('P3solexa',P3SOLEXA)]))
    p6s=g.product('p6solexa','Source retained P6solexa component from actual eluent',K.OLIGO,state=Construct([Segment('P6solexa',P6SOLEXA)]))
    elution_conditions=g.input('pcr2_elution_conditions','Source bead elution in primer-containing mix; time/temperature unassigned',K.CONDITION)
    g.react('pcr2_elute','Elute selected product directly into source PCR2 primer mix; no later fresh primer addition',S.EXTRACTION,((selected2,R.TEMPLATE),(pcr2stock,R.STRUCTURAL),(elution_conditions,R.CONDITION)),(eluate2,p3s,p6s))
    enzyme2,nt2,buf2=mastermix('pcr2_mix','Source fresh10µl2xPhusionHFmastermix')
    product2=g.product('pcr2_product','Source three-cycle final PCR2 library class; whole occupied product unassembled',K.MOLECULAR,state=Construct([Segment('final library','[source final library class]',placeholder=True,inferred=True)]))
    g.react('pcr2','Actual PCR2 three cycles; source known primer polymers retained',S.PCR,((eluate2,R.TEMPLATE),(p3s,R.STRUCTURAL),(p6s,R.STRUCTURAL),(nt2,R.STRUCTURAL),(enzyme2,R.CATALYST),(buf2,R.CONDITION)),(product2,),geometry=PCRGeometry.unavailable(product1.state,reason='Source known primers and actual3cycle amplification; whole occupied sourcePCR1 template/primer sites and full product remain unassembled.'))
    library=clean('final_cleanup',product2,'Source36µlAMPure;2x70percentwash;20µlwater')
    library=qc('library_qc',library,'Source1–2µlfinal library for HS-DNA bioanalysis')
    loaded=select('loading',library,'Selected source sequencing library contribution; volume/run cycles/pooled sample allocation unassigned',kind=K.LIBRARY)
    g.sequence(loaded,SequencingPlan(limitation='Associated manuscript prints custom P6sequencing primer, but acquired read cycles, readend count, exact final binding geometry and sample loading/pooling remain unassigned in this model.',acquisition='unspecified'))
    g.annotate(Highlight('conditional_extend','Printed function versus inferred mechanism',
        extension.hypothesis.statement+' This is one conditional first-extension path, not full PCR validation.'))
    g.annotate(DetailView('cdna_circle','Spacer-bearing cDNA topology',
        'INFERRED cDNA copying with the source-written C18 linkers retained. The circle is opened only for display; no linker is encoded as a base.'))
    g.declare_reviewed_scope(('proteinase_k','recover_rna','rt_primer','tgirt','reverse_transcribe',
                             'streptavidin','capture_hybrid','p3short','heat_manganese','release_cdna',
                             'circligase','circularize','recover_circle','p6tall','p3tall'),
                            SOURCE+'; printed reagents/operations; derived backbone/landing mechanisms explicitly inferred')
    g.declare_reviewed_scope(('adapter_oh','phosphorylate_reagent','phosphorylate','adenylate_reagent',
                             'adenylate','click_reagent','click','rna_3oh','t4_rnl1','ligate'),SOURCE)
    g.annotate(SourceDiscrepancy('click', 'Click timing across inspected source versions',
        (SourceClaim('Accepted manuscript Online Methods JATS55', '2 h at 37 °C'),
         SourceClaim('Publisher Supplementary Figure1a, PDFp1', 'click 1 hr')),
        'This source-qualified model represents the inspected accepted-manuscript procedure (2 h); publisher Figure S1a prints 1 h.',
        'These are different source versions. The final published procedural text has not been reconciled.',
        'Revision differences or figure shorthand remain possible; no publication-error claim is established.'))
    g.declare_reviewed_scope(ADDITIONAL_SOURCE_INVENTORY,SOURCE+'; actual named supplies/fractions/downstreamMethods77; C18analytical reference explicitly inferred')
    return g.require_publishable()

def selftest():
    from checks import Check
    c=Check()
    try: Adapter().adenylate()
    except ValueError: rejected=True
    else: rejected=False
    c('unphosphorylated donor cannot be adenylated',rejected,True)
    prepared=Adapter().phosphorylate().adenylate().click()
    c('source adapter end/modification identity survives preparation',
      (prepared.five_prime,prepared.internal,prepared.three_prime),('App','IR800CW conjugate','biotin'))
    c('source-positive actual finish passes normal gate',protocol_graph().publication_problems(),[])
    g=protocol_graph();ports=lambda i:tuple(p.material for p in g.steps[i].inputs)
    c('click dye covalently contributes as molecular stock',g.materials['click_reagent'].kind,MaterialKind.MOLECULAR)
    c('click dye uses structural port',next(p.role for p in g.steps['click'].inputs if p.material=='click_reagent'),InputRole.STRUCTURAL)
    c('phosphorylation ATP source contributes','phosphorylate_atp' in ports('phosphorylate'),True)
    c('adenylation ATP source contributes','adenylate_atp' in ports('adenylate'),True)
    c('preadapter half stock retained',g.terminals['preadapter_half_remaining'].disposition,Disposition.UNRESOLVED)
    c('adapter selected from actual factory stock','adapter_ligation_aliquot_selected' in ports('ligate'),True)
    c('RT does not freshly add external primer','rt_primer' in ports('reverse_transcribe'),False)
    c('RT actual nucleotide contribution','rt_nucleotides' in ports('reverse_transcribe'),True)
    c('nucleotides not added in source enzyme preincubation','rt_nucleotides' in ports('rt_preincubation'),False)
    c('RT carries enzyme from actual prior addition','retained_tgirt' in ports('reverse_transcribe'),True)
    c('RT does not add enzyme twice','tgirt' in ports('reverse_transcribe'),False)
    c('conditional reference is nonphysical',g.steps['conditional_extend'].changes_molecule,False)
    c('actual PCR1 does not consume hypothetical extension','conditional_first_strand' in ports('pcr1'),False)
    c('actual PCR2 consumes retained eluent primer','p3solexa' in ports('pcr2'),True)
    c('PCR2 oligos originate in source elution step',g.steps['pcr2_elute'].outputs,('pcr2_template','p3solexa','p6solexa'))
    c('source library QC typed',g.materials['library_qc_qc'].kind,MaterialKind.QC)
    c('run acquisition not invented',next(iter(g.sequencing.values())).acquisition,'unspecified')
    from checks import Source,have
    from storage import sources_root
    import re
    source=Source(sources_root()/Path(__file__).resolve().parents[1].name/'nmeth.3840_PMC5477425.xml.txt','AssociatedAMMethods55/71/73')
    if have(c,source,label='irCLIP whole modified source orderings'):
        text=re.sub(r'\s+','',source.path.read_text()).lower()
        for name,seq in (('adapter',ADAPTER_WRITTEN),('RTbarcode1',RT_WRITTEN),('P3short',P3SHORT),('P3tall',P3TALL),('P6tall',P6TALL),('P3solexa',P3SOLEXA),('P6solexa',P6SOLEXA),('P6sequencing',P6SEQUENCING)):
            c('whole source '+name,seq.lower() in text,True)
    uncovered=protocol_graph();uncovered.input('unreviewed_extra','unreviewed source input',MaterialKind.CONDITION)
    c('new input not silently covered by curator inventory',any('source-inventory' in x and 'unreviewed_extra' in x for x in uncovered.publication_problems()),True)
    from selftest_spacers import run
    run(c)
    c.report()



ADDITIONAL_SOURCE_INVENTORY=('phosphorylate_atp','adenylate_atp','click_conditions','phospho_precipitation_conditions','phospho_precipitation','preadapter_half','preadapter_precipitation_conditions','preadapter_precipitation','adapter_columns_conditions','adapter_columns','adapter_ligation_aliquot','proteinase_qc','rna_extraction_conditions','rna_extraction','linear_polyacrylamide','rna_precipitation_conditions','rna_precipitation','rna_recovery_qc','rna_rt_aliquot','anneal_rt','rt_nucleotides','rt_preincubation','conditional_extend','pcr1_primer_stock','pcr1_primer_components','pcr1_mix_stock','pcr1_mix_components','sybr_green','pcr1','pcr1_size1_conditions','pcr1_size1','pcr1_size2_conditions','pcr1_size2','pcr2_elution_primers','pcr2_elution_conditions','pcr2_elute','pcr2_mix_stock','pcr2_mix_components','pcr2','final_cleanup_conditions','final_cleanup','library_qc','loading','loading_selected')

def render_page():
    g=protocol_graph()
    primer=rt_primer_backbone();circle=linear_cdna_backbone().close();first=conditional_first_extension(circle)
    return assemble_page(g,(head('irCLIP'),'<div class="wrap" data-protocol-complete="true"><h1>irCLIP</h1>',
      AnnotationSlot.INTRO,info('<a href="https://doi.org/10.1038/nmeth.3840">Zarnegar et al. 2016, Nature Methods</a>. This preparation is qualified to the associated accepted manuscript PMC5477425; full final-procedure equivalence has not been verified.'),
      '<p>A preadenylated, infrared-dye/biotin DNA adapter makes crosslinked RNA complexes visible and supports streptavidin capture after TGIRT reverse transcription. Source cDNA circularization and two PCR stages produce sequencing libraries. The graph starts from actual trimmed, dephosphorylated post-IP RNA; it does not invent an upstream RBP preparation.</p>',
      info('INFERRED — RT landing and the conditional C18-stop calculation are sequence-derived hypotheses. The whole occupied cDNA, PCR populations, linker readthrough and exact acquired run geometry remain unassembled. The source describes the actual PCR finish; it does not establish a unique nucleotide-by-nucleotide product or the displayed first-extension mechanism.'),
      '<h2>Activated adapter and known oligos</h2>',panel(Adapter().phosphorylate().adenylate().click().scene().rows(),caption='Source preadenylated adapter with covalently attached IR800CW and terminal biotin; iAzideN conjugation-site chemistry unassembled; no extra nucleotide invented.'),
      panel(primer.rows(),caption='Source barcode1 RT primer: both non-nucleotide C18 linkers retained; they are not DNA bases.'),
      '<p>P3short (elution): <seq>'+P3SHORT+'</seq><br>PCR1 P3tall: <seq>'+P3TALL+'</seq><br>PCR1 P6tall: <seq>'+P6TALL+'</seq><br>PCR2 P3solexa: <seq>'+P3SOLEXA+'</seq><br>PCR2 P6solexa: <seq>'+P6SOLEXA+'</seq></p>',
      '<h2>Read setup and identifiers</h2><p>Source custom P6 sequencing primer: <seq>'+P6SEQUENCING+'</seq>. Acquired read cycles, read-end count, final primer landing and sample pooling/loading amounts are unassigned in this model. ATCACG is the source barcode1 sample identifier; it is inline in the RT primer, not a separately asserted index read. The W3/N3/N5 composite molecular-identifier role is inferred from the source degenerate-barcode/duplicate-processing context and remains labelled as such.</p>',
      info('Source-version comparison: accepted-manuscript Methods prints the IR800CW–DBCO click step as 2 h at 37 °C; publisher Supplementary Figure S1a prints click 1 hr. This graph follows the inspected manuscript, retaining both locations. A revision difference or schematic shorthand is possible; no publication error or full final-version equivalence is asserted.'),
      '<p>PCR1 is monitored in real time and stopped at the exponential-growth knee within the programmed 15 cycles. Two bead size selections precede elution into the actual PCR2 primer mix, addition of fresh master mix and three PCR2 cycles. The conditional first-copy panel below is an analytical reference, not an extra physical aliquot or PCR reaction.</p>',
      '<details><summary>Conditional C18 topology reference</summary>',panel(circle.rows(),caption='INFERRED — circle opened only for display; both source C18 nodes remain.'),panel(first.product.rows(),caption='INFERRED analytical reference only — first P6 extension under explicit stop hypothesis; not the unique actual PCR1 product.'),'</details>',
      '<details><summary>Full selected workflow and unresolved molecular checks</summary>',WorkflowSection('Activated adapter → RNA capture → circular cDNA → actual two-PCR finish'),AnnotationSlot.AFTER_WORKFLOW,'</details>',
      '<p><a href="01_irclip.html">Source notes and limitations</a></p></div>'))
