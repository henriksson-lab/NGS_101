"""Source-qualified irCLIP adapter preparation and single-strand ligation."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Construct, Segment, Scene, feature, revcomp
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import ProtocolGraph, MaterialKind, InputRole, StepKind, Terminal, Disposition, SourceDiscrepancy, SourceClaim, DetailView, Highlight
from illumina import PE_PCR_PRIMER_2
from circular import CircularizationGeometry
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
          'polymerase-stop hypothesis. Full two-primer PCR, acquired read geometry and '
          'final published procedure applicability remain unresolved.')


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
        sc.strand('adapter', [Segment('L3 DNA arm',ARM),Segment('A12 upstream','A'*12),
                  Segment('internal '+self.internal,'[modified nucleotide]',placeholder=True,length_bp=1),
                  Segment('A12 downstream','A'*12)], label='irCLIP adapter',
                  mod5=self.five_prime,mod3=self.three_prime)
        sc.note('adapter','Internal /iAzideN/ chemistry is retained; no canonical base identity is invented.')
        sc.ligation_end('adapter','L3 DNA arm',kind='single-strand',polarity='5-prime',owner='adapter')
        return sc

def protocol_graph():
    g=ProtocolGraph('irCLIP — adapter activation and RNA ligation',authority=SOURCE,
                    representation='published_instance')
    state=Adapter()
    stock=g.input('adapter_oh',ADAPTER_WRITTEN,MaterialKind.ADAPTER,state=tuple(state.scene().rows()),chemical_state=state)
    prior=stock
    for ident,label,enzyme,op,kind in (
        ('phosphorylate','OptiKinase + ATP (37 °C, 1 h)','OptiKinase',Adapter.phosphorylate,StepKind.PHOSPHORYLATION),
        ('adenylate','Mth RNA ligase + ATP (65 °C, 2 h)','Mth RNA ligase',Adapter.adenylate,StepKind.CHEMICAL_CONVERSION),
        ('click','IR800CW–DBCO conjugation (Methods:37 °C, 2 h; FigS1a:1 h)','IR800CW–DBCO',Adapter.click,StepKind.CHEMICAL_CONVERSION)):
        reagent_kind=MaterialKind.CONDITION if ident=='click' else MaterialKind.CATALYST
        reagent=g.input(ident+'_reagent',enzyme,reagent_kind)
        role=InputRole.CONDITION if ident=='click' else InputRole.CATALYST
        state=op(state)
        product=g.product(ident+'_adapter',label,MaterialKind.ADAPTER,state=tuple(state.scene().rows()),chemical_state=state)
        g.react(ident,label,kind,((prior,InputRole.TEMPLATE),(reagent,role)),(product,))
        prior=product
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
    primer=rt_primer_backbone()
    primer_input=g.input('rt_primer','manuscript barcode1 RT primer: phosphorylated DNA with two C18 links',
                         MaterialKind.OLIGO,state=tuple(primer.rows()),chemical_state=primer)
    tgirt=g.input('tgirt','TGIRT-III',MaterialKind.CATALYST)
    linear=linear_cdna_backbone()
    hybrid=g.product('rt_hybrid','INFERRED adapter-arm landing; RNA–cDNA hybrid with spacer-bearing cDNA',
                     MaterialKind.MOLECULAR,state=tuple(linear.rows()),
                     chemical_state=SpacerMixture(linear,('RNA template with biotin adapter retained',)),inferred=True)
    g.react('reverse_transcribe','TGIRT-III RT (60 °C, 2 h); exact landing/copy inferred from printed arm',
            StepKind.REVERSE_TRANSCRIPTION,((recovered,InputRole.TEMPLATE),(primer_input,InputRole.STRUCTURAL),
                                           (tgirt,InputRole.CATALYST)),(hybrid,))
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
    p6=g.input('p6tall','P6tall PCR1 primer',MaterialKind.OLIGO,state=Construct([Segment('P6tall',P6TALL)]))
    p3tall=g.input('p3tall','P3tall PCR1 primer (other copying paths unresolved)',MaterialKind.OLIGO,
                   state=Construct([Segment('P3tall',P3TALL)]))
    polymerase=g.input('phusion','Phusion HF',MaterialKind.CATALYST)
    extension=conditional_first_extension(circle)
    first=g.product('conditional_first_strand','INFERRED one P6-initiated linear extension, stopping at second C18',
                    MaterialKind.MOLECULAR,state=tuple(extension.product.rows()),chemical_state=extension,inferred=True,
                    source=SOURCE+'; conditional C18-stop hypothesis, not a printed product sequence')
    g.react('conditional_extend','INFERRED P6 first extension: cross circle closure, copy RNA-derived DNA, stop before C18',
            StepKind.EXTENSION,((clean_circle,InputRole.TEMPLATE),(p6,InputRole.STRUCTURAL),
                                (p3tall,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(first,),
            source=SOURCE+'; INFERRED conditional mechanism only; other PCR copying/product fates unresolved')
    g.mark_terminal(Terminal(first.id,Disposition.UNRESOLVED,BOUNDARY))
    g.annotate(Highlight('conditional_extend','Printed function versus inferred mechanism',
        extension.hypothesis.statement+' This is one conditional first-extension path, not full PCR validation.'))
    g.annotate(DetailView('cdna_circle','Spacer-bearing cDNA topology',
        'INFERRED cDNA copying with the source-written C18 linkers retained. The circle is opened only for display; no linker is encoded as a base.'))
    g.declare_reviewed_scope(('proteinase_k','recover_rna','rt_primer','tgirt','reverse_transcribe',
                             'streptavidin','capture_hybrid','p3short','heat_manganese','release_cdna',
                             'circligase','circularize','recover_circle','p6tall','p3tall','phusion'),
                            SOURCE+'; printed reagents/operations; derived backbone/landing mechanisms explicitly inferred')
    g.declare_reviewed_scope(('adapter_oh','phosphorylate_reagent','phosphorylate','adenylate_reagent',
                             'adenylate','click_reagent','click','rna_3oh','t4_rnl1','ligate'),SOURCE)
    g.annotate(SourceDiscrepancy('click', 'Click timing across inspected source versions',
        (SourceClaim('Accepted manuscript Online Methods JATS55', '2 h at 37 °C'),
         SourceClaim('Publisher Supplementary Figure1a, PDFp1', 'click 1 hr')),
        'This private partial model represents the inspected accepted-manuscript procedure (2 h); publisher Figure S1a prints 1 h.',
        'These are different source versions. The final published procedural text has not been reconciled.',
        'Revision differences or figure shorthand remain possible; no publication-error claim is established.'))
    return g

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
    c('spacer-dependent downstream library stays unpublished',protocol_graph().publishable,False)
    from selftest_spacers import run
    run(c)
    c.report()

if __name__ == '__main__': selftest()
