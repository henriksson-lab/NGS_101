"""Luo2018 snmC-seq2: native RP-H, actual ss Adaptase substrate, printed index pair."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Construct, Scene, Segment, feature, oligo, revcomp
from protocol_model import role
from pcr import PCRGeometry
from protocol_graph import (ProtocolGraph, MaterialKind, InputRole, StepKind,
                            SequencingPlan, Terminal, Disposition)
import seqprimers as sp

TITLE = 'snmC-seq2 — indexed single-nucleus methylomes'
NOTES = '01_snmc-seq2.html'
SOURCE = ('Luo et al., Nature Communications2018, '
          '<a href="https://doi.org/10.1038/s41467-018-06355-2">doi:10.1038/s41467-018-06355-2</a>; '
          'published PMC6147798 Methods141–153 and original MOESM1 SupplementaryMethods/Table2.')
SUMMARY = ('Converted DNA is copied from source RP-H primers and treated with ExoI/rSAP. '
           'Eight indexed reactions are pooled and recovered before Adaptase and15-cycle PCR '
           'installs dual indexes. Three final SPRI cleanups are distinct operations.')
CAVEAT = ('AD002 and Table2 index pair1 are source-printed representatives, not assigned to every '
          'experimental nucleus. Adaptase-derived span and native P7 priming geometry remain '
          'unassigned; exact printed PCR oligos do not prove a completed duplex precursor. '
          'The source10bp index definitions determine the primer subregions; no generic8bp '
          'TruSeq shell or canonical P7 boundary is substituted. Acquisition cycles remain unassigned.')
INLINE = 'CGATGT'
RP_HANDLE = 'TTCCCTACACGACGCTCTTCCGATCT'
P5_INDEX_1 = 'AATGATACGGCGACCACCGAGATCTACACACGATCAGACACTCTTTCCCTACACGACGCTCT'
P7_INDEX_1 = 'CAAGCAGAAGACGGCATACGAGATAGGCAATGGTGACTGGAGTTCAGACGTGTGCTCTT'
P5_INDEX_1_CODE = 'ACGATCAGAC'
P7_INDEX_1_CODE = 'CATTGCCTAT'
INLINE_FEATURE = feature('snmc2_inline', 'cell_barcode', 'whitelist',
                          whitelist='Luo2018 published eight RP-H barcodes', group='cell_identity',
                          note='Together with dual library-pool indexes identifies a nucleus.')
I5_FEATURE = feature('snmc2_i5', 'sample_index', 'fixed', group='cell_identity',
                     note='Published Table2 pair1 representative;10bp source code, not acquisition cycles.')
I7_FEATURE = feature('snmc2_i7', 'sample_index', 'fixed', group='cell_identity',
                     note='Published Table2 pair1 representative;10bp source code, not acquisition cycles.')


def random_parts():
    return [Segment('RP-H partial Read1 handle', RP_HANDLE),
            Segment('inline cell barcode', INLINE, 'cbc', feature=INLINE_FEATURE),
            Segment('H9 A/C/T random priming tract', 'H' * 9, placeholder=True)]


def _unique_start(sequence, word):
    if sequence.count(word) != 1:
        raise ValueError('Printed index must have one unique position in its printed primer')
    return sequence.index(word)


def library():
    i5 = _unique_start(P5_INDEX_1, P5_INDEX_1_CODE)
    i7 = _unique_start(P7_INDEX_1, revcomp(P7_INDEX_1_CODE))
    overlaps = [n for n in range(1, len(RP_HANDLE) + 1)
                if P5_INDEX_1.endswith(RP_HANDLE[:n])]
    if not overlaps:
        raise ValueError('Printed P5 primer must overlap source RP-H handle')
    overlap = max(overlaps)
    return Construct([
        Segment('P5 primer upstream region', P5_INDEX_1[:i5]),
        Segment('source i5 code', P5_INDEX_1_CODE, feature=I5_FEATURE),
        Segment('P5 primer downstream region', P5_INDEX_1[i5 + len(P5_INDEX_1_CODE):]),
        Segment('remaining RP-H handle', RP_HANDLE[overlap:]),
        *random_parts()[1:], role('copied genomic or lambda DNA insert'),
        role('additional Adaptase-dependent span outside P7 footprint; extent unassigned'),
        Segment('PCR-installed P7 annealing footprint', revcomp(P7_INDEX_1[i7 + len(P7_INDEX_1_CODE):])),
        Segment('source i7 code', P7_INDEX_1_CODE, feature=I7_FEATURE),
        Segment('P7 primer upstream region reverse complement', revcomp(P7_INDEX_1[:i7])),
    ], name='snmC-seq2 representative source-primer endpoint'), ()


FINAL_LIBRARY, SEQ_PRIMERS = library()
FINAL_CAPTION = 'Expected representative PCR endpoint from printed primer ends; native Adaptase junction and input-to-product derivation remain unvalidated.'
SEQUENCING_INTRO = 'Published NovaSeq platform. Inline RP-H barcode and dual source10bp codes are molecular regions; acquired read/index cycle allocation remains unassigned.'
SEQUENCING_ENDING = SEQUENCING_INTRO
SHOW_READ_CONTENT = False
READ_SPECS = tuple(sp.read_spec(name, 'unassigned', 'Read/index acquisition allocation unassigned; source barcode/index regions retained', 'undetermined')
                   for name in ('Read 1', 'Read 2', 'Index 1 (i7)', 'Index 2 (i5)'))
SOURCE_OLIGO_EXCEPTIONS = {'random_primer': 'Original published SupplementaryMethods prints /5SpC3/ followed by the exact handle/barcode and nine H1 mix specifications (first H1:33340033); main Results explicitly defines RP-H as A/C/T. The source guard does not expand this custom synthesis notation to H9.'}


def oligos():
    yield oligo('P5L_AD002_H', random_parts(), mods='/5SpC3/ before the printed sequence')
    yield oligo('Table2 P5 primer1', [Segment('P5 primer', P5_INDEX_1)])
    yield oligo('Table2 P7 primer1', [Segment('P7 primer', P7_INDEX_1)])


def converted_state():
    sc = Scene()
    for identity in ('nuclear genomic DNA', 'unmethylated lambda DNA spike-in'):
        for strand in (1, 2):
            sc.strand(identity + str(strand), [role(identity + ' converted strand' + str(strand) + ': physical U at converted C; protected methyl-C retained')], label=identity + '; independent original DNA strand')
    return sc


def primed_state(*, adapted=False, denatured=False):
    sc = Scene()
    for identity in ('genomic DNA-derived copy', 'lambda DNA-derived copy'):
        parts = [*random_parts(), role(identity + '; copying converted U installs A')]
        if adapted:
            parts.append(role('native Adaptase-added3′ end; sequence/extent unassigned'))
        label = identity + ('; single-strand Adaptase precursor' if adapted or denatured else '; original/copy duplex occupancy unassigned before denaturation')
        sc.strand(identity, parts, label=label, mod5='SpC3')
        sc.strand(identity + ' original', [role('original converted template beside ' + identity + '; physical U remains; later fate unassigned')], label='independent original DNA template population')
    return sc


def post_pcr_state():
 sc=Scene.duplex(list(FINAL_LIBRARY),label='expected amplified library population',unpaired=tuple(p.name+"'" for p in FINAL_LIBRARY if p.is_role_token()))
 sc.strand('residual0',[role('possible co-carried original genomic DNA; U/C chemistry/integrity/partition unassigned')],label='co-carried molecular lineage; no physical separation asserted')
 sc.strand('residual1',[role('possible co-carried original lambda DNA; U chemistry/integrity/partition unassigned')],label='co-carried molecular lineage; no physical separation asserted')
 sc.strand('residual2',[role('possible co-carried5′C3 RP-H/Adaptase genomic-copy templates; integrity/partition unassigned')],label='co-carried molecular lineage; no physical separation asserted')
 sc.strand('residual3',[role('possible co-carried5′C3 RP-H/Adaptase lambda-copy templates; integrity/partition unassigned')],label='co-carried molecular lineage; no physical separation asserted')
 return sc

POST_PCR_STATE=post_pcr_state()

def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    nucleus = g.input('nucleus', 'single nucleus sorted to384-well plate', MaterialKind.BIOLOGICAL, primary=True)
    lambda_dna = g.input('lambda', 'unmethylated lambda DNA in digestion buffer', MaterialKind.MOLECULAR)
    proteinase = g.input('proteinase', 'ProteinaseK/M-digestion buffer', MaterialKind.CATALYST)
    released = g.product('released', 'released nuclear DNA beside lambda spike-in', MaterialKind.MOLECULAR)
    g.react('release_dna', 'Digest sorted nucleus in source digestion buffer', StepKind.DEPROTEINIZATION,
            ((nucleus, InputRole.TEMPLATE), (lambda_dna, InputRole.SPIKE_IN), (proteinase, InputRole.CATALYST)), (released,))
    bs = g.input('bisulfite', '15µl CT conversion reagent;98°C8min,64°C3.5h', MaterialKind.CONDITION)
    converted = g.product('converted', 'bisulfite-converted original DNA', MaterialKind.MOLECULAR, state=converted_state())
    g.react('convert', 'Bisulfite-convert original genomic and lambda DNA', StepKind.CHEMICAL_CONVERSION,
            ((released, InputRole.TEMPLATE), (bs, InputRole.CONDITION)), (converted,))
    binding = g.input('binding', '80µl M-binding buffer; Zymo384-well binding plate', MaterialKind.SELECTION)
    cleaned = g.product('cleaned', 'desulfonated/washed DNA eluted in7µl primer solution', MaterialKind.MOLECULAR, state=converted.state)
    bs_waste = g.product('bs_waste', 'binding/desulfonation/wash non-retained fractions', MaterialKind.WASTE)
    g.react('recover_bs', 'Bind, wash, desulfonate15min, wash twice and elute converted DNA', StepKind.CLEANUP,
            ((converted, InputRole.TEMPLATE), (binding, InputRole.SELECTION)), (cleaned, bs_waste))
    g.discard(bs_waste, 'Source plate-cleanup non-retained fractions; species composition unassigned.')
    primer_state = Scene(); primer_state.strand('primer', random_parts(), label='representative RP-H AD002', mod5='SpC3')
    primer = g.input('random_primer', 'source P5L_AD002_H representative', MaterialKind.OLIGO, state=primer_state)
    heat = g.input('denature', '95°C3min then ice2min', MaterialKind.CONDITION)
    annealed = g.product('annealed', 'denatured converted DNA with RP-H', MaterialKind.MOLECULAR, state=cleaned.state)
    g.react('denature_dna', 'Denature converted DNA with indexed RP-H', StepKind.THERMAL_TREATMENT,
            ((cleaned, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL), (heat, InputRole.CONDITION)), (annealed,))
    polymerase = g.input('strand_polymerase', 'Klenow exo-minus random-priming mix', MaterialKind.CATALYST)
    primed = g.product('primed', 'RP-H-initiated DNA copies', MaterialKind.MOLECULAR, state=primed_state())
    g.react('random_prime', 'Random extension:4°C5min,25°C5min,37°C60min', StepKind.EXTENSION,
            ((annealed, InputRole.TEMPLATE), (polymerase, InputRole.CATALYST)), (primed,))
    exo = g.input('exo', 'ExonucleaseI', MaterialKind.CATALYST)
    rsap = g.input('rsap', 'rSAP', MaterialKind.CATALYST)
    treated = g.product('treated', 'treated expected RP-H copies', MaterialKind.MOLECULAR, state=primed.state)
    chemistry_products = g.product('chemistry_products', 'ExoI hydrolysis/rSAP reaction products; recovery unassigned', MaterialKind.MOLECULAR)
    g.react('exo_rsap', 'ExoI+rSAP treatment:37°C30min', StepKind.DIGESTION,
            ((primed, InputRole.TEMPLATE), (exo, InputRole.CATALYST), (rsap, InputRole.CATALYST)), (treated, chemistry_products))
    g.mark_terminal(Terminal(chemistry_products.id, Disposition.UNRESOLVED, reason='Chemical inactivation/hydrolysis is distinct from physical SPRI recovery.'))
    others = g.input('other_seven', 'seven independently inline-barcoded reactions from the two384-well plates', MaterialKind.MOLECULAR)
    pooled = g.product('pooled', 'eight reaction pool; selected AD002 lineage shown', MaterialKind.MOLECULAR, state=treated.state)
    g.react('pool_eight', 'Pool eight quadrants with distinct RP-H inline barcodes', StepKind.POOL,
            ((treated, InputRole.TEMPLATE), (others, InputRole.TEMPLATE)), (pooled,))
    spri = g.input('spri', 'SPRI beads', MaterialKind.SELECTION)
    purified = g.product('purified', 'pooled copies recovered in10µl EB', MaterialKind.MOLECULAR, state=pooled.state)
    waste = g.product('priming_waste', '0.8×SPRI non-retained fraction', MaterialKind.WASTE)
    g.react('purify_primed', 'Pool cleanup0.8×SPRI; three80%EtOH washes;10µlEB', StepKind.CLEANUP,
            ((pooled, InputRole.TEMPLATE), (spri, InputRole.SELECTION)), (purified, waste))
    g.discard(waste, 'Source SPRI non-retained fraction; exact molecular depletion unassigned.')
    adapt_heat = g.input('adapt_heat', '95°C3min then ice2min', MaterialKind.CONDITION)
    denatured = g.product('denatured', 'denatured RP-H copies for Adaptase', MaterialKind.MOLECULAR, state=primed_state(denatured=True))
    g.react('denature_copies', 'Denature recovered copies before Adaptase', StepKind.THERMAL_TREATMENT,
            ((purified, InputRole.TEMPLATE), (adapt_heat, InputRole.CONDITION)), (denatured,))
    adapter = g.input('adaptase_end', 'native Adaptase completion end; bases unassigned', MaterialKind.ADAPTER)
    adaptase = g.input('adaptase', 'Swift Adaptase module', MaterialKind.CATALYST)
    adapted = g.product('adapted', 'single-strand Adaptase-completed copies', MaterialKind.MOLECULAR, state=primed_state(adapted=True))
    g.react('adaptase_complete', 'Adaptase37°C30min then95°C2min', StepKind.ADAPTER_ATTACHMENT,
            ((denatured, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL), (adaptase, InputRole.CATALYST)), (adapted,))
    p5 = g.input('p5', 'printed Table2 P5 primer1', MaterialKind.OLIGO, state=Construct([Segment('P5 primer1', P5_INDEX_1)]))
    p7 = g.input('p7', 'printed Table2 P7 primer1', MaterialKind.OLIGO, state=Construct([Segment('P7 primer1', P7_INDEX_1)]))
    kapa = g.input('kapa', 'KAPA HiFi HotStart ReadyMix', MaterialKind.CATALYST)
    amplified = g.product('amplified', '15-cycle indexed pool; representative lineage endpoint', MaterialKind.LIBRARY, state=POST_PCR_STATE)
    g.react('index_pcr', 'Amplify pooled libraries15cycles with printed dual-index pair', StepKind.PCR,
            ((adapted, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL), (p7, InputRole.STRUCTURAL), (kapa, InputRole.CATALYST)), (amplified,),
            geometry=PCRGeometry.unavailable(adapted.state, reason='Actual blocked single-strand Adaptase product; P7 target and complete native geometry unassigned. Printed primer ends do not establish an input-site proof.'))
    other_pools = g.input('other_pools', 'other uniquely dual-indexed amplified pools for source plate consolidation', MaterialKind.MOLECULAR)
    combined = g.product('combined', 'combined amplified libraries; representative lineage shown', MaterialKind.LIBRARY, state=POST_PCR_STATE)
    g.react('combine_libraries', 'Consolidate indexed source plates during final recovery', StepKind.POOL,
            ((amplified, InputRole.TEMPLATE), (other_pools, InputRole.TEMPLATE)), (combined,))
    clean1 = g.product('clean1', 'first0.8×SPRI recovered pool', MaterialKind.LIBRARY, state=POST_PCR_STATE)
    waste1 = g.product('waste1', 'first library-cleanup non-retained fraction', MaterialKind.WASTE)
    g.react('cleanup1', 'First final0.8×SPRI cleanup', StepKind.CLEANUP,
            ((combined, InputRole.TEMPLATE), (spri, InputRole.SELECTION)), (clean1, waste1)); g.discard(waste1, 'Physical cleanup fraction.')
    clean2 = g.product('clean2', 'second0.8×SPRI recovered pool', MaterialKind.LIBRARY, state=POST_PCR_STATE)
    waste2 = g.product('waste2', 'second library-cleanup non-retained fraction', MaterialKind.WASTE)
    g.react('cleanup2', 'Second final0.8×SPRI cleanup', StepKind.CLEANUP,
            ((clean1, InputRole.TEMPLATE), (spri, InputRole.SELECTION)), (clean2, waste2)); g.discard(waste2, 'Physical cleanup fraction.')
    library_out = g.product('library', FINAL_LIBRARY.name, MaterialKind.LIBRARY, state=POST_PCR_STATE)
    waste3 = g.product('waste3', 'third library-cleanup non-retained fraction', MaterialKind.WASTE)
    g.react('cleanup3', 'Third final0.8×SPRI cleanup; Qubit QC', StepKind.CLEANUP,
            ((clean2, InputRole.TEMPLATE), (spri, InputRole.SELECTION)), (library_out, waste3)); g.discard(waste3, 'Physical cleanup fraction.')
    g.sequence(library_out, SequencingPlan(library=FINAL_LIBRARY, primers=(), limitation='NovaSeq source-stated; acquired cycle allocation and native Adaptase junction unassigned.'))
    g.declare_reviewed_scope((
        'nucleus','lambda','proteinase','release_dna','released','bisulfite','convert','converted',
        'binding','recover_bs','cleaned','random_primer','denature','denature_dna','annealed',
        'strand_polymerase','random_prime','primed','exo','rsap','exo_rsap','treated','other_seven',
        'pool_eight','pooled','spri','purify_primed','purified','adapt_heat','denature_copies','denatured',
        'adaptase_end','adaptase','adaptase_complete','adapted','p5','p7','kapa','index_pcr','amplified',
        'other_pools','combine_libraries','combined','cleanup1','clean1','cleanup2','clean2','cleanup3','library'), SOURCE)
    return g.require_publishable()
