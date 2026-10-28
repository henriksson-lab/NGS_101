"""Published native-DNA CaBagE branch; no quantitative nuclease prediction."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Construct, Scene, Segment, protocol_graph_panels
from adapter_ligation import AdapterLigationGeometry
from endprep import repair_and_dA_tail, dA_tailed_scene
from protocol_graph import InputRole as R, MaterialKind as K, StepKind as S, ProtocolGraph, SequencingPlan
from page import head, info

PAGE_FILE = 'cabage.html'
SOURCE = 'https://doi.org/10.1371/journal.pone.0241253; published Methods and Results overview'
GUIDES = {
    'c9orf72_for': 'TAACGTAGAATAGAACCCGA',
    'c9orf72_rev': 'GTTCCCTGGGTGTGAGTCAA',
}
TARGET = Construct([Segment('native region between flanking Cas9 cuts', 'X'*45, placeholder=True)],
                   name='CaBagE target (length and sequence are locus dependent)')

@dataclass(frozen=True)
class ProtectedTarget:
    """Local source-declared end occupancy, not an exonuclease simulator.

    No base-exact cut coordinate or protected footprint is implied by this class.
    """
    duplex: Construct
    bound_ends: tuple[str, str] = ('left', 'right')
    protein: str = 'guide-directed Cas9'
    pam_orientation: str = 'PAM-in'
    basis: str = 'Published Results, Fig. 1 and Methods Guide RNA design'

    def __post_init__(self):
        if self.bound_ends != ('left', 'right') or self.pam_orientation != 'PAM-in':
            raise ValueError('This local CaBagE declaration requires both ends occupied in PAM-in geometry')

    def scene(self):
        sc = Scene.duplex(list(self.duplex), label='native target')
        sc.note('top', 'Cas9 remains bound at both flanking cut sites; qualitative survival is source-declared')
        return sc

PROTECTED = ProtectedTarget(TARGET)
DA_TARGET = repair_and_dA_tail(TARGET)

def protocol_graph():
    g = ProtocolGraph('CaBagE — native-DNA branch', authority=SOURCE)
    dna = g.input('genomic_dna', 'Native genomic DNA', K.BIOLOGICAL, primary=True)
    repair = g.input('repair_mix', 'NEB FFPE Repair Mix M6630S; manufacturer protocol for user-supplied library reagents', K.CATALYST)
    repaired = g.product('repaired_dna', 'Nick-repaired genomic DNA; kit substrate effects not resolved', K.MOLECULAR)
    g.react('repair', 'Repair nicks before exonuclease challenge', S.END_REPAIR,
            ((dna,R.TEMPLATE),(repair,R.CATALYST)), (repaired,))
    guides = g.input('guides', 'Two locus-specific guides in PAM-in orientation (S1 Table)', K.OLIGO)
    cas = g.input('cas9', 'Cas9 RNP preassembled at 25°C for 10 min in CutSmart', K.CATALYST)
    target = g.product('protected_target', 'Cut target occupied by Cas9 at both flanks', K.MOLECULAR,
                       state=PROTECTED.scene(), chemical_state=PROTECTED)
    background = g.product('background', 'Other DNA fragments remain in the reaction', K.MOLECULAR)
    g.react('cas9_cut', 'Cas9 cleavage at 37°C for 15 min; retain bound complexes', S.DIGESTION,
            ((repaired,R.TEMPLATE),(guides,R.STRUCTURAL),(cas,R.CATALYST)), (target,background))
    exo = g.input('exonucleases', 'Exo I 40 U + Exo III 200 U + lambda 20 U', K.CATALYST)
    surviving = g.product('surviving_target', 'Qualitatively protected target after exonuclease challenge', K.MOLECULAR,
                          state=PROTECTED.scene(), chemical_state=PROTECTED)
    debris = g.product('digestion_products', 'Digested background and residual off-target material, still in tube', K.MOLECULAR)
    g.react('background_digestion', 'Incubate exonuclease mixture at 37°C for 2 h', S.DIGESTION,
            ((target,R.TEMPLATE),(background,R.TEMPLATE),(exo,R.CATALYST)), (surviving,debris))
    heat = g.input('heat', '80°C for 20 min: exonuclease inactivation and Cas9 dissociation', K.CONDITION)
    released = g.product('released_target', 'Exposed target duplex after heat release', K.MOLECULAR, state=TARGET)
    g.react('heat_release', 'Heat-release Cas9 and inactivate exonucleases', S.DEPROTEINIZATION,
            ((surviving,R.TEMPLATE),(heat,R.CONDITION)), (released,))
    taq = g.input('taq_datp', 'Taq DNA polymerase M0267S and dATP', K.CATALYST)
    tailed = g.product('a_tailed', 'dA-tailed target duplex', K.MOLECULAR, state=dA_tailed_scene(DA_TARGET))
    g.react('a_tail', 'A-tail at 72°C for 5 min', S.TAILING,
            ((released,R.TEMPLATE),(taq,R.CATALYST)), (tailed,))
    adapter = g.input('lsk109', 'SQK-LSK109 Adapter Mix; bases and adapter-end geometry not resolved in this pass', K.ADAPTER)
    ligase = g.input('ligase', 'Quick T4 ligase E6057 in LSK109 Ligation Buffer', K.CATALYST)
    ligated = g.product('ligated', 'LSK109-adapted target plus remaining reaction material', K.MOLECULAR)
    geometry = AdapterLigationGeometry.unavailable(
        'Published Methods names LSK109 Adapter Mix and dA-tailing; its adapter-end geometry was not established by the inspected originals. Historical kit guide remains deferred.')
    g.react('ligation', 'Ligate adapters for 10 min at room temperature', S.LIGATION,
            ((tailed,R.TEMPLATE),(adapter,R.STRUCTURAL),(ligase,R.CATALYST)), (ligated,), geometry=geometry)
    beads = g.input('beads', '0.3× AMPure XP, two 200 µL Short Fragment Buffer washes; Elution Buffer', K.SELECTION)
    library = g.product('library', 'Recovered native CaBagE library, including possible residual off-targets', K.LIBRARY)
    waste = g.product('cleanup_waste', 'Unretained background products, free adapters and reaction reagents', K.WASTE)
    g.react('cleanup', 'Bead recovery, SFB washes and elution', S.CLEANUP,
            ((ligated,R.TEMPLATE),(debris,R.TEMPLATE),(beads,R.SELECTION)), (library,waste))
    g.discard(waste, 'Removed only at the bead cleanup; no pre-cleanup physical separation is declared')
    g.sequence(library, SequencingPlan(acquired_roles=('Sequencing output',),
        read_uses={'Sequencing output':'include'},
        limitation='R9.4.1 MinION, 48 h; primer-free native-DNA reads. LSK109 adapter sequence/end geometry and per-molecule read lengths remain unresolved.'))
    g.declare_reviewed_scope((
        'genomic_dna','repair_mix','repaired_dna','repair','guides','cas9',
        'protected_target','background','cas9_cut','exonucleases','surviving_target',
        'digestion_products','background_digestion','heat','released_target','heat_release',
        'taq_datp','a_tailed','a_tail','lsk109','ligase','ligated','ligation','beads',
        'library','cleanup_waste','cleanup'), SOURCE)
    return g

def render_page():
    return '\n'.join([head('CaBagE — bounded chemistry model'), '<div class="wrap"><h1>CaBagE</h1>',
        '<p class="research-notes"><a href="01_cabage.html">Research notes</a></p>',
        info('Bounded model of the published native-DNA branch. Cas9 protection is a source-declared qualitative outcome, not computed resistance or perfect purification. Historical kit geometry remains unresolved.'),
        protocol_graph_panels(protocol_graph(), cls='long'), '</div>'])

def selftest():
    from checks import Check
    c = Check()
    g = protocol_graph()
    c('Protected target does not disappear during digestion', g.materials['surviving_target'].chemical_state, PROTECTED)
    c('Digestion products reach cleanup before discard', g._consumers['digestion_products'], ['cleanup'])
    c('No dephosphorylation imported from nCATS', any(s.kind is S.DEPHOSPHORYLATION for s in g.steps.values()), False)
    c('Adapter geometry remains unresolved', g.steps['ligation'].geometry.availability, 'unavailable')
    c('Graph publication contract for represented branch', g.publication_problems(), [])
    from storage import sources_root
    source = sources_root() / Path(__file__).resolve().parents[1].name / 'journal.pone.0241253_supp_pone.0241253.s004.xlsx.txt'
    if source.exists():
        text = source.read_text()
        for name, seq in GUIDES.items():
            c('S1 Table transcription: '+name, name+'\t'+seq in text)
    else:
        c.skip('S1 Table source transcription', 'External published workbook text unavailable')
    c.report()
