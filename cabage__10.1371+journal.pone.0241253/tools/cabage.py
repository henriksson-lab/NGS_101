"""Published native-DNA CaBagE branch; no quantitative nuclease prediction."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Construct, Scene, Segment, panel
from adapter_ligation import AdapterLigationGeometry
from endprep import repair_and_dA_tail, dA_tailed_scene
from protocol_graph import InputRole as R, MaterialKind as K, StepKind as S, ProtocolGraph, SequencingPlan, Terminal, Disposition
from page import head, info

PAGE_FILE = 'cabage.html'
SOURCE = 'https://doi.org/10.1371/journal.pone.0241253; published Methods and Results overview'
GUIDES = {
    'c9orf72_for': 'TAACGTAGAATAGAACCCGA',
    'c9orf72_rev': 'GTTCCCTGGGTGTGAGTCAA',
}
TARGET = Construct([Segment('native region between flanking Cas9 cuts', 'X'*45, placeholder=True, inferred=True, note='Display width only; not a forty-five-base target or a genomic coordinate claim.')],
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
    g = ProtocolGraph('CaBagE — native-DNA branch', authority=SOURCE, representation='published_instance')
    dna = g.input('genomic_dna', 'Native genomic DNA', K.BIOLOGICAL, primary=True)
    repair = g.input('repair_mix', 'NEB FFPE Repair Mix M6630S; manufacturer protocol for user-supplied library reagents', K.CATALYST)
    repaired = g.product('repaired_dna', 'Nick-repaired genomic DNA; kit substrate effects not resolved', K.MOLECULAR)
    g.react('repair', 'Repair nicks before exonuclease challenge', S.END_REPAIR,
            ((dna,R.TEMPLATE),(repair,R.CATALYST)), (repaired,))
    guides = g.input('guides', 'C9 reference guide pair; source spacer sequences known, full crRNA/tracrRNA constructs unassembled', K.OLIGO, chemical_state={'source_spacers':GUIDES,'whole_guide_polymers':'unassembled; no scaffold/end-modification import'})
    cas = g.input('cas9', 'Cas9 protein supply for RNP preassembly', K.MOLECULAR)
    rnp = g.product('rnp', 'Source Cas9/guide RNP collection; target pairing/occupancy quantitative details unassembled', K.CATALYST)
    preassembly = g.input('rnp_conditions', '25°C10min in1xCutSmart; source23.5uL preassembly', K.CONDITION)
    g.react('assemble_rnp', 'Preassemble the source guide/Cas9 RNPs', S.BINDING,
            ((cas,R.TEMPLATE),(guides,R.STRUCTURAL),(preassembly,R.CONDITION)), (rnp,))
    target = g.product('protected_target', 'Cut target occupied by Cas9 at both flanks', K.MOLECULAR,
                       state=PROTECTED.scene(), chemical_state=PROTECTED)
    background = g.product('background', 'Other DNA fragments remain in the reaction', K.MOLECULAR)
    g.react('cas9_cut', 'INFERRED target geometry — Cas9 cleavage at 37°C for 15 min; retain bound complexes', S.DIGESTION,
            ((repaired,R.TEMPLATE),(rnp,R.CATALYST)), (target,background))
    exo = g.input('exonucleases', 'Exo I 40 U + Exo III 200 U + lambda 20 U', K.CATALYST)
    surviving = g.product('surviving_target', 'Qualitatively protected target after exonuclease challenge', K.MOLECULAR,
                          state=PROTECTED.scene(), chemical_state=PROTECTED)
    debris = g.product('digestion_products', 'Digested background and residual off-target material, still in tube', K.MOLECULAR)
    g.react('background_digestion', 'INFERRED target geometry — Incubate exonuclease mixture at 37°C for 2 h', S.DIGESTION,
            ((target,R.TEMPLATE),(background,R.TEMPLATE),(exo,R.CATALYST)), (surviving,debris))
    heat = g.input('heat', '80°C for 20 min: exonuclease inactivation and Cas9 dissociation', K.CONDITION)
    released = g.product('released_target', 'Exposed target duplex after heat release', K.MOLECULAR, state=TARGET)
    protein_remnant = g.product('protein_remnant', 'Released/inactivated protein material remains in tube; composition and carryover unassigned', K.MOLECULAR)
    g.react('heat_release', 'INFERRED target geometry — Heat-release Cas9 and inactivate exonucleases', S.DEPROTEINIZATION,
            ((surviving,R.TEMPLATE),(heat,R.CONDITION)), (released,protein_remnant))
    taq = g.input('taq_datp', 'Taq DNA polymerase M0267S,1uL', K.CATALYST)
    datp = g.input('datp', 'Source1uL10mm dATP; printed concentration wording retained', K.MOLECULAR)
    tailed = g.product('a_tailed', 'dA-tailed target duplex', K.MOLECULAR, state=dA_tailed_scene(DA_TARGET))
    g.react('a_tail', 'INFERRED target geometry — A-tail at 72°C for 5 min', S.TAILING,
            ((released,R.TEMPLATE),(taq,R.CATALYST),(datp,R.STRUCTURAL)), (tailed,))
    adapter_stock = g.input('lsk109', 'Source5uL SQK-LSK109 Adapter Mix; whole adapter and motor packaging geometry unassembled', K.MOLECULAR)
    adapter = g.product('adapter_component', 'Adapter component of the actual source supply; whole polymer/end geometry unknown', K.ADAPTER)
    motor = g.product('motor_component', 'SourceFig1 motor-protein class associated with adapter supply; identity/occupancy unassigned', K.CATALYST)
    g.react('adapter_components', 'Expose source motor/adapter inventory', S.MAGIC,
            ((adapter_stock,R.TEMPLATE),),(adapter,motor),changes_molecule=False,
            reason='Figure1 depicts motor protein plus adapter and Methods names the source Adapter Mix. Inventory views of the same supplied complex, not additional free motor or a claimed historical kit recipe; molecular packaging and occupancy unassembled.')
    ligase = g.input('ligase', 'Quick T4 ligase E6057 in LSK109 Ligation Buffer', K.CATALYST)
    ligated = g.product('ligated', 'LSK109-adapted target plus remaining reaction material', K.MOLECULAR)
    geometry = AdapterLigationGeometry('unavailable', 'unavailable', (dA_tailed_scene(DA_TARGET),),
        'INFERRED genomic context — source-known dA-tailed substrate reference only; whole LSK109 adapter and occupied junction remain unassembled. Historical kit guide is deferred, not confirmed to omit geometry.')
    g.react('ligation', 'Ligate adapters for 10 min at room temperature', S.LIGATION,
            ((tailed,R.TEMPLATE),(adapter,R.STRUCTURAL),(ligase,R.CATALYST),(motor,R.CATALYST)), (ligated,), geometry=geometry)
    beads = g.input('beads', '0.3× AMPure XP, two 200 µL Short Fragment Buffer washes; Elution Buffer', K.SELECTION)
    library = g.product('library', 'Native CaBagE eluate16.6uL; possible residual off-targets', K.LIBRARY)
    waste = g.product('cleanup_waste', 'Unretained background products, free adapters and reaction reagents', K.WASTE)
    g.react('cleanup', 'Bead recovery, SFB washes and elution', S.CLEANUP,
            ((ligated,R.TEMPLATE),(debris,R.TEMPLATE),(protein_remnant,R.TEMPLATE),(beads,R.SELECTION)), (library,waste))
    g.discard(waste, 'Removed only at the bead cleanup; no pre-cleanup physical separation is declared')
    retained = g.product('retained_library', 'Source15.8uL library retained for loading', K.LIBRARY)
    remainder = g.product('eluate_remainder', 'Untransferred0.8uL eluate; subtraction-derived amount', K.LIBRARY)
    g.react('select_library', 'Retain15.8uL of16.6uL eluate', S.SPLIT, ((library,R.TEMPLATE),),(retained,remainder),changes_molecule=False)
    g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'Source16.6uL elution/15.8uL retained; remainder fate not specified by selected scope.'))
    load_supply = g.input('load_supply', 'Actual source loading mixture:26uLSequencingBuffer,9.5uLLoadingBeads,0.5uLSequencingTether; packaging/whole molecular composition unassembled', K.MOLECULAR)
    priming = g.input('flowcell_priming', 'SourceR9.4.1 QC>800pores;800uLFlushBuffer then70uLeachSequencingBuffer/water/FlushBuffer priming; no upstream nCATS import', K.CONDITION)
    loaded = g.product('loaded_library', 'Selected native library with source loading supplies', K.LIBRARY)
    g.react('prepare_loading', 'Prepare source native-library loading mixture', S.MAGIC,
            ((retained,R.TEMPLATE),(load_supply,R.STRUCTURAL),(priming,R.CONDITION)),(loaded,),changes_molecule=False,
            reason='Source-selected15.8uL retained library contribution; kit motor/tether molecular architecture and read lengths remain unassembled. No claim that all sample DNA was recovered or sequenced.')
    g.sequence(loaded, SequencingPlan(acquired_roles=('Sequencing output',),
        read_uses={'Sequencing output':'include'},
        limitation='R9.4.1 MinION, 48 h; primer-free native-DNA reads. LSK109 adapter sequence/end geometry and per-molecule read lengths remain unresolved.'))
    g.declare_reviewed_scope((
        'rnp_conditions','assemble_rnp','protein_remnant','datp','adapter_components','select_library','eluate_remainder','load_supply','flowcell_priming','prepare_loading','loaded_library','genomic_dna','repair_mix','repaired_dna','repair','guides','cas9',
        'protected_target','background','cas9_cut','exonucleases','surviving_target',
        'digestion_products','background_digestion','heat','released_target','heat_release',
        'taq_datp','a_tailed','a_tail','lsk109','ligase','ligated','ligation','beads',
        'library','cleanup_waste','cleanup'), SOURCE)
    return g.require_publishable()

def render_page():
    from protocol_views import assemble_page,WorkflowSection,OligoSection,OligoEntry
    g=protocol_graph()
    return assemble_page(g,(head('CaBagE — Cas9-protected native DNA'),'<div class="wrap" data-protocol-complete="true"><h1>CaBagE</h1>',
        '<p class="research-notes"><a href="01_cabage.html">Research notes</a></p>',info('<a href="https://doi.org/10.1371/journal.pone.0241253">Published CaBagE article and S1 guide table</a>'),
        '<h2>Protect the target while digesting its background</h2>',
        '<p>CaBagE places Cas9 at both flanks of a native genomic region. Exonucleases digest exposed background DNA while target-associated Cas9 provides qualitative protection. Heating releases the target; dA-tailing, native adapter ligation and bead cleanup prepare it for nanopore sequencing without PCR.</p>',
        info('INFERRED — target drawing width is schematic, not a45-base target or a mapped cut footprint. The source demonstrates qualitative protection, not perfect survival, purity, specificity or yield. Exact genomic sequence, guide scaffold/modifications, LSK109 whole adapter and occupied junction remain unassembled.'),
        '<h2>Known guide spacers, incomplete guide constructs</h2>',
        OligoSection('C9 reference spacers as printed in S1 Table',tuple(OligoEntry(k,(Segment(k,v),)) for k,v in GUIDES.items()),note_html='These twenty-base source spacer strings are not full orderable crRNA/tracrRNA or proof of exact genomic pairing. No whole guide scaffold or terminal modification is borrowed.'),
        '<h2>Source stocks and readout</h2>',
        '<p>RNPs are preassembled before genomic cleavage. Background digestion products and heat-released protein material remain until cleanup. The source uses Short Fragment Buffer washes: no uninspected S1Figure comparison is used to claim superiority. Native elution yields16.6uL, of which15.8uL is retained; the derived0.8uL remainder stays explicit. Loading supplies and flow-cell priming follow the source Methods.</p>',
        '<p>R9.4.1 MinION sequencing runs48hours. This view does not assign a fixed read length, read-primer site or sequence-exact motor/adapter architecture. Cas9 cleavage and exonuclease protection are local source declarations; no nCATS dephosphorylation step or complete historical kit recipe is inherited.</p>',
        WorkflowSection('Full native-DNA preparation and molecular limits'),'</div>'))

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
    c('dATP is a molecular contribution distinct from Taq',any(p.material=='datp' and p.role is R.STRUCTURAL for p in g.steps['a_tail'].inputs),True)
    c('preassembled RNP is used without duplicate guide addition',[p.material for p in g.steps['cas9_cut'].inputs],['repaired_dna','rnp'])
    c('protein material is not silently deleted by heating',g._consumers['protein_remnant'],['cleanup'])
    c('untransferred eluate survives loading',g.terminals['eluate_remainder'].disposition,Disposition.UNRESOLVED)
    c('sequencing uses the selected loading preparation',tuple(g.sequencing),('loaded_library',))
    c.report()
