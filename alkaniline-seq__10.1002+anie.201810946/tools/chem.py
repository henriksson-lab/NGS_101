"""Source-owned AlkAniline cleavage selection; historical kit geometry stays explicit."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal)
from protocol_model import molecule, role

TITLE = 'AlkAniline-Seq — positive selection of modification-derived RNA cleavage'
NOTES = '01_alkaniline_seq.html'
SOURCE = 'Marchand et al., Angewandte Chemie (2018), <a href="https://doi.org/10.1002/anie.201810946">doi:10.1002/anie.201810946</a>, Figure 1 and main-text cleavage/ligation explanation.'
SUMMARY = 'Alkaline treatment creates abasic sites at susceptible RNA modifications; background dephosphorylation precedes aniline cleavage, whose new 5′ phosphates select RNA beginning at N+1 for library conversion.'
CAVEAT = ('The native graph follows the 2018 article and Figure 1. Its detailed Supporting Information and historical kit implementation are not reconstructed here. '
          'The exact oligos listed below belong to a later NEB manual v9.0 (July2025), shown as comparator inventory only. '
          'Neither their current formulation nor a current index member is substituted for the 2018 experiment. '
          'The chemistry also responds to dihydrouridine and other abasic-site-forming lesions; the signal is not unique to one modification.')
FINAL_LIBRARY = Construct([role('historical Illumina left library end; not reconstructed'),
                           role('cDNA/read beginning at N+1 after aniline-sensitive position'),
                           role('historical Illumina right library end; not reconstructed')], name='AlkAniline source-limited library endpoint')
FINAL_CAPTION = 'The informative read 5′ coordinate is one nucleotide downstream of the susceptible position.'
SEQ_PRIMERS = ()
READ_LENGTHS = None
SEQUENCING_INTRO = 'The 2018 article places the beginning of Read 1 at N+1; its exact historical primer and acquisition recipe remain outside this reconstruction.'
SEQUENCING_ENDING = SEQUENCING_INTRO
SEQUENCING_UNAVAILABLE = ('Read 1 start semantics are source-positive. The earlier curation note reportsSR50, but the primary passages re-read here do not establish acquired cycles or exact primer sites. The later kit comparator does not supply an experimental run.')

def _rna(label, parts, *, mod5='', mod3=''):
    scene = Scene()
    scene.strand('RNA', list(parts), label=label, mod5=mod5, mod3=mod3)
    return scene

def oligos():
    return ("<b>NEB manual v9.0 (July2025) comparator only; these oligos are not assigned to the 2018 native graph.</b>", 'NEBNext E7300/E7330 3′ SR adapter — 5′-rAppAGATCGGAAGAGCACACGTCT-NH2-3′', 'NEBNext E7300/E7330 5′ SR adapter — 5′-rGrUrUrCrArGrArGrUrUrCrUrArCrArGrUrCrCrGrArCrGrArUrC-3′', 'NEBNext E7300/E7330 SR RT primer — 5′-AGACGTGTGCTCTTCCGATCT-3′', 'NEBNext E7300/E7330 SR PCR primer — 5′-AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCG*A-3′', "NEBNext E7300/E7330 index 1 PCR primer (kit example; experiment's index choice unstated) — 5′-CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC*T-3′")


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    start = g.input('input', 'total or enriched RNA', MaterialKind.BIOLOGICAL,
        state=_rna('input RNA', [role('RNA containing labile modified positions')]), primary=True)
    alkali = g.input('alkaline_condition', 'bicarbonate pH 9.2 at 96 °C (Figure 1a)', MaterialKind.CONDITION)
    fragmented = g.product('alkali_rna', 'fragmented RNA with abasic intermediates', MaterialKind.MOLECULAR,
        state=_rna('alkali-treated RNA', [role('RNA fragments containing abasic positions')]))
    g.react('hydrolyse', 'Partially hydrolyse RNA and generate modification-derived abasic sites', StepKind.CHEMICAL_CONVERSION,
        ((start, InputRole.TEMPLATE), (alkali, InputRole.CONDITION)), (fragmented,))
    phosphatase = g.input('phosphatase', 'alkaline phosphatase (Figure 1a; exact historical preparation unassigned)', MaterialKind.CATALYST)
    reset = g.product('dephosphorylated', 'RNA fragments reset to 5′-OH/3′-OH ends', MaterialKind.MOLECULAR,
        state=_rna('background-dephosphorylated RNA', [role('RNA with abasic positions')], mod5='OH', mod3='OH'))
    g.react('dephosphorylate', 'Remove pre-existing 5′ and 3′ phosphates', StepKind.DEPHOSPHORYLATION,
        ((fragmented, InputRole.TEMPLATE), (phosphatase, InputRole.CATALYST)), (reset,))
    aniline = g.input('aniline_condition', 'aniline abasic-site cleavage condition', MaterialKind.CONDITION)
    downstream = g.product('cleaved', 'downstream cleavage fragments with new 5′ phosphate at N+1', MaterialKind.MOLECULAR,
        state=_rna('ligation-competent downstream RNA', [role('RNA starting at originalN+1')], mod5='P', mod3='OH'))
    other = g.product('other_cleavage_products', 'upstream cleavage products and fragments lacking a new 5′ phosphate', MaterialKind.MOLECULAR,
        state=_rna('other RNA outcomes', [role('upstream or uncleaved RNA; separate molecules')], mod5='OH'))
    g.react('aniline', 'Cleave abasic sites and expose the N+1 fragment 5′ phosphate', StepKind.CHEMICAL_CONVERSION,
        ((reset, InputRole.TEMPLATE), (aniline, InputRole.CONDITION)), (downstream, other))
    g.mark_terminal(Terminal(other.id, Disposition.UNRESOLVED,
        'Fragments without the new 5′ phosphate do not meet the selective 5′-ligation substrate rule. This does not demonstrate physical discard or resolve their later recovery.'))
    adapter5 = g.input('adapter5', 'historical 5′ sequencing adapter (ordered identity unassigned)', MaterialKind.ADAPTER)
    adapter3 = g.input('adapter3', 'historical 3′ sequencing adapter (ordered identity unassigned)', MaterialKind.ADAPTER)
    ligases = g.input('adapter_ligases', 'historical adapter-ligation enzyme preparation(s); not assigned from the 2025 comparator', MaterialKind.CATALYST)
    ligated = g.product('selected_rna', 'adapter-ligated RNA beginning at the N+1 cleavage coordinate', MaterialKind.MOLECULAR,
        state=_rna('selected adapter-bearing RNA', [role('historical 5′ adapter'), role('RNA beginning at N+1'), role('historical 3′ adapter')]))
    ligated.state.junction('RNA', 'historical 5′ adapter', 'RNA beginning at N+1', '5′ adapter ligation')
    ligated.state.junction('RNA', 'RNA beginning at N+1', 'historical 3′ adapter', '3′ adapter ligation')
    g.react('ligate', 'Install sequencing adapters through selective recognition of the new 5′ phosphate', StepKind.LIGATION,
        ((downstream, InputRole.TEMPLATE), (adapter5, InputRole.STRUCTURAL), (adapter3, InputRole.STRUCTURAL),
         (ligases, InputRole.CATALYST)), (ligated,), geometry=AdapterLigationGeometry.unavailable(
            'The 2018 article/Figure 1 establishes adapter ligation after aniline and the new 5′-phosphate substrate requirement. Exact historical adapter end sequences and individual ligation ordering have not been reconstructed; current2025 NEB oligos remain comparators.'))
    rt_primer = g.input('rt_primer', 'historical library RT primer; sequence/chemistry unassigned', MaterialKind.OLIGO)
    rt = g.input('rt_enzyme', 'historical reverse-transcription preparation; identity unassigned', MaterialKind.CATALYST)
    cdna = g.product('cdna', 'cleavage-selected first-strand cDNA precursor', MaterialKind.MOLECULAR,
        state=molecule('cDNA strand entering amplification', (role('historical RT-derived end'),
            role('cleavage-selected cDNA'), role('opposite adapter-derived cDNA end'))))
    g.react('rt', 'Convert the selected RNA to a cDNA precursor', StepKind.REVERSE_TRANSCRIPTION,
        ((ligated, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (cdna,),
        reason='The 2018 main text identifies RNA inputs and dsDNA amplicons; this is the RNA-to-cDNA phase of that conversion. Its reagent identities and finer intermediates are not imported from a later kit.')
    forward = g.input('library_forward', 'historical library forward primer; identity unassigned', MaterialKind.OLIGO)
    reverse = g.input('library_reverse', 'historical library opposite-end primer; identity unassigned', MaterialKind.OLIGO)
    pol = g.input('library_polymerase', 'historical library PCR preparation; identity unassigned', MaterialKind.CATALYST)
    library = g.product('library', 'AlkAniline dsDNA library', MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    g.react('library_pcr', 'Amplify cleavage-selected cDNA into the source dsDNA library', StepKind.PCR,
        ((cdna, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL), (reverse, InputRole.STRUCTURAL), (pol, InputRole.CATALYST)),
        (library,), geometry=PCRGeometry.unavailable(cdna.state, reason=
            'The actual precursor is the cDNA strand shown. Historical primer sequences, copied adapter ends and initial complementary-strand formation are not reconstructed; neither a current NEB primer nor the output duplex can supply this missing model state.'))
    g.sequence(library, SequencingPlan(library=FINAL_LIBRARY, acquired_roles=('Read 1',),
        limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(('input','alkaline_condition','hydrolyse','phosphatase','dephosphorylate',
        'aniline_condition','aniline','adapter5','adapter3','adapter_ligases','ligate','rt_primer','rt_enzyme','rt',
        'library_forward','library_reverse','library_polymerase','library_pcr','library'), SOURCE)
    return g.require_publishable()
