"""Source-owned CUTAC, with an explicitly conditional in-house adapter branch."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, oligo, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal
from protocol_model import molecule, role

TITLE = 'CUTAC — low-salt tethered tagmentation of adjacent accessible DNA'
NOTES = '01_cutac.html'
SOURCE = '<a href="https://doi.org/10.7554/eLife.63274">Henikoff et al., eLife (2020)</a>; author <a href="https://doi.org/10.17504/protocols.io.bmbfk2jn">CUTAC Version2, September2020</a> and <a href="https://doi.org/10.17504/protocols.io.8yrhxv6">in-house pA–Tn5 loading</a> originals.'
SUMMARY = 'Antibodies tether pA(G)–Tn5 near transcription-associated nucleosomes; low-salt activation favors nearby accessible DNA. SDS release and Triton neutralization precede gap filling and PCR in the same tube.'
CAVEAT = ('INFERRED — this molecular branch combines the explicitly cited in-house pA–Tn5 adapter pool with the permitted barcoded Buenrostro PCR family. '
          'The paper also uses commercial pAG–Tn5; its loaded oligos are not assigned here, and the in-house pool is not evidence of their identity. '
          'The detailed protocol permits universal OR barcoded i5; the endpoint below is conditional on the barcoded-i5 option, with selected index members unassigned. '
          'Its explicit instruction is “Do not use Nextera or NEBNext primers.” Canonical adapter sequence names do not override that source-specific primer choice. '
          'Fragment-end assemblies, same-end possibilities and their frequencies/fates remain inferred or unassigned; no source-wide absence claim is made.')

SOURCE_OLIGO_EXCEPTIONS = {
    'loaded_adapter_pool': 'Original PDF table columns interleave bases with numeric fields; whole source strands including [Phos] are checked against decoded anonymous API29425 table cells by local selftest.py. The generic guard currently reads text twins, not public API JSON.'
}

# Author-linked 8yrhxv6 original table assigns these bases to the in-house pool.
ME_A = nx.S5 + nx.ME
ME_B = nx.S7 + nx.ME
ME_REVERSE = nx.ME_RC
I5_FEATURE = feature('sample_index_i5', 'sample_index', 'unknown')
INDEX_FEATURE = feature('sample_index_i7', 'sample_index', 'unknown')
I5 = Segment('candidate i5; selected member unassigned', 'N' * 8, 'cbc', placeholder=True, feature=I5_FEATURE)
P5_POOL = Construct((Segment('P5', il.P5, 'p5'), I5, Segment('Buenrostro i5 annealing arm', nx.S5 + nx.ME[:10], 'r1')), name='Buenrostro barcoded i5 family; optional source branch')
I7 = Segment('candidate i7 index; selected member unassigned', 'N' * 8, 'cbc',
             placeholder=True, feature=INDEX_FEATURE)
# Source-allowed Buenrostro2015 i7 family: P7 + eight bases + S7 + first7 ME bases.
P7_POOL = Construct((Segment('P7', il.P7, 'p7'), I7,
                     Segment('Buenrostro i7 annealing arm', nx.S7 + nx.ME[:7], 'r2')),
                    name='source-allowed Buenrostro2015 i7 family; member unassigned')

def paired_end_construct(left, right, *, name):
    """Conditional fragment assembly; adapter identities are explicit local choices."""
    return Construct((Segment('left entry arm', left[:-len(nx.ME)], 's5', inferred=True),
                      Segment('left ME', nx.ME, 'me', inferred=True),
                      role('target-proximal DNA insert'),
                      Segment('right ME reverse complement', nx.ME_RC, 'me', inferred=True),
                      Segment('right entry reverse complement', revcomp(right[:-len(nx.ME)]), 's7', inferred=True)), name=name)

PREAMP_TEMPLATE = paired_end_construct(ME_A, ME_B, name='INFERRED — A/B gap-filled fragment')
FINAL_LIBRARY = Construct((Segment('P5', il.P5, 'p5', inferred=True),
    Segment('candidate i5', 'N' * 8, 'cbc', placeholder=True, feature=I5_FEATURE, inferred=True),
    *PREAMP_TEMPLATE.segments,
    Segment('candidate i7 reverse complement', 'N' * 8, 'cbc', placeholder=True, feature=INDEX_FEATURE, inferred=True),
    Segment('P7 reverse complement', il.P7_RC, 'p7', inferred=True)),
    name='INFERRED — in-house A/B product with barcoded Buenrostro i5/i7')
FINAL_CAPTION = 'INFERRED — conditional in-house A/B branch and permitted barcoded-i5/Buenrostro-i7 family. The source also allows universal i5 and commercial pAG–Tn5; no sample-specific primer assignment or commercial oligo identity is implied.'
SEQ_PRIMERS = ()
SEQUENCING_ENDING = 'The original article reports paired-end sequencing on an Illumina HiSeq2500 Turbo flow cell; acquired read/index cycles are unassigned in this model.'
SEQUENCING_UNAVAILABLE = 'Exact loaded run-primer oligos and index orientation are not declared; a source PCR arm does not establish the loaded sequencing primer.'
READ_SPECS = (
    sp.read_spec('Read 1', 'run-specific', 'accessible genomic DNA', 'include'),
    sp.read_spec('Index 1 (i7)', '8 in represented custom family', 'sample index; selected member unassigned', 'demultiplex'),
    sp.read_spec('Index 2 (i5)', '8 in represented barcoded option', 'sample index in this conditional branch; universal-i5 alternative exists', 'demultiplex'),
    sp.read_spec('Read 2', 'run-specific', 'opposite end of accessible genomic DNA', 'include'),
)

def pcr_geometry():
    return PCRGeometry.exact(PREAMP_TEMPLATE,
        PCRPrimer('Buenrostro i5 family', P5_POOL.top(), anneal_nt=len(nx.S5) + 10),
        PCRPrimer('Buenrostro i7 family', P7_POOL.top(), anneal_nt=len(nx.S7) + 7),
        statement='INFERRED — exact inward binding of the source custom PCR arms on the conditional in-house A/B gap-filled template; selected index members unassigned.')

def adapter_pool_scene():
    sc = Scene()
    for name, seq, tag in (('A', ME_A, 's5'), ('B', ME_B, 's7')):
        core = name + ' mosaic end'
        sc.strand(name, [Segment(name + ' single-strand entry', seq[:-len(nx.ME)], tag),
                         Segment(core, nx.ME, 'me')], label='source ME-' + name)
        sc.anneal(name + '_rev', [Segment(name + ' ME reverse', ME_REVERSE, 'me')],
                  to=name, pair=(name + ' ME reverse', core), label='source phosphorylated MErev', mod5='/5Phos/')
    return sc

def chromatin_scene(label, note=''):
    sc = molecule(label, (role('chromatin DNA'),))
    if note: sc.note('molecule', note)
    return sc

def same_end_scene():
    sc = Scene()
    for key, seq in (('AA', ME_A), ('BB', ME_B)):
        state = paired_end_construct(seq, seq, name='INFERRED — ' + key)
        sc.strand(key, list(state), label='INFERRED — ' + key + ' fragment')
    return sc

def oligos():
    return (oligo('in-house ME-A', [Segment('ME-A', ME_A)]),
            oligo('in-house ME-B', [Segment('ME-B', ME_B)]),
            oligo('in-house Tn5MErev', [Segment('MErev', ME_REVERSE)], mods='/5Phos/'))

def protocol_graph():
    graph = ProtocolGraph(TITLE, authority='Henikoff2020 Methods; bmbfk2jn Version2; in-house loading8yrhxv6; custom Buenrostro2015 PCR family')
    cells = graph.input('cells', 'prepared, ConA-bound nuclei', MaterialKind.BIOLOGICAL,
                         state=chromatin_scene('chromatin in prepared cells'), primary=True)
    primary = graph.input('primary_antibody', 'target-specific primary antibody', MaterialKind.SELECTION)
    secondary = graph.input('secondary_antibody', 'appropriate secondary antibody for Protein A recruitment', MaterialKind.SELECTION)
    adapters = graph.input('loaded_adapter_pool', 'in-house preannealed ME-A/ME-B pool; commercial identity not assigned', MaterialKind.ADAPTER,
                            state=adapter_pool_scene(), source='Author-linked 8yrhxv6 original loading table; in-house branch only')
    tn5 = graph.input('pa_tn5', 'Protein A–Tn5 adapter complex', MaterialKind.CATALYST)
    targeted = graph.product('targeted_cells', 'washed antibody-tethered pA–Tn5 cells', MaterialKind.MOLECULAR,
        state=chromatin_scene('antibody-targeted chromatin', 'antibodies and pA–Tn5 interact noncovalently; dimer adapter pairing unassigned'))
    unbound = graph.product('unbound_reagents', 'unbound antibody and pA–Tn5', MaterialKind.WASTE,
                            state=molecule('wash material', (role('unbound targeting reagents'),)))
    graph.react('target_chromatin', 'Bind antibodies and adapter-loaded pA–Tn5; wash', StepKind.AFFINITY_CAPTURE,
                ((cells, InputRole.TEMPLATE), (primary, InputRole.SELECTION), (secondary, InputRole.SELECTION),
                 (adapters, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)), (targeted, unbound))
    graph.discard(unbound, 'The CUTAC detailed method washes away unbound antibodies and pA–Tn5 before activation.')
    magnesium = graph.input('magnesium', 'source low-salt magnesium activation option', MaterialKind.CONDITION)
    tagged = graph.product('gapped_ab', 'INFERRED — target-proximal A/B gapped product', MaterialKind.MOLECULAR,
        state=molecule('gapped DNA; exact strand ends unassembled', (role('A/B-bearing chromatin fragment with transposition gaps'),)),
        inferred=True, source='Derived A/B end assignment from the source adapter pool; gap filling explicitly reported in CUTAC Methods')
    same = graph.product('same_end_possibilities', 'INFERRED — A/A and B/B end possibilities; fates unassigned', MaterialKind.MOLECULAR,
        state=same_end_scene(), inferred=True, source='End-identity possibilities of the declared pool; not measured product frequencies or absolute amplification claims')
    graph.react('tagment', 'INFERRED — fragment-end classes after source-reported low-salt CUTAC tagmentation', StepKind.TAGMENTATION,
                ((targeted, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST), (magnesium, InputRole.CONDITION)), (tagged, same))
    graph.mark_terminal(Terminal(same.id, Disposition.UNRESOLVED,
        'Same-end products are shown explicitly. Their abundance and physical fate are unassigned; absence of both distinct library ends is not proof that single-primer amplification is impossible.'))
    taps = graph.input('taps_edta_wash', '50 µL 10 mM TAPS pH 8.5, 0.2 mM EDTA bead wash', MaterialKind.CONDITION)
    washed = graph.product('washed_beads', 'bead-bound tagmented chromatin after TAPS wash', MaterialKind.MOLECULAR, state=tagged.state)
    taps_waste = graph.product('tagmentation_liquid', 'removed tagmentation solution and TAPS wash', MaterialKind.WASTE,
                               state=molecule('withdrawn liquid', (role('tagmentation solution and wash; composition not characterized'),)))
    graph.react('taps_wash', 'Chill, magnetize, remove liquid and wash beads in TAPS/EDTA', StepKind.CLEANUP,
                ((tagged, InputRole.TEMPLATE), (taps, InputRole.CONDITION)), (washed, taps_waste))
    graph.discard(taps_waste, 'Article Methods and bmbfk2jn steps35–36 withdraw the tagmentation liquid and the TAPS wash.')
    release = graph.input('sds_release', '5 µL 0.1% SDS + 10 µL TAPS at 58 °C, then 15 µL 0.67% Triton X-100 neutralization', MaterialKind.CONDITION)
    recovered = graph.product('released_dna', 'released tagmented DNA in neutralized reaction', MaterialKind.MOLECULAR, state=tagged.state)
    graph.react('release_dna', 'SDS release at 58 °C then Triton neutralization before direct PCR', StepKind.DEPROTEINIZATION,
                ((washed, InputRole.TEMPLATE), (release, InputRole.CONDITION)), (recovered,))
    mix = graph.input('pcr_mix', 'NEBNext 2X PCR Master Mix (non-hot-start)', MaterialKind.CATALYST)
    filled = graph.product('gap_filled', 'INFERRED — gap-filled A/B template', MaterialKind.MOLECULAR,
                            state=PREAMP_TEMPLATE, inferred=True, source='Derived A/B fragment assembly; CUTAC detailed protocol reports gap filling before denaturation')
    graph.react('gap_fill', 'INFERRED — A/B assembly after source-reported gap filling (bmbfk2jn: 58 °C 5 min, 72 °C 5 min)', StepKind.EXTENSION,
                ((recovered, InputRole.TEMPLATE), (mix, InputRole.CATALYST)), (filled,))
    i5 = graph.input('barcoded_i5', 'Buenrostro barcoded i5 option; member unassigned', MaterialKind.OLIGO,
                     state=P5_POOL, source='bmbfk2jn PCR universal OR barcoded i5; represented barcoded branch uses original Buenrostro family')
    i7 = graph.input('barcoded_i7', 'source-allowed Buenrostro i7 family; member unassigned', MaterialKind.OLIGO, state=P7_POOL,
                     source='Buenrostro2015 MOESM36 indexed Ad2 family; CUTAC Version2 explicitly requires custom family, excludes Nextera/NEBNext primers')
    library = graph.product('library', FINAL_LIBRARY.name, MaterialKind.LIBRARY, state=FINAL_LIBRARY,
                            inferred=True, source='Conditional in-house A/B and custom barcoded-i5/i7 endpoint')
    graph.react('pcr', 'INFERRED — 12-cycle PCR (98 °C 30 s; 98 °C 10 s / 60 °C 10 s) of the represented A/B branch with custom barcoded i5/i7 option', StepKind.PCR,
                ((filled, InputRole.TEMPLATE), (i5, InputRole.STRUCTURAL), (i7, InputRole.STRUCTURAL), (mix, InputRole.CATALYST)), (library,),
                geometry=pcr_geometry())
    clean = graph.input('spri_cleanup', '65 µL (1.3 vol) SPRI bead slurry, 80% ethanol washes, elution in 20 µL 1 mM Tris-HCl pH 8, 0.1 mM EDTA', MaterialKind.CONDITION)
    ready = graph.product('sequencing_library', FINAL_LIBRARY.name, MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True)
    wash = graph.product('pcr_cleanup_wash', 'post-PCR cleanup wash material', MaterialKind.WASTE)
    graph.react('cleanup_library', 'Clean up the amplified CUTAC library', StepKind.CLEANUP,
                ((library, InputRole.TEMPLATE), (clean, InputRole.CONDITION)), (ready, wash))
    graph.discard(wash, 'Source post-PCR cleanup removes wash material; exact molecular recovery is not quantified here.')
    graph.sequence(ready, SequencingPlan(library=FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE))
    graph.declare_reviewed_scope(
        ('target_chromatin', 'tagment', 'release_dna', 'gap_fill', 'pcr', 'cleanup_library', 'cells', 'primary_antibody', 'secondary_antibody', 'loaded_adapter_pool', 'pa_tn5', 'targeted_cells', 'unbound_reagents', 'magnesium', 'gapped_ab', 'same_end_possibilities', 'taps_wash', 'washed_beads', 'tagmentation_liquid', 'taps_edta_wash', 'sds_release', 'released_dna', 'pcr_mix', 'gap_filled', 'barcoded_i5', 'barcoded_i7', 'library', 'spri_cleanup', 'sequencing_library', 'pcr_cleanup_wash'), graph.authority)
    return graph.require_publishable()
