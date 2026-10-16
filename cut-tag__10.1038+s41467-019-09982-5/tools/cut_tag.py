"""Source-owned original2019 bulk CUT&Tag; no shared endpoint recipe."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, oligo, revcomp
from pcr import PCRGeometry
from protocol_graph import Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal
from protocol_model import molecule, role

TITLE = 'CUT&Tag — original 2019 bulk antibody-tethered tagmentation'
NOTES = '01_cut-tag.html'
SOURCE = '<a href="https://doi.org/10.1038/s41467-019-09982-5">Kaya-Okur et al., <i>Nature Communications</i> (2019)</a>; exact loaded oligos: cited <a href="https://doi.org/10.1101/gr.177881.114">Picelli et al. (2014)</a>; historical author protocol <a href="https://doi.org/10.17504/protocols.io.wnufdew">Version 1 (March 2019)</a>.'
SUMMARY = 'Antibodies recruit adapter-loaded Protein A–Tn5 to chromatin. Magnesium activates targeted tagmentation; DNA extraction precedes the gap-fill extension and library PCR.'
CAVEAT = ('This page represents the original bulk method, whose PCR uses a universal i5 primer and a uniquely barcoded i7 primer. '
          'The exact universal i5 oligo and choice between the source-allowed Nextera/Buenrostro primer families are not assigned here. '
          'The paper’s separate single-cell method uses 72 i5 and 72 i7 indices. '
          'INFERRED — fragment-end assemblies and possible same-end products below are derived from the two loaded adapter identities; '
          'dimer composition, product frequencies and same-end product fates are not established. '
          'The model follows the published article’s AMPure DNA extraction. Its cited March2019 protocol v1 instead uses phenol/chloroform and ethanol precipitation; these are distinct documented routes, not a claimed publication error.')

# Picelli2014 Methods explicitly assign these canonical bases to the loaded oligos.
ME_A = nx.S5 + nx.ME
ME_B = nx.S7 + nx.ME
ME_REVERSE = nx.ME_RC
INDEX_FEATURE = feature('sample_index_i7', 'sample_index', 'unknown')
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
FINAL_LIBRARY = Construct((role('universal i5-derived 5-prime region; configuration unassigned'),
    *PREAMP_TEMPLATE.segments,
    Segment('candidate i7 reverse complement', 'N' * 8, 'cbc', placeholder=True,
            feature=INDEX_FEATURE, inferred=True),
    Segment('P7 reverse complement', il.P7_RC, 'p7', inferred=True)),
    name='INFERRED — bulk CUT&Tag endpoint with source-allowed i7 family')
FINAL_CAPTION = 'INFERRED — conditional A/B-derived endpoint with the cited Buenrostro i7 family. The exact universal i5-derived region, selected i7 member and loaded run primers are unassigned; no dual-variable-index bulk run is implied.'
SEQ_PRIMERS = ()
SEQUENCING_ENDING = 'The published article specifies paired-end Illumina sequencing and describes economical 25×25-base pooling. It does not assign one acquired run configuration to every experiment.'
SEQUENCING_UNAVAILABLE = 'Exact loaded run-primer oligos and universal-i5 architecture are not declared in this model; PCR-site and run-primer identities are kept separate.'
READ_SPECS = (
    sp.read_spec('Read 1', '25-base economical example; run-specific', 'target-proximal genomic DNA', 'include'),
    sp.read_spec('Index 1 (i7)', '8 in represented candidate family', 'sample barcode; selected member and sample assignment unassigned', 'demultiplex'),
    sp.read_spec('Index 2 (i5)', 'not assigned', 'universal i5 primer used across the represented bulk samples; exact configuration and acquired cycles unassigned', 'undetermined'),
    sp.read_spec('Read 2', '25-base economical example; run-specific', 'opposite end of target-proximal genomic DNA', 'include'),
)

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
    return (oligo('Picelli ME-A', [Segment('ME-A', ME_A)]),
            oligo('Picelli ME-B', [Segment('ME-B', ME_B)]),
            oligo('Picelli Tn5MErev', [Segment('MErev', ME_REVERSE)], mods='/5Phos/'))

def protocol_graph():
    graph = ProtocolGraph(TITLE, authority='Kaya-Okur2019 bulk Methods; Picelli2014 oligos; protocols.io wnufdew v1')
    cells = graph.input('cells', 'permeabilized, ConA-bound cells (bulk entry)', MaterialKind.BIOLOGICAL,
                         state=chromatin_scene('chromatin in prepared cells'), primary=True)
    primary = graph.input('primary_antibody', 'target-specific primary antibody', MaterialKind.SELECTION)
    secondary = graph.input('secondary_antibody', 'appropriate secondary antibody for Protein A recruitment', MaterialKind.SELECTION)
    adapters = graph.input('loaded_adapter_pool', 'equimolar preannealed ME-A/ME-B adapter pool', MaterialKind.ADAPTER,
                            state=adapter_pool_scene(), source='Picelli2014 Methods110; Kaya-Okur2019 transposome preparation')
    tn5 = graph.input('pa_tn5', 'Protein A–Tn5 adapter complex', MaterialKind.CATALYST)
    targeted = graph.product('targeted_cells', 'washed antibody-tethered pA–Tn5 cells', MaterialKind.MOLECULAR,
        state=chromatin_scene('antibody-targeted chromatin', 'antibodies and pA–Tn5 interact noncovalently; dimer adapter pairing unassigned'))
    unbound = graph.product('unbound_reagents', 'unbound antibody and pA–Tn5', MaterialKind.WASTE,
                            state=molecule('wash material', (role('unbound targeting reagents'),)))
    graph.react('target_chromatin', 'Bind antibodies and adapter-loaded pA–Tn5; wash', StepKind.AFFINITY_CAPTURE,
                ((cells, InputRole.TEMPLATE), (primary, InputRole.SELECTION), (secondary, InputRole.SELECTION),
                 (adapters, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)), (targeted, unbound))
    graph.discard(unbound, 'The published bulk method washes away unbound antibodies and pA–Tn5 before activation.')
    magnesium = graph.input('magnesium', 'magnesium activation condition', MaterialKind.CONDITION)
    tagged = graph.product('gapped_ab', 'INFERRED — target-proximal A/B gapped product', MaterialKind.MOLECULAR,
        state=molecule('gapped DNA; exact strand ends unassembled', (role('A/B-bearing chromatin fragment with transposition gaps'),)),
        inferred=True, source='Derived A/B end assignment from the source adapter pool; gap filling explicitly reported in bulk Methods')
    same = graph.product('same_end_possibilities', 'INFERRED — A/A and B/B end possibilities; fates unassigned', MaterialKind.MOLECULAR,
        state=same_end_scene(), inferred=True, source='End-identity possibilities of the declared pool; not measured product frequencies or absolute amplification claims')
    graph.react('tagment', 'INFERRED — fragment-end classes after source-reported targeted tagmentation', StepKind.TAGMENTATION,
                ((targeted, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST), (magnesium, InputRole.CONDITION)), (tagged, same))
    graph.mark_terminal(Terminal(same.id, Disposition.UNRESOLVED,
        'Same-end products are shown explicitly. Their abundance and physical fate are unassigned; absence of both distinct library ends is not proof that single-primer amplification is impossible.'))
    extraction = graph.input('dna_extraction', 'EDTA/SDS/Proteinase K and AMPure DNA extraction (published bulk route)', MaterialKind.CONDITION)
    recovered = graph.product('recovered_dna', 'extracted tagmented DNA', MaterialKind.MOLECULAR, state=tagged.state)
    extraction_waste = graph.product('extraction_waste', 'DNA-extraction wash material', MaterialKind.WASTE,
                                     state=molecule('discarded wash fractions', (role('material not recovered with DNA'),)))
    graph.react('extract_dna', 'Stop tagmentation and extract DNA using the published article route', StepKind.CLEANUP,
                ((tagged, InputRole.TEMPLATE), (extraction, InputRole.CONDITION)), (recovered, extraction_waste))
    graph.discard(extraction_waste, 'The published article withdraws AMPure supernatant and ethanol wash fractions.')
    mix = graph.input('pcr_mix', 'NEBNext HiFi 2X PCR Master Mix', MaterialKind.CATALYST)
    filled = graph.product('gap_filled', 'INFERRED — gap-filled A/B template', MaterialKind.MOLECULAR,
                            state=PREAMP_TEMPLATE, inferred=True, source='Derived A/B fragment assembly; original bulk protocol reports gap filling before denaturation')
    graph.react('gap_fill', 'INFERRED — A/B assembly after source-reported gap filling', StepKind.EXTENSION,
                ((recovered, InputRole.TEMPLATE), (mix, InputRole.CATALYST)), (filled,))
    i5 = graph.input('universal_i5', 'source universal i5 primer; exact oligo unassigned', MaterialKind.OLIGO,
                     state=molecule('universal i5 PCR primer', (role('exact universal i5 oligo unassigned'),)))
    i7 = graph.input('barcoded_i7', 'source-allowed Buenrostro i7 family; member unassigned', MaterialKind.OLIGO, state=P7_POOL,
                     source='Buenrostro2015 MOESM36 indexed Ad2 family; author v1PCR footnote permits this family')
    library = graph.product('library', FINAL_LIBRARY.name, MaterialKind.LIBRARY, state=FINAL_LIBRARY,
                            inferred=True, source='Conditional endpoint; source-allowed i7 family with universal-i5 architecture unassigned')
    graph.react('pcr', 'INFERRED — amplify the represented A/B branch with universal i5 and barcoded i7', StepKind.PCR,
                ((filled, InputRole.TEMPLATE), (i5, InputRole.STRUCTURAL), (i7, InputRole.STRUCTURAL), (mix, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(PREAMP_TEMPLATE,
                    reason='The exact universal i5 oligo and choice between the source-allowed Nextera/Buenrostro primer variants are not declared; selected i7 is unassigned.'))
    graph.sequence(library, SequencingPlan(library=FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE))
    graph.declare_reviewed_scope(
        ('target_chromatin', 'tagment', 'extract_dna', 'gap_fill', 'pcr', 'cells', 'primary_antibody', 'secondary_antibody', 'loaded_adapter_pool', 'pa_tn5', 'targeted_cells', 'unbound_reagents', 'magnesium', 'gapped_ab', 'same_end_possibilities', 'dna_extraction', 'recovered_dna', 'extraction_waste', 'pcr_mix', 'gap_filled', 'universal_i5', 'barcoded_i7', 'library'), graph.authority)
    return graph.require_publishable()
