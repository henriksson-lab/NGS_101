"""Provisional PhIP-Seq assay branch from the inspected author supplement.

The final-published detailed procedure has not been reconciled. This module is
private at notes status and does not reconstruct the upstream T7 library cloning.
"""
from dataclasses import dataclass
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
import illumina as il
import seqprimers as sp
from chemdraw import (Construct, Scene, Segment, feature, oligo, panel, revcomp,
                      workflow_from_protocol_graph, workflow_panel)
from page import head, info
from pcr import PCRGeometry
from protocol_graph import (ProtocolGraph, InputRole as R, MaterialKind as K,
                            StepKind as S, SequencingPlan, SourceItem)
from protocol_model import role

PAGE_FILE = 'phip-seq.html'
SOURCE = 'Larman2011 DOI10.1038/nbt.1856; associated accepted manuscript and NIHMS284271 author supplement pp14–15; final procedural applicability unverified'
READ_PRIMER = 'AGGTGTGATGCTCGGGGATCCAGGAATTCC'
BARCODES = ('ACT', 'GCC', 'CCT', 'TCT', 'GAT', 'GGT', 'GTT', 'GCG')
FIXED_AFTER_BARCODE = 'GCGC'
P5_PREFIX = il.P5[:20]
# This source oligo is not the complete canonical modern P7 primer.
REVERSE_PRIMER = il.P7[:20] + 'ACACTGAACCCCTCAAGACCCGTTTA'
INDEX = feature('phip_inline_sample', 'sample_index', 'whitelist',
                whitelist='NIHMS284271 author supplementp14 HsORF-FL-mmBC1–8-F',
                note='Eight printed three-base members; actual sample assignment unassigned')
GENOME = Construct((role('T7FNS2 phage DNA carrying a peptide-coding insert'),), name='clone-dependent phage genome')


def genome_scene():
    sc = Scene()
    sc.strand('top', list(GENOME), label='duplex phage DNA; exact strands unassembled')
    return sc


def forward_primer(barcode):
    if barcode not in BARCODES:
        raise ValueError('Select one of the eight printed author-supplement barcodes')
    return P5_PREFIX + READ_PRIMER + barcode + FIXED_AFTER_BARCODE


@dataclass(frozen=True)
class DisplayParticle:
    """A source-declared genotype/display association, not a binding predictor."""
    phase: str

    def __post_init__(self):
        if self.phase not in {'input pool', 'antibody incubation', 'bead-selected pool'}:
            raise ValueError('Unknown source-assay phase')

    def scene(self):
        sc = genome_scene()
        sc.note('top', 'Cognate peptide displayed as a T7 10B capsid fusion; protein is not covalently joined to this DNA strand')
        sc.note('top', self.phase + '; clone identities, binding affinities and background carryover are unassigned')
        return sc


def library():
    """Inferred primer-derived order, with no assigned study barcode member."""
    return Construct((
        Segment('printed P5 prefix', P5_PREFIX, 'p5', inferred=True),
        Segment('custom sequencing-primer site', READ_PRIMER, 'r1', inferred=True),
        Segment('three variable sample bases', 'N' * 3, 'cbc', placeholder=True, feature=INDEX, inferred=True),
        Segment('four fixed bases in seven-cycle parsing window', FIXED_AFTER_BARCODE, inferred=True),
        role('peptide-library insert and unassembled phage flanks', inferred=True),
        Segment('reverse-primer reverse complement', revcomp(REVERSE_PRIMER), 'p7', inferred=True),
    ), name='INFERRED — author-supplement PhIP PCR endpoint')


PRIMERS = (sp.custom('Read 1', 'mmBC-FL_seq_prim', READ_PRIMER, SOURCE),)


def protocol_graph():
    g = ProtocolGraph('PhIP-Seq — provisional author-supplement assay branch', authority=SOURCE)
    particles = g.input('particles', 'Prebuilt T7-Pep display library', K.BIOLOGICAL,
                        state=DisplayParticle('input pool').scene(), primary=True)
    antibodies = g.input('antibodies', 'Patient antibodies; SAPK4 antibody positive control', K.BIOLOGICAL)
    incubation = g.input('incubation', 'BSA-blocked tubes; 24 h at4C in M9LB', K.CONDITION)
    incubated = g.product('incubated', 'Phage/antibody incubation mixture; binding is not assumed universal', K.MOLECULAR,
                         state=DisplayParticle('antibody incubation').scene())
    g.react('incubate', 'Incubate display library with antibody sample', S.BINDING,
            ((particles, R.TEMPLATE), (antibodies, R.ANALYTE), (incubation, R.CONDITION)), (incubated,))
    beads = g.input('beads', '1:1 Protein A/Protein G magnetic bead slurry;4 h at4C; six IP-buffer washes', K.SELECTION)
    retained = g.product('retained', 'Bead-selected phage population, including possible nonspecific background', K.MOLECULAR,
                         state=DisplayParticle('bead-selected pool').scene())
    washes = g.product('washes', 'Unretained particles and soluble material in actual IP washes', K.WASTE)
    g.react('capture', 'Capture antibody-associated material and wash beads six times', S.AFFINITY_CAPTURE,
            ((incubated, R.TEMPLATE), (beads, R.SELECTION)), (retained, washes))
    g.discard(washes, 'Physically removed supernatants/washes; no perfect specificity or yield asserted')
    heat = g.input('heat', 'Resuspend beads in30uL water;90C10min; transfer aqueous template for PCR', K.CONDITION)
    released = g.product('released', 'Released phage DNA in transferred aqueous material; copurified material unassigned', K.MOLECULAR,
                         state=genome_scene())
    g.react('release', 'Heat-disrupt phage and recover aqueous DNA template', S.DEPROTEINIZATION,
            ((retained, R.TEMPLATE), (heat, R.CONDITION)), (released,))
    forward = g.input('forward', 'One selected member of HsORF-FL-mmBC1–8-F; sample assignment unassigned', K.OLIGO,
                      state=Construct((Segment('P5 prefix', P5_PREFIX), Segment('read site', READ_PRIMER),
                                       Segment('sample bases', 'N'*3, placeholder=True, feature=INDEX),
                                       Segment('fixed tail', FIXED_AFTER_BARCODE))))
    reverse = g.input('reverse', 'P7-T7Down exactly as printed', K.OLIGO,
                      state=Construct((Segment('P7-T7Down', REVERSE_PRIMER),)))
    polymerase = g.input('polymerase', 'TaKaRa HS Ex RR006A', K.CATALYST)
    cycling = g.input('cycling', 'Author PDFp15 table:98C10s/56C15s/72C25s;40cycles, optional45; final72C7min', K.CONDITION)
    amplified = g.product('amplified', 'INFERRED — inline-barcoded primer-derived PCR product', K.MOLECULAR,
                          state=library(), inferred=True)
    g.react('pcr', 'INFERRED — PCR with inline sample barcode', S.PCR,
            ((released, R.TEMPLATE), (forward, R.STRUCTURAL), (reverse, R.STRUCTURAL),
             (polymerase, R.CATALYST), (cycling, R.CONDITION)), (amplified,),
            geometry=PCRGeometry.unavailable(GENOME, reason='Exact phage template flanks and mismatch-introducing primer geometry are not reconstructed. Primer-derived endpoint is conditional; final-published procedure applicability also remains unverified.'))
    gel = g.input('gel', 'Individual gel purification of PCR products', K.SELECTION)
    purified = g.product('purified', 'Purified individual sample library', K.LIBRARY, state=library(), inferred=True)
    rejected = g.product('gel_rejected', 'Unselected gel bands and unrecovered material', K.WASTE)
    g.react('gel_cleanup', 'Gel-purify each sample before multiplex pooling', S.CLEANUP,
            ((amplified, R.TEMPLATE), (gel, R.SELECTION)), (purified, rejected))
    g.discard(rejected, 'Removed during gel purification; exact band limits unassigned')
    other = g.input('other_libraries', 'Other individually gel-purified barcoded PhIP samples', K.LIBRARY)
    pool_condition = g.input('pool_mass', '500ng of each barcoded sample; author supplementp15', K.CONDITION)
    pooled = g.product('pooled', 'Multiplexed PhIP libraries', K.LIBRARY, state=library(), inferred=True)
    g.react('pool', 'Pool individually purified sample libraries by mass', S.POOL,
            ((purified, R.TEMPLATE), (other, R.TEMPLATE), (pool_condition, R.CONDITION)), (pooled,))
    g.sequence(pooled, SequencingPlan(library=library(), primers=PRIMERS, acquired_roles=('Read 1',),
        limitation='Author supplement names custom forward sequencing primer and seven leading sample-parsing cycles, not the total acquired read length or an index-read mask. Three variable bases plus four fixed GCGC bases are represented; final publication applicability and exact PCR product remain unverified.'))
    g.declare_reviewed_scope(('particles','antibodies','incubation','incubated','incubate','beads','retained','washes','capture','heat','released','release','forward','reverse','polymerase','cycling','amplified','pcr','gel','purified','gel_rejected','gel_cleanup','other_libraries','pool_mass','pooled','pool'), SOURCE)
    g.add_source_item(SourceItem('upstream_display_library', 'Synthetic oligo pools, T7FNS2 cloning, packaging and propagation',
                                'Associated manuscript Figure1 and author supplementp13; exact vector/oligo assemblies not modelled by this assay branch'))
    return g


def oligos():
    rows = [oligo('HsORF-FL-mmBC'+str(i+1)+'-F', [Segment('P5 prefix',P5_PREFIX),
            Segment('custom read site',READ_PRIMER),
            Segment('sample barcode',bc,'cbc',feature=feature('phip_printed_'+str(i+1),'sample_index','fixed')),
            Segment('fixed GCGC',FIXED_AFTER_BARCODE)]) for i,bc in enumerate(BARCODES)]
    rows.extend((oligo('P7-T7Down',[Segment('source polymer',REVERSE_PRIMER)]),
                 oligo('mmBC-FL_seq_prim',[Segment('custom read primer',READ_PRIMER)])))
    return rows


def render_page():
    graph = protocol_graph()
    endpoint = Scene()
    endpoint.strand('library', list(library()), label='conditional primer-derived strand; insert length unassigned')
    # A private partial view, with no complete-protocol marker. The actual
    # publication gate remains unsatisfied, as checked in selftest().
    return '\n'.join((head('PhIP-Seq — provisional source model'), '<main class="wrap"><h1>PhIP-Seq</h1>',
        '<p class="research-notes"><a href="01_phip-seq.html">Research notes</a></p>',
        info('INFERRED — Provisional model of the inspected author-supplement assay branch. Final published detailed procedure, full T7 library construction and exact PCR geometry remain unverified. Binding and display association are source-declared, not simulated.'),
        '<h2>Author-supplement oligos</h2>', *oligos(),
        workflow_panel(workflow_from_protocol_graph(graph, 'pooled'), cls='long'),
        panel(endpoint.rows(), caption='INFERRED — primer-derived order; three variable sample bases and four fixed bases in the seven-cycle parsing window.'),
        sp.section(library(), PRIMERS, required_roles=('Read 1',)), '</main>'))


def selftest():
    import re
    from checks import Check, Source, have
    from storage import sources_root
    c = Check()
    source = Source(sources_root()/HERE.parent.name/'nbt.1856_supp_NIHMS284271-supplement-1.pdf.txt', 'author supplement pp14–15; provisional applicability')
    if have(c,source,label='PhIP author-supplement transcription'):
        text = re.sub(r'\s+', '', source.path.read_text())
        for i,bc in enumerate(BARCODES):c('printed forward barcode '+str(i+1), forward_primer(bc) in text, True)
        c('printed reverse primer', REVERSE_PRIMER in text, True)
        c('printed custom sequencing primer', READ_PRIMER in text, True)
    c('custom sequencing primer lands on conditional endpoint', sp.verify(library(), PRIMERS, required_roles=('Read 1',)), [])
    c('upstream construction remains a real source/model obligation', bool(protocol_graph().publication_problems()), True)
    c.report()
