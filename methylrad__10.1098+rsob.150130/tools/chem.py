"""Source-local MethylRAD standard-library preparation, with unassembled products."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
import illumina
from chemdraw import Scene, Segment, revcomp, feature
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import (ProtocolGraph, MaterialKind as K, InputRole as R,
                            StepKind as S, Terminal, Disposition as D,
                            SequencingPlan, SourceItem, SourceDiscrepancy, SourceClaim)

AUTHORITY = ('Wang et al., published Open Biology2015 DOI10.1098/rsob.150130; '
             'PMC4680569 Methods5.1–5.3/XML368–382 and overview2.1/Fig1; '
             'published rsob150130supp1.pdf TableS5 actualp6. '
             'Selected standard NNNN library; full occupied genomic products unassembled.')
CORE = illumina.TRUSEQ_READ1[-13:]
ANTI_CORE = revcomp(CORE)
P1 = illumina.TRUSEQ_READ1[:23]
P2 = illumina.TRUSEQ_READ2[:24]
P3 = illumina.P5 + P1[4:]
INDEX_TAIL = illumina.TRUSEQ_READ2[:21]
INDEX_WIDTH = 6
INDEX_ORDER = illumina.P7 + 'X' * INDEX_WIDTH + INDEX_TAIL
READ_LENGTH = 36
POLARITY_TEXT = 'Published overview2.1: arbitrary four-base 5′ overhangs; Figure1 caption: four-base 3′ overhangs.'
# Independent whole-ordering fixtures transcribed from rendered TableS5, not a kit recipe.
SOURCE_ORDERS = {
    'adap1': 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT',
    'adap2': 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT',
    'anti_std': 'NNNNAGATCGGAAGAGC(AminoC6)',
    'anti_rtr1': 'TNNNAGATCGGAAGAGC(AminoC6)',
    'anti_rtr2': 'CNNNAGATCGGAAGAGC(AminoC6)',
    'p1': 'ACACTCTTTCCCTACACGACGCT',
    'p2': 'GTGACTGGAGTTCAGACGTGTGCT',
    'p3': 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCT',
    'index': 'CAAGCAGAAGACGGCATACGAGATXXXXXXGTGACTGGAGTTCAGACGTGT',
}

def adapter_scene(which, *, selective=False):
    """Computed terminal pairing of the independently source-written oligos.

    The declared pairing is a narrow relation, not a preferred-fold prediction or
    evidence that every source reagent molecule was annealed. No 5′ phosphate is
    added to an order that does not print one.
    """
    if which not in (1, 2):
        raise ValueError('adapter must be 1 or 2')
    sense = illumina.TRUSEQ_READ1 if which == 1 else illumina.TRUSEQ_READ2
    end = ('TNNN' if which == 1 else 'CNNN') if selective else 'NNNN'
    sc = Scene()
    sc.strand('sense', [Segment('outer arm', sense[:-len(CORE)]), Segment('paired core', CORE)])
    sc.anneal('antisense', [Segment('cohesive synthesis pattern', end),
                           Segment('paired core complement', ANTI_CORE)],
              to='sense', pair=('paired core complement', 'paired core'), mod3='AminoC6')
    sc.overhang('antisense', 'cohesive synthesis pattern', polarity='5-prime',
                text='Source ordering: four-base 5′ end; not a UMI')
    return sc

def ligation_geometry():
    return AdapterLigationGeometry.known(
        'cohesive', adapter_scene(1), adapter_scene(2),
        statement=('TableS5 source adapter orders and overview2.1 support the four-base '
                   '5′ adapter ends. Figure1 caption instead says3′; both are recorded. '
                   'Genomic junction identities, terminal phosphorylation and complete '
                   'sealed products remain unassembled.'))

def primer_scene(name):
    sc = Scene()
    if name == 'index':
        sc.strand(name, [Segment('P7', illumina.P7),
                        Segment('sample index', 'X' * INDEX_WIDTH, placeholder=True,
                                feature=feature('methylrad.sample', 'sample_index', 'unknown',
                                                note='Six source X positions; study-specific value unassigned.')),
                        Segment('source primer binding tail', INDEX_TAIL)])
    else:
        sc.strand(name, [Segment(name, {'p1': P1, 'p2': P2, 'p3': P3}[name])])
    return sc

def protocol_graph():
    g = ProtocolGraph('MethylRAD standard FspEI library: partial preparation',
                      authority=AUTHORITY, representation='published_instance')
    x = g.input('genomic_dna', 'Actual Col-0 aerial-tissue CTAB-extracted genomic DNA stock',
                K.BIOLOGICAL, primary=True, chemical_state={
                    'scope': 'Prepared DNA starting boundary; plant growth/CTAB recipe not reconstructed.',
                    'input': 'Source1–200ng range; chosen aliquot quantity unassigned.',
                    'modification': 'FspEI-dependent modified-cytosine enrichment; no separate5mC/5hmC resolution.'})
    def stock(key, label, kind=K.CATALYST, role=R.CATALYST, **kw):
        return g.input(key, label, kind, **kw), role
    def condition(key, label):
        return stock(key, label, K.CONDITION, R.CONDITION)
    def operate(key, label, kind, extra=(), *, reason=None, geometry=None):
        nonlocal x
        y = g.product(key + '_product', label + '; occupied substrates/products and yield unassembled', K.MOLECULAR)
        kw = {}
        if kind == S.MAGIC:
            kw['reason'] = reason or ('Positive published operation; complete occupied molecular products, '
                                      'historical formulation and recovery remain unassembled.')
        if kind == S.GENERIC:
            kw['changes_molecule'] = False
        if geometry is not None:
            kw['geometry'] = geometry
        g.react(key, label, kind, ((x, R.TEMPLATE),) + tuple(extra), (y,), **kw)
        x = y
    def split(key, selected, remainder, *, qc=False):
        nonlocal x
        y = g.product(key + '_selected', selected, K.MOLECULAR)
        r = g.product(key + '_remainder', remainder, K.MOLECULAR)
        g.react(key, selected, S.SPLIT, ((x, R.TEMPLATE),), (y, r))
        if qc:
            g.mark_terminal(Terminal(y.id, D.QC, 'Actual analytical material excluded from production.'))
            x = r
        else:
            g.mark_terminal(Terminal(r.id, D.UNRESOLVED, remainder))
            x = y
    def cleanup(key, label, extra=()):
        nonlocal x
        y = g.product(key + '_retained', label + ': actual retained fraction; recovery/carryover unassigned', K.MOLECULAR)
        waste = g.product(key + '_waste', 'Actual separated workup fractions', K.WASTE)
        g.react(key, label, S.CLEANUP, ((x, R.TEMPLATE), condition(key + '_condition', label)) + tuple(extra), (y, waste))
        g.discard(waste, 'Source workup fractions excluded without a claim of perfect molecular selection.')
        x = y
    def pcr(key, label, first, second, programme):
        # Methods names separate enzyme/NT/buffer additions, not an invented premix.
        supplies = (stock(key + '_enzyme', 'Actual0.4U Phusion high-fidelity DNA polymerase'),
                    stock(key + '_nt', 'Actual0.3mM dNTP preparation in20ul PCR', K.MOLECULAR, R.STRUCTURAL),
                    condition(key + '_buffer', 'Actual1x Phusion HF buffer; formulation unassembled'),
                    stock(key + '_forward', first + ': source0.2uM primer', K.OLIGO, R.STRUCTURAL, state=primer_scene(first)),
                    stock(key + '_reverse', second + ': source0.2uM primer', K.OLIGO, R.STRUCTURAL, state=primer_scene(second)),
                    condition(key + '_programme', programme))
        operate(key, label, S.MAGIC, supplies,
                reason='Source-positive PCR and whole primer envelopes known; occupied genomic templates, exact junctions/ends and whole products unassembled. No source-external library constructor.')

    split('dna_allocation', 'Actual DNA aliquot for15ul digest; mass within source1–200ng range unassigned',
          'Untransferred extracted genomic DNA stock')
    operate('digest', 'FspEI digestion37C4h in15ul', S.MAGIC,
            (stock('fspei', 'Actual4U FspEI NEB'), condition('digest_conditions', 'Source37C4h; complete buffer/activator allocation unassembled')),
            reason=('SourceFspEI cutsmodifiedDNA N12/N16; eligible symmetric sites yield short tags alongside other genomic classes. '
                    'Genomic site occupancy/distribution, activator/buffer and complete products unassembled; not a digestion simulation.'))
    split('digest_qc', 'Actual5ul digest analytical1percent agarose gel aliquot',
          'Actual10ul digest remains for production after5of15ul QC removal', qc=True)
    master = g.input('ligation_master', 'Actual10ul source ligation master added to remaining digest', K.MOLECULAR,
                     chemical_state={'stock': 'One actual added master; component quantities not extra fresh additions.',
                                     'concentrations': 'Source0.2uM each adapter/1mM ATP/800U ligase; master-versus-final concentration convention unassigned.'})
    a1 = g.product('adap1', 'Standard Adap-1 duplex components from actual ligation master', K.OLIGO, state=adapter_scene(1))
    a2 = g.product('adap2', 'Standard Adap-2 duplex components from actual ligation master', K.OLIGO, state=adapter_scene(2))
    atp = g.product('ligation_atp', 'ATP component from actual10ul ligation master', K.MOLECULAR)
    ligase = g.product('ligase', 'Actual800U T4 DNA ligase NEB in same10ul master', K.CATALYST)
    g.react('master_components', 'Copresent source master components; not physical separation or extra doses',
            S.GENERIC, ((master, R.TEMPLATE),), (a1, a2, atp, ligase), changes_molecule=False)
    operate('ligate', 'Source cohesive adapter ligation4C6–12h', S.LIGATION,
            ((a1, R.STRUCTURAL), (a2, R.STRUCTURAL), (atp, R.STRUCTURAL), (ligase, R.CATALYST),
             condition('ligation_temperature', '4C6–12h; annealing/end phosphorylation not added by assumption')),
            geometry=ligation_geometry())
    split('pcr1_allocation', 'Actual7ul ligation mixture selected for20ul PCR1',
          'Untransferred ligation mixture; calculated13ul if source15−5+10ul allocations are exhaustive')
    programme = '98C5s/60C20s/72C10s cycles; final72C5min'
    pcr('pcr1', 'First amplification source16–22cycles', 'p1', 'p2', programme + ';16–22cycles')
    split('gel_loading', 'Actual PCR product gel-loading allocation, quantity unassigned', 'Untransferred PCR1 stock')
    cleanup('gel', '8percent PAGE approximately100bp target-band excision',
            (stock('page_reagents', 'Actual PAGE workup reagent allocation; formulation unassembled', K.SELECTION, R.SELECTION),))
    cleanup('diffusion', 'Source nuclease-free-water diffusion6–12h4C',
            (stock('diffusion_water', 'Actual nuclease-free water allocation, volume unassigned', K.MOLECULAR, R.STRUCTURAL),))
    split('index_template', 'Actual25ng gel-recovered DNA for20ul indexPCR', 'Untransferred gel-recovered DNA stock')
    pcr('index_pcr', 'Second amplification/index introduction source5–7cycles', 'p3', 'index', programme + ';5–7cycles')
    cleanup('column', 'QIAquick PCR purification',
            (stock('column_reagents', 'Actual source QIAquick purification preparation; historical components unassembled', K.SELECTION, R.SELECTION),))
    split('sample_contribution', 'Actual selected indexed sample contribution to sequencing pool, quantity unassigned', 'Untransferred indexed sample library')
    other = g.input('other_indexed_samples', 'Other actual indexed sample contributions perFig1; amounts/bodies unassembled', K.MOLECULAR)
    pool = g.product('sequencing_pool', 'Actual pooled indexed sample libraries', K.MOLECULAR)
    g.react('pool', 'Figure1 physical pooling after sample indexing', S.POOL, ((x, R.TEMPLATE), (other, R.TEMPLATE)), (pool,))
    x = pool
    split('load_allocation', 'Actual pooled loading aliquot, quantity unassigned', 'Untransferred prepared sequencing pool')
    lib = g.product('library', 'Source standard MethylRAD sequencing library; genomic products and end combinations unassembled', K.LIBRARY)
    g.react('load', 'Source HiSeq2000 single-end36 sequencing', S.GENERIC, ((x, R.TEMPLATE),), (lib,), changes_molecule=False)
    g.sequence(lib, SequencingPlan(read_lengths={'R1': READ_LENGTH}, limitation='SourceSE36 only; index primer has6X positions but index acquisition cycles/orientation and complete primer landing unassembled.'))
    g.add_source_item(SourceItem('main', 'PMC4680569 Methods5.1–5.3/overview2.1/Fig1', 'Source-positive preparation; no missing-information claim.'))
    g.add_source_item(SourceItem('s5', 'Published rsob150130supp1.pdf TableS5 actualp6', 'Nine unique whole ordering envelopes including alternative RTR pair; selected graph standardNNNN only.'))
    g.annotate(SourceDiscrepancy('ligate', 'Published overhang polarity wording',
        (SourceClaim('Published overview2.1/XML70', 'Four-base5prime overhangs'),
         SourceClaim('Published Figure1 caption/XML74/PDFp3', 'Four-base3prime overhangs')),
        'Adapter scenes retain the source-written5primeNNNN and computedterminal13bp pairing; genomic cut products remain unassembled.',
        'RenderedTableS5 placesNNNN at the antisense5prime end; overview2.1 also specifies5prime.',
        'The inspected published caption conflicts with the main wording; complete enzymatic product/end chemistry is not established by this adapter relation.'))
    return g

def selftest():
    from checks import Check, Source, have
    from storage import sources_root
    c = Check()
    actual = {'adap1': illumina.TRUSEQ_READ1, 'adap2': illumina.TRUSEQ_READ2,
              'anti_std': 'NNNN' + ANTI_CORE + '(AminoC6)',
              'anti_rtr1': 'TNNN' + ANTI_CORE + '(AminoC6)',
              'anti_rtr2': 'CNNN' + ANTI_CORE + '(AminoC6)',
              'p1': P1, 'p2': P2, 'p3': P3, 'index': INDEX_ORDER}
    source = sources_root() / Path(__file__).resolve().parents[1].name / 'rsob.150130_supp_rsob150130supp1.pdf.txt'
    if have(c, Source(source, 'https://pmc.ncbi.nlm.nih.gov/articles/PMC4680569/')):
        text = source.read_text()
        for name, written in SOURCE_ORDERS.items():
            c('TableS5 whole transcription ' + name, actual[name], written)
            c('Actual source contains whole ordering ' + name, written in text, True)
    g = protocol_graph()
    c('Analytical5ul does not rejoin production', g.steps['ligate'].inputs[0].material, 'digest_qc_remainder')
    c('Ligation master allocated once', tuple(x.material for x in g.steps['master_components'].inputs), ('ligation_master',))
    c('FirstPCR uses real7ul allocation', g.steps['pcr1'].inputs[0].material, 'pcr1_allocation_selected')
    c('IndexPCR uses separate25ng gel product', g.steps['index_pcr'].inputs[0].material, 'index_template_selected')
    c('Only selected pooled aliquot enters sequencer', g.steps['load'].inputs[0].material, 'load_allocation_selected')
    c('Index X positions are sample identifier, not UMI', primer_scene('index').strands['index'].segs[1].feature.role, 'sample_index')
    c('Index acquisition not inferred from primer width', g.sequencing['library'].read_lengths, {'R1': READ_LENGTH})
    c('Full library remains unassembled', g.materials['library'].state, None)
    c('Partial preparation cannot be published as complete', bool(g.publication_problems()), True)
    c.report()

if __name__ == '__main__':
    selftest()
