"""Published LM-Seq, selected 100ng/overnight/SE branch; private partial model.

Known physical oligos do not establish occupied cDNA or full PCR/read geometry.
"""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'lib'))
from chemdraw import Scene, Segment, revcomp, feature
from illumina import TRUSEQ_SMALL_RNA_RA3, INDEX2_PRIMER_RC, NEBNEXT_UNIVERSAL_PRIMER, P7
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (ProtocolGraph, MaterialKind as K, InputRole as R,
                            StepKind as S, Terminal, Disposition, SequencingPlan, SourceItem)

AUTHORITY = 'Published PMC4381617 main55/73/83–95/Fig1p2; published SI Methods1–8/Table2; two-original closed source pass.'
RT_TAIL = revcomp(TRUSEQ_SMALL_RNA_RA3)
INDEX_EXTENSION = 'GTGACTGGAGTT'
ADAPTER_LENGTH = len(INDEX2_PRIMER_RC) + 1
EXAMPLE_INDEX = 'TTGGTTAGGT'  # Published Table2 ZH001; not selected study sample.
SOURCE_OLIGO_EXCEPTIONS = {
    'hzg885': 'Published main Methods91 prints full /5phos/AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTddC. The shared scraper strips terminalddC entirely and reports only33nt body; physical donor retains appendedC as34th nucleotide. Local selftest checks the complete modified source ordering and physical terminalTC/length/blocker relation. Narrow text-parser exception only, not a source correction or chemistry exemption.'
}

@dataclass(frozen=True)
class PhysicalDonor:
    body: str = INDEX2_PRIMER_RC
    terminal_base: str = 'C'
    five_prime: str = 'P'
    three_prime: str = 'ddC'

    def __post_init__(self):
        if (self.body != INDEX2_PRIMER_RC or self.terminal_base != 'C'
                or self.five_prime != 'P' or self.three_prime != 'ddC'):
            raise ValueError('This constructor represents only source HZG885 with appended terminal ddC')

    @property
    def sequence(self):
        return self.body + self.terminal_base

    def accept_at_three_oh(self):
        raise ValueError('Physical HZG885 terminal ddC has no 3primeOH; it cannot be this ligation acceptor')

def rt_primer_scene():
    sc = Scene()
    sc.strand('HZG883', [Segment('RT tail', RT_TAIL),
                         Segment('random priming hexamer', 'NNNNNN',
                                 note='Priming population, not a molecule identifier or UMI')])
    return sc

def adapter_scene():
    donor = PhysicalDonor()
    sc = Scene()
    sc.strand('HZG885', [Segment('donor body', donor.body),
                         Segment('appended terminal ddC', donor.terminal_base,
                                 note='Additional nucleotide with dideoxyribose; not copied as a sugar modification')],
              mod5='P', mod3='ddC; no 3primeOH')
    sc.ligation_end('HZG885', 'donor body', kind='single-strand', polarity='5-prime', owner='adapter', text='Source donor 5prime phosphate')
    return sc

def intended_acceptor_scene():
    """Inferred end interpretation by default; never a source-known occupied cDNA."""
    sc = Scene()
    sc.strand('intended cDNA', [Segment('unassembled cDNA', '[unassembled cDNA]',
                                      placeholder=True, inferred=True)], mod3='INFERRED 3primeOH')
    sc.ligation_end('intended cDNA', 'unassembled cDNA', kind='single-strand',
                    polarity='3-prime', owner='substrate', text='INFERRED intended cDNA 3primeOH acceptor')
    return sc

def adapter_geometry():
    return AdapterLigationGeometry.known('single-strand', adapter_scene(), intended_acceptor_scene(),
        statement='Narrow source-known HZG885 donor5P/terminalddC; intended cDNA3OH is derived and marked inferred. Published Fig1caption says5prime cDNA; actual occupied substrate/whole junction remain unassembled.')

def index_primer_scene(index='NNNNNNNNNN'):
    index = index.upper()
    if len(index) != 10 or set(index) - set('ACGTN'):
        raise ValueError('LM-Seq sample index must have ten DNA positions')
    sc = Scene()
    sc.strand('index primer', [Segment('P7', P7),
        Segment('sample index', index, placeholder='N' in index, feature=feature('lmseq.sample_index', 'sample_index',
                'unknown' if 'N' in index else 'fixed', note='Primer ordering; acquired index orientation unassembled')),
        Segment('source extension', INDEX_EXTENSION), Segment('RT tail', RT_TAIL)])
    return sc

def lc056_scene():
    sc = Scene(); sc.strand('LC056', [Segment('source LC056', NEBNEXT_UNIVERSAL_PRIMER)]); return sc

def unavailable_scene():
    sc = Scene(); sc.strand('unassembled', [Segment('whole copied products', '[unassembled products]', placeholder=True)]); return sc

def sequence_of(sc):
    return ''.join(s.top for strand in sc.strands.values() for s in strand.segs)

def protocol_graph():
    g = ProtocolGraph('LM-Seq selected 100ng overnight ligation and SE51 library', authority=AUTHORITY,
                      representation='published_instance')
    x = g.input('total_rna', 'Source100ng purified totalRNA; actual transcripts/sequence population unassigned', K.MOLECULAR, primary=True)

    def reagent(key, label, kind=K.CONDITION, role=R.CONDITION, **kw):
        return g.input(key, label, kind, **kw), role

    def op(key, label, kind, inputs=(), product='Source reaction mixture; whole molecules/off-products/occupancy and recovery unassembled', **kw):
        nonlocal x
        y = g.product(key + '_product', product, K.MOLECULAR)
        if kind == S.GENERIC: kw['changes_molecule'] = False
        g.react(key, label, kind, ((x, R.TEMPLATE),) + tuple(inputs), (y,), **kw); x = y

    def cleanup(key, label):
        nonlocal x
        kept = g.product(key + '_kept', 'Source physically retained fraction; copresent/off-product recovery unassigned', K.MOLECULAR)
        waste = g.product(key + '_waste', 'Actual discarded workup fractions; imperfect recovery and molecular carryover unassigned', K.WASTE)
        g.react(key, label, S.CLEANUP, ((x, R.TEMPLATE), reagent(key + '_buffer', label)), (kept, waste))
        g.discard(waste, 'Physical workup, not ideal class-specific removal'); x = kept

    def transfer(key, label):
        nonlocal x
        used = g.product(key + '_used', label, K.MOLECULAR)
        remainder = g.product(key + '_remainder', 'Actual untransferred eluate/reaction stock; subsequent fate unassigned', K.MOLECULAR)
        g.react(key, label, S.SPLIT, ((x, R.TEMPLATE),), (used, remainder))
        g.mark_terminal(Terminal(remainder.id, Disposition.UNRESOLVED, 'Untransferred remainder not automatically discarded')); x = used

    def qc(key, label):
        nonlocal x
        stock = g.product(key + '_stock', 'Actual retained production stock after QC aliquot; amount unassigned', K.MOLECULAR)
        aliquot = g.product(key + '_aliquot', label, K.MOLECULAR)
        g.react(key, label, S.SPLIT, ((x, R.TEMPLATE),), (stock, aliquot))
        g.mark_terminal(Terminal(aliquot.id, Disposition.QC, label + '; assay preparations/volumes unassembled')); x = stock

    op('polya_bind', 'Source15ulNEBNextpolyA beads prewashed twice100ul2xbindingbuffer/resuspended50ul;50ulRNA+50ulbeads/one postbinding wash', S.MAGIC,
       (reagent('polya_beads', 'Actual prepared15ulNEBNextoligo-dTbead aliquot:prewash twice100ul2xbindingbuffer/resuspend50ul2xbindingbuffer; factory bead manufacture opaque', K.MOLECULAR, R.STRUCTURAL),),
       reason='Source polyA capture positive; historical kit binding/incubation and full populations unassembled.')
    cleanup('polya_wash', 'Source magnetic mRNA selection and one postbinding wash; bead-bound fraction retained')
    primer = g.input('hzg883', 'Source0.75ul20uMHZG883 in8.25ulEPF stock', K.OLIGO, state=rt_primer_scene())
    op('fragment_elute', 'SourceEPF5ulwater/2.5ul5xbuffer/.75ulprimer;85C7min/coolice', S.MAGIC,
       ((primer, R.STRUCTURAL), reagent('epf_buffer', 'Source8.25ulEPF stock conditions')),
       reason='Source heat fragmentation/primer carryover positive; exact RNA cleavage ends/sites and occupied priming unassembled.')
    # Remove magnetic beads, preserving actual remaining supernatant separately.
    sup = g.product('epf_supernatant', 'Actual fragmentedRNA/previously added HZG883 in supernatant; amounts unassigned', K.MOLECULAR)
    bead = g.product('epf_beads', 'Magnet-retained beads and possible residual bound molecules; fate unassigned', K.MOLECULAR)
    g.react('epf_magnet', 'Source magnet separation after cooling', S.SPLIT, ((x, R.TEMPLATE),), (sup, bead))
    g.mark_terminal(Terminal(bead.id, Disposition.UNRESOLVED, 'Beads physically excluded from RT; molecular loss/fate unassigned')); x = sup
    transfer('epf_transfer', 'Source actual6.5ul supernatant transfer from8.25ulEPF')
    op('reverse_transcribe', 'Source10ulRT:6.5ulcarryover+1water+.25DTT+1dNTP+.25inhibitor+1SmartScribe;23C10min/42C30min', S.MAGIC,
       (reagent('smartscribe', 'Source1ul SmartScribe RT100U', K.CATALYST, R.CATALYST),
        reagent('rt_nts', 'Source1ul10mM-each dNTP', K.MOLECULAR, R.STRUCTURAL),
        reagent('rt_additives', 'Source water/DTT/RNase inhibitor; no fresh primer addition')),
       reason='Random priming and RT are source-positive; actual RNA/primer occupancy/full copied products unassembled.',
       product='RT mixture: cDNA/RNA hybrids, free RNA, retained HZG883 and unreacted/copied classes, not pure cDNA')
    op('remove_rna', 'Source+.5ulRNaseH+.5ulRNaseA→11ul;37C15min/95C5min', S.MAGIC,
       (reagent('rnases', 'Source RNaseA/H catalysts', K.CATALYST, R.CATALYST), reagent('rnase_conditions', 'Source11ulRNA digestion/heating conditions')),
       reason='Source RNA hydrolysis positive; exact full cDNA/substrate/end/population states unassembled.')
    cleanup('cdna_clean', 'Source20ulAMPure on11ul;15minRT/wash twice200ul80percentEtOH/dry/elute12ulEB')
    transfer('cdna_transfer', 'Source10ul transfer from12ul eluate')
    donor = g.input('hzg885', 'Source1.5ul20uM physical HZG885 with appendedddC', K.ADAPTER,
        state=adapter_scene(), chemical_state={'five_prime': 'P', 'three_prime': 'ddC',
           'additional_terminal_modified_base': (len(INDEX2_PRIMER_RC), 'ddC'),
           'source_ordering': '/5phos/' + INDEX2_PRIMER_RC + 'ddC'})
    op('ligate', 'Source30ul:10cDNA/3buffer/9.5PEG50/3ATP10mM/1.5donor/3RNA1;22C12–16h', S.LIGATION,
       ((donor, R.STRUCTURAL), reagent('rna_ligase1', 'Source3ul T4RNA ligase1 stock10000Uml', K.CATALYST, R.CATALYST),
        reagent('ligation_atp', 'Source3ul10mM ATP', K.MOLECULAR, R.STRUCTURAL),
        reagent('ligation_conditions', 'Source3ul10xbuffer/9.5ul50percentPEG8000')),
       geometry=adapter_geometry(), product='Intended cDNA3OH-to-donor5P ligation population (derived direction), unligated cDNA/donor/free primers/off-products remain; whole junction unassembled')
    op('stop_ligation', 'Source+20ul10mMTris8/15mMEDTA→50ul', S.GENERIC,
       (reagent('stop_buffer', 'Source20ulTris/EDTA stop'),))
    cleanup('ligation_clean1', 'Source90ulAMPure on50ul;15minRT/twoEtOH washes/elute22ul')
    transfer('ligation_transfer1', 'Source20ul transfer from22ul eluate')
    cleanup('ligation_clean2', 'Source36ulAMPure on20ul;15minRT/twoEtOH washes/elute22ul')
    transfer('ligation_transfer2', 'Source20ul transfer from22ul eluate')
    f = g.input('lc056', 'Source.5ul20uMLC056', K.OLIGO, state=lc056_scene())
    r = g.input('index_primer', 'Source.5ul20uM10nt sample-index primer; actualTable2choice unassigned', K.OLIGO, state=index_primer_scene())
    premix = g.input('failsafe_premix', 'Actual source30ul FailSafePCR PremixE prepared stock', K.MOLECULAR)
    premix_nts = g.product('pcr_nts', 'Nucleotide component of the same actual30ulPremixE; separate dose unassigned', K.MOLECULAR)
    premix_buffer = g.product('pcr_premix_buffer', 'Buffer/condition component of the same actual30ulPremixE', K.CONDITION)
    g.react('premix_components', 'Represent components of source30ulPremixE, not fresh independent additions', S.GENERIC,
            ((premix, R.TEMPLATE),), (premix_nts, premix_buffer), changes_molecule=False)
    op('pcr', 'Source60.1ulFailSafePCR;95C1min/18cycles95C10s55C10s68C3min/68C7min', S.PCR,
       ((f, R.STRUCTURAL), (r, R.STRUCTURAL),
        reagent('failsafe', 'Source.6ulFailSafePCR enzyme', K.CATALYST, R.CATALYST),
        (premix_nts, R.STRUCTURAL), (premix_buffer, R.CONDITION),
        reagent('pcr_water', 'Actual source8.5ulwater addition;20ul transferred material already supplies template')),
       geometry=PCRGeometry.unavailable(unavailable_scene(), reason='Source full primers known; complete occupied template/landing/strand products unassembled. Physical donor ddC is not copied as a sugar modification.'),
       product='PCR copiedDNA plus remaining native donor/cDNA/primers/off-products; physical3ddC applies only to source original donor, not newly made DNA')
    cleanup('pcr_clean', 'Source60ulAMPure on60.1ulPCR (nominal1x);elute22ul')
    transfer('pcr_transfer', 'Source20ul transfer from22ul eluate')
    qc('single_qc', 'Source separate Qubit/PicoGreen concentration aliquot')
    other = g.input('other_indexed_libraries', 'Other actually prepared indexed library contributions; chosen samples/molecules/amounts unassigned', K.MOLECULAR)
    op('pool', 'Source physically pool indexed libraries; selectedless-than48 branch', S.POOL, ((other, R.TEMPLATE),))
    cleanup('pool_clean', 'Source nominal1xAMPure pool cleanup; selectedless-than48:elute14ul')
    transfer('pool_transfer', 'Source12ul transfer from14ul eluate;48/96 alternative27/25 outside branch')
    qc('pool_qc', 'Source separate pooled QubitHS/BioAnalyzerHS aliquots; volumes unassigned')
    transfer('loading', 'Actual sequencer-loading fraction of pooled stock; volume unassigned')
    lib = g.product('library', 'Selected LM-Seq pooled indexed Illumina library; whole molecules unassembled', K.LIBRARY)
    g.react('load', 'Source HiSeq2500 SE51 plus10base index read', S.GENERIC, ((x, R.TEMPLATE),), (lib,), changes_molecule=False)
    g.sequence(lib, SequencingPlan(limitation='Published SE51/index10 and SmallRNA index-primer name; exact primer ordering/landing/acquired index orientation unassembled. N6 is random primer, not UMI; analysis filters are computational, not physical cleanup.'))
    g.add_source_item(SourceItem('ligation_direction', 'Published mainFig1caption PDFp2 says donor ligates to5prime ss-cDNA; main55/91 donor5P+terminalddC/T4RNA1 support derived cDNA3OH acceptor. Both claims retained, no silent repair.', 'Model uses inferred intended cDNA3OH-to-donor5P junction; full actual product and resolution of caption discrepancy remain unassembled.'))
    g.add_source_item(SourceItem('excluded_branches', '100ng/SmartScribe/overnight/<48pool/SE selected; proposedPE/rRNAdepletion/lowinput and short-ligation comparisons not inherited.', 'Two originals only; no vendor/predecessor/revision/correction original newly opened.'))
    return g

def selftest():
    from checks import Check, Source, have
    from storage import sources_root
    c = Check(); g = protocol_graph()
    root = sources_root() / HERE.parent.name
    main = root / 'srep09570_PMC4381617.xml.txt'
    si = root / 'srep09570_supp_srep09570-s1.pdf.txt'
    if have(c, Source(main, 'https://doi.org/10.1038/srep09570')):
        t = ''.join(main.read_text().split()).upper()
        for name, ordering in [('HZG883', sequence_of(rt_primer_scene())),
                               ('HZG885-modified', '/5phos/' + INDEX2_PRIMER_RC + 'ddC'),
                               ('LC056', sequence_of(lc056_scene())),
                               ('Index-template', sequence_of(index_primer_scene()))]:
            c('Published whole source ordering ' + name, ordering.upper() in t, True)
    if have(c, Source(si, 'https://doi.org/10.1038/srep09570')):
        t = ''.join(si.read_text().split()).upper()
        c('PublishedTable2 ZH001 ordering', sequence_of(index_primer_scene(EXAMPLE_INDEX)) in t, True)
    donor = PhysicalDonor()
    c('AppendedddC adds physical nucleotide', donor.sequence[-2:], 'TC')
    c('Whole physical donor length', len(sequence_of(adapter_scene())), 34)
    try: donor.accept_at_three_oh(); rejected = False
    except ValueError: rejected = True
    c('Source blocked donor cannot be selected3OH acceptor', rejected, True)
    c('IntendedcDNA constructor inferred by default', intended_acceptor_scene().strands['intended cDNA'].segs[0].inferred, True)
    c('Random hexamer has no identifier feature', rt_primer_scene().strands['HZG883'].segs[-1].feature, None)
    c('Sample index has structured meaning', index_primer_scene().strands['index primer'].segs[1].feature.role, 'sample_index')
    try: index_primer_scene('ACGT'); rejected = False
    except ValueError: rejected = True
    c('Reject malformed sampleindex', rejected, True)
    c('RT consumes actual EPF transfer', g.steps['reverse_transcribe'].inputs[0].material, 'epf_transfer_used')
    c('No fresh RT primer addition', any(p.material == 'hzg883' for p in g.steps['reverse_transcribe'].inputs), False)
    c('Ligation follows actual10ul transfer', g.steps['ligate'].inputs[0].material, 'cdna_transfer_used')
    c('Ligation known narrow donorend', g.steps['ligate'].geometry.availability, 'known')
    c('PCR consumes secondcleanup transfer', g.steps['pcr'].inputs[0].material, 'ligation_transfer2_used')
    c('Pool uses postQC stock', g.steps['pool'].inputs[0].material, 'single_qc_stock')
    c('QC is not sequencing production', g.terminals['pool_qc_aliquot'].disposition, Disposition.QC)
    c('Loading remainder not discarded', g.terminals['loading_remainder'].disposition, Disposition.UNRESOLVED)
    c('Whole PCR copied products unassembled', g.materials['pcr_product'].state, None)
    c('Private incomplete graph', bool(g.publication_problems()), True)
    c.report()

if __name__ == '__main__': selftest()
