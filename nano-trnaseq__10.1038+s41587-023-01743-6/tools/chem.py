"""Source-owned Nano-tRNAseq padding/RNA002 preparation with source-qualified public views."""
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'lib'))
from chemdraw import Scene, Segment, revcomp
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import (ProtocolGraph, MaterialKind as K, InputRole as R,
                            StepKind as S, Terminal, Disposition, SequencingPlan, SourceItem, SourceCoverage)

AUTHORITY = ('Published PMC10791586 Methods202/206, Results87 and actualFig1bPDFp3; '
             'published MOESM3 TableS25 original sheet25 B4:B9. Two originals, no newer kit inheritance.')
RAW = {
    'rna5': 'rCrCrUrArArGrArGrCrArArGrArArGrArArGrCrCrUrGrGrN',
    'rna3': '/5Phos/rGrGrCrUrUrCrUrUrCrUrUrGrCrUrCrUrUrArGrGrArArArArArArArArArAAAA',
    'rta_a': '/5Phos/GGCTTCTTCTTGCTCTTAGGTAGTAGGTTC',
    'rta_b': 'GAGGCGAGCGGTCAATTTTCCTAAGAGCAAGAAGAAGCCTTTTTTTTTT',
}
COMPARATORS = {
    'rna_only_3': '/5Phos/rGrGrCrUrUrCrUrUrCrUrUrGrCrUrCrUrUrArGrGrArArArArArArArArArA',
    'rt_test': '/5Phos/ACTTGCCTGTCGCTCTATCTTCTTTTTTTTTTTTTTTTTTTTVN',
}


def residues(raw):
    """Keep r-prefix sugar identity; an rA followed by AAAA means one RNA plus three DNA."""
    text = raw.removeprefix('/5Phos/')
    tokens = re.findall(r'r[ACGUN]|[ACGTNV]', text)
    if ''.join(tokens) != text:
        raise ValueError('Unrepresented source oligo notation')
    return tuple(('RNA', t[1]) if t.startswith('r') else ('DNA', t) for t in tokens)


def bases(which):
    return ''.join(base.replace('U', 'T') for sugar, base in residues(RAW[which]))


def sugar_count(which, sugar):
    return sum(s == sugar for s, base in residues(RAW[which]))


def descriptor(which):
    return {'source_modified_ordering': RAW[which], 'residues': residues(RAW[which]),
            'five_prime': 'P' if RAW[which].startswith('/5Phos/') else 'unassigned',
            'scope': 'Source oligo chemistry only; exact occupied joining product is unassembled.'}


def oligo_scene(which):
    sc = Scene()
    parts = []
    for i, (sugar, base) in enumerate(residues(RAW[which])):
        if parts and parts[-1][0] == sugar:
            parts[-1][1] += base.replace('U', 'T')
        else:
            parts.append([sugar, base.replace('U', 'T')])
    sc.strand('oligo', [Segment(sugar + ' block ' + str(i), seq,
                              note='Source ' + sugar + ' sugar; T displays source U when RNA; RNA5 terminal N is source RNA discriminator-binding class, not UMI.')
                        for i, (sugar, seq) in enumerate(parts)],
              mod5='P' if RAW[which].startswith('/5Phos/') else '')
    return sc


def padding_scene():
    """Compute the published adapter duplex; no complete tRNA or sealed product supplied."""
    a, b = bases('rna5'), bases('rna3')
    # Source Figure1b shows the 5' ends of the adapters paired. Determine that core.
    candidates = [n for n in range(1, min(len(a), len(b)) + 1)
                  if a[:n] == revcomp(b[:n])]
    if not candidates:
        raise ValueError('Published padding oligos lack their proposed paired core')
    n = max(candidates)
    sc = Scene()
    sc.strand('RNA5', [Segment('paired RNA5 core', a[:n], note='RNA sugar'),
                        Segment('NCCA binding arm', a[n:],
                                note='Source mature-tRNA/discriminator binding class; not a UMI')])
    # The RNA:DNA boundary comes from the full modified ordering, not the drawing.
    rna_count = sum(sugar == 'RNA' for sugar, base in residues(RAW['rna3']))
    sc.anneal('RNA3', [Segment('paired RNA3 core', b[:n], note='RNA sugar'),
                        Segment('RNA A tail', b[n:rna_count], note='RNA sugar'),
                        Segment('DNA A tail', b[rna_count:], note='DNA sugar')], mod5='P',
              to='RNA5', pair=('paired RNA3 core', 'paired RNA5 core'))
    sc.ligation_end('RNA3', 'paired RNA3 core', kind='splinted', polarity='5-prime',
                    text='Source 5P RNA donor; complete tRNA acceptor geometry unassembled')
    return sc


def rta_scene():
    a, b = bases('rta_a'), bases('rta_b')
    # Locate the longest exact antiparallel core once, retaining both complete strands.
    candidates = [(n, j) for n in range(1, len(a) + 1)
                  for j in range(len(b) - n + 1) if a[:n] == revcomp(b[j:j+n])]
    n, j = max(candidates)
    sc = Scene()
    sc.strand('RTA A', [Segment('RTA A core', a[:n]), Segment('RTA A tail', a[n:])], mod5='P')
    # Fig1b shows the divergent RTA arms. They are not an additional paired stem.
    sc.anneal('RTA B', [Segment('RTA B 5 arm', b[:j], note='Unpaired Y arm in source Fig1b'),
                        Segment('RTA B core', b[j:j+n]), Segment('RTA B T tail', b[j+n:])],
              to='RTA A', pair=('RTA B core', 'RTA A core'), unpaired=('RTA B 5 arm',))
    sc.ligation_end('RTA A', 'RTA A core', kind='splinted', polarity='5-prime',
                    text='Known source RTA-A 5P donor; RNA polyA register unassembled')
    return sc


def rta_core_scene():
    """A cropped paired-stem view; plain-text duplexes cannot show divergent Y arms."""
    full = rta_scene()
    a = next(s.top for s in full.strands['RTA A'].segs if s.name == 'RTA A core')
    b = next(s.top for s in full.strands['RTA B'].segs if s.name == 'RTA B core')
    sc = Scene()
    sc.strand('A stem fragment', [Segment('A paired stem', a, note='Cropped fragment, tail omitted')])
    sc.anneal('B stem fragment', [Segment('B paired stem', b, note='Internal fragment, both tails omitted')],
              to='A stem fragment', pair=('B paired stem', 'A paired stem'))
    return sc


def protocol_graph():
    g = ProtocolGraph('Nano-tRNAseq selected commercial total yeast tRNA branch',
                      authority=AUTHORITY, representation='published_instance')
    x = g.input('trna', 'Source commercial S.cerevisiae total tRNA AM7119; native heterogeneous molecules',
                K.MOLECULAR, primary=True)

    def inp(key, label, kind=K.CONDITION, role=R.CONDITION):
        return g.input(key, label, kind), role

    def op(key, label, ports=(), product='Source mixture; occupied molecules, yields and off-products unassembled',
           kind=S.MAGIC, **kw):
        nonlocal x
        y = g.product(key + '_product', product, K.MOLECULAR)
        if kind == S.MAGIC:
            kw.setdefault('reason', 'Positive source operation; whole molecules and joining/product geometry unassembled.')
        if kind == S.GENERIC:
            kw['changes_molecule'] = False
        g.react(key, label, kind, ((x, R.TEMPLATE),) + tuple(ports), (y,), **kw)
        x = y

    def clean(key, label):
        nonlocal x
        kept = g.product(key + '_kept', 'Actual retained production fraction', K.MOLECULAR)
        waste = g.product(key + '_waste', 'Actual discarded workup/wash fractions', K.WASTE)
        g.react(key, label, S.CLEANUP, ((x, R.TEMPLATE), inp(key + '_conditions', label)), (kept, waste))
        g.discard(waste, 'Source physical workup; no perfect specificity or yield claimed')
        x = kept

    def aliquot(key, label, qc=False):
        nonlocal x
        used = g.product(key + '_used', label, K.QC if qc else K.MOLECULAR)
        rest = g.product(key + '_rest', 'Other actual stock fraction; amount/fate unassigned', K.MOLECULAR)
        g.react(key, label, S.SPLIT, ((x, R.TEMPLATE),), (used, rest))
        if qc:
            g.mark_terminal(Terminal(used.id, Disposition.QC, label))
            x = rest
        else:
            g.mark_terminal(Terminal(rest.id, Disposition.UNRESOLVED, 'Untransferred stock; amount/fate unassigned'))
            x = used

    def oligo(key):
        return g.input(key, 'Source TableS25 ' + key, K.OLIGO,
                       state=oligo_scene(key), chemical_state=descriptor(key))

    op('deacylate', 'Source 10ul tRNA water /95ul100mMTris9;37C30min',
       (inp('deacylation_buffer', 'Actual100mMTris9/37C30min'),),
       product='Deacylated native tRNA class plus possible remaining charged/off-products; no selected genome sequence')
    clean('deacylation_cleanup', 'Source ZymoR1016 retain>=17nt;1.3xethanol after RNAprepBuffer')
    aliquot('deacylation_qc', 'Actual TapeStationRNAHS aliquot, amount unassigned', qc=True)
    p5, p3 = oligo('rna5'), oligo('rna3')
    pad = g.product('padding_stock', 'Actual preannealed RNA padding stock; unbound species unassigned',
                    K.ADAPTER, state=padding_scene(), chemical_state={'oligos': [descriptor('rna5'), descriptor('rna3')]})
    g.react('padding_anneal', 'Source1:1/50ngul;Tris7.5/NaCl/RNasin;75C15s→25C at.1C/s', S.HYBRIDIZATION,
            ((p5, R.TEMPLATE), (p3, R.STRUCTURAL),
             inp('padding_conditions', '10mMTris7.5/50mMNaCl; prepared stock allocation unassigned'),
             inp('padding_anneal_rnasin','Actual1ul RNasin protein contribution to padding annealing',K.MOLECULAR,R.STRUCTURAL)), (pad,))
    pad_used = g.product('padding_used', 'Actual padding contribution, tRNA:padding1.2:1 using source90nt average', K.ADAPTER, state=padding_scene(), chemical_state={'oligos':[descriptor('rna5'),descriptor('rna3')]})
    pad_rest = g.product('padding_rest', 'Unused prepared padding stock', K.ADAPTER, state=padding_scene(), chemical_state={'oligos':[descriptor('rna5'),descriptor('rna3')]})
    g.react('padding_allocation', 'Take actual padding contribution; stock total/remainder unassigned', S.SPLIT,
            ((pad, R.TEMPLATE),), (pad_used, pad_rest))
    g.mark_terminal(Terminal(pad_rest.id, Disposition.UNRESOLVED, 'Unused adapter stock; fate unassigned'))
    op('padding_ligate', 'Source50ul/20percentPEG8000/1xRnl2buffer/4ul6mgmlRnl2/1ulRNasin;RT2h',
       ((pad_used, R.STRUCTURAL), inp('rnl2', 'Actual4ul6mgml in-houseT4RNAligase2', K.CATALYST, R.CATALYST),
        inp('rnl2_conditions', 'Actual20percent PEG8000/1xRnl2buffer; stock packaging unassigned'),
        inp('padding_ligation_rnasin','Actual fresh1ul RNasin protein contribution to padding ligation',K.MOLECULAR,R.STRUCTURAL)),
       product='Source dual-end padded native tRNA class plus free/partly ligated/off-products; whole junctions unassembled',
       geometry=AdapterLigationGeometry.known('splinted', padding_scene(),
                  statement='Known adapter duplex and RNA3 5P donor. RNA5 5end and whole occupied tRNA ends/junctions not established by this local view.'))
    clean('padding_cleanup', 'Source2xRNACleanXP/RT15min/70percentEtOH; water elution10min')
    aliquot('padded_qc', 'Actual RNAHSQubit aliquot; volume unassigned', qc=True)
    aliquot('padded_input', 'Actual200ng padded tRNA selected; remaining stock not silently consumed')
    a, b = oligo('rta_a'), oligo('rta_b')
    rta = g.product('rta_stock', 'Source preannealed RTA A/B; free components unassigned', K.ADAPTER, state=rta_scene())
    g.react('rta_anneal', 'RTA A/B annealed as padding recipe; no RNA004 substitution', S.HYBRIDIZATION,
            ((a, R.TEMPLATE), (b, R.STRUCTURAL), inp('rta_conditions', 'Source same annealing conditions; amounts unassigned')), (rta,))
    ru = g.product('rta_used', 'Actual8.6pmol RTA to source~4.3pmol/200ng RNA', K.ADAPTER, state=rta_scene(), chemical_state={'oligos':[descriptor('rta_a'),descriptor('rta_b')]})
    rr = g.product('rta_rest', 'Unused preannealed RTA stock', K.ADAPTER, state=rta_scene(), chemical_state={'oligos':[descriptor('rta_a'),descriptor('rta_b')]})
    g.react('rta_allocation', 'ActualRNA:RTA1:2 contribution; unused prepared stock separate', S.SPLIT,
            ((rta, R.TEMPLATE),), (ru, rr))
    g.mark_terminal(Terminal(rr.id, Disposition.UNRESOLVED, 'Unused RTA stock; amount/fate unassigned'))
    op('rta_ligate', 'Source15ul/1xQuickbuffer/1.5ulT4DNA ligase/.5ulRNasin;RT30min',
       ((ru, R.STRUCTURAL), inp('rta_ligase', 'Actual fresh1.5ulT4DNA ligase M0202M', K.CATALYST, R.CATALYST),
        inp('rta_buffer', 'Actual1x QuickLigation buffer'),
        inp('rta_ligation_rnasin','Actual fresh0.5ul RNasin protein contribution to RTA ligation',K.MOLECULAR,R.STRUCTURAL)),
       product='RTA joined to padded nativeRNA plus copresent/unjoined species; full RNA/polyA/RTA register unassembled',
       geometry=AdapterLigationGeometry.known('splinted', rta_scene(),
                  statement='Source RTA A/B whole duplex/5P donor known. PolyA13 versus T10 binding register and complete recipient/product unassembled.'))
    master = g.input('rt_master', 'ONE actual25ulMaxima master:13water+2dNTP+8buffer+2enzyme', K.MOLECULAR)
    enzyme = g.product('rt_enzyme', 'Actual2ulMaximaHMinus from same25ulmaster', K.CATALYST)
    nts = g.product('rt_nts', 'Actual2ul10mMdNTP from same25ulmaster', K.MOLECULAR)
    conditions = g.product('rt_conditions', 'Actual13ulwater/8ul5xMaximabuffer from same25ulmaster', K.CONDITION)
    g.react('rt_components', 'Component inventory of one actual prepared RTmaster; no fresh duplicate additions', S.GENERIC,
            ((master, R.TEMPLATE),), (enzyme, nts, conditions), changes_molecule=False)
    op('protective_rt', 'Add entire25ulmaster DIRECTLY to15ulRTAreaction;60C1h/85C5min/4C',
       ((enzyme, R.CATALYST), (nts, R.STRUCTURAL), (conditions, R.CONDITION)),
       product='NativeRNA with protective cDNA/RT-truncated/uncopied/off-product classes; RNA remains sequencing substrate')
    clean('rt_cleanup', 'Source2xRNACleanXP as prior ligation cleanup; elution unassigned')
    rmx = g.input('rmx', 'Actual6ul historicalRNA002 RMX motor-adapter oligo class; whole strands/mods unassembled', K.ADAPTER)
    op('rmx_ligate', 'Source40ul/1xQuickbuffer/3ulT4DNA ligase/6ulRMX;RT30min',
       ((rmx, R.STRUCTURAL), inp('rmx_ligase', 'Actual fresh3ulT4DNA ligase M0202M; distinct post-cleanup addition', K.CATALYST, R.CATALYST),
        inp('rmx_buffer', 'Actual QuickLigation buffer')),
       product='Source motor-adapted nativeRNA/cDNA-containing library mixture; whole RMX geometry unassembled',
       geometry=AdapterLigationGeometry.unavailable('Known padding/RTA donor geometry retained separately; historical RMX recipient/strand ends and full motor joining geometry unassembled.'))
    clean('rmx_cleanup', 'Source2xRNACleanXP/RT10min;two150ulWSB washes/elute20ulELB10min')
    op('loading_mix', 'Source add17.5ulwater+37.5ulvortexedRRB; keepice',
       (inp('loading_supplies', 'Actual17.5ulwater/37.5ulRRB'),), kind=S.GENERIC)
    aliquot('loading_allocation', 'Actual MinION loading allocation; input volume unassigned')
    lib = g.product('native_library', 'Nano-tRNAseq native RNA loading material; full molecules unassembled', K.LIBRARY)
    g.react('load', 'Source FLO-MIN106QC/priming/loading, historicalSQK-RNA002', S.GENERIC,
            ((x, R.TEMPLATE),), (lib,), changes_molecule=False)
    g.sequence(lib, SequencingPlan(limitation='NativeRNA readout; protective cDNA is not read substrate. Whole tRNA/RMX/occupied primer landing/modified-base readout remains unassembled. No PCR or RNA004 inheritance.'))
    g.add_source_item(SourceItem('selected_branch', 'Methods202/206 commercialtotalyeasttRNA/nativeRNA002 route; Fig1b and TableS25.',
                               'Yeast culture/extraction, synthetic-tRNA factory, alternative A/B/RNA-only strategies and LCMS gel/nucleoside digestion are separate branches, not missing source information.'))
    g.add_source_item(SourceItem('computational_settings', 'Main250:1sadapter/2sstrand; ExtendedData4caption300:2sadapter/1sstrand.',
                               'Conflicting published descriptions retained in note; per-run configuration unresolved. Not physical size selection or a proven actual error.'))
    return _publication_coverage(g)


def selftest():
    from checks import Check, Source, have
    from storage import sources_root
    c = Check()
    p = sources_root() / HERE.parent.name / 's41587-023-01743-6_supp_41587_2023_1743_MOESM3_ESM.xlsx.txt'
    if have(c, Source(p, 'https://doi.org/10.1038/s41587-023-01743-6')):
        table = p.read_text().split('## sheet Table S25\n')[1].split('## sheet ')[0]
        for key, raw in {**RAW, **COMPARATORS}.items():
            c('Whole modified TableS25 ordering ' + key, raw in table, True)
    # Independent source wording distinguishes sugars that raw A-letter counting conflates.
    main = p.parent / 's41587-023-01743-6_PMC10791586.xml.txt'
    if have(c, Source(main, 'https://doi.org/10.1038/s41587-023-01743-6')):
        c('Source says three terminal DNA residues', 'three DNA bases' in main.read_text(), True)
        c('Typed sugar count agrees with source prose', sum(s == 'DNA' for s, b in residues(RAW['rna3'])), 3)
    padding_scene()
    rta_scene()
    rta_core_scene()
    g = protocol_graph()
    c('RNA production proceeds from selected200ng stock', g.steps['rta_ligate'].inputs[0].material, 'padded_input_used')
    c('RT retains entire RTA reaction', g.steps['protective_rt'].inputs[0].material, 'rta_ligate_product')
    c('RT enzyme supplied by the same prepared master', g.producer_of('rt_enzyme').id, 'rt_components')
    c('Motor ligation follows RT cleanup', g.steps['rmx_ligate'].inputs[0].material, 'rt_cleanup_kept')
    c('QC not sequenced', g.terminals['padded_qc_used'].disposition, Disposition.QC)
    c('Unused RTA not silently discarded', g.terminals['rta_rest'].disposition, Disposition.UNRESOLVED)
    c('No whole native RNA product invented', g.materials['native_library'].state, None)
    c('Selected source-qualified publication gate',g.publication_problems(),[])
    c('Whole RNA5 has24RNA residues including source terminalN',sugar_count('rna5','RNA'),24)
    c('RNA5 terminal rN retains RNA sugar',sugar_count('rna5','DNA'),0)
    c('Used and remaining padding preserve same source stock state',g.materials['padding_used'].state is not None and g.materials['padding_rest'].chemical_state==g.materials['padding_used'].chemical_state,True)
    c('Analytical allocations are typed QC',g.materials['padded_qc_used'].kind,K.QC)
    c('Distinct source RNasin additions consumed by their own stages',all(key in {i.material for i in g.steps[step].inputs} and g.materials[key].kind==K.MOLECULAR for key,step in [('padding_anneal_rnasin','padding_anneal'),('padding_ligation_rnasin','padding_ligate'),('rta_ligation_rnasin','rta_ligate')]),True)
    c.report()


PAGE_FILE='nano-trnaseq.html'
SOURCE_TARGETS=(
 ('PublishedmainMethods202: selectedcommercialtotalyeasttRNA/deacylation/Zymo>=17nt and analyticalTapeStation. No upstreamculture/syntheticRNAfactory imported.',
  ('trna','deacylation_buffer','deacylate','deacylation_cleanup_conditions','deacylation_cleanup','deacylation_qc')),
 ('Publishedmain206/Fig1b and actualTableS25B4/B6: source5RNA/3RNA-DNApadding, anneal/allocation/Rnl2/cleanup/Qubit/200nginput. Whole occupied tRNA products unassembled.',
  ('rna5','rna3','padding_anneal_rnasin','padding_ligation_rnasin','padding_conditions','padding_anneal','padding_allocation','rnl2','rnl2_conditions','padding_ligate','padding_cleanup_conditions','padding_cleanup','padded_qc','padded_input')),
 ('Samepublishedmain206/TableS25B7:B8: actualRTAoligos/stock/allocation,15ulRTAreaction thenONE25ulRTmaster addedwithoutcleanup, subsequentRNAClean. Knownoligostocks retained, whole recipient/polyAregister unknown.',
  ('rta_a','rta_b','rta_ligation_rnasin','rta_conditions','rta_anneal','rta_allocation','rta_ligase','rta_buffer','rta_ligate','rt_master','rt_components','protective_rt','rt_cleanup_conditions','rt_cleanup')),
 ('Samepublishedmain206: historicalRNA002 RMX/actual fresh3ulligase/cleanup20ulELB/water+RRB/loading. RMXwholebody/mods and readlanding unassembled; not a newer kitrecipe.',
  ('rmx','rmx_ligase','rmx_buffer','rmx_ligate','rmx_cleanup_conditions','rmx_cleanup','loading_supplies','loading_mix','loading_allocation','load')),
)

def _publication_coverage(g):
    g.cover_source_item(SourceCoverage('selected_branch','modelled','native_library'))
    g.cover_source_item(SourceCoverage('computational_settings','modelled','native_library'))
    for source,targets in SOURCE_TARGETS:
        for target in targets:
            item='reviewed_'+target
            g.add_source_item(SourceItem(item,'Selected source scope: '+target,source))
            g.cover_source_item(SourceCoverage(item,'modelled',target))
    return g.require_publishable()


def render_page():
    from page import head,caveat,table
    from chemdraw import panel
    from protocol_views import assemble_page,WorkflowSection
    g=protocol_graph()
    return assemble_page(g,[head('Nano-tRNAseq — padded native tRNA sequencing'),
      '<main class="wrap" data-protocol-complete="true"><h1>Nano-tRNAseq</h1><p>RNA splint adapters pad native tRNAs before RNA002 motor-adapter preparation. Reverse transcription adds protective cDNA; the RNA remains the nanopore sequencing substrate.</p><p><a href="https://doi.org/10.1038/s41587-023-01743-6">Defining published article</a> · <a href="01_nano-trnaseq.html">Research note and source scope</a></p>',
      caveat('<strong>Selected published preparation.</strong> This page follows the commercial total yeast tRNA branch in Methods202/206 and the original TableS25 oligos. It does not assemble an occupied biological tRNA, complete ligation junctions or the historical RNA002 RMX/motor complex. Known local duplexes are reagent views, not proof of a unique product or perfect ligation.'),
      '<h2>Preparation</h2><ol><li>Deacylate commercial total yeast tRNA and recover RNA with the source ≥17-nt cleanup; allocate a separate TapeStation assay.</li><li>Anneal the RNA5 and mixed RNA3:DNA splint adapters, allocate prepared stock and ligate with T4 RNA ligase2. Clean up, quantify an analytical aliquot and select 200 ng padded RNA.</li><li>Ligate the preannealed RTA duplex in 15 µl. Add the entire 25 µl Maxima RT master directly to that reaction, then clean up.</li><li>Add historical RNA002 RMX with a fresh 3 µl T4 DNA ligase contribution, recover in 20 µl ELB and add water/RRB. Keep stock and loading allocations distinct.</li></ol>',
      '<h2>Known molecular parts</h2>',table(('Part','Source-supported role and limit'),(('RNA5','24 RNA residues, including the terminal source N discriminator-binding class. It is not a UMI.'),('RNA3:DNA','30 RNA residues followed by three DNA A residues. Source written rA followed by AAAA includes one RNA A plus three DNA A, not four DNA residues.'),('RTA A/B','Known whole DNA oligos and paired stem with divergent arms. Poly(A)/T-tail binding register and occupied tRNA junction remain unassembled.'),('RMX','Actual historical motor-adapter inventory; whole strands, motor and chemical ends remain unresolved. No RNA004 component is substituted.'))),
      panel(padding_scene().rows(),caption='Published padding reagent duplex; source sugar blocks retained, complete tRNA acceptor not drawn'),
      panel(rta_scene().rows(),caption='Published RTA reagent duplex with divergent Y arms; source-known local pairing only'),
      '<h2>Read setup</h2><p>Native RNA on MinION FLO-MIN106 using historical SQK-RNA002. There is no source PCR step in this branch. Protective cDNA is a companion, not the RNA read substrate. Exact motor/landing geometry and modified-base calling behavior remain unassembled; no fixed Illumina cycle layout applies.</p>',
      '<h2>Source limits and differing descriptions</h2><p>Methods250 describes custom read classification as 1 second adapter and 2 seconds strand. Extended Data Fig4 caption300 instead describes 2 seconds adapter and 1 second strand. Both published descriptions are retained; per-run configuration remains unresolved. These are computational read-classification thresholds, not physical size-selection operations.</p><p>The TableS25 terminal rN is an RNA discriminator-binding base class, not a UMI. Whole raw orderings are preserved in the research note. Synthetic tRNA manufacture, yeast culture, alternative adapter strategies and LC–MS assays are separate source branches. No absence claim or whole-kit equivalence is made.</p>',
      '<details><summary>Full workflow, stock fractions and unresolved chemistry</summary>',WorkflowSection('Published selected preparation'),'</details></main>'])


if __name__ == '__main__':
    selftest()
