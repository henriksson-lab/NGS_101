#!/usr/bin/env python3
"""Meaningful boundary/identity tests for the offline sequence component export."""
from pathlib import Path
from types import SimpleNamespace
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'lib'))
from checks import Check
from chemdraw import Construct, Scene, Segment, feature, revcomp
from protocol_graph import ProtocolGraph, MaterialKind
import illumina as il
import seqprimers as sp
from pcr import PCRGeometry, PCRPrimer
from build_sequence_index import public_url
from sequence_catalogue import combine, extract_graph, informative, index_axis, template


def main():
    check = Check()
    seq = il.TRUSEQ_READ1
    check('short sequences excluded', informative(seq[:11]), False)
    check('homopolymers excluded', informative('T' * 30), False)
    check('simple repeats excluded', informative('ACGT' * 8), False)
    check('canonical primer retained', informative(seq), True)
    check('ambiguous key excluded', informative(seq + 'N'), False)
    graph = ProtocolGraph('fixture', authority='fixture source')
    idx = feature('fixture.index', 'sample_index', 'fixed')
    umi = feature('fixture.umi', 'umi', 'random')
    indexed = [Segment('P5', il.P5), Segment('index', 'AACCGTTA', feature=idx),
               Segment('read1', il.TRUSEQ_READ1)]
    graph.input('primer', 'named source primer', MaterialKind.OLIGO, state=Construct(indexed))
    graph.input('rc', 'reverse alias', MaterialKind.OLIGO,
                state=Construct([Segment('reverse', revcomp(il.P5))]))
    graph.input('analyte', 'example analyte', MaterialKind.BIOLOGICAL,
                state=Construct([Segment('arbitrary example', 'AGTTCCGATAGCTAGCTAGTTCCA')]))
    graph.input('ambiguous', 'ambiguous polymer', MaterialKind.OLIGO,
                state=Construct([Segment('bounded ambiguity', il.P5 + 'NNN' + seq)]))
    graph.input('unknown', 'unassigned polymer', MaterialKind.OLIGO,
                state=Construct([Segment('left', il.P5), Segment('unknown', '[insert]', placeholder=True),
                                 Segment('right', seq)]))
    scene = Scene()
    scene.strand('first', [Segment('one', il.P5)])
    scene.strand('second', [Segment('two', seq)])
    graph.input('independent', 'two polymers', MaterialKind.ADAPTER, state=scene)
    extracted = extract_graph(graph, 'fixture')
    catalogue = combine([{'id': 'fixture', 'title': 'Fixture', 'url': 'fixture.html'}], [extracted])
    p5 = next(c for c in catalogue['components'] if c['sequence'] == min(il.P5, revcomp(il.P5)))
    check('reverse complement deduplication preserves orientation',
          set(o['orientation'] for o in p5['occurrences']), {'+', '-'})
    check('reverse alias retained', 'reverse' in p5['aliases'], True)
    keys = {item['sequence'] for item in extracted['occurrences']}
    check('known runs around within-segment ambiguity retained',
          {item['sequence'] for item in extracted['occurrences'] if 'material:ambiguous/' in item['occurrence']['locator']
           and not item['occurrence'].get('canonical_reference')},
          {il.P5, seq})
    check('unknown boundary never joined', il.P5 + seq in keys, False)
    check('independent polymers never joined',
          any('material:independent' in item['occurrence']['locator'] and item['sequence'] == il.P5 + seq
              for item in extracted['occurrences']), False)
    check('arbitrary analyte excluded',
          any('material:analyte' in item['occurrence']['locator'] for item in extracted['occurrences']), False)
    made = extracted['templates'][0]
    check('i5 role established by canonical geometry', made['index_spans'][0]['role'], 'i5')
    check('i5 forward axis', made['index_spans'][0]['to_forward'], 'forward')
    check('i7 oligo has reverse forward-library axis',
          index_axis(il.P7 + 'AACCGTTA' + il.TRUSEQ_READ2, len(il.P7), len(il.P7) + 8)['to_forward'],
          'reverse_complement')
    reversed_template = revcomp(il.P7 + 'AACCGTTA' + il.TRUSEQ_READ2)
    check('reverse i7 polymer transforms coordinates and axis',
          index_axis(reversed_template, len(il.TRUSEQ_READ2), len(il.TRUSEQ_READ2) + 8)['to_forward'], 'forward')
    symbolic = template([Segment('P5', il.P5),
                         Segment('index', '[index]', placeholder=True, length_bp=8, feature=idx),
                         Segment('read1', seq)], 'fixture', 'symbolic', 'input')
    check('bounded placeholder uses Ns', symbolic['sequence'][len(il.P5):len(il.P5)+8], 'N'*8)
    check('unknown placeholder width rejects template',
          template([Segment('index', '[index]', placeholder=True, feature=idx)],
                   'fixture', 'unknown', 'input'), None)
    check('UMI alone is never a sample index template',
          template([Segment('umi', 'NNNNNNNN', placeholder=True, feature=umi)],
                   'fixture', 'UMI', 'input'), None)
    check('unresolved index axis stays unresolved',
          template([Segment('index', 'AACCGTTA', feature=idx)], 'fixture', 'unresolved', 'input')['index_spans'],
          [{'start': 0, 'end': 8, 'feature_id': idx.id}])
    # Derive actual bases reported by the independently implemented sequencing
    # primer locator, then compare the export's index orientation with those reads.
    i5, raw_i7 = 'ACGTTCGA', 'GATCTGCA'
    fwd = il.P5 + i5 + il.TRUSEQ_READ1
    rev = il.P7 + raw_i7 + il.TRUSEQ_READ2
    library = Construct([Segment('fwd', fwd), Segment('insert', 'GGCCTACAGATTCGACGTCA'),
                         Segment('rev', revcomp(rev))])
    i1_read = sp.locate(library, sp.TRUSEQ['I1'], n=8).reads
    i2_rc_read = sp.locate(library, sp.TRUSEQ['I2'], n=8).reads
    def converted(raw, axis):
        return revcomp(raw) if axis['to_forward'] == 'reverse_complement' else raw
    check('P7 index export agrees with bases actually reported by Index 1 primer',
          converted(raw_i7, index_axis(rev, len(il.P7), len(il.P7)+8)), i1_read)
    check('P5 forward index export is reverse complement of RC-workflow Index 2 read',
          converted(i5, index_axis(fwd, len(il.P5), len(il.P5)+8)), revcomp(i2_rc_read))
    graph.input('mixed_evidence', 'mixed source and inferred', MaterialKind.OLIGO,
                state=Construct([Segment('known left', il.P5),
                                 Segment('inferred middle', il.P7, inferred=True),
                                 Segment('known right', seq)]))
    graph.input('inferred_pcr', 'inferred PCR primer', MaterialKind.OLIGO,
                state=Construct([Segment('inferred primer', fwd, inferred=True)]))
    geometry = PCRGeometry.exact(library, PCRPrimer('forward', fwd), PCRPrimer('reverse', rev),
                                 statement='fixture geometry')
    graph.steps['fixture_pcr'] = SimpleNamespace(id='fixture_pcr', geometry=geometry, source='fixture PCR')
    graph.sequencing['fixture_library'] = SimpleNamespace(primers=(sp.TRUSEQ['I1'],))
    evidence = extract_graph(graph, 'fixture')
    mixed = [item for item in evidence['occurrences']
             if 'material:mixed_evidence/' in item['occurrence']['locator']]
    check('inferred boundaries prevent joined source-known runs',
          {item['sequence'] for item in mixed if not item['occurrence'].get('canonical_reference')},
          {il.P5, il.P7, seq})
    check('known component retains source-known status next to inferred component',
          all(not item['occurrence']['inferred'] for item in mixed if item['sequence'] == seq), True)
    check('PCR primer inherits matching explicit inferred reagent declaration',
          next(item['occurrence']['inferred'] for item in evidence['occurrences']
               if item['occurrence']['locator'] == 'step:fixture_pcr/forward/fixed-run:0'), True)
    check('PCR occurrence preserves geometry scope',
          next(item['occurrence']['geometry_availability'] for item in evidence['occurrences']
               if item['occurrence']['locator'] == 'step:fixture_pcr/reverse/fixed-run:0'), 'exact')
    reference = next(item['occurrence'] for item in evidence['occurrences']
                     if item['occurrence']['locator'].startswith('sequencing:'))
    check('sequencing primer is reference rather than loaded-primer assertion',
          (reference['role'], reference['reference_only'], reference['reference_role']),
          ('sequencing reference', True, 'Index 1 (i7)'))
    check('standalone URL discovery resolves wrapper output',
          public_url('smart-seq__10.1038+nbt.2282'), 'smart-seq__10.1038+nbt.2282/smart-seq.html')
    check('standalone URL discovery resolves wrapper-free output',
          public_url('smart-seq2__10.1038+nmeth.2639'), 'smart-seq2__10.1038+nmeth.2639/smart-seq2.html')
    undivided = ProtocolGraph('canonical containment', authority='fixture source')
    undivided.input('whole', 'undivided source primer', MaterialKind.OLIGO,
                    state=Construct([Segment('whole primer', fwd, inferred=True)]))
    cut = len(il.P5) // 2
    undivided.input('broken', 'broken P5', MaterialKind.OLIGO,
                    state=Construct([Segment('left half', il.P5[:cut]),
                                     Segment('unknown', 'NNNN', placeholder=True),
                                     Segment('right half', il.P5[cut:])]))
    contained = extract_graph(undivided, 'fixture')['occurrences']
    contained_p5 = [item for item in contained
                    if item['occurrence'].get('canonical_reference') == 'illumina.P5']
    check('canonical P5 extracted from undivided primer', bool(contained_p5), True)
    check('canonical component keeps parent inference/source/locator',
          all(item['occurrence']['inferred'] and item['occurrence']['source'] == 'fixture source'
              and item['occurrence']['parent_locator'].startswith('material:whole/')
              for item in contained_p5), True)
    check('canonical identity has no implied experimental primer role',
          {item['occurrence']['role'] for item in contained_p5}, {'canonical sequence part'})
    check('canonical parts never bridge unknown sequence',
          any(item['occurrence'].get('canonical_reference') == 'illumina.P5'
              and 'material:broken/' in item['occurrence']['locator'] for item in contained), False)
    check('canonical parts preserve sequencing reference qualification',
          all(item['occurrence'].get('reference_only') is True for item in evidence['occurrences']
              if item['occurrence'].get('canonical_reference')
              and item['occurrence']['locator'].startswith('sequencing:')), True)
    check('deterministic JSON-ready output', json.dumps(catalogue, sort_keys=True),
          json.dumps(combine([{'id': 'fixture', 'title': 'Fixture', 'url': 'fixture.html'}], [extracted]), sort_keys=True))
    check.report()


if __name__ == '__main__':
    main()
