"""Read-only export of declared oligos for the browser's exact sequence search.

Coordinates are zero-based, half-open and always follow each polymer 5′ to 3′.
Names and drawing colours never establish molecular roles. This exports declarations,
not independent source verification or a prediction of chemical interactions.
"""
from __future__ import annotations

from dataclasses import asdict, fields, is_dataclass
from functools import lru_cache
import hashlib
import json
import math
import re

from chemdraw import Construct, Scene, Segment, revcomp
from protocol_graph import MaterialKind

MIN_LENGTH = 12


def _identifier(prefix, text):
    return prefix + hashlib.sha256(text.encode()).hexdigest()[:20]


def informative(sequence, minimum=MIN_LENGTH):
    """Exclude short, ambiguous, homopolymeric and simple-repeat search keys."""
    if len(sequence) < minimum or set(sequence) - set('ACGT'):
        return False
    entropy = -sum((sequence.count(b) / len(sequence)) * math.log2(sequence.count(b) / len(sequence))
                   for b in set(sequence))
    if entropy < 1.2:
        return False
    return not any(all(base == sequence[i % width] for i, base in enumerate(sequence))
                   for width in range(1, 5))


@lru_cache(maxsize=1)
def canonical_parts():
    """Curated modules' public fixed-DNA constants; never arbitrary k-mers.

    Components are associated with a protocol only through an exact contained
    sequence in one accepted concrete run. Module identity supplies the reference
    name, while the parent declaration continues to supply protocol provenance.
    """
    import illumina
    import nextera
    return tuple((module.__name__ + '.' + name, sequence)
                 for module in (illumina, nextera)
                 for name, sequence in sorted(vars(module).items())
                 if name.isupper() and not name.startswith('_')
                 and isinstance(sequence, str) and informative(sequence))


def polymers(state, path='state', seen=None):
    """Walk model containers, preserving independent polymers; never use rendered rows."""
    seen = set() if seen is None else seen
    if state is None or isinstance(state, (str, int, float, bool)) or id(state) in seen:
        return
    seen.add(id(state))
    if isinstance(state, Segment):
        yield path, [state]
    elif isinstance(state, Construct):
        yield path, list(state)
    elif isinstance(state, Scene):
        for name, strand in state.strands.items():
            yield path + '/strand:' + name, list(strand.segs)
        # Components are already represented by Scene.strands; no joining or re-walk.
    elif isinstance(state, dict):
        for key, value in sorted(state.items(), key=lambda x: str(x[0])):
            yield from polymers(value, path + '/' + str(key), seen)
    elif isinstance(state, (tuple, list)):
        for i, value in enumerate(state):
            yield from polymers(value, path + '/' + str(i), seen)
    elif is_dataclass(state):
        for field in fields(state):
            yield from polymers(getattr(state, field.name), path + '/' + field.name, seen)


def index_axis(sequence, start, end):
    """Canonical contiguous arm geometry is required for an exportable index axis.

    Forward means the index bases in a P5-oriented complete library. A P7-side
    primer therefore carries the reverse complement of the forward i7 value.
    """
    import illumina as il
    import nextera as nx
    matches = []
    for reverse in (False, True):
        seq = revcomp(sequence) if reverse else sequence
        a, b = (len(seq) - end, len(seq) - start) if reverse else (start, end)
        for role, outer, inners in (('i5', il.P5, (il.TRUSEQ_READ1, nx.S5)),
                                     ('i7', il.P7, (il.TRUSEQ_READ2, nx.S7))):
            if seq[:a].endswith(outer) and any(seq[b:].startswith(inner) for inner in inners):
                rc = (role == 'i7') != reverse
                matches.append({'role': role, 'to_forward': 'reverse_complement' if rc else 'forward',
                                'basis': 'Contiguous canonical ' + ('P5/read-1' if role == 'i5' else 'P7/read-2') + ' arms flanking the structured sample index.'})
    unique = {(m['role'], m['to_forward']) for m in matches}
    return matches[0] if len(unique) == 1 else {}


def template(segments, protocol, label, locator, source='', inferred=False):
    """Export only bounded polymers with structured sample indices.

    A schematic placeholder width is not molecular length. Only explicit length_bp
    or an ambiguity-string feature with a declared position count establishes it.
    """
    sequence, parts, indices = '', [], []
    for segment in segments:
        text = segment.top.upper()
        if segment.placeholder:
            width = segment.length_bp
            if width is None and segment.feature and text and set(text) <= set('ACGTRYSWKMBDHVN'):
                width = len(text)
            if width is None:
                return None
            text = 'N' * width
        elif not text or set(text) - set('ACGTN'):
            return None
        start = len(sequence)
        sequence += text
        part = {'name': segment.name, 'start': start, 'end': len(sequence),
                'placeholder': segment.placeholder, 'inferred': segment.inferred}
        if segment.feature:
            part['feature'] = asdict(segment.feature)
            if segment.feature.role == 'sample_index':
                indices.append({'start': start, 'end': len(sequence), 'feature_id': segment.feature.id})
        parts.append(part)
    if not indices:
        return None
    for index in indices:
        index.update(index_axis(sequence, index['start'], index['end']))
    return {'id': _identifier('t_', protocol + '\0' + locator + '\0' + sequence),
            'protocol': protocol, 'label': label, 'source': source, 'locator': locator,
            'inferred': inferred or any(s.inferred for s in segments),
            'sequence': sequence, 'segments': parts, 'index_spans': indices}


def extract_graph(graph, protocol):
    """Extract explicit reagent polymers and declared PCR/sequencing primers.

    Arbitrary analyte examples and computed library products are not searched as
    whole reference components. Their source-identifiable reagent components enter
    through explicit oligo declarations instead.
    """
    occurrences, templates, unrepresented = [], [], []
    declared_sequences = {}

    def add(sequence, name, locator, role, source, inferred=False, feature=None, **metadata):
        sequence = sequence.upper()
        if not informative(sequence):
            return
        occurrence = {'protocol': protocol, 'name': name or role, 'locator': locator,
                      'role': role, 'source': source, 'inferred': bool(inferred)}
        occurrence.update(metadata)
        if feature:
            occurrence['feature'] = asdict(feature)
        occurrences.append({'sequence': sequence, 'occurrence': occurrence})
        for reference, canonical in canonical_parts():
            # Explicit references can name opposite strands of one component.
            # Keep their aliases/provenance; combine() deduplicates the bases.
            for query in dict.fromkeys((canonical, revcomp(canonical))):
                start = sequence.find(query)
                while start >= 0:
                    part = dict(occurrence, name=reference, role='canonical sequence part',
                                locator=locator + '/canonical:' + reference + '/bases:' + str(start),
                                parent_locator=locator, canonical_reference=reference,
                                match_basis='Exact contained sequence; identity alone does not establish primer or reagent usage.')
                    occurrences.append({'sequence': query, 'occurrence': part})
                    start = sequence.find(query, start + 1)

    for material in graph.materials.values():
        if material.kind not in (MaterialKind.OLIGO, MaterialKind.ADAPTER):
            continue
        declared_polymers = list(polymers(material.state))
        if not declared_polymers:
            unrepresented.append(material.id)
        for path, segments in declared_polymers:
            locator = 'material:' + material.id + '/' + path
            if all(not s.placeholder and set(s.top.upper()) <= set('ACGT') for s in segments):
                full = ''.join(s.top.upper() for s in segments)
                declared_sequences.setdefault(full, []).append(
                    material.inferred or any(s.inferred for s in segments))
            made = template(segments, protocol, material.label, locator, material.source, material.inferred)
            if made:
                templates.append(made)
            run, run_inferred, run_names = '', material.inferred, []

            def flush():
                if run:
                    add(run, material.label + (' — ' + ', '.join(run_names) if run_names else ''),
                        locator + '/fixed-run:' + str(run_start), material.kind.value,
                        material.source, run_inferred)

            run_start = 0
            for i, segment in enumerate(segments):
                text = segment.top.upper()
                if segment.placeholder:
                    flush()
                    run, run_names, run_inferred = '', [], material.inferred
                    run_start = i + 1
                    continue
                inferred = material.inferred or segment.inferred
                if run and run_inferred != inferred:
                    flush()
                    run, run_names, run_inferred = '', [], inferred
                    run_start = i
                chunks = list(re.finditer('[ACGT]+', text))
                if not chunks or chunks[0].start() > 0:
                    flush()
                    run, run_names, run_inferred = '', [], material.inferred
                    run_start = i
                for chunk in chunks:
                    add(chunk.group(), segment.name or material.label,
                        locator + '/segment:' + str(i) + '/bases:' + str(chunk.start()),
                        segment.feature.role if segment.feature else 'fixed segment', material.source,
                        material.inferred or segment.inferred, segment.feature)
                    run += chunk.group()
                    run_names.append(segment.name)
                    run_inferred = run_inferred or segment.inferred
                    if chunk.end() < len(text):
                        flush()
                        run, run_names, run_inferred = '', [], material.inferred
                        run_start = i
            flush()
    for step in graph.steps.values():
        geometry = step.geometry
        # PCRPrimer currently has no segment/feature structure; never invent one.
        for side in ('forward', 'reverse'):
            primer = getattr(geometry, side, None)
            sequence = getattr(primer, 'sequence', None)
            if isinstance(sequence, str):
                matching = declared_sequences.get(sequence.upper(), [])
                inferred = (getattr(primer, 'inferred', False) or getattr(geometry, 'inferred', False)
                            or bool(matching and all(matching)))
                metadata = {}
                if hasattr(geometry, 'availability'):
                    metadata['geometry_availability'] = geometry.availability
                    metadata['geometry_statement'] = geometry.statement
                for i, fixed in enumerate(re.split('[^ACGT]+', sequence.upper())):
                    add(fixed, primer.name, 'step:' + step.id + '/' + side + '/fixed-run:' + str(i),
                        'PCR primer', step.source, inferred, **metadata)
    for material_id, plan in graph.sequencing.items():
        for i, primer in enumerate(plan.primers):
            add(primer.seq, primer.name, 'sequencing:' + material_id + '/primer:' + str(i),
                'sequencing reference', primer.source, reference_only=True,
                reference_role=primer.role, note=primer.note)
    return {'occurrences': occurrences, 'templates': templates,
            'unrepresented_oligos': sorted(unrepresented)}


def combine(protocols, extractions, report=()):
    groups, templates = {}, []
    for extracted in extractions:
        templates.extend(extracted['templates'])
        for item in extracted['occurrences']:
            seq = item['sequence']
            canonical = min(seq, revcomp(seq))
            entry = groups.setdefault(canonical, {'id': _identifier('s_', canonical), 'sequence': canonical,
                                                  'aliases': set(), 'occurrences': {}, 'protocols': set()})
            occurrence = dict(item['occurrence'], orientation='+' if seq == canonical else '-')
            entry['aliases'].add(occurrence['name'])
            entry['protocols'].add(occurrence['protocol'])
            entry['occurrences'][json.dumps(occurrence, sort_keys=True)] = occurrence
    components = []
    for seq, entry in sorted(groups.items()):
        components.append(dict(entry, aliases=sorted(entry['aliases']), protocols=sorted(entry['protocols']),
                               occurrences=[entry['occurrences'][key] for key in sorted(entry['occurrences'])]))
    templates = sorted({t['id']: t for t in templates}.values(), key=lambda t: (t['protocol'], t['id']))
    return {'schema_version': 1, 'min_length': MIN_LENGTH,
            'protocols': sorted(protocols, key=lambda p: p['id']), 'components': components,
            'templates': templates, 'report': sorted(report, key=lambda p: p['protocol']),
            'stats': {'protocols': len(protocols), 'components': len(components), 'templates': len(templates),
                      'occurrences': sum(len(c['occurrences']) for c in components)}}
