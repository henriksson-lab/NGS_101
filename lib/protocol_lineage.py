"""Curated protocol ancestry and independent sequence-component similarity.

Neither relation is an input to a protocol's molecular constructor or chemistry checks.
"""
from __future__ import annotations

from collections import defaultdict
from hashlib import sha256
import html
import json
import math
from pathlib import Path
from urllib.parse import quote, urlparse

BASES = str.maketrans('ACGT', 'TGCA')
BASES_SET = set('ACGT')
BASES_MIN = 16
SIMILARITY_MIN = .15
SIMILARITY_NEIGHBOURS = 3
RELATIONS = {'derived_from', 'adapted_component'}
EVIDENCE_BASES = {'documented', 'inferred'}


def validate_registry(data: dict, known: set[str]) -> list[dict]:
    if not isinstance(data, dict) or data.get('schema_version') != 1 or not isinstance(data.get('relations'), list):
        raise ValueError('Unsupported lineage registry')
    edges, seen, children = [], set(), defaultdict(list)
    for entry in data['relations']:
        if not isinstance(entry, dict):
            raise ValueError('A lineage relation must be an object')
        parent, child = entry.get('parent'), entry.get('child')
        if parent not in known or child not in known:
            raise ValueError(f'Unknown protocol in lineage: {parent} → {child}')
        if parent == child:
            raise ValueError('A protocol cannot be its own ancestor')
        if (parent, child) in seen:
            raise ValueError(f'Duplicate lineage pair: {parent} → {child}')
        if entry.get('basis') not in EVIDENCE_BASES or entry.get('relation') not in RELATIONS:
            raise ValueError('Invalid lineage relation or evidence basis')
        for key in ('scope', 'rationale', 'limits', 'reviewed'):
            if not isinstance(entry.get(key), str) or not entry[key].strip():
                raise ValueError(f'Lineage needs {key}')
        evidence = entry.get('evidence')
        if not isinstance(evidence, list) or not evidence:
            raise ValueError('Lineage needs identified evidence, including for inference')
        for source in evidence:
            if not isinstance(source, dict) or not isinstance(source.get('url'), str) or not isinstance(source.get('locator'), str):
                raise ValueError('Lineage evidence must identify a source and locator')
            url = urlparse(source.get('url', ''))
            if url.scheme != 'https' or not url.netloc or not source.get('locator', '').strip():
                raise ValueError('Lineage evidence needs an HTTPS source and inspected locator')
        seen.add((parent, child))
        children[parent].append(child)
        edge = dict(entry)
        edge['id'] = 'r_' + sha256((parent + '\0' + child).encode()).hexdigest()[:16]
        edges.append(edge)
    # Kahn traversal avoids recursion limits for a future 1,000-protocol chain.
    indegree = {node: 0 for edge in edges for node in (edge['parent'], edge['child'])}
    for edge in edges:
        indegree[edge['child']] += 1
    queue = [node for node, degree in indegree.items() if degree == 0]
    for node in queue:
        for child in children[node]:
            indegree[child] -= 1
            if indegree[child] == 0:
                queue.append(child)
    if len(queue) != len(indegree):
        raise ValueError('Curated protocol ancestry contains a cycle')
    return edges


def sequence_signatures(catalogue: dict, eligible: set[str]) -> tuple[dict[str, set[str]], dict[str, dict]]:
    """Deduplicate exact, non-inferred components on either strand.

    Components include primers, adapters and fixed known sequence parts. We do not
    call a shared platform handle a shared protocol or a chemistry equivalence.
    """
    signatures = defaultdict(set)
    metadata = {}
    for component in catalogue.get('components', []):
        sequence = component.get('sequence', '').upper()
        if (len(sequence) < BASES_MIN or set(sequence) - BASES_SET
                or len(set(sequence)) < 3):
            continue
        key = min(sequence, sequence.translate(BASES)[::-1])
        occurrences = [o for o in component.get('occurrences', [])
                       if o.get('protocol') in eligible and not o.get('inferred', True)]
        if not occurrences:
            continue
        entry = metadata.setdefault(key, {'sequence': key, 'names': set(), 'protocols': set()})
        for occurrence in occurrences:
            pid = occurrence['protocol']
            signatures[pid].add(key)
            entry['protocols'].add(pid)
            entry['names'].add(occurrence.get('name') or 'Known sequence component')
    return dict(signatures), metadata


def sequence_similarity(catalogue: dict, eligible: set[str], *, neighbours=SIMILARITY_NEIGHBOURS,
                        minimum=SIMILARITY_MIN) -> tuple[list[dict], dict]:
    signatures, components = sequence_signatures(catalogue, eligible)
    count = len(signatures)
    if not count:
        return [], {'protocols_with_components': 0, 'components': 0}
    weights = {seq: math.log((count + 1) / (len(info['protocols']) + 1))
               for seq, info in components.items()}
    totals = {pid: sum(weights[s] for s in seqs) for pid, seqs in signatures.items()}
    intersections = defaultdict(float)
    shared = defaultdict(list)
    for sequence, info in components.items():
        if weights[sequence] <= 0:
            continue
        members = sorted(info['protocols'])
        for i, first in enumerate(members):
            for second in members[i+1:]:
                pair = (first, second)
                intersections[pair] += weights[sequence]
                shared[pair].append(sequence)
    candidates, per_protocol = [], defaultdict(list)
    for pair, intersection in intersections.items():
        if len(shared[pair]) < 2:
            continue
        union = totals[pair[0]] + totals[pair[1]] - intersection
        score = intersection / union if union else 0
        if score < minimum:
            continue
        edge = {'first': pair[0], 'second': pair[1], 'score': round(score, 6),
                'shared_count': len(shared[pair]),
                'shared': [{'sequence': seq, 'names': sorted(components[seq]['names'])[:3],
                            'protocol_count': len(components[seq]['protocols'])}
                           for seq in sorted(shared[pair], key=lambda s: (-weights[s], s))[:5]]}
        candidates.append(edge)
        for pid in pair:
            per_protocol[pid].append(edge)
    selected = set()
    for pid, edges in per_protocol.items():
        for edge in sorted(edges, key=lambda e: (-e['score'], e['first'], e['second']))[:neighbours]:
            selected.add((edge['first'], edge['second']))
    return [e for e in sorted(candidates, key=lambda e: (-e['score'], e['first'], e['second']))
            if (e['first'], e['second']) in selected], {
                'protocols_with_components': count, 'components': len(components),
                'metric': 'IDF-weighted Jaccard of exact known sequence components (strand-canonical)',
                'weight': 'ln((N+1)/(df+1)); ubiquitous components have zero weight',
                'minimum_length': BASES_MIN, 'minimum_shared_components': 2,
                'minimum_score': minimum, 'nearest_neighbours': neighbours,
                'scope': 'Similarity of declared primers/adapters/fixed parts only; not ancestry, experimental similarity or complete chemistry equivalence. Missing sequence coverage is unknown.'}


def build_payload(rows: list[dict], eligible: list[dict], registry: dict, components: dict) -> dict:
    known = {r['dir']: r for r in rows}
    relations = validate_registry(registry, set(known))
    supplied = {r['id']: r for r in eligible}
    nodes = []
    for pid, item in sorted(supplied.items()):
        row = known.get(pid)
        if not row:
            raise ValueError('Unknown eligible lineage protocol: ' + pid)
        if (row.get('section') != 'published' or row.get('status') not in {'documented', 'construct', 'assay'}
                or row.get('representation') == 'illustrative_example'):
            continue
        if not item.get('url'):
            continue
        nodes.append({'id': pid, 'title': row['protocol'], 'url': item['url'],
                      'modality': row['modality'], 'status': row['status']})
    public = {n['id'] for n in nodes}
    visible = [r for r in relations if r['parent'] in public and r['child'] in public]
    similarities, metric = sequence_similarity(components, public)
    return {'schema_version': 1, 'nodes': nodes, 'relations': visible,
            'similarities': similarities, 'similarity_method': metric,
            'scope': registry.get('scope', ''),
            'stats': {'protocols': len(nodes), 'curated_relations': len(visible),
                      'documented': sum(r['basis'] == 'documented' for r in visible),
                      'inferred': sum(r['basis'] == 'inferred' for r in visible),
                      'curated_protocols': len({r[k] for r in visible for k in ('parent', 'child')}),
                      'similarities': len(similarities)}}


def local_view(protocol: str, payload: dict, *, prefix='../') -> str:
    nodes = {n['id']: n for n in payload.get('nodes', [])}
    if protocol not in nodes:
        return ''
    edges = [r for r in payload.get('relations', []) if protocol in (r['parent'], r['child'])]
    parents = sorted({r['parent'] for r in edges if r['child'] == protocol}, key=lambda p: nodes[p]['title'])
    children = sorted({r['child'] for r in edges if r['parent'] == protocol}, key=lambda p: nodes[p]['title'])
    height = max(104, 54 * max(len(parents), len(children), 1) + 24)
    escape = lambda value: html.escape(str(value), quote=True)
    link = prefix + 'lineage.html?protocol=' + quote(protocol)
    parts = ['<details class="protocol-lineage-local"><summary>Protocol relationships</summary>',
             '<style>.protocol-lineage-local{margin:16px 0}.protocol-lineage-local summary{cursor:pointer;font-weight:600}.protocol-lineage-local .lineage-local-scroll{overflow-x:auto}.protocol-lineage-local svg{max-width:none;display:block}.protocol-lineage-local p{font-size:13px}.protocol-lineage-local svg text{font:12px sans-serif}</style>']
    if edges:
        parts.append(f'<div class="lineage-local-scroll"><svg xmlns="http://www.w3.org/2000/svg" width="760" height="{height}" viewBox="0 0 760 {height}" role="img" aria-label="Local protocol ancestry">')
        positions = {protocol: (290, height / 2 - 20)}
        for members, x in ((parents, 10), (children, 570)):
            for i, pid in enumerate(members):
                positions[pid] = (x, (height - len(members) * 54) / 2 + i * 54)
        for edge in edges:
            a, b = positions[edge['parent']], positions[edge['child']]
            color = '#287d89' if edge['basis'] == 'documented' else '#a67415'
            dash = '' if edge['basis'] == 'documented' else ' stroke-dasharray="5 4"'
            x1, y1, x2, y2 = a[0] + 180, a[1] + 20, b[0], b[1] + 20
            parts.append(f'<path d="M{x1},{y1} C{x1+35},{y1} {x2-35},{y2} {x2},{y2}" fill="none" stroke="{color}" stroke-width="1.7"{dash}><title>{escape(edge["basis"]+": "+edge["scope"])}</title></path><path d="M{x2-7},{y2-4} L{x2},{y2} L{x2-7},{y2+4}" fill="none" stroke="{color}"/>')
        for pid, (x, y) in positions.items():
            label = nodes[pid]['title']
            short = label if len(label) <= 26 else label[:24] + '…'
            url = link if pid == protocol else prefix + nodes[pid]['url']
            parts.append(f'<a href="{escape(url)}"><rect x="{x}" y="{y}" width="180" height="40" rx="6" fill="{"#e8f3f5" if pid == protocol else "#f5f7fa"}" stroke="#ccd7df"/><text x="{x+90}" y="{y+24}" text-anchor="middle" fill="#21364b">{escape(short)}</text><title>{escape(label)}</title></a>')
        parts.append('</svg></div><p>Solid: source-documented · dashed: inferred. Links describe development or component scope, not identical chemistry.</p>')
    else:
        parts.append('<p>No curated ancestry links yet. This does not establish independent origins.</p>')
    parts.append(f'<p><a href="{escape(link)}">Explore relationships and sequence similarity</a></p></details>')
    return ''.join(parts)


def enrich_page(page: str, protocol: str, payload: dict) -> str:
    marker = '<details class="protocol-lineage-local">'
    if marker in page:
        raise ValueError('Protocol lineage enrichment already present')
    view = local_view(protocol, payload)
    if not view:
        return page
    position = page.find('</h1>')
    if position < 0:
        raise ValueError('Protocol page has no h1 for lineage enrichment')
    return page[:position+5] + view + page[position+5:]
