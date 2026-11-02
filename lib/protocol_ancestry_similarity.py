"""Explainable earlier-method suggestions; scores never establish ancestry.

Use curated chemistry facets and exact non-inferred sequence components in separate
IDF-weighted Jaccard channels. Common platform handles contribute little or nothing.
Missing coverage is unknown. Publication years only constrain candidate discovery;
original date/version identity must be checked before a guessed edge is curated.
"""
from collections import Counter
from math import log
import re
from protocol_lineage import sequence_signatures

CHEMISTRY_FACETS = {'starting_material', 'fragmentation', 'adapter_installation',
                    'amplification', 'topology', 'strand_handling', 'selection',
                    'conversion', 'identifiers', 'partitioning'}
UNKNOWN = {'not yet curated', 'unassigned in model', 'not reported', 'not applicable',
           'none', 'no index', 'no barcode', 'no UMI'}


def publication_years(rows, metadata):
    """Use only exact defining DOI metadata, not another cited/grouped paper."""
    known = {row['dir']: row for row in rows}
    dates = {}
    for row in metadata:
        pid = row.get('our_dir')
        if pid not in known:
            continue
        doi = known[pid].get('doi', '').strip().lower()
        if doi and doi == row.get('doi', '').strip().lower() and row.get('year', '').isdigit():
            dates.setdefault(pid, set()).add(int(row['year']))
    return {pid: next(iter(values)) for pid, values in dates.items() if len(values) == 1}


def chemistry_features(properties):
    return {pid: {key + ':' + value for key, values in facets.items()
                  if key in CHEMISTRY_FACETS for value in values if value not in UNKNOWN}
            for pid, facets in properties.items()}


def _weights(features):
    active = [set(value) for value in features.values() if value]
    counts = Counter(token for values in active for token in values)
    return {token: log((len(active) + 1) / (count + 1)) for token, count in counts.items()}


def _channel(first, second, weights):
    if not first or not second:
        return None, []
    denominator = sum(weights[token] for token in first | second)
    shared = sorted((token for token in first & second if weights[token] > 0), key=lambda token: (-weights[token], token))
    return (sum(weights[token] for token in shared) / denominator if denominator else 0), shared


def cited_protocols(protocol, rows, original_text):
    """Locate known defining DOIs in a selected original; citations are not descent.

    This returns line pointers, never source facts or reaction inheritance. It
    does not treat missing DOI strings as evidence that a method was not cited.
    """
    result = []
    for row in rows:
        doi = row.get('doi', '').strip()
        if row['dir'] == protocol or not doi:
            continue
        pattern = re.compile(r'(?<![A-Za-z0-9])' + re.escape(doi) + r'(?![A-Za-z0-9_/+\-]|\.[A-Za-z0-9])', re.I)
        lines = sorted({original_text.count('\n', 0, match.start()) + 1 for match in pattern.finditer(original_text)})
        if lines:
            result.append({'protocol': row['dir'], 'title': row['protocol'], 'doi': doi,
                           'lines': lines, 'scope': 'Defining DOI occurs in selected source; inspect context. Citation, inspiration, branch identity and ancestry remain distinct.'})
    return sorted(result, key=lambda item: (item['lines'][0], item['protocol']))


def rank_candidates(protocol, rows, properties, years, components=None, *, limit=10):
    """Suggest earlier papers using chemistry and sequence evidence, never names.

    Equal-year and unknown-date candidates are excluded rather than guessing their
    order. An actual same-year predecessor can still be curated after date review.
    """
    known = {row['dir']: row for row in rows}
    if protocol not in known:
        raise ValueError('Unknown protocol: ' + protocol)
    year = years.get(protocol)
    features = chemistry_features(properties)
    sequence_eligible = {pid for pid, row in known.items() if row.get('representation') not in {'illustrative_example', 'own_design'}}
    sequences, _ = sequence_signatures(components or {'components': []}, sequence_eligible)
    facet_weights, sequence_weights = _weights(features), _weights(sequences)
    result = []
    if isinstance(year, int):
        for pid, row in known.items():
            earlier = years.get(pid)
            if pid == protocol or not isinstance(earlier, int) or earlier >= year:
                continue
            if row.get('section') != 'published' or row.get('representation') == 'own_design':
                continue
            facet_score, shared_features = _channel(features.get(protocol, set()), features.get(pid, set()), facet_weights)
            sequence_score, shared_sequences = _channel(sequences.get(protocol, set()), sequences.get(pid, set()), sequence_weights)
            if len(shared_features) < 2 and len(shared_sequences) < 2:
                continue
            channels = [(score, weight) for score, weight in [(facet_score, .75), (sequence_score, .25)] if score is not None]
            if not channels:
                continue
            # Rank by the lower bound, not a renormalization that rewards a
            # missing sequence channel over a measured partial match. Unknown
            # channels remain an interval rather than an asserted mismatch.
            score = sum(value * weight for value, weight in channels)
            upper = score + 1 - sum(weight for _, weight in channels)
            if score <= 0:
                continue
            result.append({'protocol': pid, 'title': row['protocol'], 'doi': row.get('doi', ''),
                           'year': earlier, 'method_concept_only': row.get('representation') == 'illustrative_example', 'score': round(score, 6), 'score_upper_bound': round(upper, 6),
                           'chemistry_score': round(facet_score, 6) if facet_score is not None else None,
                           'sequence_score': round(sequence_score, 6) if sequence_score is not None else None,
                           'shared_chemistry': shared_features, 'shared_sequences': shared_sequences[:5],
                           'coverage': {'chemistry_features': len(features.get(pid, [])), 'sequence_components': len(sequences.get(pid, []))}})
    result.sort(key=lambda item: (-item['score'], -item['year'], item['protocol']))
    return {'protocol': protocol, 'year': year,
            'status': 'ranked_candidates' if isinstance(year, int) else 'unknown_target_publication_year',
            'coverage': {'chemistry_features': len(features.get(protocol, [])), 'sequence_components': len(sequences.get(protocol, []))},
            'scope': 'Discovery suggestions only. Scores do not establish citation, inspiration, descent, source identity or chemical compatibility. Earlier years are cached bibliographic metadata, not verified protocol-version chronology. Unknown channels remain unknown. Same-year candidates require independent exact-date review. A published method represented by an illustrative page is a method-concept candidate only; invented molecular sequences never contribute.',
            'metric': 'IDF-weighted Jaccard in chemistry facets (75%) and exact known sequence components (25%); ranked by the lower bound with unknown-channel upper bounds; at least two informative shared chemistry features or two informative sequence components.',
            'candidates': result[:limit]}
