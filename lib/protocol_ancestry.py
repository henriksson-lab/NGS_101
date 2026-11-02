"""Review coverage for protocol ancestry, independent of molecular validation."""
from collections import Counter
from urllib.parse import urlparse
from protocol_lineage import validate_registry, validate_guess_chronology, RELATIONS, EVIDENCE_BASES

REVIEW_STATES = {'unreviewed', 'existing_curated', 'reviewed', 'source_blocked', 'not_applicable'}
FORMAL_STATES = {'unassessed', 'internal', 'external', 'internal_and_external', 'cannot_formally_assign', 'not_applicable'}
INFERENCE_STATES = {'not_assessed', 'supported', 'no_defensible_candidate', 'not_needed', 'blocked', 'guessed_candidate'}


def _text(value, label):
    if not isinstance(value, str) or not value.strip():
        raise ValueError('Ancestry review needs ' + label)


def _url(value, *, require_https=False):
    _text(value, 'an original URL')
    parsed = urlparse(value)
    if parsed.scheme not in ({'https'} if require_https else {'http', 'https'}) or not parsed.netloc or parsed.username or parsed.password:
        raise ValueError('Ancestry evidence needs an HTTPS original URL')


def validate_reviews(data, known, relations):
    if data.get('schema_version') != 1 or not isinstance(data.get('reviews'), dict):
        raise ValueError('Unsupported ancestry review registry')
    reviews = data['reviews']
    for pid, review in reviews.items():
        if pid not in known:
            raise ValueError('Unknown protocol in ancestry reviews: ' + pid)
        if review.get('status') not in REVIEW_STATES or review.get('formal_assignment') not in FORMAL_STATES:
            raise ValueError('Invalid ancestry review/formal assignment: ' + pid)
        if review.get('inference_status') not in INFERENCE_STATES:
            raise ValueError('Invalid likely-predecessor assessment: ' + pid)
        for field in ('reason', 'limits', 'inference_reason'):
            _text(review.get(field), field + ': ' + pid)
        for source in review.get('sources', []):
            _url(source.get('url'))
            for field in ('locator', 'identity', 'availability'):
                _text(source.get(field), 'source ' + field + ': ' + pid)
        for blocker in review.get('blockers', []):
            _url(blocker.get('url'))
            _text(blocker.get('reason'), 'actual source blocker: ' + pid)
        if review['status'] == 'source_blocked' and not review.get('blockers'):
            raise ValueError('Source-blocked ancestry needs an actual source target: ' + pid)
        if review['status'] in {'reviewed', 'existing_curated'} and not review.get('sources'):
            raise ValueError('Reviewed ancestry needs inspected or prior-curated evidence: ' + pid)
        if review['status'] != 'unreviewed':
            _text(review.get('reviewed'), 'review date: ' + pid)
        if review['status'] == 'reviewed' and review['formal_assignment'] == 'unassessed':
            raise ValueError('Reviewed ancestry needs an explicit formal-attribution outcome: ' + pid)
        if review['status'] in {'reviewed', 'source_blocked'} and review['inference_status'] == 'not_assessed':
            raise ValueError('Reviewed ancestry must assess a likely predecessor: ' + pid)
        external = review.get('external_predecessors', [])
        external_seen = set()
        for parent in external:
            identity = (parent.get('url'), parent.get('title', '').casefold())
            if identity in external_seen:
                raise ValueError('Duplicate external predecessor: ' + pid)
            external_seen.add(identity)
            _text(parent.get('title'), 'external predecessor identity: ' + pid)
            _url(parent.get('url'), require_https=True)
            if parent.get('relation') not in RELATIONS or parent.get('basis') not in EVIDENCE_BASES:
                raise ValueError('Invalid external predecessor relation: ' + pid)
            validate_guess_chronology(parent)
            for field in ('scope', 'rationale', 'limits', 'reviewed'):
                _text(parent.get(field), 'external predecessor ' + field + ': ' + pid)
            if not parent.get('evidence'):
                raise ValueError('External predecessor needs located evidence: ' + pid)
            for evidence in parent['evidence']:
                _url(evidence.get('url'), require_https=True)
                _text(evidence.get('locator'), 'external predecessor evidence locator: ' + pid)
        internal = [edge for edge in relations if edge['child'] == pid]
        documented_internal = any(edge['basis'] == 'documented' for edge in internal)
        documented_external = any(edge['basis'] == 'documented' for edge in external)
        formal = review['formal_assignment']
        if formal in {'internal', 'internal_and_external'} and not documented_internal:
            raise ValueError('Formal internal assignment has no documented registry edge: ' + pid)
        if formal in {'external', 'internal_and_external'} and not documented_external:
            raise ValueError('Formal external assignment has no documented external predecessor: ' + pid)
        if formal == 'cannot_formally_assign' and (documented_internal or documented_external):
            raise ValueError('A scoped formal predecessor is already documented: ' + pid)
        if review['inference_status'] == 'supported' and not any(edge['basis'] == 'inferred' for edge in [*internal, *external]):
            raise ValueError('Supported inference has no explicitly inferred edge: ' + pid)
        if review['inference_status'] == 'guessed_candidate' and not any(edge['basis'] == 'guessed' for edge in [*internal, *external]):
            raise ValueError('Guessed candidate needs an explicitly guessed relationship: ' + pid)
    return reviews


def ancestry_audit(rows, registry, review_registry=None, *, visible_ids=None):
    """Return per-protocol coverage; missing parents never mean independent origin.

    Internal parents are counted against the complete registry, even if they are
    hidden from the public graph. External predecessors are recorded separately.
    """
    known = {row['dir'] for row in rows}
    relations = validate_registry(registry, known)
    reviews = validate_reviews(review_registry, known, relations) if review_registry is not None else {}
    report = []
    for row in sorted(rows, key=lambda item: item['dir']):
        pid = row['dir']
        incoming = [edge for edge in relations if edge['child'] == pid]
        outgoing = [edge for edge in relations if edge['parent'] == pid]
        review = reviews.get(pid, {})
        external = review.get('external_predecessors', [])
        report.append({'id': pid, 'title': row['protocol'], 'catalogue_status': row['status'],
                       'review_status': review.get('status', 'unreviewed'),
                       'formal_assignment': review.get('formal_assignment', 'unassessed'),
                       'inference_status': review.get('inference_status', 'not_assessed'),
                       'internal_parents': [edge['parent'] for edge in incoming],
                       'documented_parents': [edge['parent'] for edge in incoming if edge['basis'] == 'documented'],
                       'inferred_parents': [edge['parent'] for edge in incoming if edge['basis'] == 'inferred'],
                       'guessed_parents': [edge['parent'] for edge in incoming if edge['basis'] == 'guessed'],
                       'external_predecessors': external,
                       'without_internal_parent': not any(edge['basis'] != 'guessed' for edge in incoming),
                       'without_recorded_candidate': not incoming and not external,
                       'without_formal_parent': not any(edge['basis'] == 'documented' for edge in [*incoming, *external]),
                       'without_visible_parent': not (any(edge['parent'] in visible_ids for edge in incoming) or bool(external) and pid in visible_ids) if visible_ids is not None else None,
                       'unlinked': not incoming and not outgoing and not external,
                       'reason': review.get('reason', 'Ancestry has not yet been reviewed.'),
                       'inference_reason': review.get('inference_reason', ''),
                       'has_source_blocker': review.get('status') == 'source_blocked' or bool(review.get('blockers')),
                       'blockers': review.get('blockers', []),
                       'sources': review.get('sources', []), 'limits': review.get('limits', ''), 'reviewed': review.get('reviewed', '')})
    return {'scope': 'Curated predecessor coverage, not proof of independent origin or molecular compatibility.',
            'protocols': report,
            'stats': {'protocols': len(report),
                      'without_internal_parent': sum(item['without_internal_parent'] for item in report),
                      'without_formal_parent': sum(item['without_formal_parent'] for item in report),
                      'unlinked': sum(item['unlinked'] for item in report),
                      'review_states': dict(Counter(item['review_status'] for item in report)),
                      'formal_assignments': dict(Counter(item['formal_assignment'] for item in report))}}


def protocols_without_predecessors(rows, registry, review_registry=None, *, include_inferred=True, include_guessed=False):
    """List protocols with no curated in-catalogue parent (externals stay labelled)."""
    items = ancestry_audit(rows, registry, review_registry)['protocols']
    return [item for item in items if not (item['documented_parents'] + (item['inferred_parents'] if include_inferred else []) + (item['guessed_parents'] if include_guessed else []))]
