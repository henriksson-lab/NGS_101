"""Exact two-strand pairing discovery shared by Python and the browser prototype.

Candidates describe possible complementary geometry, never binding probability,
preferred folds, polymerase activity or covalent closure. Inputs are actual ssDNA
strands in molecular 5′→3′ order; a target's unprovided opposite strand is not added.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import re
import xml.etree.ElementTree as ET

from chemdraw import Scene, Segment, panel_svg, revcomp
from topology import PadlockCaptureTopology

MAX_BASES = 3000
MAX_CANDIDATES = 250
MAX_WORK = 500_000


def normalize(value: str) -> str:
    """Plain DNA or a single FASTA record. N is an unknown-base barrier; RNA and other ambiguity symbols are rejected."""
    if not isinstance(value, str):
        raise ValueError('Enter a DNA sequence.')
    lines = value.strip().splitlines()
    if lines and lines[0].lstrip().startswith('>'):
        lines = lines[1:]
    if any(line.lstrip().startswith('>') for line in lines):
        raise ValueError('Use one FASTA record per input.')
    sequence = re.sub(r'\s+', '', ''.join(lines)).upper()
    if not sequence or set(sequence) - set('ACGTN'):
        raise ValueError('Use A/C/G/T (N for unknown unpaired bases), written 5′ to 3′.')
    if len(sequence) > MAX_BASES:
        raise ValueError(f'Each input is limited to {MAX_BASES:,} bases.')
    return sequence


@dataclass(frozen=True)
class LinearCandidate:
    first_start: int
    second_start: int
    length: int
    kind: str = 'linear'

    def scene(self, first: str, second: str) -> Scene:
        a, b, k = self.first_start, self.second_start, self.length
        if k < 1 or a < 0 or b < 0 or a + k > len(first) or b + k > len(second):
            raise ValueError('Pairing coordinates are outside the supplied strands.')
        if set(first[a:a + k] + second[b:b + k]) - set('ACGT') or first[a:a + k] != revcomp(second[b:b + k]):
            raise ValueError('The selected spans are not exact complements.')
        def parts(seq, start):
            return ([Segment('5′ overhang', seq[:start])] if start else []) + [
                Segment('paired span', seq[start:start + k], 'cbc')] + (
                [Segment('3′ overhang', seq[start + k:])] if start + k < len(seq) else [])
        scene = Scene()
        scene.strand('sequence 1', parts(first, a))
        scene.anneal('sequence 2', parts(second, b), to='sequence 1',
                    pair=('paired span', 'paired span'), unpaired=('5′ overhang', '3′ overhang'))
        return scene


@dataclass(frozen=True)
class PadlockCandidate:
    target: int
    extension_start: int
    extension_length: int
    ligation_start: int
    ligation_length: int
    kind: str = 'padlock'

    @property
    def length(self):
        return self.extension_length + self.ligation_length

    def topology(self, first: str, second: str) -> PadlockCaptureTopology:
        target, probe = (first, second) if self.target == 1 else (second, first)
        if self.target not in (1, 2):
            raise ValueError('Target must be sequence 1 or 2.')
        e, l = self.extension_length, self.ligation_length
        if e < 1 or l < 1 or e + l >= len(probe):
            raise ValueError('Two terminal arms must leave a non-empty padlock backbone.')
        if set(probe[:l] + probe[-e:]) - set('ACGT'):
            raise ValueError('Padlock arms must contain known A/C/G/T bases only.')
        return PadlockCaptureTopology(
            'sequence ' + str(3 - self.target) + ' padlock', target,
            probe[:l], probe[l:len(probe) - e], probe[-e:],
            ligation_site=self.ligation_start, extension_site=self.extension_start)


class _Budget:
    def __init__(self):
        self.used = 0
        self.truncated = False

    def spend(self, amount=1):
        self.used += amount
        if self.used > MAX_WORK:
            self.truncated = True
            return False
        return True


def _sites(sequence, motif):
    start = 0
    while True:
        start = sequence.find(motif, start)
        if start < 0:
            return
        yield start
        start += 1


def _linear(first, second, minimum, budget):
    reverse = revcomp(second)
    # Seeds locate maximal exact runs. Interior seeds of one run are skipped.
    for a in range(len(first) - minimum + 1):
        seed = first[a:a + minimum]
        if 'N' in seed:
            continue
        for b in _sites(reverse, seed):
            if not budget.spend():
                return
            if a and b and first[a - 1] in 'ACGT' and first[a - 1] == reverse[b - 1]:
                continue
            k = minimum
            while a + k < len(first) and b + k < len(reverse) and first[a + k] in 'ACGT' and first[a + k] == reverse[b + k]:
                if not budget.spend():
                    return
                k += 1
            yield LinearCandidate(a, len(second) - b - k, k)


def _padlocks(target, probe, target_number, minimum, budget):
    reverse = revcomp(probe)
    if len(probe) < 2 * minimum + 1 or 'N' in probe[:minimum] + probe[-minimum:]:
        return
    extensions, ligations = [], []
    for start in _sites(target, reverse[:minimum]):
        length = minimum
        while start + length < len(target) and length < len(reverse) and target[start + length] in 'ACGT' and target[start + length] == reverse[length]:
            if not budget.spend():
                return
            length += 1
        extensions.append((start, length))
    for start in _sites(target, reverse[-minimum:]):
        end, length = start + minimum, minimum
        while end - length > 0 and length < len(reverse) and target[end - length - 1] in 'ACGT' and target[end - length - 1] == reverse[-length - 1]:
            if not budget.spend():
                return
            length += 1
        ligations.append((end, length))
    for ext_start, ext_max in extensions:
        for lig_end, lig_max in ligations:
            if not budget.spend():
                return
            # Prefix/suffix arms must be disjoint in both probe and target. Keep
            # Pareto-maximal arm lengths; enumerate alternative valid boundaries.
            capacity = len(probe) - 1
            if ext_start < lig_end:
                capacity = min(capacity, lig_end - ext_start)
            if capacity < 2 * minimum:
                continue
            upper = min(ext_max, capacity - minimum)
            lower = max(minimum, min(ext_max, capacity - lig_max))
            for e in range(lower, upper + 1):
                l = min(lig_max, capacity - e)
                if l < minimum or not (ext_start + e <= lig_end - l or lig_end <= ext_start):
                    continue
                if not budget.spend():
                    return
                yield PadlockCandidate(target_number, ext_start, e, lig_end - l, l)


def search(first, second, mode='both', minimum=8, limit=MAX_CANDIDATES):
    """Enumerate bounded maximal exact candidates; truncation is always explicit.

    Padlock mode tries each supplied strand as target. Maximal arms are geometric
    candidates; shorter nested versions are not separately listed. When the two arms
    compete for sequence, alternative maximal boundaries are preserved.
    """
    first, second = normalize(first), normalize(second)
    if mode not in ('linear', 'padlock', 'both'):
        raise ValueError('Choose linear, padlock or both.')
    if isinstance(minimum, bool) or not isinstance(minimum, int) or not 2 <= minimum <= 100:
        raise ValueError('Minimum paired length must be an integer from 2 to 100.')
    if isinstance(limit, bool) or not isinstance(limit, int) or not 1 <= limit <= MAX_CANDIDATES:
        raise ValueError(f'Candidate limit must be from 1 to {MAX_CANDIDATES}.')
    candidates, budget = [], _Budget()
    families = []
    if mode in ('linear', 'both'):
        families.append(_linear(first, second, minimum, budget))
    if mode in ('padlock', 'both'):
        families.extend((_padlocks(first, second, 1, minimum, budget),
                         _padlocks(second, first, 2, minimum, budget)))
    # Round-robin preserves both modes/target assignments under a result cap.
    while families and len(candidates) <= limit:
        remaining = []
        for family in families:
            try:
                candidates.append(next(family))
                remaining.append(family)
            except StopIteration:
                pass
            if len(candidates) > limit:
                break
        families = remaining
    truncated = budget.truncated or len(candidates) > limit
    candidates = sorted(candidates[:limit], key=lambda c: (-c.length, c.kind, repr(c)))
    return {'first': first, 'second': second, 'candidates': [asdict(c) for c in candidates],
            'truncated': truncated, 'work': budget.used, 'minimum': minimum, 'mode': mode}


def _standalone(svg):
    """Keep shared SVG semantics and embed their styles for a portable export."""
    if not svg.startswith('<svg'):
        svg = svg[svg.index('<svg'):svg.index('</svg>') + 6]
    element = ET.fromstring(svg)
    element.set('xmlns', 'http://www.w3.org/2000/svg')
    style = ET.Element('style')
    style.text = ('svg{--ink:#202831;--surface-2:#eef2f6;--c-cbc:#26767a;'
                  '--c-p5:#536caa;--c-p7:#9c5776;--c-r2:#737c38;--c-r1:#697fb3;}'
                  'text{font-family:monospace;fill:#202831;}'
                  '.chem-strand-box{fill:#e4eef0;stroke:#68858c;stroke-width:.6;}'
                  '.chem-svg.long{font-size:11.8px;}.chem-svg.small{font-size:13.4px;}'
                  '.fold-pair{stroke:#879a9f;stroke-width:1;}'
                  '.fold-base{font-family:monospace;font-size:11px;}'
                  '.fold-label,.fold-end{font-family:sans-serif;font-size:11px;}'
                  '.chem-cbc{fill:#26767a;}.chem-unp{fill:#66717b;}')
    element.insert(0, style)
    return ET.tostring(element, encoding='unicode')


def render(first, second, candidate):
    """Validate the selected arrangement again, then reuse the shared SVG renderer."""
    first, second = normalize(first), normalize(second)
    fields = dict(candidate)
    kind = fields.pop('kind', None)
    if kind == 'linear':
        model = LinearCandidate(**fields)
        svg = panel_svg(model.scene(first, second).rows())
        summary = f'{model.length} exact base pairs'
        coordinates = {
            'first': [model.first_start, model.first_start + model.length],
            'second': [model.second_start, model.second_start + model.length],
            'first_overhangs': [model.first_start, len(first) - model.first_start - model.length],
            'second_overhangs': [model.second_start, len(second) - model.second_start - model.length]}
    elif kind == 'padlock':
        model = PadlockCandidate(**fields)
        topology = model.topology(first, second)
        svg = topology.folded_svg()
        e_start, e_end = model.extension_start, model.extension_start + model.extension_length
        l_start, l_end = model.ligation_start, model.ligation_start + model.ligation_length
        gap = [e_end, l_start] if e_end <= l_start else [l_end, e_start]
        probe_length = len(second if model.target == 1 else first)
        coordinates = {'target_sequence': model.target, 'probe_sequence': 3 - model.target,
                       'target_extension_arm': [e_start, e_end], 'target_ligation_arm': [l_start, l_end],
                       'probe_extension_arm': [probe_length - model.extension_length, probe_length],
                       'probe_ligation_arm': [0, model.ligation_length],
                       'probe_backbone': [model.ligation_length, probe_length - model.extension_length],
                       'target_gap': gap}
        summary = (f'{model.ligation_length} + {model.extension_length} paired arm bases; '
                   f'{len(topology.backbone)} backbone bases; {gap[1] - gap[0]} target gap bases')
    else:
        raise ValueError('Unknown candidate type.')
    standalone = _standalone(svg)
    return {'svg': standalone, 'summary': summary, 'coordinates': coordinates,
            'id': hashlib.sha256((first + second + repr(candidate)).encode()).hexdigest()[:16]}
