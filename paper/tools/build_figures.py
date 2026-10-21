#!/usr/bin/env python3
"""Generate editable paper SVGs with Python's standard library only."""
from __future__ import annotations

import ast
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'paper' / 'figures'
NS = 'http://www.w3.org/2000/svg'
ET.register_namespace('', NS)
INK = '#172b3a'
MUTED = '#516373'
BLUE = '#285c88'
TEAL = '#176c68'
BORDER = '#bdcbd5'
PALE = '#edf4f8'
GREEN = '#e7f3ed'
AMBER = '#fff3dc'
GREY = '#eff1f4'
RED = '#fae8e8'


def enums(path, name):
    tree = ast.parse((ROOT / path).read_text())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == name)
    return [(n.targets[0].id, ast.literal_eval(n.value)) for n in cls.body
            if isinstance(n, ast.Assign)]


KINDS = enums('lib/protocol_graph.py', 'StepKind')
RELATIONS = enums('lib/protocol_graph.py', 'MolecularRelation')
STATUSES = enums('lib/chemistry_validation.py', 'Status')
RELATION_TEXT = {
    'CIRCLE_CLOSURE': ('Circle closure', 'Linear substrate and circle preserve strand order'),
    'RESTRICTION_FILL': ('Restriction-end fill', 'New tracts follow the declared cut overhang'),
    'ROLLING_CIRCLE_COPY': ('Rolling-circle copying', 'Repeats and primer landing derive from the circle'),
    'CONTACT_JUNCTION': ('Contact junction', 'Product contains the derived restriction junction'),
    'HAIRPIN_STEM': ('Hairpin stem', 'Declared intramolecular stem can base-pair'),
    'STRAND_ORDER': ('Strand-order preservation', 'Declared covalent sequence remains unchanged'),
    'END_COMBINATIONS': ('End combinations', 'End classes match declared downstream primer sites'),
}
assert set(RELATION_TEXT) == {name for name, _ in RELATIONS}, 'Update figure for new relations'
assert [name for name, _ in STATUSES] == ['ESTABLISHED', 'CONDITIONAL', 'UNRESOLVED', 'CONTRADICTORY']


class SVG:
    def __init__(self, name, height, title, description):
        self.name = name
        self.height = height
        self.root = ET.Element(f'{{{NS}}}svg', {
            'width': '180mm', 'height': f'{height / 8:g}mm',
            'viewBox': f'0 0 1440 {height}', 'role': 'img',
            'aria-labelledby': 'title description',
        })
        ET.SubElement(self.root, 'title', id='title').text = title
        ET.SubElement(self.root, 'desc', id='description').text = description
        defs = ET.SubElement(self.root, 'defs')
        marker = ET.SubElement(defs, 'marker', {
            'id': 'arrow', 'viewBox': '0 0 10 10', 'refX': '9', 'refY': '5',
            'markerWidth': '8', 'markerHeight': '8', 'orient': 'auto-start-reverse',
        })
        ET.SubElement(marker, 'path', d='M 0 0 L 10 5 L 0 10 z', fill=MUTED)
        self.add('rect', x=0, y=0, width=1440, height=height, fill='white')

    def add(self, tag, **attrs):
        return ET.SubElement(self.root, tag, {
            key.replace('_', '-'): str(value) for key, value in attrs.items()
        })

    def text(self, x, y, lines, *, size=22, weight=400, color=INK,
             leading=31, anchor='start'):
        if isinstance(lines, str):
            lines = [lines]
        node = self.add('text', x=x, y=y, font_family='Arial, Helvetica, sans-serif',
                        font_size=size, font_weight=weight, fill=color,
                        text_anchor=anchor)
        for i, line in enumerate(lines):
            ET.SubElement(node, 'tspan', x=str(x), dy=str(0 if i == 0 else leading)).text = line
        return node

    def box(self, x, y, w, h, *, fill='white', stroke=BORDER, radius=12, dash=False):
        attrs = dict(x=x, y=y, width=w, height=h, rx=radius,
                     fill=fill, stroke=stroke, stroke_width=1.6)
        if dash:
            attrs['stroke_dasharray'] = '7 5'
        return self.add('rect', **attrs)

    def arrow(self, points, *, dash=False):
        attrs = {'d': 'M ' + ' L '.join(f'{x} {y}' for x, y in points),
                 'fill': 'none', 'stroke': MUTED, 'stroke_width': 2,
                 'marker_end': 'url(#arrow)'}
        if dash:
            attrs['stroke_dasharray'] = '6 5'
        return self.add('path', **attrs)

    def panel(self, letter, title, y):
        self.text(48, y, letter, size=26, weight=700, color=BLUE)
        self.text(87, y, title, size=25, weight=700)

    def save(self):
        _indent_svg(self.root)
        p = OUT / self.name
        p.write_bytes(ET.tostring(self.root, encoding='utf-8', xml_declaration=True))
        print(p.relative_to(ROOT))


def definition_figure():
    s = SVG('figure_1_protocol_definition.svg', 1080,
            'Protocol-owned definitions form a typed provenance hypergraph',
            'A: source-qualified materials, steps, outcomes and evidence. '
            'B: an illustrative material-step graph with explicit reagent roles '
            'and a physical cleanup branch. C: views derived from the declarations.')
    s.text(48, 53, 'Protocol definitions', size=36, weight=700)
    s.text(48, 88, 'Source-specific Python objects; shared computation and rendering', color=MUTED)
    s.panel('A', 'What each protocol declares', 139)
    cards = [
        ('Materials', ['Material kind and identity', 'Sequence or symbolic state', 'Features and modifications', 'Inference and provenance']),
        ('Steps', ['Operation label (StepKind)', 'Inputs with explicit roles', 'Products and geometry', 'Optional molecular relations']),
        ('Outcomes', ['Retained and residual classes', 'Physical discard, QC or gaps', 'Sequencing plan / limitation', 'Read and identifier layout']),
        ('Evidence', ['Defining source and revision', 'Manual source inventory', 'Review scope and exclusions', 'Discrepancies and limitations']),
    ]
    for i, (title, body) in enumerate(cards):
        x = 48 + i * 342
        s.box(x, 165, 318, 218, fill=PALE)
        s.text(x + 20, 204, title, size=26, weight=700, color=BLUE)
        s.text(x + 20, 249, body, size=21, leading=32)
    s.panel('B', 'Materials and steps form a provenance hypergraph', 440)
    # Reagents are material nodes; roles are carried on the input ports.
    s.box(285, 491, 225, 82, fill=PALE, radius=28)
    s.text(397, 524, ['Primers', 'structural input'], anchor='middle', size=21, leading=27)
    s.box(545, 491, 250, 82, fill=PALE, radius=28)
    s.text(670, 524, ['Catalyst + conditions', 'role-labelled inputs'], anchor='middle', size=21, leading=27)
    s.box(855, 491, 225, 82, fill=PALE, radius=28)
    s.text(967, 524, ['Selection reagent', 'selection input'], anchor='middle', size=21, leading=27)
    s.arrow([(397, 573), (397, 630)])
    s.arrow([(670, 573), (670, 599), (463, 599), (463, 630)])
    s.arrow([(967, 573), (967, 630)])
    s.arrow([(258, 671), (300, 671)])
    s.arrow([(510, 671), (560, 671)])
    s.arrow([(810, 671), (855, 671)])
    s.arrow([(1080, 671), (1130, 671)])
    s.arrow([(967, 714), (967, 784)])
    s.box(48, 620, 210, 102, fill=PALE, radius=28)
    s.text(153, 662, ['Native input', 'molecular state'], anchor='middle', size=22, leading=29)
    s.box(300, 630, 210, 84, fill='#e4efef', stroke=TEAL)
    s.text(405, 664, ['PCR', 'declared geometry'], anchor='middle', size=22, leading=29)
    s.box(560, 620, 250, 102, fill=PALE, radius=28)
    s.text(685, 649, ['Mixture: daughters,', 'possible originals', 'and residuals'], anchor='middle', size=21, leading=27)
    s.box(855, 630, 225, 84, fill='#e4efef', stroke=TEAL)
    s.text(967, 664, ['Physical cleanup', 'declared recovery'], anchor='middle', size=22, leading=29)
    s.box(1130, 620, 262, 102, fill=PALE, radius=28)
    s.text(1261, 662, ['Retained mixture', 'sequencing plan'], anchor='middle', size=22, leading=29)
    s.box(855, 784, 225, 82, fill=GREY, radius=28)
    s.text(967, 817, ['Removed fraction', 'discard terminal'], anchor='middle', size=21, leading=27)
    s.box(48, 784, 747, 82, fill='#f7f9fa', stroke='#dce4e9')
    s.text(68, 817, ['Illustrative topology; no experimental recipe is specified.',
                    'Independent polymers remain distinct within a mixture.'], size=21, leading=28, color=MUTED)
    s.text(1130, 804, ['Rounded nodes: materials', 'Rectangles: operations', 'Arrows: input / output ports'], size=20, leading=28, color=MUTED)
    s.panel('C', 'Shared readers use the same declarations', 931)
    for x, title, body in [
        (48, 'Schematic rendering', 'States, reactions and evidence notices'),
        (504, 'Sequencing interpretation', 'Primer sites, read cycles and identifiers'),
        (960, 'Validation and audit', 'Graph contracts and scoped chemistry claims'),
    ]:
        s.box(x, 958, 432, 86, fill='white')
        s.text(x + 18, 992, title, size=23, weight=700, color=BLUE)
        s.text(x + 18, 1023, body, size=19, color=MUTED)
    s.save()


def validation_figure():
    s = SVG('figure_2_validation_architecture.svg', 1160,
            'Broad operation vocabulary, explicitly declared narrow molecular checks',
            'Operation labels do not automatically map to molecular relations. '
            'The author supplies states, geometry and optional relation tags. '
            'Seven relation checks and four typed-state check families yield scoped '
            'claims with four statuses; source fidelity is assessed separately.')
    s.text(48, 53, 'Validation of declared molecular claims', size=36, weight=700)
    s.text(48, 88, 'A reusable check kernel acts on explicit states and relations', color=MUTED)
    s.panel('A', 'Operation vocabulary', 145)
    s.text(390, 145, 'B', size=26, weight=700, color=BLUE)
    s.text(428, 145, 'Explicit declarations', size=25, weight=700)
    s.text(750, 145, 'C', size=26, weight=700, color=BLUE)
    s.text(788, 145, 'Narrow relation kernel', size=25, weight=700)
    s.box(48, 170, 282, 486, fill=PALE)
    s.text(68, 211, f'{len(KINDS)} operation labels', size=25, weight=700, color=BLUE)
    s.text(68, 243, 'StepKind (examples)', size=20, color=MUTED)
    examples = ['PCR', 'Reverse transcription', 'Ligation', 'Circularization', 'Tagmentation',
                'Fragmentation', 'End repair', 'Chemical conversion', 'Folding',
                'Cleanup', 'Sampling, split and pool']
    s.text(68, 284, examples, size=22, leading=31)
    s.box(390, 223, 300, 375, fill='#e4efef', stroke=TEAL)
    s.text(410, 263, ['Protocol supplies'], size=25, weight=700, color=TEAL)
    s.text(410, 307, ['Input / output states', 'Actual template geometry', 'Optional relation tags'], size=22, leading=37)
    s.add('line', x1=410, y1=411, x2=670, y2=411, stroke=BORDER)
    s.text(410, 448, ['No automatic', 'StepKind → relation', 'lookup is performed.'], size=23, weight=700, leading=32)
    s.text(410, 562, ['States determine what', 'can actually be checked.'], size=20, color=MUTED, leading=27)
    s.arrow([(330, 405), (390, 405)])
    s.arrow([(690, 405), (750, 405)])
    s.box(750, 170, 642, 486, fill='white')
    s.text(770, 211, f'{len(RELATIONS)} explicit MolecularRelation types', size=25, weight=700, color=BLUE)
    for i, (name, _) in enumerate(RELATIONS):
        title, detail = RELATION_TEXT[name]
        y = 254 + i * 55
        s.text(770, y, title, size=22, weight=700)
        s.text(770, y + 25, detail, size=19, color=MUTED)
    s.text(48, 692, ['Step labels enforce graph-shape contracts.',
                    'Generic and unmodelled labels are included.'], size=19, color=MUTED, leading=27)
    s.arrow([(540, 598), (540, 728)])
    s.box(390, 728, 1002, 145, fill=PALE)
    s.text(410, 766, 'Additional checks from typed states and geometry', size=24, weight=700, color=BLUE)
    for x, name, data in [
        (410, 'Base pairing', 'Scene columns'),
        (656, 'PCR primer sites', 'PCRGeometry'),
        (902, 'Restriction cuts', 'Digest'),
        (1148, 'View agreement', 'RCA state vs drawing'),
    ]:
        s.text(x, 811, name, size=22, weight=700)
        s.text(x, 842, data, size=19, color=MUTED)
    s.panel('D', 'Each checked question produces its own claim status', 926)
    status_info = [
        ('Established', ['Local relation supported', 'by exact declared data'], GREEN, TEAL),
        ('Conditional', ['Compatible under uncertain', 'or inferred bases'], AMBER, '#855f13'),
        ('Unresolved', ['Missing relation, states', 'or required geometry'], GREY, MUTED),
        ('Contradictory', ['Declared states or geometry', 'disagree with the check'], RED, '#953838'),
    ]
    for i, (name, body, fill, color) in enumerate(status_info):
        x = 48 + i * 342
        s.box(x, 949, 318, 119, fill=fill)
        s.text(x + 18, 988, name, size=25, weight=700, color=color)
        s.text(x + 18, 1022, body, size=21, color=INK, leading=28)
    s.text(48, 1111, 'A successful local claim does not establish the whole reaction or protocol.', size=23, weight=700)
    s.text(48, 1143, 'Source fidelity, graph completeness and experimental performance are separate questions.', size=21, color=MUTED)
    s.save()


def _indent_svg(element, level=0):
    """Pretty-print structure without inserting whitespace into molecular text."""
    local = element.tag.rsplit('}', 1)[-1]
    if local in {'text', 'tspan', 'style'} or element.get('{http://www.w3.org/XML/1998/namespace}space') == 'preserve':
        return
    children = list(element)
    if not children:
        return
    if not element.text or not element.text.strip():
        element.text = '\n' + '  ' * (level + 1)
    for i, child in enumerate(children):
        _indent_svg(child, level + 1)
        if not child.tail or not child.tail.strip():
            child.tail = '\n' + '  ' * (level if i == len(children) - 1 else level + 1)


def _embed_model(figure, markup, x, y, width, height, viewbox=None):
    """Embed the actual model SVG, preserving its text and semantic titles."""
    root = ET.fromstring(markup)
    if root.tag.rsplit('}', 1)[-1] != 'svg':
        root = next(e for e in root.iter() if e.tag.rsplit('}', 1)[-1] == 'svg')
    palette = {
        '--ink': INK, '--ink-muted': MUTED, '--surface': 'white',
        '--surface-2': '#f2f6f9', '--c-cbc': TEAL, '--c-r1': BLUE,
        '--c-r2': TEAL, '--c-s5': '#7755a2', '--note-rule': '#9a6b28',
    }
    for e in root.iter():
        if not e.tag.startswith('{'):
            e.tag = f'{{{NS}}}' + e.tag
        for key, value in list(e.attrib.items()):
            for name, color in palette.items():
                value = value.replace(f'var({name})', color)
            e.set(key, value)
    # The paper supplies consistent legends outside cropped molecular viewports.
    if 'fold-svg' in root.get('class', ''):
        for e in list(root):
            if e.get('class') == 'fold-label' or (e.tag.rsplit('}', 1)[-1] == 'line' and e.get('stroke-width') == '4'):
                root.remove(e)
    root.set('x', str(x)); root.set('y', str(y))
    root.set('width', str(width)); root.set('height', str(height))
    if viewbox is not None:
        root.set('viewBox', ' '.join(str(v) for v in viewbox))
    style = ET.Element(f'{{{NS}}}style')
    style.text = f'''
.chem-svg {{font-family: "DejaVu Sans Mono", "Liberation Mono", monospace; fill:{INK};}}
.chem-svg.small {{font-size:13.4px;}}
.chem-strand-box {{fill:#f2f6f9;stroke:{MUTED};stroke-width:1;}}
.chem-binding-span {{fill:none;stroke:{BLUE};stroke-width:1.5;}}
.chem-binding-label {{fill:{MUTED};font-family:Arial,Helvetica,sans-serif;font-size:.92em;}}
.chem-p5 {{fill:{BLUE};}} .chem-umi {{fill:#7755a2;}} .chem-r2 {{fill:{TEAL};}}
.chem-inferred {{text-decoration:underline;text-decoration-style:dotted;text-decoration-color:#9a6b28;}}
.fold-svg .fold-pair {{stroke:{MUTED};stroke-width:1;opacity:.55;}}
.fold-svg .fold-base {{font-size:11px;font-weight:600;}}
.fold-svg .fold-end {{fill:{MUTED};font-size:11px;}}
.fold-svg .fold-label {{fill:{INK};font-family:Arial,Helvetica,sans-serif;font-size:12px;}}
.fold-svg .fold-ligation {{stroke:#9a6b28;stroke-width:3;}}
'''
    root.insert(0, style)
    figure.root.append(root)


def annotation_figure():
    from figure_examples import annotation_example
    from chemdraw import panel_svg
    construct, scene = annotation_example()
    s = SVG('figure_3_sequence_annotation.svg', 1120,
            'Sequence annotations drive semantic molecular rendering',
            'Sequence segments carry molecular content, identifier features, styling and inference. '
            'Scene assembly computes pairing coordinates and typed visuals. The actual shared '
            'renderer generates selectable SVG text and geometry for a synthetic example.')
    s.text(48, 53, 'From annotated sequence to SVG', size=36, weight=700)
    s.text(48, 88, 'Molecular content, identifier meaning and presentation are separate declarations', color=MUTED)
    s.panel('A', 'Annotate named sequence segments', 145)
    cards = [
        ('Molecular content', ['Segment: name + sequence', 'Placeholder / inferred flags', 'Known or unassigned length', 'Strand-end modifications']),
        ('Identifier semantics', ['MolecularFeature: stable ID', 'Role: UMI, barcode, index…', 'Encoding and group / part', 'Independent of colour tags']),
        ('Presentation', ['Segment.tag: display colour', 'Marks anchor to segments', 'Polarity follows strand order', 'Tooltips derive from content']),
    ]
    for i, (title, body) in enumerate(cards):
        x = 48 + i * 456
        s.box(x, 170, 432, 206, fill=PALE)
        s.text(x + 20, 211, title, size=26, weight=700, color=BLUE)
        s.text(x + 20, 254, body, size=22, leading=31)
    s.panel('B', 'Assemble once; derive geometry and text', 433)
    for x, title, body in [
        (48, 'Construct + Scene', ['Order segments; anneal strands', 'Pair named regions; check columns']),
        (540, 'Semantic row metadata', ['StrandVisual, SpanVisual,', 'ArrowVisual and CommentVisual']),
        (1032, 'SVG primitives', ['Polygons and paths', 'Selectable text and tspan regions']),
    ]:
        s.box(x, 460, 360, 127, fill='#e4efef' if x == 540 else 'white')
        s.text(x + 18, 500, title, size=25, weight=700, color=TEAL if x == 540 else BLUE)
        s.text(x + 18, 540, body, size=19, leading=29)
    s.arrow([(408, 523), (540, 523)])
    s.arrow([(900, 523), (1032, 523)])
    s.panel('C', 'Actual shared-renderer output from synthetic annotations', 644)
    s.box(48, 671, 1344, 250, fill='white')
    _embed_model(s, panel_svg(scene.rows(), cls='small'), 65, 690, 860, 205)
    s.text(980, 718, 'UMI feature', size=25, weight=700, color='#7755a2')
    s.text(980, 756, ['id: paper_umi', 'role: umi', 'encoding: random'], size=22, leading=31)
    s.text(980, 865, ['Colour is presentation;', 'the feature carries meaning.'], size=20, leading=27, color=MUTED)
    s.panel('D', 'Coordinates and hover information are computed', 973)
    s.box(48, 998, 660, 87, fill=PALE)
    s.text(68, 1032, 'Span location: derived from the named UMI segment', size=22)
    s.text(68, 1064, 'Complement and alignment: computed, not hand-spaced', size=21, color=MUTED)
    s.box(732, 998, 660, 87, fill=PALE)
    s.text(752, 1032, 'Example hover title from Segment.hover_text()', size=22)
    s.text(752, 1064, construct.segments[0].hover_text(), size=21, color=MUTED)
    s.save()


def _fold_view(layout, width=340, height=145):
    points = [point for strand in layout.strands.values() for point in strand.points]
    cx = (min(x for x, _ in points) + max(x for x, _ in points)) / 2
    cy = (min(y for _, y in points) + max(y for _, y in points)) / 2
    return (round(cx - width / 2, 2), round(cy - height / 2, 2), width, height)


def folding_figure():
    from figure_examples import hairpin_example, padlock_example
    padlock, hairpin = padlock_example(), hairpin_example()
    pad_layout, hair_layout = padlock.folded_layout(), hairpin.folded_layout()
    s = SVG('figure_4_constraint_checked_folding.svg', 1110,
            'Constraint-checked geometric layouts for padlocks and hairpins',
            'Declared strand identity and pairing constraints feed topology-specific placement '
            'and flexible routing. FoldLayout checks pairing, collinearity and crossings. '
            'Actual synthetic padlock and hairpin drawings illustrate this schematic layout; '
            'it does not predict free energy or a preferred molecular conformation.')
    s.text(48, 53, 'Constraint-checked folding views', size=36, weight=700)
    s.text(48, 88, 'Sequence-resolved geometric layout of a declared topology', color=MUTED)
    s.panel('A', 'From molecular constraints to a drawing', 145)
    stages = [
        ('Declare topology', ['Strand identities + sequences', 'Paired spans / arm sites', 'Straight hints and bonds']),
        ('Place and route', ['Topology-specific placement', 'Flexible spacer coordinates', 'Candidate routes where used']),
        ('Validate layout', ['Antiparallel pairing + bounds', 'Required collinearity', 'No unintended crossing']),
        ('Render SVG', ['One point per nucleotide', 'Strand paths and pair rungs', 'Polarity and junction marks']),
    ]
    for i, (title, body) in enumerate(stages):
        x = 48 + i * 342
        s.box(x, 170, 318, 191, fill=PALE)
        s.text(x + 18, 210, title, size=25, weight=700, color=BLUE)
        s.text(x + 18, 252, body, size=20, leading=32)
        if i < 3:
            s.arrow([(x + 318, 331), (x + 342, 331)])
    s.text(48, 405, ['Some complex topologies compare deterministic routes using length and bend costs.',
                    'The examples below use topology-specific routes, rather than a route search.'], size=22, leading=31, color=MUTED)
    s.panel('B', 'Actual renderer output: synthetic examples', 496)
    s.box(48, 522, 660, 350, fill='white')
    s.box(732, 522, 660, 350, fill='white')
    s.text(68, 563, 'Padlock capture', size=27, weight=700, color=TEAL)
    s.text(68, 596, 'Two arm sites on one target; flexible backbone', size=21, color=MUTED)
    _embed_model(s, pad_layout.svg(620, 330, 'Synthetic padlock capture'),
                 68, 611, 620, 245, _fold_view(pad_layout))
    s.text(752, 563, 'Hairpin', size=27, weight=700, color=BLUE)
    s.text(752, 596, 'Two declared stem spans; an unpaired loop', size=21, color=MUTED)
    _embed_model(s, hair_layout.svg(520, 280, 'Synthetic hairpin'),
                 752, 611, 620, 245, _fold_view(hair_layout))
    for x, color, label in [(80, INK, 'target'), (264, TEAL, 'padlock probe'), (764, BLUE, 'hairpin')]:
        s.add('line', x1=x, y1=850, x2=x + 25, y2=850, stroke=color, stroke_width=4)
        s.text(x + 34, 856, label, size=19, color=MUTED)
    s.panel('C', 'What this establishes', 927)
    s.box(48, 953, 660, 116, fill=GREEN)
    s.text(68, 989, 'Geometric and sequence consistency', size=25, weight=700, color=TEAL)
    s.text(68, 1024, ['Declared pairing, strand identity and layout constraints',
                     'are checked for these supplied molecular objects.'], size=22, leading=29)
    s.box(732, 953, 660, 116, fill=AMBER)
    s.text(752, 989, 'Thermodynamic prediction is outside this model', size=24, weight=700, color='#855f13')
    s.text(752, 1024, ['No free-energy minimization or preferred-fold search;',
                     'a schematic fold does not establish reaction yield.'], size=22, leading=29)
    s.save()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    definition_figure()
    validation_figure()
    annotation_figure()
    folding_figure()
    for filename, rows in [('operation_vocabulary.tsv', KINDS), ('molecular_relations.tsv', RELATIONS)]:
        with (OUT / filename).open('w', newline='') as f:
            writer = csv.writer(f, delimiter='\t', lineterminator='\n')
            writer.writerow(['identifier', 'label'])
            writer.writerows(rows)
    sources = ['lib/protocol_graph.py', 'lib/chemistry_validation.py', 'lib/chemdraw.py',
               'lib/seqprimers.py', 'lib/pcr.py', 'lib/topology.py', 'lib/page.py',
               'paper/tools/figure_examples.py', 'DESIGN.md']
    snapshot = {
        'schema_version': 1,
        'commit': subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip(),
        'counts': {'operation_labels': len(KINDS), 'relation_types': len(RELATIONS), 'claim_statuses': len(STATUSES)},
        'source_sha256': {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in sources},
        'scope': 'Implementation architecture, not a measurement of validation coverage or protocol correctness.',
    }
    (OUT / 'implementation_snapshot.json').write_text(json.dumps(snapshot, indent=2) + '\n')


if __name__ == '__main__':
    main()
