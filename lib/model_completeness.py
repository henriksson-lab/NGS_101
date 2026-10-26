"""Computed scope of a declared molecular model, separate from chemistry/source proof.

The report travels with a generated graph page. The index reads it without importing
protocol code or inspecting prose for words such as 'missing'. Variable inserts and
barcode placeholders are legitimate model inputs, not automatically missing detail.
"""
from __future__ import annotations

from dataclasses import fields, is_dataclass
from html.parser import HTMLParser
import json

from chemdraw import Construct, Row, Scene, Segment
from protocol_graph import Disposition, MaterialKind, StepKind


LEVELS = {
    "none": "No schematic",
    "unassessed": "Model not assessed",
    "unresolved": "Unresolved steps or material fates",
    "limited": "Molecular or readout limits",
    "inferred": "Model with inferred regions",
    "declared": "No declared model gaps",
}
SIGNALS = ("unmodelled_steps", "unresolved_outputs", "limited_geometry",
           "limited_readouts", "unrepresented_molecules", "inferred_objects")
LIMITS = SIGNALS[:-1]


def _objects(value, seen=None):
    """Walk model data only; never call a protocol's drawing/chemistry methods."""
    if seen is None:
        seen = set()
    if value is None or isinstance(value, (str, bytes, int, float, bool)):
        return
    if id(value) in seen:
        return
    seen.add(id(value))
    yield value
    if isinstance(value, Construct):
        children = value.segments
    elif isinstance(value, Scene):
        children = value.strands.values()
    elif is_dataclass(value) and not isinstance(value, type):
        children = (getattr(value, field.name) for field in fields(value))
    elif isinstance(value, dict):
        children = value.values()
    elif isinstance(value, (tuple, list, set, frozenset)):
        children = value
    else:
        children = ()
    for child in children:
        yield from _objects(child, seen)


def assess(graph) -> dict:
    """Summarize explicit gaps/inference; do not grade reaction-check coverage."""
    values = [*graph.materials.values(), *(s.geometry for s in graph.steps.values())]
    objects = tuple(_objects(values))
    molecular = {MaterialKind.MOLECULAR, MaterialKind.LIBRARY,
                 MaterialKind.OLIGO, MaterialKind.ADAPTER}
    signals = {
        "unmodelled_steps": sum(s.kind is StepKind.MAGIC for s in graph.steps.values()),
        "unresolved_outputs": sum(t.disposition is Disposition.UNRESOLVED
                                  for t in graph.terminals.values()),
        "limited_geometry": sum(getattr(o, "availability", None) in {"unavailable", "regional"}
                                for o in objects),
        "limited_readouts": sum(bool(p.limitation.strip()) for p in graph.sequencing.values()),
        "unrepresented_molecules": sum(
            not any(isinstance(o, Segment) for o in _objects(
                m.chemical_state if m.chemical_state is not None else m.state))
            for m in graph.materials.values() if m.kind in molecular),
        "inferred_objects": sum(
            getattr(o, "inferred", False) is True
            or isinstance(o, Row) and any(chunk[2] for chunk in o.chunks)
            for o in objects),
    }
    return {"version": 1, "signals": signals}


def level(report: dict) -> str:
    """Reject incomplete/old reports conservatively instead of assigning a high level."""
    if (not isinstance(report, dict) or type(report.get("version")) is not int
            or report["version"] != 1):
        return "unassessed"
    signals = report.get("signals")
    if (not isinstance(signals, dict) or set(signals) != set(SIGNALS)
            or any(type(n) is not int or n < 0 for n in signals.values())):
        return "unassessed"
    if signals["unmodelled_steps"] or signals["unresolved_outputs"]:
        return "unresolved"
    if any(signals[key] for key in LIMITS):
        return "limited"
    return "inferred" if signals["inferred_objects"] else "declared"


def markup(graph) -> str:
    return ('<script type="application/json" data-model-completeness="1">'
            + json.dumps(assess(graph), separators=(",", ":")) + '</script>')


class _Reports(HTMLParser):
    def __init__(self):
        super().__init__()
        self.reports = []
        self.current = None

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "script" and "data-model-completeness" in attrs:
            if (attrs.get("data-model-completeness") != "1"
                    or attrs.get("type") != "application/json"):
                self.reports.append(None)
            else:
                self.current = []

    def handle_data(self, data):
        if self.current is not None:
            self.current.append(data)

    def handle_endtag(self, tag):
        if tag == "script" and self.current is not None:
            try:
                report = json.loads("".join(self.current))
            except (ValueError, TypeError):
                report = None
            self.reports.append(report)
            self.current = None


def page_level(page_html: str) -> str:
    parser = _Reports()
    parser.feed(page_html)
    if parser.current is not None or not parser.reports:
        return "unassessed"
    levels = [level(report) for report in parser.reports]
    # Every represented graph must supply a current, readable report. Once assessed,
    # combine multiple graph views conservatively rather than selecting the best one.
    for candidate in ("unassessed", "unresolved", "limited", "inferred", "declared"):
        if candidate in levels:
            return candidate
    return "unassessed"
