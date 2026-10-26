"""Behavior checks for the index's computed model-scope classification."""
from dataclasses import replace
import json
from pathlib import Path
import sys
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from chemdraw import Construct, Scene, Segment
from model_completeness import assess, level, markup, page_level, SIGNALS
from protocol_graph import (ProtocolGraph, Disposition, MaterialKind as K,
                            InputRole as R, StepKind as S, SequencingPlan, Terminal)


def fixture(inferred=False, variable=False):
    g = ProtocolGraph("fixture", authority="fixture source", required_output=Disposition.CONSTRUCT)
    con = Construct((Segment("reference", "NNNN" if variable else "ACGT",
                             inferred=inferred, placeholder=variable),))
    a = g.input("input", "input reference", K.MOLECULAR, primary=True,
                state=Scene.duplex(tuple(con)))
    b = g.product("output", "output reference", K.MOLECULAR, state=con)
    waste = g.product("waste", "cleanup waste", K.WASTE)
    g.react("recover", "recover reference", S.CLEANUP, [(a, R.TEMPLATE)], [b, waste])
    g.discard(waste, "fixture cleanup waste")
    g.constructed(b, "fixture construct")
    g.declare_reviewed_scope(("input", "recover", "output"), "fixture source scope")
    return g.require_publishable()


def run(check):
    check.section("computed model scope is separate from source and reaction validation")
    g = fixture()
    check("an explicitly represented reference has no declared gaps", level(assess(g)), "declared")
    check("a deliberate variable insert is not missing model detail",
          level(assess(fixture(variable=True))), "declared")
    check("inferred sequence remains distinct from declared model limits",
          level(assess(fixture(inferred=True))), "inferred")
    g.steps["recover"] = replace(g.steps["recover"], kind=S.MAGIC, reason="unmodelled chemistry")
    check("an explicit unmodelled transition lowers inferred/known references",
          level(assess(g)), "unresolved")
    g = fixture()
    g.terminals["output"] = Terminal("output", Disposition.UNRESOLVED, "unassigned fate")
    check("unresolved material fate is an explicit workflow limit", level(assess(g)), "unresolved")
    g = fixture()
    g.materials["output"] = replace(g.materials["output"], state=None)
    check("a missing molecular state cannot receive the highest level", level(assess(g)), "limited")
    g = fixture()
    g.steps["recover"] = replace(g.steps["recover"], geometry=SimpleNamespace(availability="unavailable"))
    check("unavailable geometry is an explicit model limit", level(assess(g)), "limited")
    g = fixture()
    # A limitation can explain primer-free acquisition as well as an unresolved end.
    # Conservatively classify it as a declared limit without calling it a primer gap.
    g.sequencing["output"] = SequencingPlan(limitation="Primer-free nanopore acquisition")
    check("a declared readout limitation is retained without interpreting its prose",
          level(assess(g)), "limited")
    g = fixture()
    parsed = json.loads(markup(g).split(">", 1)[1].split("</script>")[0])
    check("generated metadata carries no asserted source-review or validation grade",
          set(parsed), {"version", "signals"})
    check("index recovers the computed level from a generated report", page_level(markup(g)), "declared")
    check("an older page without a report stays unassessed", page_level('<div data-protocol-complete="true">'), "unassessed")
    check("malformed report stays unassessed", page_level('<script type="application/json" data-model-completeness="1">oops</script>'), "unassessed")
    check("a truncated report stays unassessed", page_level(markup(g)[:-9]), "unassessed")
    report = assess(g)
    for name, broken in (
            ("old version", report | {"version": 0}),
            ("boolean version", report | {"version": True}),
            ("missing signal", {"version": 1, "signals": {}}),
            ("negative signal", report | {"signals": dict.fromkeys(SIGNALS, -1)}),
            ("boolean counts", report | {"signals": dict.fromkeys(SIGNALS, False)})):
        check(name + " stays unassessed", level(broken), "unassessed")
    check("a stale outer report marker stays unassessed",
          page_level(markup(g).replace('data-model-completeness="1"', 'data-model-completeness="0"')), "unassessed")
    check("multiple graph reports combine conservatively",
          page_level(markup(g) + markup(fixture(inferred=True))), "inferred")
    check("one unassessed graph view prevents an optimistic aggregate",
          page_level(markup(g) + '<script data-model-completeness="0"></script>'), "unassessed")


if __name__ == "__main__":
    from checks import Check
    c = Check()
    run(c)
    c.report()
