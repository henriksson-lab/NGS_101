"""Bounded workflow conversion for staged molecular scenes."""
from chemdraw import MolecularState, Scene, Workflow


def _workflow(initial: Scene, steps) -> Workflow:
    out = Workflow(MolecularState("Starting material", tuple(initial.rows())))
    for action, state, note in steps:
        out.react(action, state.rows(), name=action, note=note)
    return out
