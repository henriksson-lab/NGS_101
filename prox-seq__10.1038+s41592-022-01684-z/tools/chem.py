from adapter_ligation import AdapterLigationGeometry
from chemdraw import Scene, Segment, Workflow
from rna_cell_assays import get

SPEC = get("prox-seq")
_base_workflow = SPEC.workflow


def connector_geometry():
    sc = Scene(); sc.strand("connector", [Segment("connector A half", "X" * 10, placeholder=True),
                                           Segment("connector B half", "X" * 10, placeholder=True)],
                                      label="universal connector")
    sc.anneal("probe ends", [Segment("probe-B 5-prime ligating arm", "X" * 10, placeholder=True),
                              Segment("probe-A 3-prime ligating arm", "X" * 10, placeholder=True)],
              to="connector", pair=("probe-B 5-prime ligating arm", "connector B half"),
              label="adjacent antibody-probe ends", mod5="p",
              unpaired=("probe-B 5-prime ligating arm", "probe-A 3-prime ligating arm",
                        "connector A half", "connector B half"))
    sc.ligation_end("probe ends", "probe-B 5-prime ligating arm", kind="splinted",
                    polarity="5-prime", owner="adapter", text="probe 5′ nick partner")
    return AdapterLigationGeometry.known(
        "splinted", sc,
        statement="The defining paper states that 20 bases of the 119-base reporter hybridize to the universal connector during ligation; accessible text does not print those bases, so their lengths and nick geometry are shown without invented sequence.")


def workflow():
    old = _base_workflow()
    wf = Workflow(old.initial)
    for reaction in old.reactions:
        if "universal connector" in reaction.action:
            wf.ligate_adapter(reaction.action, connector_geometry(), reaction.after.rows,
                              name=reaction.after.name, note=reaction.note)
        else:
            wf.react(reaction.action, reaction.after.rows,
                     name=reaction.after.name, note=reaction.note)
    return wf


SPEC.workflow = workflow
