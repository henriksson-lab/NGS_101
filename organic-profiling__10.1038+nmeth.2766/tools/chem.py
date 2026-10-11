from protein_dna_mapping import get
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Scene, Segment, Workflow

globals().update(vars(get("organic-profiling")))
_base_workflow = workflow


def adapter_geometry():
    adapter = Scene(); adapter.strand("universal adapter", [
        Segment("published universal-adapter body", "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC"),
        Segment("terminal phosphorothioate dT", "T")], label="Universal Adapter")
    adapter.ligation_end("universal adapter", "terminal phosphorothioate dT",
                         kind="cohesive", polarity="3-prime", owner="adapter",
                         text="published 3′-dT cohesive end")
    insert = Scene(); insert.strand("footprint", [Segment("native footprint", "X" * 28, placeholder=True),
                                                   Segment("footprint 3-prime dA", "A")],
                                        label="end-prepared footprint", mod5="p")
    insert.ligation_end("footprint", "footprint 3-prime dA", kind="cohesive",
                        polarity="3-prime", owner="substrate", text="3′-dA")
    return AdapterLigationGeometry.known(
        "cohesive", adapter, insert,
        statement="The authors' detailed protocol prints the Universal Adapter with a terminal phosphorothioate-protected dT and separately specifies one-base dA tailing of the footprint.")


def workflow():
    old = _base_workflow(); wf = Workflow(old.initial)
    for reaction in old.reactions:
        if reaction.action == "Ligate indexed adapters and PCR-enrich":
            wf.ligate_adapter(reaction.action, adapter_geometry(), reaction.after.rows,
                              name=reaction.after.name, note=reaction.note)
        else:
            wf.react(reaction.action, reaction.after.rows,
                     name=reaction.after.name, note=reaction.note)
    return wf
