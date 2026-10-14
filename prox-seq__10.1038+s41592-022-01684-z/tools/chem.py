from adapter_ligation import AdapterLigationGeometry
from chemdraw import Scene, Segment, Workflow
from rna_cell_assays import get
from proximity_topology import ProximityLigationTopology

SPEC = get("prox-seq")
SPEC.CAVEAT = ("The schematic combines the invariant PLA reporter across published Drop-seq, "
               "Smart-seq2 and 10x implementations. Supplementary Tables 7–10 print "
               "branch-specific oligos and custom sequencing primers; the schematic "
               "does not choose one platform-specific library end.")
SPEC.SEQUENCING_UNAVAILABLE = (
    "The paper publishes custom sequencing primers for its droplet and plate branches, "
    "but a single primer-binding map cannot describe all Prox-seq back ends at once.")
SPEC.oligos = lambda: (
    "<b>Universal ligation connector</b> — 5′-TTTCACGACACGACACGATTTAGGTC-3′",
    "<b>Droplet Read 1</b> — 5′-GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC-3′",
    "<b>Droplet Read 2</b> — 5′-AGATGCATCGTCTCGTGGGCTCGG-3′",
    "<b>Droplet i7</b> — 5′-CCGAGCCCACGAGACGATGCATCT-3′",
    "<b>Plate Read 1</b> — 5′-TATTTGGCATCGTCTCGTGGGCTCGG-3′",
    "<b>Plate i5</b> — 5′-CCGAGCCCACGAGACGATGCCAAATA-3′",
    "<b>Plate i7</b> — 5′-ACTCTGCGTTGATACCACTGCTTCGGC-3′",
)
SPEC.TOPOLOGY_PANELS = ProximityLigationTopology().panels()
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
    published = Scene()
    published.strand("published connector", [Segment(
        "full connector oligo", "TTTCACGACACGACACGATTTAGGTC")],
        label="source-written 26-base connector")
    return AdapterLigationGeometry.known(
        "splinted", sc, published,
        statement="The defining paper prints the 26-base connector and says 20 reporter bases "
                  "hybridize across the ligation junction. It does not assign the 20 paired "
                  "bases to exact positions on each antibody-probe end; the paired span "
                  "therefore remains role-level while the complete connector is shown separately.")


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
