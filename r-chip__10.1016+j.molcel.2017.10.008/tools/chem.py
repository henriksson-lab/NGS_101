from chemdraw import Workflow
from hybrid_base_protocols import get, rchip_adapter_geometry
SPEC = get("r-chip")
SPEC.CAVEAT = ("The defining paper prints the tailed N9 primer; the authors’ detailed "
               "2019 protocol also prints both adapter strands and the PCR primers. "
               "The final adapter body remains a role-level segment pending reconstruction "
               "of its one-ended amplification product.")
SPEC.SEQUENCING_UNAVAILABLE = (
    "The library-preparation PCR primers are published; a separate instrument "
    "sequencing-primer synthesis oligo is not specified for this run.")
SPEC.oligos = lambda: (
    "<b>Tailed N9 primer</b> — 5′-/invddt/CAAGCAGAAGACGGCATACGAGNNNNNNNNN-3′",
    "<b>Adapter oligo A</b> — 5′-/Phos/GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT-3′",
    "<b>Adapter oligo B</b> — 5′-AGACGTGTGCTCTTCCGATCT-3′",
    "<b>PCR primer</b> — 5′-CAAGCAGAAGACGGCATACGAG-3′",
    "<b>Barcoded PCR primer</b> — 5′-AATGATACGGCGACCACCGAGATCTACACNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT-3′",
)


def adapter_geometry():
    return rchip_adapter_geometry()


_base_workflow = SPEC.workflow


def workflow():
    old = _base_workflow()
    wf = Workflow(old.initial)
    for reaction in old.reactions:
        if "Ligate one Illumina adapter" in reaction.action:
            wf.ligate_adapter(reaction.action, adapter_geometry(), reaction.after.rows,
                              name=reaction.after.name, note=reaction.note)
        else:
            wf.react(reaction.action, reaction.after.rows,
                     name=reaction.after.name, note=reaction.note)
    return wf


SPEC.workflow = workflow
