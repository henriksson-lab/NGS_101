"""Source-scoped, role-level molecular path of Chen et al. 2022 superRCA."""
from __future__ import annotations

from chemdraw import Construct, MolecularState, Scene, Workflow
from circular import CircularizationGeometry
from pcr import PCRGeometry
from protocol_model import closed_circle, molecule, role


DOI = "10.1038/s41467-022-31397-y"
TITLE = "superRCA — two-stage rolling-circle mutation detection"
SOURCE = ('<a href="https://doi.org/10.1038/s41467-022-31397-y">'
          'Chen et al., <i>Nature Communications</i> (2022)</a>, Fig. 1 and Methods.')

# The article establishes roles and order, while exact oligos remain unassigned in this bounded pass.
# The original supplementary material remains deferred. Role tokens carry no invented bases or lengths. These constructors are
# intentionally distinct from the exact-sequence circular.rolling_circle() path.
TARGET = role("target region containing the candidate allele")
GENOMIC = Construct((role("genomic flank upstream"), TARGET,
                     role("genomic flank downstream")), "patient genomic DNA")
AMPLICON = (role("forward-primer-defined amplicon end"), TARGET,
            role("reverse-primer-defined amplicon end"))
CIRCULARIZATION = CircularizationGeometry.known(
    "splinted single-strand ligation",
    statement="Fig. 1B and Methods: a ligation-template oligo anneals to both ends of one amplicon strand; Ampligase seals those ends. The main text does not resolve their terminal modification.")
PADLOCK_CLOSURE = CircularizationGeometry.known(
    "splinted single-strand ligation",
    statement="Fig. 1D–E: matched padlock ends are joined by Ampligase while annealed to the repeated target sequence; no gap-fill step is reported.")
PCR = PCRGeometry.unavailable(
    GENOMIC, reason="Targeted forward and reverse primers are used, but their sequences and exact landing sites remain unassigned pending the original supplementary material.")


def cluster(*, labelled: bool = False) -> Scene:
    """Keep primary and secondary strands separate: they form an associated object."""
    sc = Scene()
    sc.strand("primary", [role("primary target concatemer")],
              label="one starting amplicon circle → primary RCA concatemer")
    sc.strand("secondary", [role("allele-specific secondary RCA concatemers")],
              label="many matched padlock circles → secondary RCA concatemers")
    if labelled:
        sc.strand("fluor", [role("hybridized allele-specific fluorescent probes")],
                  label="fluorophore probes on secondary repeats")
    sc.footer("Separate DNA strands remain associated in one superRCA object; no covalent fusion is asserted.",
              "secondary")
    return sc


def workflow() -> Workflow:
    w = Workflow(MolecularState("patient genomic DNA",
                                tuple(molecule("patient genomic DNA", tuple(GENOMIC)).rows())),
                 track_graph=False)
    w.pcr("High-fidelity target PCR: 15 cycles", PCR,
          molecule("target amplicon", AMPLICON).rows(), name="target amplicon",
          note="SuperFi PCR enriches the patient-specific target before circle formation.")
    w.react("Anneal ligation template across the amplicon ends",
            molecule("ligation-ready amplicon strand", (
                role("amplicon 5-prime end"), TARGET,
                role("amplicon 3-prime end"))).rows(),
            name="splinted amplicon ends", note=CIRCULARIZATION.statement)
    w.react("Ampligase seals the amplicon strand into a circle",
            closed_circle("target amplicon circle", AMPLICON,
                          closure_label="amplicon end-to-end ligation").rows(),
            name="single-stranded target circle",
            note="One target-bearing amplicon strand is circularized; its exact end bases are unresolved.")
    w.react("Phi29 copies the target circle from the ligation-template primer",
            molecule("primary RCA product", (
                role("target-circle complement, repeat 1"),
                role("target-circle complement, repeat 2"),
                role("target-circle complement, repeat 3"))).rows(),
            name="primary target concatemer",
            note="Three conceptual repeats are drawn. The physical product has hundreds of repeats.")
    w.react("Anneal mutant and wild-type padlocks to repeated allele sites",
            molecule("allele-probed primary concatemer", (
                role("primary target concatemer"),
                role("mutant-specific padlocks at matching copies"),
                role("wild-type-specific padlocks at matching copies"))).rows(),
            name="genotyping padlocks on target repeats",
            note="Each padlock's ligation junction discriminates the target allele.")
    w.react("Ampligase circularizes matched padlocks",
            closed_circle("allele-specific padlock circle", (
                role("allele-specific padlock ligation arm"),
                role("padlock backbone and detection tag"),
                role("opposite padlock arm")),
                closure_label="allele-specific padlock ligation").rows(),
            name="padlock circles around primary RCA product",
            note=PADLOCK_CLOSURE.statement + " Exact arms and backbone remain unassigned pending original supplementary material.")
    w.react("Added primers and Phi29 perform secondary RCA",
            cluster().rows(),
            name="superRCA object",
            note="Many secondary products remain associated with one primary product.")
    w.react("Hybridize fluorophore probes and count objects by flow cytometry",
            cluster(labelled=True).rows(),
            name="cytometry-ready superRCA object",
            note="Fluorescence identifies the allele; no sequencing library is produced.")
    return w
