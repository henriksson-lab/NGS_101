"""MALBAC constructs stated in Zong et al. Science 2012.

The accessible paper specifies the 27-nt common handle and N8 random foot, but the exact
handle bases are confined to an unavailable supplement. The schematic therefore models
those regions honestly as length-preserving placeholders rather than promoting a
secondary-source sequence to a published fact.
"""
from __future__ import annotations

from chemdraw import Construct, Segment, complement_segments
from pcr import PCRGeometry
from protocol_graph import (InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind)
from topology import HairpinTopology

HANDLE_NT = 27
RANDOM_NT = 8
PREAMPLIFICATION_CYCLES = 5
MELT_C = 94
ANNEAL_C = 0
EXTEND_C = 65
LOOP_C = 58


def seg(name: str, seq: str, **kw) -> Segment:
    return Segment(name=name, top=seq, **kw)


def handle(name: str = "common handle") -> Segment:
    return seg(name, "H" * HANDLE_NT, tag="tso", placeholder=True, inferred=True,
               note="27-nt handle length is published; exact bases are unavailable")


def random_foot(name: str = "N8") -> Segment:
    return seg(name, "N" * RANDOM_NT, placeholder=True)


def random_primer() -> Construct:
    return Construct([handle(), random_foot()], name="MALBAC random primer")


def pcr_primer() -> Construct:
    return Construct([handle()], name="MALBAC PCR primer")


def semi_amplicon(insert_nt: int = 28) -> Construct:
    return Construct([handle(), random_foot(),
                      seg("genomic copy", "X" * insert_nt, placeholder=True)],
                     name="MALBAC semi-amplicon")


def full_amplicon(insert_nt: int = 28) -> Construct:
    left = [handle("left handle"), random_foot("left N8"),
            seg("genomic copy", "X" * insert_nt, placeholder=True)]
    return Construct([*left, *complement_segments([handle("right handle")])],
                     name="MALBAC full amplicon")


def loop_topology(insert_nt: int = 28) -> HairpinTopology:
    """The end-paired full amplicon, without inventing unpublished handle bases."""
    loop_nt = RANDOM_NT + insert_nt
    return HairpinTopology(
        "looped MALBAC full amplicon",
        "N" * (HANDLE_NT + loop_nt + HANDLE_NT),
        first_stem_start=0, second_stem_start=HANDLE_NT + loop_nt,
        stem_length=HANDLE_NT, colour="var(--c-tso)")


def unresolved_sequencing_library() -> Construct:
    return Construct([
        seg("unidentified left platform arm","[PLATFORM ARM]",placeholder=True,
            inferred=True,note="downstream library-preparation kit not reported"),
        seg("MALBAC-amplified genomic insert","X"*36,placeholder=True),
        seg("unidentified right platform arm","[PLATFORM ARM]",placeholder=True,
            inferred=True,note="downstream library-preparation kit not reported")],
        name="MALBAC downstream sequencing-library boundary")


def protocol_graph() -> ProtocolGraph:
    source="Defining source: https://doi.org/10.1126/science.1229164"
    graph=ProtocolGraph("MALBAC",authority=source)
    genome=graph.input("genome","heat-denatured single-cell genomic DNA",
                       MaterialKind.BIOLOGICAL,primary=True)
    primer=graph.input("random_primer","27-nt-handle plus N8 random primer",
                       MaterialKind.OLIGO,state=random_primer(),inferred=True,
                       source=source+" (handle bases unavailable; length and N8 architecture published)")
    polymerase=graph.input("strand_displacing_polymerase","strand-displacing DNA polymerase",
                           MaterialKind.CATALYST)
    semi=graph.product("semi_amplicon",semi_amplicon().name,MaterialKind.MOLECULAR,
                       state=semi_amplicon(),inferred=True,
                       source=source+" (handle bases unavailable)")
    graph.react("first_extension","Random-prime and displace a first-generation semi-amplicon",
                StepKind.EXTENSION,
                ((genome,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),
                 (polymerase,InputRole.CATALYST)),(semi,))
    full=graph.product("full_amplicon",full_amplicon().name,MaterialKind.MOLECULAR,
                       state=full_amplicon(),inferred=True,
                       source=source+" (handle bases unavailable)")
    graph.react("second_extension","Copy a semi-amplicon to install complementary handles",
                StepKind.EXTENSION,
                ((semi,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),
                 (polymerase,InputRole.CATALYST)),(full,))
    loop_temp=graph.input("loop_temperature","58-degree intramolecular annealing condition",
                          MaterialKind.CONDITION)
    looped=graph.product("looped_amplicon","looped full amplicon withdrawn from random priming",
                         MaterialKind.MOLECULAR,state=full_amplicon(),inferred=True,
                         source=source+" (handle bases unavailable)")
    graph.react("loop","Anneal complementary handles intramolecularly after each cycle",
                StepKind.FOLDING,
                ((full,InputRole.TEMPLATE),(loop_temp,InputRole.CONDITION)),(looped,))
    pcr_primer_input=graph.input("pcr_primer","common 27-nt MALBAC PCR primer",
                                 MaterialKind.OLIGO,state=pcr_primer(),inferred=True,
                                 source=source+" (exact bases unavailable)")
    pcr_pol=graph.input("pcr_polymerase","PCR polymerase",MaterialKind.CATALYST)
    amplified=graph.product("amplified_genome","PCR-amplified MALBAC genomic DNA",
                            MaterialKind.MOLECULAR,state=full_amplicon(),inferred=True,
                            source=source+" (handle bases unavailable)")
    graph.react("pcr","Amplify full amplicons with the common handle primer",StepKind.PCR,
                ((looped,InputRole.TEMPLATE),(pcr_primer_input,InputRole.STRUCTURAL),
                 (pcr_pol,InputRole.CATALYST)),(amplified,),
                geometry=PCRGeometry.unavailable(
                    full_amplicon(),reason="The paper gives the common-primer length but not its bases, so exact landing geometry cannot be asserted.",
                    symmetric_species=True))
    library_kit=graph.input("library_kit","unspecified downstream sequencing-library kit",
                            MaterialKind.ADAPTER,inferred=True,
                            source=source+" (kit identity and adapter chemistry not reported)")
    unresolved=unresolved_sequencing_library()
    library=graph.product("library",unresolved.name,MaterialKind.LIBRARY,state=unresolved,
                          inferred=True,source=source+" (downstream kit unspecified)")
    graph.react("downstream_library_prep",
                "INFERRED — convert amplified DNA with the unreported sequencing-library kit",
                StepKind.ADAPTER_ATTACHMENT,
                ((amplified,InputRole.TEMPLATE),(library_kit,InputRole.STRUCTURAL)),(library,),
                source=source+" (the paper sequences MALBAC products but does not identify the kit)")
    graph.sequence(library,SequencingPlan(
        acquired_roles=("Read 1","Read 2"),
        read_uses={"Read 1":"include","Read 2":"include"},
        limitation="The defining paper leaves the downstream sequencing-library kit and primer sites unspecified."))
    graph.declare_reviewed_scope((
        "genome","random_primer","strand_displacing_polymerase","first_extension",
        "second_extension","loop_temperature","loop","pcr_primer","pcr_polymerase",
        "pcr","library_kit","downstream_library_prep","library"),source)
    return graph.require_publishable()

