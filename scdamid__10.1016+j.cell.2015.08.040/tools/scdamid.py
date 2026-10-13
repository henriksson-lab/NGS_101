"""Molecular model for single-cell DamID (Kind et al., Cell 2015)."""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)

# The accessible defining paper names these oligos but does not print their sequences.
# Exact bases are from the upstream scg_lib_structs reconstruction and are therefore
# created only through secondary(), which makes inferred styling unavoidable.
ADRT = "CTAATACGACTCACTATAGGGCAGCGTGGTCGCGGCCGAGGA"
ADRB = "TCCTCGGCCG"
PREAMP_STUB = "GTGGTCGCGGCCGAGGA"
RANDOM_NT = 4
INDEX_NT = 6
RUN_ROLES = ("Read 1", "Index 1 (i7)")


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def secondary(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    """A source-bounded segment that can never render as established sequence."""
    return seg(name, top, tag, inferred=True, **kw)


def damid_adaptor_scene() -> Scene:
    paired = ADRT[-len(ADRB):]
    sc = Scene()
    sc.strand("AdRt", [secondary("5' tail", ADRT[:-len(ADRB)]),
                       secondary("paired end", paired)], label="AdRt")
    sc.anneal("AdRb", [secondary("paired end reverse complement", ADRB)],
              to="AdRt", pair=("paired end reverse complement", "paired end"),
              label="AdRb")
    return sc


def dpni_fragment(insert_nt: int = 28) -> Construct:
    """Blunt genomic fragment selected between two Dam-methylated GATC sites."""
    return Construct([
        seg("left GATC", "GATC"),
        seg("Dam-marked genomic DNA", "X" * insert_nt, placeholder=True),
        seg("right GATC", "GATC"),
    ], name="DpnI-selected genomic fragment")


def preamp_product(insert_nt: int = 28) -> Construct:
    """One strand after the single DamID pre-amplification primer acts at both ends."""
    return Construct([
        secondary("N4, left", "N" * RANDOM_NT, placeholder=True),
        secondary("DamID stub, left", PREAMP_STUB),
        secondary("restored GATC, left", "TC"),
        seg("Dam-marked genomic DNA", "X" * insert_nt, placeholder=True),
        secondary("restored GATC, right", "GA"),
        secondary("DamID stub, right", revcomp(PREAMP_STUB)),
        secondary("N4, right", "N" * RANDOM_NT, placeholder=True),
    ], name="scDamID pre-amplification product")


SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"])


def final_library(insert_nt: int = 28) -> Construct:
    """Single-index Illumina library; every outer base is secondary-source here."""
    lib = Construct([
        secondary("P5", il.P5, "p5"),
        secondary("Read 1 remainder", il.TRUSEQ_READ1[4:], "r1"),
        *preamp_product(insert_nt),
        secondary("Read 2 reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
        secondary("i7", "I" * INDEX_NT, "cbc", placeholder=True,
                  feature=feature("cell_i7", "cell_barcode", "unknown")),
        secondary("P7 reverse complement", il.P7_RC, "p7"),
    ], name="scDamID sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS, required_roles=RUN_ROLES)
    if problems:
        raise ValueError("invalid scDamID final library: " + "; ".join(problems))
    return lib


def preamp_duplex() -> Scene:
    return Scene.duplex(preamp_product().segments)


def _validate() -> None:
    # This is construction logic: the published adaptor can only be exposed if its short
    # strand actually pairs, and the final run must retain both used primer sites.
    damid_adaptor_scene().rows()
    Scene.duplex(dpni_fragment().segments).rows()
    preamp_duplex().rows()
    final_library()


_validate()


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1016/j.cell.2015.08.040"
    secondary_source = source + "; exact oligo sequence reconstructed by upstream scg_lib_structs"
    graph = ProtocolGraph("scDamID", authority=source)
    genome = graph.input("dam_marked_genome", "Dam-methylated single-cell genomic DNA",
                         MaterialKind.BIOLOGICAL, primary=True)
    dpni = graph.input("dpni", "DpnI", MaterialKind.CATALYST)
    fragments = graph.product("dpni_fragments", "DpnI-selected GATC fragments",
                              MaterialKind.MOLECULAR, state=dpni_fragment())
    graph.react("dpni_digest", "Digest Dam-methylated GATC sites with DpnI",
                StepKind.DIGESTION,
                ((genome, InputRole.TEMPLATE), (dpni, InputRole.CATALYST)), (fragments,))
    adaptor = graph.input("damid_adaptor", "DamID AdRt/AdRb adaptor",
                          MaterialKind.ADAPTER, state=damid_adaptor_scene(), inferred=True,
                          source=secondary_source)
    ligase = graph.input("ligase", "DNA ligase", MaterialKind.CATALYST)
    adapted = graph.product("adapted_fragments", "DamID-adaptor-ligated fragments",
                            MaterialKind.MOLECULAR, inferred=True, source=secondary_source)
    graph.react("adaptor_ligation", "Ligate the DamID adaptor to both blunt ends",
                StepKind.LIGATION,
                ((fragments, InputRole.TEMPLATE), (adaptor, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The defining paper establishes adaptor ligation but does not print the adaptor terminal bases."))
    preamp_primer = graph.input("preamp_primer", "single N4-tailed pre-amplification primer",
                                MaterialKind.OLIGO, inferred=True, source=secondary_source)
    polymerase = graph.input("polymerase", "pre-amplification and library PCR polymerase",
                             MaterialKind.CATALYST)
    preamp = graph.product("preamp_product", "single-primer DamID pre-amplification product",
                           MaterialKind.MOLECULAR, state=preamp_product(), inferred=True,
                           source=secondary_source)
    graph.react("preamp_pcr", "Amplify with the single N4-tailed primer", StepKind.PCR,
                ((adapted, InputRole.TEMPLATE), (preamp_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (preamp,),
                geometry=PCRGeometry.unavailable(
                    preamp_product(), symmetric_species=True,
                    reason="The exact symmetric primer is secondary-source reconstruction."))
    end_repair = graph.input("end_repair_enzymes", "end-repair enzymes",
                             MaterialKind.CATALYST)
    repaired = graph.product("repaired_product", "end-repaired DamID amplicon",
                             MaterialKind.MOLECULAR, inferred=True, source=secondary_source)
    graph.react("end_repair", "End-repair the DamID amplicon", StepKind.END_REPAIR,
                ((preamp, InputRole.TEMPLATE), (end_repair, InputRole.CATALYST)), (repaired,))
    tailing = graph.input("tailing_polymerase", "dA-tailing polymerase", MaterialKind.CATALYST)
    tailed = graph.product("da_tailed_product", "dA-tailed DamID amplicon",
                           MaterialKind.MOLECULAR, inferred=True, source=secondary_source)
    graph.react("da_tailing", "Add 3′ dA overhangs", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    indexed_adapter = graph.input("indexed_adapter", "indexed Illumina adaptor",
                                  MaterialKind.ADAPTER, inferred=True, source=secondary_source)
    adapter_ligated = graph.product("adapter_ligated", "indexed-adaptor-ligated library precursor",
                                    MaterialKind.MOLECULAR, inferred=True,
                                    source=secondary_source)
    graph.react("illumina_ligation", "Ligate the indexed Illumina adaptor", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (indexed_adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapter_ligated,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The accessible defining paper does not print the indexed adapter terminal bases."))
    p5 = graph.input("p5_primer", "P5 library PCR primer", MaterialKind.OLIGO,
                     inferred=True, source=secondary_source)
    p7 = graph.input("p7_primer", "indexed P7 library PCR primer", MaterialKind.OLIGO,
                     inferred=True, source=secondary_source)
    final = final_library()
    library = graph.product("library", final.name, MaterialKind.LIBRARY, state=final,
                            inferred=True, source=secondary_source)
    graph.react("library_pcr", "PCR-complete the indexed Illumina library", StepKind.PCR,
                ((adapter_ligated, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The completed library is secondary-source reconstruction."))
    graph.sequence(library, SequencingPlan(
        library=final, primers=SEQ_PRIMERS, acquired_roles=RUN_ROLES,
        read_lengths={"Read 1": 51, "Index 1 (i7)": INDEX_NT},
        read_uses={"Read 1": "include", "Index 1 (i7)": "demultiplex"},
        limitation="Exact outer library sequence is secondary-source reconstruction because the accessible defining paper omits its oligo table."))
    graph.declare_reviewed_scope((
        "dam_marked_genome", "dpni", "dpni_digest", "damid_adaptor", "ligase",
        "adaptor_ligation", "preamp_primer", "polymerase", "preamp_pcr",
        "end_repair_enzymes", "end_repair", "tailing_polymerase", "da_tailing",
        "indexed_adapter", "illumina_ligation", "p5_primer", "p7_primer",
        "library_pcr", "library"), source)
    return graph.require_publishable()
