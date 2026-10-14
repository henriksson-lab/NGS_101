"""Molecular model for single-cell DamID (Kind et al., Cell 2015).

Provenance of the bases drawn here, none of which comes from the scg_lib_structs page
(inspiration only):

* The defining paper's own oligo figure (S1E) and indexing-primer table (S2) are PMC
  supplements behind the download gate, so the paper itself gives no sequence.  What it
  does state is that detecting Dam methylation "is similar to that of conventional DamID
  (Vogel et al., 2007) and uses largely the same reagents", and that the pre-amplification
  primer carries an "NNNN: random nucleotide tail added at the 5' end" producing a 51-nt
  read (Fig. S1E legend, both in the accessible text).
* The conventional-DamID oligos it reuses -- AdRt, AdRb and AdR_PCR -- are published, and
  are held in this repository's fetched sources as an oligo table that prints them
  verbatim and credits Vogel et al. 2007 (EpiDamID, Mol Cell 2022, supplement mmc4 and
  the article XML).  They enter the model through ``predecessor()``, which forces inferred
  styling: the reuse is stated, the identity of each oligo in *this* paper is not.
* The outer Illumina arms are canonical ``illumina`` constants, also inferred, because
  Table S2 is unavailable and the paper names no kit.
"""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)

# Conventional-DamID oligos, as published and as reprinted with their Vogel et al. 2007
# attribution in the EpiDamID oligo table.  The paper under study says it reuses these
# reagents but prints nothing itself, so every segment built from them is inferred.
ADRT = "CTAATACGACTCACTATAGGGCAGCGTGGTCGCGGCCGAGGA"   # "AdRt", 42 nt
ADRB = "TCCTCGGCCGCG"                                  # "AdRb", 12 nt
ADR_PCR = "GGTCGCGGCCGAGGATC"                          # "AdR_PCR", 17 nt
# The PCR primer spans the adaptor/GATC junction: its last two bases are the restored
# half-site a DpnI end supplies, so it extends only from a ligated fragment.
JUNCTION = ADR_PCR[-2:]
PREAMP_STUB = ADR_PCR[:-len(JUNCTION)]
RANDOM_NT = 4          # "NNNN ... added at the 5' end of the primer" (Fig. S1E legend)
INDEX_NT = 6           # drawing width; Table S2 is unavailable
RUN_ROLES = ("Read 1", "Index 1 (i7)")


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def predecessor(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    """A segment whose bases come from the reused predecessor protocol, never this paper.

    Always inferred: the defining paper states the reuse but its own oligo figure is
    unavailable, so no base here may render as established for scDamID.
    """
    return seg(name, top, tag, inferred=True, **kw)


def damid_adaptor_scene() -> Scene:
    paired = ADRT[-len(ADRB):]
    sc = Scene()
    sc.strand("AdRt", [predecessor("5' tail", ADRT[:-len(ADRB)]),
                       predecessor("paired end", paired)], label="AdRt")
    sc.anneal("AdRb", [predecessor("paired end reverse complement", ADRB)],
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
        predecessor("N4, left", "N" * RANDOM_NT, placeholder=True),
        predecessor("DamID stub, left", PREAMP_STUB),
        predecessor("restored GATC, left", JUNCTION),
        seg("Dam-marked genomic DNA", "X" * insert_nt, placeholder=True),
        predecessor("restored GATC, right", revcomp(JUNCTION)),
        predecessor("DamID stub, right", revcomp(PREAMP_STUB)),
        predecessor("N4, right", "N" * RANDOM_NT, placeholder=True),
    ], name="scDamID pre-amplification product")


SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"])


def final_library(insert_nt: int = 28) -> Construct:
    """Single-index Illumina library; no outer base is established for this paper."""
    lib = Construct([
        predecessor("P5", il.P5, "p5"),
        predecessor("Read 1 remainder", il.TRUSEQ_READ1[4:], "r1"),
        *preamp_product(insert_nt),
        predecessor("Read 2 reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
        predecessor("i7", "I" * INDEX_NT, "cbc", placeholder=True,
                    feature=feature("cell_i7", "cell_barcode", "unknown")),
        predecessor("P7 reverse complement", il.P7_RC, "p7"),
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
    if revcomp(ADRB) != ADRT[-len(ADRB):]:
        raise ValueError("AdRb is no longer the reverse complement of the AdRt 3' end")
    if ADR_PCR != ADRT[-(len(ADR_PCR) - len(JUNCTION)):] + JUNCTION:
        raise ValueError("AdR_PCR no longer spans the AdRt 3' end and the restored GATC")
    damid_adaptor_scene().rows()
    Scene.duplex(dpni_fragment().segments).rows()
    preamp_duplex().rows()
    final_library()


_validate()


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1016/j.cell.2015.08.040"
    secondary_source = (source + "; Figure S1E and Table S2 are behind the PMC download "
                        "gate, so the drawn bases are the conventional-DamID oligos the "
                        "paper says it reuses (Vogel et al., Nat Protoc 2007) and canonical "
                        "Illumina arms, all inferred for this paper")
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
                    reason="The paper's own primer figure is unavailable; the symmetric primer drawn is the reused conventional-DamID AdR_PCR with the stated NNNN tail."))
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
                    final, reason="The paper names no library-prep kit or PCR primer, and Table S2 is unavailable."))
    graph.sequence(library, SequencingPlan(
        library=final, primers=SEQ_PRIMERS, acquired_roles=RUN_ROLES,
        read_lengths={"Read 1": 51, "Index 1 (i7)": INDEX_NT},
        read_uses={"Read 1": "include", "Index 1 (i7)": "demultiplex"},
        limitation="The 51-nt single-end read is stated by the Figure S1E legend; the index width and every adapter base are inferred, because Figure S1E itself and Table S2 are behind the PMC download gate."))
    graph.declare_reviewed_scope((
        "dam_marked_genome", "dpni", "dpni_digest", "damid_adaptor", "ligase",
        "adaptor_ligation", "preamp_primer", "polymerase", "preamp_pcr",
        "end_repair_enzymes", "end_repair", "tailing_polymerase", "da_tailing",
        "indexed_adapter", "illumina_ligation", "p5_primer", "p7_primer",
        "library_pcr", "library"), source)
    return graph.require_publishable()
