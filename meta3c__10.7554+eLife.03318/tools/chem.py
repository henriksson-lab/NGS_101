"""meta3C restriction-ligation library from mixed microbial populations."""
from multibranch import illumina_branch
from protocol_extensions import proximity_product
from chemdraw import Row
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role
TITLE="meta3C — metagenomic chromosome-conformation capture"
NOTES="01_meta3c.html"
SOURCE='Defining source: <a href="https://doi.org/10.7554/eLife.03318">Marbouty et al., <i>eLife</i> (2014)</a>.'
SUMMARY="Crosslinked cells from an intact microbial mixture are restriction-digested and diluted-ligated; ordinary shotgun library construction preserves both contact junctions and unligated genomic fragments without biotin enrichment."
MM70="GTANNNNNNAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG"
MM71="ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNTACT"
MM76="CAGNNNNNNAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG"
MM77="ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNCTGT"
def oligos():
    return [f"<b>{name}</b> (Table 1, 5′→3′): <code>{sequence}</code>"
            for name, sequence in (("MM70, yeast/river", MM70), ("MM71, yeast/river", MM71),
                                   ("MM76, bacteria", MM76), ("MM77, bacteria", MM77))]
FINAL_LIBRARY,SEQ_PRIMERS=illumina_branch("meta3C","microbial 3C fragment", inferred_adapters=True)
FINAL_CAPTION="Paired-end reads include ordinary genomic fragments and restriction-junction contact products."
SEQUENCING_INTRO="Declared Illumina primers bind the finished meta3C library."
def sections():
    return [("Crosslink the mixed microbial community",[Row(chunks=[("species A cells + species B cells + … → formaldehyde-fixed chromosomes", None, False)])],"Crosslinking occurs before lysis so contacts retain cellular provenance."),("Digest fixed chromatin",[Row(chunks=[("bacteria: HpaII C^CGG; yeast: DpnII G^ATC", "me", False)])],"The enzyme is chosen for the organisms and crosslinking conditions."),("Dilute and ligate proximal restriction ends",proximity_product(left="microbial locus A",right="microbial locus B",bridge="restriction junction").rows(),"Ligation produces within-cell contact chimeras; no biotin fill-in or junction pull-down is used."),("Reverse crosslinks, shear and make a paired-end library",[Row(chunks=[("3C pool → Covaris shear → custom PE adapters → size select → PCR", "r1", False)])],"Keeping non-junction fragments also supplies shotgun coverage for metagenome assembly."),]


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    sample = graph.input(
        "microbial_community", "mixed microbial community", MaterialKind.BIOLOGICAL,
        state=molecule("mixed microbial cells", (role("species-resolved chromosomes"),)),
        primary=True)
    formaldehyde = graph.input("formaldehyde", "formaldehyde cross-linking condition",
                               MaterialKind.CONDITION)
    fixed = graph.product("fixed_chromatin", "cross-linked microbial chromatin",
                          MaterialKind.MOLECULAR, state=sample.state)
    graph.react("crosslink", "Cross-link the intact microbial community",
                StepKind.CROSS_LINKING,
                ((sample, InputRole.TEMPLATE), (formaldehyde, InputRole.CONDITION)),
                (fixed,))

    restriction_enzyme = graph.input(
        "restriction_enzyme", "organism-appropriate HpaII or DpnII restriction enzyme",
        MaterialKind.CATALYST)
    digested = graph.product("digested_chromatin", "restriction-digested fixed chromatin",
                             MaterialKind.MOLECULAR, state=fixed.state)
    graph.react("digest", "Digest fixed chromatin", StepKind.DIGESTION,
                ((fixed, InputRole.TEMPLATE),
                 (restriction_enzyme, InputRole.CATALYST)), (digested,))

    ligase = graph.input("proximity_ligase", "dilute-ligation DNA ligase",
                         MaterialKind.CATALYST)
    pool = graph.product(
        "three_c_pool", "3C pool containing contact junctions and ordinary genomic fragments",
        MaterialKind.MOLECULAR,
        state=molecule("meta3C molecular pool", (
            role("within-cell restriction-junction contact products"),
            role("unligated genomic restriction fragments retained for shotgun coverage"))))
    graph.react("proximity_ligation", "Dilute-ligate proximal restriction ends",
                StepKind.PROXIMITY_LIGATION,
                ((digested, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (pool,))

    reversal = graph.input("crosslink_reversal", "cross-link reversal condition",
                           MaterialKind.CONDITION)
    purified = graph.product("purified_three_c", "de-cross-linked meta3C DNA pool",
                             MaterialKind.MOLECULAR, state=pool.state)
    graph.react("reverse_crosslinks", "Reverse cross-links and recover DNA",
                StepKind.DEPROTEINIZATION,
                ((pool, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)), (purified,))
    covaris = graph.input("covaris", "Covaris shearing condition", MaterialKind.CONDITION)
    sheared = graph.product("sheared_three_c", "sheared meta3C DNA pool",
                            MaterialKind.MOLECULAR, state=purified.state)
    graph.react("shear", "Shear the recovered 3C pool", StepKind.FRAGMENTATION,
                ((purified, InputRole.TEMPLATE), (covaris, InputRole.CONDITION)), (sheared,))

    adapter = graph.input(
        "paired_end_adapter", "custom paired-end sequencing adapter (Table 1: MM70/MM71 or MM76/MM77)",
        MaterialKind.ADAPTER,
        source=SOURCE + " (both adapter oligos printed in Table 1; ligation-end geometry not specified)")
    adapter_ligase = graph.input("adapter_ligase", "library-adapter DNA ligase",
                                 MaterialKind.CATALYST)
    ligated = graph.product("adapter_ligated", "adapter-ligated meta3C fragments",
                            MaterialKind.MOLECULAR, state=sheared.state)
    graph.react("ligate_adapters", "Ligate paired-end sequencing adapters",
                StepKind.LIGATION,
                ((sheared, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (adapter_ligase, InputRole.CATALYST)), (ligated,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The defining paper names custom paired-end adapters but does not print their ligation-end geometry."))
    size_matrix = graph.input("size_selection", "library size-selection matrix",
                              MaterialKind.SELECTION)
    selected = graph.product("selected_fragments", "size-selected adapter-ligated fragments",
                             MaterialKind.MOLECULAR, state=ligated.state)
    rejected = graph.product("rejected_sizes", "out-of-range library fragments",
                             MaterialKind.WASTE)
    graph.react("size_select", "Select the library fragment-size range",
                StepKind.SELECTION,
                ((ligated, InputRole.TEMPLATE), (size_matrix, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "falls outside the retained library-size range")

    forward = graph.input("pcr_forward", "custom paired-end PCR forward primer",
                          MaterialKind.OLIGO, inferred=True,
                          source=SOURCE + " (primer sequence not printed)")
    reverse = graph.input("pcr_reverse", "custom paired-end PCR reverse primer",
                          MaterialKind.OLIGO, inferred=True,
                          source=SOURCE + " (primer sequence not printed)")
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(
        "sequencing_library", "meta3C paired-end sequencing library",
        MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True,
        source=SOURCE + " (canonical compatible shell inferred from paired-end sequencing)")
    graph.react("library_pcr", "PCR-amplify the paired-end library", StepKind.PCR,
                ((selected, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,
                    reason="The defining paper does not print the custom paired-end PCR primers."))
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=tuple(SEQ_PRIMERS),
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "microbial_community", "formaldehyde", "crosslink", "restriction_enzyme",
        "digest", "proximity_ligase", "proximity_ligation", "crosslink_reversal",
        "reverse_crosslinks", "covaris", "shear", "paired_end_adapter",
        "adapter_ligase", "ligate_adapters", "size_selection", "size_select",
        "pcr_forward", "pcr_reverse", "pcr_polymerase", "library_pcr",
        "sequencing_library"), SOURCE)
    return graph.require_publishable()
