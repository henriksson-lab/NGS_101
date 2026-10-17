"""meta3C restriction-ligation library from mixed microbial populations."""
from protocol_extensions import proximity_product
from chemdraw import Row, Scene, Segment, Construct, feature
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)
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
# Endpoint geometry is not yet reconstructed from the historical kit original.
# Do not assign canonical index reads to the published inline-tag library.
FINAL_LIBRARY = Construct([role("custom PE end / inline tag / random bases"),
                           role("microbial 3C fragment"),
                           role("opposite custom PE end / inline tag / random bases")],
                          name="source-limited meta3C library")
SEQ_PRIMERS = ()
CAVEAT = "The graph selects the bacterial HpaII/MM76–MM77 route. Table 1 also lists yeast/river adapters; river chromatin uses HaeIII. Exact historical kit primer and adapter-end geometry remain unverified, so final ends are role tokens."
SEQUENCING_ENDING = "HiSeq2000 paired-end (Methods 2 × 100; Results gives 2 × 104 for the bacterial library); inline tags demultiplex libraries and six random bases on each adapter identify PCR duplicates."
SEQUENCING_UNAVAILABLE = "Historical PE-930-1001 primer sites and custom-adapter ligation geometry have not been reconstructed."
READ_CYCLES = "100 (Methods) or 104 (Results, bacterial library)"
READ_SPECS = tuple(sp.read_spec(r, READ_CYCLES,
    "Table 1 places six random bases and a short fixed library tag at the adapter end adjacent to the insert; "
    "the paper says reads begin with the tag. Exact cycle boundaries are unresolved because the kit primer site is not reconstructed.",
    "extract_then_include") for r in ("Read 1", "Read 2"))

def bacterial_adapters():
    """Table 1 bacterial ordering strands, split only where the paper assigns meaning."""
    umi = feature("adapter_n6", "umi", "random",
                  note="Methods: six Ns on each adapter collapse PCR duplicates")
    tag = feature("meta3c_library_tag", "inline_barcode", "fixed",
                  note="Methods: small tag in the custom adapters demultiplexes libraries")
    sc = Scene()
    # Independent ordering strands: this does not assert annealed geometry.
    sc.strand("MM76", [Segment("MM76 tag", MM76[:3], feature=tag),
                       Segment("MM76 N6", MM76[3:9], feature=umi),
                       Segment("MM76 PE adapter", MM76[9:])],
              label="MM76 (ordering strand)")
    sc.strand("MM77", [Segment("MM77 PE adapter", MM77[:33]),
                       Segment("MM77 N6", MM77[33:39], feature=umi),
                       Segment("MM77 tag", MM77[39:], feature=tag)],
              label="MM77 (ordering strand)")
    return sc
FINAL_CAPTION="Paired-end reads include ordinary genomic fragments and restriction-junction contact products."
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
        "restriction_enzyme", "HpaII (selected bacterial demonstration; yeast DpnII, river HaeIII)",
        MaterialKind.CATALYST)
    digested = graph.product("digested_chromatin", "restriction-digested fixed chromatin",
                             MaterialKind.MOLECULAR, state=fixed.state)
    graph.react("digest", "Digest fixed chromatin", StepKind.DIGESTION,
                ((fixed, InputRole.TEMPLATE),
                 (restriction_enzyme, InputRole.CATALYST)), (digested,))

    ligase = graph.input("proximity_ligase", "T4 DNA ligase for dilute proximity ligation",
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

    reversal = graph.input("crosslink_reversal", "proteinase K / EDTA cross-link reversal, followed by precipitation, RNase and phenol–chloroform recovery",
                           MaterialKind.CONDITION)
    purified = graph.product("purified_three_c", "de-cross-linked meta3C DNA pool",
                             MaterialKind.MOLECULAR, state=pool.state)
    graph.react("reverse_crosslinks", "Reverse cross-links and recover DNA",
                StepKind.DEPROTEINIZATION,
                ((pool, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)), (purified,))
    covaris = graph.input("covaris", "Covaris S220 shearing of 5 µg 3C DNA", MaterialKind.CONDITION)
    sheared = graph.product("sheared_three_c", "sheared meta3C DNA pool",
                            MaterialKind.MOLECULAR, state=purified.state)
    graph.react("shear", "Shear the recovered 3C pool", StepKind.FRAGMENTATION,
                ((purified, InputRole.TEMPLATE), (covaris, InputRole.CONDITION)), (sheared,))

    qiaquick = graph.input("qiaquick", "QIAquick column recovery", MaterialKind.SELECTION)
    recovered = graph.product("column_recovered", "column-purified sheared DNA",
                              MaterialKind.MOLECULAR, state=sheared.state)
    column_flow = graph.product("column_flowthrough", "QIAquick flow-through and washes",
                                MaterialKind.WASTE)
    graph.react("column_cleanup", "Recover sheared DNA on QIAquick columns", StepKind.CLEANUP,
                ((sheared, InputRole.TEMPLATE), (qiaquick, InputRole.SELECTION)),
                (recovered, column_flow))
    graph.discard(column_flow, "unbound material removed by the column")
    kit = graph.input("library_kit", "Illumina PE-930-1001 processing reagents", MaterialKind.CONDITION)
    prepared = graph.product("kit_prepared", "kit-processed DNA before custom-adapter ligation",
                             MaterialKind.MOLECULAR)
    graph.react("kit_processing", "Process DNA with PE-930-1001 before ligation", StepKind.MAGIC,
                ((recovered, InputRole.TEMPLATE), (kit, InputRole.CONDITION)), (prepared,),
                reason="Historical PE-930-1001 kit original remains uninspected; terminal preparation is not reconstructed.")
    adapter = graph.input(
        "paired_end_adapter", "bacterial custom paired-end adapter pair MM76/MM77 (Table 1)",
        MaterialKind.ADAPTER, state=bacterial_adapters(),
        source=SOURCE + " (both adapter oligos printed in Table 1; ligation-end geometry not specified)")
    adapter_ligase = graph.input("adapter_ligase", "T4 DNA ligase for custom PE adapter ligation",
                                 MaterialKind.CATALYST)
    ligated = graph.product("adapter_ligated", "adapter-ligated meta3C fragments",
                            MaterialKind.MOLECULAR, state=None)
    graph.react("ligate_adapters", "Ligate custom PE adapters (3 h RT), then 65 °C 20 min",
                StepKind.LIGATION,
                ((prepared, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (adapter_ligase, InputRole.CATALYST)), (ligated,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The defining paper names custom paired-end adapters but does not print their ligation-end geometry."))
    size_matrix = graph.input("size_selection", "PippinPrep 400–800 bp size selection",
                              MaterialKind.SELECTION)
    selected = graph.product("selected_fragments", "size-selected adapter-ligated fragments",
                             MaterialKind.MOLECULAR, state=ligated.state)
    rejected = graph.product("rejected_sizes", "out-of-range library fragments",
                             MaterialKind.WASTE)
    graph.react("size_select", "PippinPrep: retain 400–800 bp fragments",
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
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase (not named in the paper)",
                             MaterialKind.CATALYST)
    library = graph.product(
        "amplified_library", "PCR-amplified meta3C library",
        MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True,
        source=SOURCE + " (role-only endpoint; historical kit geometry unresolved)")
    graph.react("library_pcr", "Test PCRs set the cycle number; eight large-scale PCRs amplify the library", StepKind.PCR,
                ((selected, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    selected.state,
                    reason="The defining paper does not print the custom paired-end PCR primers."))
    minelute = graph.input("minelute", "Qiagen MinElute column recovery", MaterialKind.SELECTION)
    clean = graph.product("sequencing_library", "MinElute-purified meta3C sequencing library",
                          MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True)
    pcr_flow = graph.product("minelute_flowthrough", "MinElute flow-through and washes",
                             MaterialKind.WASTE)
    graph.react("final_cleanup", "Recover PCR products on MinElute columns", StepKind.CLEANUP,
                ((library, InputRole.TEMPLATE), (minelute, InputRole.SELECTION)),
                (clean, pcr_flow))
    graph.discard(pcr_flow, "unbound PCR-reaction material removed by the column")
    graph.sequence(clean, SequencingPlan(limitation=SEQUENCING_ENDING))
    graph.annotate(SourceDiscrepancy(
        "sequencing_library", "Published read lengths differ for the bacterial library",
        (SourceClaim("Results, first meta3C paragraph (PMC4381813 JATS line 129)",
                     "bacterial mix sequenced on HiSeq2000, paired end 2 × 104 bp"),
         SourceClaim("Materials and methods, Illumina sequencing (PMC4381813 JATS line 393)",
                     "libraries sequenced on HiSeq2000, PE 2 × 100")),
        "The model records both values and does not assign one acquired cycle count.",
        "Both statements are in the published eLife article (identical in the eLife v1 XML).",
        "The article does not say whether 104 includes extra cycles or the Methods value is generic; the river library reports 91 useful bases."))
    graph.declare_reviewed_scope((
        "microbial_community", "formaldehyde", "crosslink", "restriction_enzyme",
        "digest", "proximity_ligase", "proximity_ligation", "crosslink_reversal",
        "reverse_crosslinks", "covaris", "shear", "paired_end_adapter",
        "adapter_ligase", "ligate_adapters", "size_selection", "size_select",
        "pcr_forward", "pcr_reverse", "pcr_polymerase", "library_pcr",
        "sequencing_library", "qiaquick", "column_cleanup", "library_kit", "kit_processing",
        "minelute", "final_cleanup"), SOURCE)
    return graph.require_publishable()
