"""Ramani et al. 2017 sci-Hi-C chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, MolecularState, Row, Scene, Segment, Workflow, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from restriction import DPNII, Digest
from single_cell_hic import ContactWorkflow
WORKFLOW = ContactWorkflow("sci-Hi-C", DPNII, "combinatorial indexing", "bridge-adaptor", "PCR", barcoding_rounds=2)
TITLE = "sci-Hi-C — two-level combinatorial indexing of contacts"
NOTES = "01_sci-hi-c.html"
CITATION = 'Defining source: <a href="https://doi.org/10.1038/nmeth.4155">Ramani et al., <i>Nature Methods</i> (2017)</a>.'
SUMMARY = "A first well barcode enters through indexed bridge adaptors before proximity ligation; a second well barcode enters through indexed library adaptors after redistribution."
LIBRARY_CAVEAT = ("<b>Published processing-order conflict.</b> Ramani et al. 2017 main-text overview says pooled dual-barcoded material is captured with streptavidin beads and then restriction-digested. Supplementary Methods steps 8–12 instead place MseI + HaeIII digestion, end repair and dA-tailing in the second-round wells <i>before</i> indexed Y-adaptor ligation, pooling and biotin pull-down. This page's workflow follows the detailed Supplementary Methods; the overview order is not silently treated as equivalent.")
STEPS = (("Digest in 96 wells", "Distribute fixed nuclei and digest chromatin with DpnII."),
         ("Barcode round 1", "Ligate biotinylated indexed bridge adaptors to restriction ends in each first-round well."),
         ("Pool, ligate, redistribute", "Pool nuclei, proximity-ligate bridged ends, then sort a small number of nuclei into second-round wells."),
         ("Barcode round 2 and PCR", "Reverse-crosslink, digest with MseI/HaeIII, end-repair and dA-tail in the second wells; add indexed Illumina Y adaptors, enrich biotin-bearing molecules and PCR-amplify."))
READOUT = "The first eight bases of each genomic read carry BC2; the two BC1 copies flank the internal bridge. The barcode pair assigns the contact to a nucleus."
JUNCTION_CAPTION = "The contact junction contains the first-round indexed bridge (BC1); ** marks each ligation boundary."
LIBRARY_CAPTION = "BC1 remains internal. BC2 is an 8-bp inline barcode at both read starts, followed by the source-printed CGT fixed bases; sci-Hi-C does not use an index read for BC2."
JUNCTIONS = (("locus A", "DpnII end A", "restriction end"),
             ("DpnII end A", "round-1 barcode A", "adapter ligation"),
             ("round-1 barcode B", "DpnII end B", "adapter ligation"))
BC1_EXAMPLE = "GAGCTCGA"  # A1 bridge-adaptor barcode, Supplementary Spreadsheet 1
BRIDGE_CORE = BC1_EXAMPLE + "G"  # pairs to MKVII_R_1; final G begins the EcoRI site
BRIDGE_ECORI_OVERHANG = "AATT"
BRIDGE_F_A1 = "GATC" + BRIDGE_CORE + BRIDGE_ECORI_OVERHANG
BRIDGE_R_A1 = revcomp(BRIDGE_CORE)

def contact_product():
    return Construct([Segment("locus A", "X"*22, placeholder=True), Segment("DpnII end A", "GATC", "me"),
        Segment("round-1 barcode A", BC1_EXAMPLE, "cbc",
                feature=feature("cell_bc1_a", "cell_barcode", "combinatorial", group="cell_id", part="round 1, end A")), Segment("bridge EcoRI site", "GAATTC", "me"),
        Segment("round-1 barcode B", revcomp(BC1_EXAMPLE), "cbc",
                feature=feature("cell_bc1_b", "cell_barcode", "combinatorial", group="cell_id", part="round 1, end B")), Segment("DpnII end B", "GATC", "me"),
        Segment("locus B", "X"*22, placeholder=True)], name="sci-Hi-C bridged contact")


def bridge_adapter_geometry():
    """MKVII_F_1/MKVII_R_1: DpnII-compatible end, BC1 and EcoRI overhang."""
    sc = Scene()
    sc.strand("MKVII_F_1", [Segment("5′ GATC", "GATC", "me"),
                             Segment("BC1 + EcoRI G", BRIDGE_CORE, "cbc"),
                             Segment("3′ AATT", BRIDGE_ECORI_OVERHANG, "me")],
              mod3="internal biotin-dT in AATT")
    sc.anneal("MKVII_R_1", [Segment("paired reverse", revcomp(BRIDGE_CORE), "cbc")],
              to="MKVII_F_1", pair=("paired reverse", "BC1 + EcoRI G"))
    sc.ligation_end("MKVII_F_1", "5′ GATC", kind="cohesive", polarity="5-prime",
                    owner="adapter", text="5′-GATC DpnII-compatible end")
    sc.overhang("MKVII_F_1", "3′ AATT", polarity="3-prime",
                text="AATT proximity-ligation overhang")
    return AdapterLigationGeometry.known(
        "cohesive", sc,
        statement="Ramani et al. Supplementary Spreadsheet 1 prints MKVII_F_1 and MKVII_R_1, including the 5′ GATC end, eight-base barcode and internal biotin-dT.")


def y_adapter_geometry():
    """SeqShort_BC1 HI/LO: forked TruSeq adapter with a one-base 3′ dT end."""
    hi = il.TRUSEQ_READ1 + BC2_EXAMPLE + "CGT"
    # The printed LO arm joins ``...AAG`` directly to ``ATCG...``: it omits the first G
    # of the standalone 34-nt Index-1 primer constant.
    lo = "CGCAAGGTCAAG" + il.INDEX1_PRIMER[1:]
    paired_hi = hi[-24:-1]
    paired_lo = lo[:23]
    if paired_hi != revcomp(paired_lo):
        raise ValueError("sci-Hi-C Y-adapter HI/LO core does not pair")
    sc = Scene()
    sc.strand("SeqShort_BC1_HI", [Segment("unpaired Read-1 arm", hi[:-24], "r1"),
                                  Segment("paired HI core", paired_hi, "cbc"),
                                  Segment("3′ dT", hi[-1], "me")], mod3="phosphorothioate")
    sc.anneal("SeqShort_BC1_LO paired end", [Segment("paired LO core", paired_lo, "cbc")],
              to="SeqShort_BC1_HI", pair=("paired LO core", "paired HI core"), mod5="phosphate")
    sc.ligation_end("SeqShort_BC1_HI", "3′ dT", kind="cohesive", polarity="3-prime",
                    owner="adapter", text="3′-dT pairs with a dA-tailed insert")
    lo_sc = Scene()
    lo_sc.strand("complete SeqShort_BC1_LO",
                 [Segment("paired LO core", paired_lo, "cbc"),
                  Segment("unpaired Index-1 arm", lo[23:], "r2")], mod5="phosphate")
    lo_sc.note("complete SeqShort_BC1_LO",
               "The forked unpaired arm is shown separately so its bases remain readable.")
    return AdapterLigationGeometry.known(
        "cohesive", sc, lo_sc,
        statement="Ramani et al. Supplementary Spreadsheet 1 prints the complete SeqShort_BC1 HI/LO pair and its modifications; BC1 is representative of 96 Y adapters.")

def adapter_workflow():
    start=Scene.duplex([Segment("DpnII fragment", "X"*20+"GATC", placeholder=True)],label="DpnII-digested chromatin")
    w=Workflow(MolecularState("DpnII restriction ends",tuple(start.rows())))
    bridge=bridge_adapter_geometry()
    w.ligate_adapter("Ligate first-round bridge adapters",bridge,Scene.duplex(list(contact_product()),label="bridged contact").rows(),name="bridge-adapter contact",note="BC1 remains internal after proximity ligation.")
    yend=y_adapter_geometry()
    w.ligate_adapter("Ligate second-round indexed Y adapters",yend,Scene.duplex(list(final_library()),label="sci-Hi-C library").rows(),name="BC2 library",note="BC2 begins each genomic read.")
    return w
BC2_EXAMPLE = "TGACCTTG"  # A1 adapter, Supplementary Data oligo table
Y_HI_A1 = il.TRUSEQ_READ1 + BC2_EXAMPLE + "CGT"
Y_LO_A1 = "CGCAAGGTCAAG" + il.INDEX1_PRIMER[1:]
SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["R2"])
REQUIRED_ROLES = ("Read 1", "Read 2")
def final_library():
    read_prefix = BC2_EXAMPLE + "CGT"
    lib = Construct([
        Segment("P5", il.P5, "p5"),
        Segment("Read 1 arm", il.TRUSEQ_READ1[4:], "r1"),
        Segment("BC2 A1", BC2_EXAMPLE, "cbc",
                feature=feature("cell_bc2_a", "cell_barcode", "combinatorial", group="cell_id", part="round 2, end A")), Segment("fixed read prefix", "CGT", "me"),
        *list(contact_product()),
        Segment("opposite fixed prefix and BC2", revcomp(read_prefix), "cbc",
                feature=feature("cell_bc2_b_with_prefix", "cell_barcode", "combinatorial", group="cell_id", part="round 2, end B", note="segment also contains the fixed CGT prefix")),
        Segment("Read 2 site", revcomp(il.TRUSEQ_READ2), "r2"),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="sci-Hi-C two-barcode library")
    errors = sp.verify(lib, SEQ_PRIMERS, required_roles=REQUIRED_ROLES)
    if errors: raise ValueError("invalid sci-Hi-C library: " + "; ".join(errors))
    return lib


def protocol_graph() -> ProtocolGraph:
    """The two published barcode rounds and their retained contact path."""
    source = "Ramani et al. 2017; doi:10.1038/nmeth.4155"
    graph = ProtocolGraph(TITLE, authority=source)
    nuclei = graph.input("nuclei", "fixed nuclei containing chromosomal DNA",
                         MaterialKind.BIOLOGICAL, primary=True)
    dpnii = graph.input("dpnii", "DpnII restriction enzyme", MaterialKind.CATALYST)
    digested = graph.product("digested", "DpnII-digested nuclei",
                             MaterialKind.MOLECULAR,
                             chemical_state=Digest(DPNII, "X" * 18, "X" * 18))
    graph.react("digest", "Digest chromatin with DpnII", StepKind.DIGESTION,
                ((nuclei, InputRole.TEMPLATE), (dpnii, InputRole.CATALYST)),
                (digested,))

    first_plate = graph.input("first_plate", "first 96-well barcode plate",
                              MaterialKind.CONDITION)
    round1 = graph.product("round1_a1", "DpnII fragments in representative A1 well",
                           MaterialKind.MOLECULAR, state=digested.state)
    other_round1 = graph.product("round1_other", "fragments in the other 95 first wells",
                                 MaterialKind.MOLECULAR, state=digested.state)
    graph.react("partition_round1", "Distribute intact nuclei into first-round wells",
                StepKind.PARTITION,
                ((digested, InputRole.TEMPLATE), (first_plate, InputRole.CONDITION)),
                (round1, other_round1), changes_molecule=False)
    bridge = graph.input("bridge_adapter", "A1 indexed biotinylated bridge duplex",
                         MaterialKind.ADAPTER, state=(BRIDGE_F_A1, BRIDGE_R_A1),
                         source=source + " Supplementary Data, MKVII_F_1/MKVII_R_1")
    bridge_ligase = graph.input("bridge_ligase", "bridge-adaptor DNA ligase",
                                MaterialKind.CATALYST)
    bridged = graph.product("bridged", "first-round indexed restriction fragments",
                            MaterialKind.MOLECULAR)
    graph.react("ligate_bridge", "Ligate biotinylated BC1 bridge adaptors",
                StepKind.LIGATION,
                ((round1, InputRole.TEMPLATE), (bridge, InputRole.STRUCTURAL),
                 (bridge_ligase, InputRole.CATALYST)), (bridged,),
                geometry=bridge_adapter_geometry())
    other_bridges = graph.input("other_bridges", "the other 95 indexed bridge duplexes",
                                MaterialKind.ADAPTER, inferred=True,
                                source=source + " Supplementary Data lists all 96 members")
    other_bridged = graph.product("other_bridged", "other first-round indexed fragments",
                                  MaterialKind.MOLECULAR)
    graph.react("ligate_other_bridges", "Ligate the other BC1 bridge adaptors",
                StepKind.LIGATION,
                ((other_round1, InputRole.TEMPLATE),
                 (other_bridges, InputRole.STRUCTURAL),
                 (bridge_ligase, InputRole.CATALYST)), (other_bridged,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The A1 duplex is drawn as the representative; the other 95 indexed duplexes are not assembled in this model."))

    pool_condition = graph.input("pool_condition", "pool intact first-round nuclei",
                                 MaterialKind.CONDITION)
    pooled = graph.product("pooled_round1", "pooled BC1-tagged nuclei",
                           MaterialKind.MOLECULAR, state=bridged.state)
    graph.react("pool_round1", "Pool first-round wells before contact ligation",
                StepKind.POOL,
                ((bridged, InputRole.TEMPLATE),
                 (other_bridged, InputRole.TEMPLATE),
                 (pool_condition, InputRole.CONDITION)), (pooled,),
                changes_molecule=False)
    contact_ligase = graph.input("contact_ligase", "proximity DNA ligase",
                                 MaterialKind.CATALYST)
    contacts = graph.product("contacts", "BC1-flanked proximity-ligation products",
                             MaterialKind.MOLECULAR, state=contact_product())
    graph.react("proximity_ligation", "Ligate spatially adjacent bridged fragments",
                StepKind.PROXIMITY_LIGATION,
                ((pooled, InputRole.TEMPLATE),
                 (contact_ligase, InputRole.CATALYST)), (contacts,))

    second_plate = graph.input("second_plate", "dilute into second 96-well plate",
                               MaterialKind.CONDITION)
    round2 = graph.product("round2_wells", "contact molecules in second-round wells",
                           MaterialKind.MOLECULAR, state=contacts.state)
    other_round2 = graph.product("round2_other", "contacts in the other second-round wells",
                                 MaterialKind.MOLECULAR, state=contacts.state)
    graph.react("partition_round2", "Redistribute at most 25 nuclei per second well",
                StepKind.PARTITION,
                ((contacts, InputRole.TEMPLATE), (second_plate, InputRole.CONDITION)),
                (round2, other_round2), changes_molecule=False)
    crosslink_reversal = graph.input("crosslink_reversal", "SDS and proteinase K reverse-crosslinking mix",
                                     MaterialKind.CONDITION)
    ampure = graph.input("ampure", "AMPure bead cleanup", MaterialKind.SELECTION)
    msei = graph.input("msei", "MseI four-cutter", MaterialKind.CATALYST)
    haeiii = graph.input("haeiii", "HaeIII four-cutter", MaterialKind.CATALYST)
    repair = graph.input("end_repair_klenow", "Klenow end-repair polymerase",
                         MaterialKind.CATALYST)
    datp = graph.input("datp", "dATP tailing substrate", MaterialKind.CONDITION)
    exo_minus = graph.input("tailing_klenow_exo_minus", "Klenow exo- dA-tailing polymerase",
                            MaterialKind.CATALYST)
    tailed_wells = []
    for suffix, well in (("a1", round2), ("other", other_round2)):
        deproteinized = graph.product("deproteinized_" + suffix,
                                      "reverse-crosslinked contact DNA in " + suffix,
                                      MaterialKind.MOLECULAR)
        graph.react("reverse_crosslink_" + suffix, "Reverse crosslinks after second-well dilution",
                    StepKind.DEPROTEINIZATION,
                    ((well, InputRole.TEMPLATE),
                     (crosslink_reversal, InputRole.CONDITION)), (deproteinized,))
        purified = graph.product("purified_" + suffix, "AMPure-cleaned contact DNA in " + suffix,
                                 MaterialKind.MOLECULAR)
        cleanup_waste = graph.product("cleanup_waste_" + suffix,
                                      "unretained cleanup material in " + suffix,
                                      MaterialKind.WASTE)
        graph.react("ampure_cleanup_" + suffix, "Clean contact DNA on AMPure beads",
                    StepKind.CLEANUP,
                    ((deproteinized, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)),
                    (purified, cleanup_waste))
        graph.discard(cleanup_waste, "unretained material leaves the well")
        digested_again = graph.product("four_cutter_product_" + suffix,
                                       "MseI/HaeIII-digested contact DNA in " + suffix,
                                       MaterialKind.MOLECULAR)
        graph.react("four_cutter_digest_" + suffix,
                    "Digest second-well contact DNA with MseI and HaeIII",
                    StepKind.DIGESTION,
                    ((purified, InputRole.TEMPLATE), (msei, InputRole.CATALYST),
                     (haeiii, InputRole.CATALYST)), (digested_again,))
        repaired = graph.product("repaired_" + suffix, "end-repaired contact DNA in " + suffix,
                                 MaterialKind.MOLECULAR)
        graph.react("end_repair_" + suffix, "End-repair four-cutter fragments",
                    StepKind.END_REPAIR,
                    ((digested_again, InputRole.TEMPLATE),
                     (repair, InputRole.CATALYST)), (repaired,))
        tailed = graph.product("da_tailed_" + suffix, "dA-tailed contact DNA in " + suffix,
                               MaterialKind.MOLECULAR)
        graph.react("da_tail_" + suffix, "Add dA to repaired contact ends",
                    StepKind.TAILING,
                    ((repaired, InputRole.TEMPLATE), (datp, InputRole.CONDITION),
                     (exo_minus, InputRole.CATALYST)), (tailed,))
        tailed_wells.append(tailed)
    y_adapter = graph.input("y_adapter", "A1 barcoded SeqShort HI/LO Y adapter",
                            MaterialKind.ADAPTER, state=(Y_HI_A1, Y_LO_A1),
                            source=source + " Supplementary Data, SeqShort_BC1 HI/LO")
    y_ligase = graph.input("y_ligase", "Y-adaptor DNA ligase",
                           MaterialKind.CATALYST)
    adapter_ligated = graph.product("adapter_ligated", "BC2-tagged contact DNA",
                                    MaterialKind.MOLECULAR)
    graph.react("ligate_y_adapter", "Ligate second-round indexed Y adaptors",
                StepKind.LIGATION,
                ((tailed_wells[0], InputRole.TEMPLATE), (y_adapter, InputRole.STRUCTURAL),
                 (y_ligase, InputRole.CATALYST)), (adapter_ligated,),
                geometry=y_adapter_geometry())
    other_y_adapters = graph.input("other_y_adapters", "the other 95 indexed Y duplexes",
                                   MaterialKind.ADAPTER, inferred=True,
                                   source=source + " Supplementary Data lists all 96 members")
    other_ligated = graph.product("other_adapter_ligated", "other BC2-tagged contacts",
                                  MaterialKind.MOLECULAR)
    graph.react("ligate_other_y_adapters", "Ligate the other BC2 Y adaptors",
                StepKind.LIGATION,
                ((tailed_wells[1], InputRole.TEMPLATE),
                 (other_y_adapters, InputRole.STRUCTURAL),
                 (y_ligase, InputRole.CATALYST)), (other_ligated,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The A1 Y duplex is drawn as the representative; the other 95 indexed duplexes are not assembled in this model."))

    pool2_condition = graph.input("pool2_condition", "pool second-round wells",
                                  MaterialKind.CONDITION)
    pooled2 = graph.product("pooled_round2", "pooled dual-barcoded contact DNA",
                            MaterialKind.MOLECULAR)
    graph.react("pool_round2", "Pool second-round wells", StepKind.POOL,
                ((adapter_ligated, InputRole.TEMPLATE),
                 (other_ligated, InputRole.TEMPLATE),
                 (pool2_condition, InputRole.CONDITION)), (pooled2,),
                changes_molecule=False)
    beads = graph.input("streptavidin", "streptavidin capture beads",
                        MaterialKind.SELECTION)
    captured = graph.product("captured", "biotin-enriched contact DNA",
                             MaterialKind.MOLECULAR)
    rejected = graph.product("uncaptured", "unbound non-contact fragments",
                             MaterialKind.WASTE)
    graph.react("capture", "Capture biotinylated bridge junctions",
                StepKind.AFFINITY_CAPTURE,
                ((pooled2, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "lacks a retained biotinylated contact junction")

    forward = graph.input("pcr_forward", "library-PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input("pcr_reverse", "library-PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "library-PCR polymerase",
                             MaterialKind.CATALYST)
    final = final_library()
    library = graph.product("library", final.name, MaterialKind.LIBRARY,
                            state=final)
    graph.react("library_pcr", "Amplify the two-barcode contact library",
                StepKind.PCR,
                ((captured, InputRole.TEMPLATE),
                 (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The post-capture pre-PCR contact ends and primer pair are not printed as one exact molecule."))
    graph.sequence(library, SequencingPlan(
        library=final, primers=SEQ_PRIMERS, acquired_roles=REQUIRED_ROLES))
    graph.declare_reviewed_scope((
        "nuclei", "dpnii", "digest", "first_plate", "partition_round1",
        "bridge_adapter", "bridge_ligase", "ligate_bridge", "other_bridges",
        "ligate_other_bridges", "pool_condition",
        "pool_round1", "contact_ligase", "proximity_ligation", "second_plate",
        "partition_round2", "y_adapter", "y_ligase", "ligate_y_adapter",
        "other_y_adapters", "ligate_other_y_adapters",
        "pool2_condition", "pool_round2", "streptavidin", "capture",
        "crosslink_reversal", "ampure", "msei", "haeiii", "end_repair_klenow",
        "datp", "tailing_klenow_exo_minus",
        *(name + suffix for suffix in ("a1", "other") for name in (
            "reverse_crosslink_", "ampure_cleanup_", "four_cutter_digest_",
            "end_repair_", "da_tail_")),
        "pcr_forward",
        "pcr_reverse", "pcr_polymerase", "library_pcr", "library",
    ), source)
    return graph.require_publishable()
