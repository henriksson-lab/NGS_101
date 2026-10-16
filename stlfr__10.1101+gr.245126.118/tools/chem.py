"""Source-owned single-transposon / 3-prime branch stLFR declaration."""
from protocol_specs import Definition, get
from pcr import PCRGeometry
from chemdraw import Construct, Segment
from protocol_graph import (AnnotationSlot, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind)
import stlfr_chem as S

SOURCE = "Wang2019 Methods and Cheng2018 doi:10.1038/protex.2018.116, original steps61–81/workbook"
DEFINITION = Definition(
    "stLFR — single-tube long-fragment reads", "10.1101/gr.245126.118",
    "Wang et al., Genome Research (2019)",
    "Long DNA is captured on barcoded beads. This page selects the single-Transposon1T route with 3′ branch ligation used for stLFR-1–3; it does not combine that route with the two-transposon alternative.",
    "dnb", ("three-set linked-read barcode", "genomic subfragment"),
    "The original detailed protocol reports paired-end 100-base genomic reads and a 42-position barcode/common-adapter string appended to Read 1. The current author repository describes a different paired-end input convention (barcode appended to Read 2) and supports both 42- and 54-position variants; these source versions are not treated as the same run layout. The workbook separately lists two-transposon and3′-branch sequencing oligos. Exact primer binding and the physical DNB circle are not represented here.",
    "DNA",
    "The ordered workbook supplies common, transposon, branch, sequencing and alternative Illumina oligos. Alternative-table entries are not all active inputs to this selected BGI route. Bead occupancy is statistical: multiple long molecules may share a barcode; no physical rejection of ambiguous assignments is asserted. The physical circle-closing state remains a model question, not confirmed publication absence.")


def roles(name):
    return Construct([Segment(name, "[" + name + "]", placeholder=True)], name=name)


FINAL_LIBRARY = Construct([
    Segment("BGI readout flank", "[BGI readout flank]", placeholder=True),
    *S.barcode_read().segments,
    Segment("genomic subfragment", "[genomic subfragment]", placeholder=True),
    Segment("circle assembly", "[physical circle not constructed]", placeholder=True)],
    name="BGI readout roles; not a constructed physical circle")
SPEC = get("stlfr", definition=DEFINITION, final_override=FINAL_LIBRARY,
           oligos_fn=S.oligos,
           sequencing_unavailable="Published sequencing oligos are listed above; exact binding on a source-matched physical circle is not represented by this model.")


def protocol_graph():
    g = ProtocolGraph(DEFINITION.title, authority=SOURCE)
    dna = g.input("dna", "Long genomic DNA", MaterialKind.BIOLOGICAL,
                  state=roles("long genomic DNA"), primary=True)
    transposome = g.input("transposome", "Transposon1T / uracil-bearing TransposonB transposome", MaterialKind.ADAPTER)
    tn5 = g.input("tn5", "Tn5 activity in the loaded transposome", MaterialKind.CATALYST)
    tagged = g.product("tagged", "Capture-transposon-tagged long DNA", MaterialKind.MOLECULAR,
                       state=roles("capture-tagged long DNA"))
    g.react("tagment", "Insert the single capture-transposon species", StepKind.TAGMENTATION,
            ((dna, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)), (tagged,))
    bead = g.input("barcoded_bead", "Preassembled three-set barcoded capture beads", MaterialKind.SELECTION,
                   state=S.assembled_bead_oligo())
    bridge = g.input("bridge", "RNA capture bridge", MaterialKind.OLIGO,
                     state=Construct([Segment("RNA bridge", S.BRIDGE_RNA)]))
    captured = g.product("captured", "Bead-associated tagged DNA; multiple occupancy possible", MaterialKind.MOLECULAR,
                         state=roles("bead-associated tagged DNA"))
    g.react("capture", "Hybridize tagged long DNA with barcoded beads", StepKind.HYBRIDIZATION,
            ((tagged, InputRole.TEMPLATE), (bead, InputRole.SELECTION), (bridge, InputRole.STRUCTURAL)), (captured,))
    ligase = g.input("ligase", "T4 DNA ligase", MaterialKind.CATALYST)
    barcoded = g.product("barcoded", "Bead-barcode ligated DNA", MaterialKind.MOLECULAR,
                         state=roles("bead-barcode ligated DNA"))
    g.react("ligate_barcode", "Transfer bead barcode at the RNA-splinted nick", StepKind.LIGATION,
            # One existing capture complex supplies both ligating partners;
            # this does not introduce a second bead-barcode reagent.
            ((captured, InputRole.TEMPLATE), (captured, InputRole.STRUCTURAL),
             (ligase, InputRole.CATALYST)), (barcoded,),
            geometry=S.barcode_transfer_geometry())
    exo1 = g.input("exo1", "Exonuclease I", MaterialKind.CATALYST)
    exo3 = g.input("exo3", "Exonuclease III", MaterialKind.CATALYST)
    digested = g.product("adapter_digested", "Adapter-digested bead-associated DNA", MaterialKind.MOLECULAR,
                         state=roles("adapter-digested bead-associated DNA"))
    adapter_waste = g.product("adapter_waste", "Digested excess adapter oligos", MaterialKind.WASTE)
    g.react("digest_adapters", "Digest adapter oligos with Exonuclease I and III", StepKind.DIGESTION,
            ((barcoded, InputRole.TEMPLATE), (exo1, InputRole.CATALYST), (exo3, InputRole.CATALYST)),
            (digested, adapter_waste))
    g.discard(adapter_waste, "Source-reported adapter digestion and bead washing.")
    sds = g.input("sds", "SDS deproteinization condition", MaterialKind.CONDITION)
    washed = g.product("washed", "Washed DNA after Tn5 removal", MaterialKind.MOLECULAR,
                       state=roles("washed bead-associated DNA"))
    wash_waste = g.product("wash_waste", "Removed Tn5 and unbound material", MaterialKind.WASTE)
    g.react("remove_tn5", "Remove Tn5 with SDS and bead washes", StepKind.DEPROTEINIZATION,
            ((digested, InputRole.TEMPLATE), (sds, InputRole.CONDITION)), (washed, wash_waste))
    g.discard(wash_waste, "Material removed during the reported SDS/wash steps; no barcode-collision rejection is implied.")
    udg = g.input("udg", "UDG", MaterialKind.CATALYST)
    ape1 = g.input("ape1", "APE1", MaterialKind.CATALYST)
    processed = g.product("capture_digested", "Capture-oligo-digested bead-associated DNA", MaterialKind.MOLECULAR,
                          state=roles("capture-oligo-digested bead-associated DNA"))
    capture_waste = g.product("capture_waste", "Digested excess capture oligos", MaterialKind.WASTE)
    g.react("digest_capture", "Digest excess capture oligos with UDG and APE1", StepKind.DIGESTION,
            ((washed, InputRole.TEMPLATE), (udg, InputRole.CATALYST), (ape1, InputRole.CATALYST)),
            (processed, capture_waste))
    g.discard(capture_waste, "Source-reported capture-oligo digestion and washing; this selected route does not add the dual-transposon Exonuclease I treatment here.")
    ssb = g.input("ssb", "ET SSB binding protein", MaterialKind.MOLECULAR)
    bound = g.product("ssb_bound", "SSB-bound bead-associated DNA", MaterialKind.MOLECULAR,
                      state=roles("SSB-bound bead-associated DNA"))
    g.react("bind_ssb", "Bind SSB before 3′ branch ligation", StepKind.BINDING,
            ((processed, InputRole.TEMPLATE), (ssb, InputRole.STRUCTURAL)), (bound,))
    branch = g.input("branch_adapter", "Phosphorylated BranchT /3′-ddG BranchB adapter", MaterialKind.ADAPTER,
                     state=Construct([Segment("BranchT", S.BRANCH_T)]))
    branch_product = g.product("branch_ligated", "3′ branch-adapter-ligated DNA", MaterialKind.MOLECULAR,
                              state=roles("3′ branch-adapter-ligated DNA"))
    g.react("branch_ligation", "Ligate adapter 5′ phosphate to recessed genomic 3′ hydroxyl", StepKind.LIGATION,
            ((bound, InputRole.TEMPLATE), (branch, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)),
            (branch_product,), geometry=S.branch_adapter_geometry())
    pcr1 = g.input("pcr1", "PCR1; workbook common-primer variant", MaterialKind.OLIGO,
                   state=Construct([Segment("PCR1", S.PCR1)]))
    pcr2 = g.input("pcr2", "PCR2; branch-adapter side", MaterialKind.OLIGO,
                   state=Construct([Segment("PCR2", S.PCR2)]))
    polymerase = g.input("polymerase", "PfuCx PCR enzyme", MaterialKind.CATALYST)
    amplified = g.product("amplified", "Amplified co-barcoded genomic library", MaterialKind.MOLECULAR,
                          state=roles("amplified co-barcoded genomic library"))
    g.react("pcr", "Amplify with source PCR1/PCR2", StepKind.PCR,
            ((branch_product, InputRole.TEMPLATE), (pcr1, InputRole.STRUCTURAL),
             (pcr2, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (amplified,),
            geometry=PCRGeometry.unavailable(roles("3′ branch-ligated DNA"),
                reason="Both ordered primers are source-known. Their full physical branch-product binding state is not constructed here; this is a model limitation, not a publication-absence claim."))
    library = g.product("library", "Reported BGI sequencing endpoint; circle model pending", MaterialKind.LIBRARY,
                        state=FINAL_LIBRARY)
    g.react("circle_and_dnb", "Reported circularization and DNB formation", StepKind.MAGIC,
            ((amplified, InputRole.TEMPLATE),), (library,),
            reason="The paper explicitly reports circularization/DNB formation and cites Drmanac 2010. A source-matched physical circle and complete reaction geometry are not constructed here; the named predecessor/source scope remains under review. This is a model gap, not a publication-absence claim.")
    g.sequence(library, SequencingPlan(limitation=DEFINITION.sequencing))
    g.declare_reviewed_scope(("dna", "transposome", "tn5", "tagment", "barcoded_bead", "bridge", "capture",
        "ligase", "ligate_barcode", "exo1", "exo3", "digest_adapters", "sds", "remove_tn5", "udg", "ape1",
        "digest_capture", "ssb", "bind_ssb", "branch_adapter", "branch_ligation", "pcr1", "pcr2", "polymerase",
        "pcr", "circle_and_dnb", "library"), SOURCE)
    g.annotate(SourceDiscrepancy("barcoded_bead", "Published common-adapter variants",
        (SourceClaim("Original Protocol Exchange Supplementary Figure 1",
                     "/52-Bio/AAAAAAAACCCGTAGCCATGTCGTTCTG; 8A linker and a different common-primer region."),
         SourceClaim("Original Supplement_0.xlsx, Common adapter sequences, BeadCommonT/PCR1",
                     "/52-Bio/AAAAAAAAAATGTGAGCCAAGGAGTTG; 10A linker, with matching PCR1 TGTGAGCCAAGGAGTTG.")),
        "The model uses the ordered workbook variant and its matching PCR1.",
        "The workbook supplies explicit oligo order sequences and a self-consistent primer pair; the protocol also allows changing the common-primer sequence for different sequencing technologies.",
        "The figure and workbook do not explain which physical preparations used each variant. This is an unresolved source variant, not a demonstrated author error."))
    return g.require_publishable()


SPEC.protocol_graph = protocol_graph
