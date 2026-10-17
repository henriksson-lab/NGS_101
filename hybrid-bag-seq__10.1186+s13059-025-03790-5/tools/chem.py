"""Local source and reaction model for hybrid-bag-seq."""
from chemdraw import Construct, Segment, Workflow, feature
import seqprimers as sp
from illumina import P5
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry
from protocol_model import role
from rare_cell_and_affinity import _scene, _workflow
from types import SimpleNamespace

TRANSCRIPT = role("transcript-derived insert")
GENOMIC = role("genomic insert")
BAG_DNA_CAPTURE_WRITTEN = "/5ACryd//iSp18/TGTGTTGGGTGTGTTTGG NBBBB KKKKKKKGKKKKKKKKNN"
BAG_RNA_CAPTURE_WRITTEN = "/5ACryd//iSp18/TGTGTTGGGTGTGTTTGG DDDDDNNN T(30)"
HYBRID_TSO_WRITTEN = "AAGGCTCTTGACATGAAGGATAGrGrGrG"
CUSTOM_READ1_PRIMER = "CGCTGTCCAAGGCTCTTGACACGAAGGATAG"
P5_ADAPTER3 = "AATGATACGGCGACCACCGAGATCTACACCGCTGTCCAAGGCTCTTGACACGAAGGATAG"
N701_PRIMER = "CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGG"
ROUND3_TG_PRIMER = "TGTGTTGGGTGTGTTTGG"
ROUND3_ADAPTER3_EXAMPLE = "AAGGCTCTTGACACGAAGGATAGTACCGTNNGACCGACTCGCATTACCCTAT"
BAG_BARCODE = feature("bag_barcode", "cell_barcode", "combinatorial",
                      group="bag_barcode")
BAG_VARIETAL = feature("varietal_tag", "umi", "combinatorial",
                       group="varietal_tag")

ADAPTER1_COMMON = "ACGCAGAGTCCGCTC"          # Adapter1-BC1-vt1-X 5' arm; Adapter2 3' arm
ADAPTER2_COMMON = "GACCGACTCGCATTACCCTAT"    # Adapter2 5' arm; Adapter3 3' arm
ADAPTER3_COMMON = "AAGGCTCTTGACACGAAGGATAG"  # Adapter3 5' arm
P5_EXTENSION = "CGCTGTCC"                    # P5-Adapter3 / Read1-Adapter3 spacer


# The common arms are the shared overlaps of the printed oligo families; refuse to
# assemble if a transcription change breaks those overlaps.
if (CUSTOM_READ1_PRIMER != P5_EXTENSION + ADAPTER3_COMMON
        or P5_ADAPTER3 != P5 + CUSTOM_READ1_PRIMER
        or not ROUND3_ADAPTER3_EXAMPLE.startswith(ADAPTER3_COMMON)
        or not ROUND3_ADAPTER3_EXAMPLE.endswith(ADAPTER2_COMMON)):
    raise ValueError("hybrid BAG adapter arms do not overlap as printed")


def _ident(name, n, feat, tag):
    return Segment(name, "N" * n, tag, placeholder=True, feature=feat)


def bag_final_library() -> Construct:
    """Read-1-side order derived from the author-protocol oligo polarities.

    Adapter3's 3' arm equals Adapter2's 5' arm, and Adapter2's 3' arm equals
    Adapter1's 5' arm, so the P5-Adapter3 strand reads Adapter3, BC3, Adapter2,
    BC2/vt2, Adapter1, BC1/vt1 and then CATG. The article places the genomic
    four varietal bases two bases beyond CATG. The Tn5 end is unassigned here.
    """
    return Construct([
        Segment("P5", P5, "p5"),
        Segment("P5-Adapter3 spacer", P5_EXTENSION, "r1"),
        Segment("Adapter3 common arm", ADAPTER3_COMMON, "r1"),
        _ident("BAG barcode round 3", 6, BAG_BARCODE, "cbc"),
        Segment("round-3 2N stagger (3N family also printed)", "NN", placeholder=True),
        Segment("Adapter2 common arm", ADAPTER2_COMMON, "linker"),
        _ident("BAG barcode round 2", 6, BAG_BARCODE, "cbc"),
        _ident("varietal tag round 2", 4, BAG_VARIETAL, "umi"),
        Segment("Adapter1 common arm", ADAPTER1_COMMON, "linker"),
        _ident("BAG barcode round 1", 5, BAG_BARCODE, "cbc"),
        _ident("varietal tag round 1", 4, BAG_VARIETAL, "umi"),
        Segment("CATG NlaIII site", "CATG"),
        role("two template bases after CATG", length_bp=2),
        role("template-derived varietal bases", length_bp=4, molecular_feature=BAG_VARIETAL),
        role("remaining DNA- or RNA-derived insert"),
        role("Tn5 end, N7 index and P7 (Nextera adapter not printed in hybrid sources)"),
    ], name="hybrid BAG-seq library")

def bag_workflow() -> Workflow:
    captured = _scene("hybrid BAG", role("polyacrylamide gel matrix"),
                      role("5′-Acrydite random-TG DNA primer"),
                      role("5′-Acrydite poly(dT) RNA primer"),
                      role("single nucleus"))
    copies = _scene("gel-bound copies", role("Acrydite anchor"), GENOMIC,
                    role("Acrydite anchor"), TRANSCRIPT, role("TSO-derived NlaIII site"))
    cut = _scene("NlaIII-cleaved copies", role("Acrydite anchor"),
                 role("DNA- or RNA-derived copy"), role("CATG sticky end"))
    tagged = _scene("three-round tagged molecules", role("UP1"),
                    role("BAG barcode 1 + varietal tag 1"),
                    role("BAG barcode 2 + varietal tag 2"),
                    role("BAG barcode 3"), role("DNA- or RNA-derived copy"))
    final = bag_final_library()
    return _workflow(_scene("single nucleus", role("single sorted nucleus")), (
        ("Co-encapsulate at 85 °C and polymerize the hybrid BAG", captured,
         "Both 5′-Acrydite primers copolymerize into the gel; the lower denaturation temperature preserves RNA."),
        ("Reverse-transcribe RNA and Klenow-extend genomic DNA", copies,
         "Both copy types remain covalently gel-bound; the TSO supplies a double-stranded NlaIII site to cDNA."),
        ("Treat with Exonuclease I and RNase H", _scene("enzyme-treated copies", role("Acrydite anchor"),
                                                         role("DNA- or RNA-derived copy")),
         "Author protocol steps 17–19 (37 °C, 1.5 h, then washes); the protocol does not state the molecular targets."),
        ("Digest the bound copies with NlaIII (two sequential reactions)", cut,
         "Cleavage exposes the sticky end used by the first split-pool extension."),
        ("Ligation-extension tag, NaOH denaturation, Bst-extension tag, then split PCR", tagged,
         "Round 1 Adapter1 (Klenow exo− + Quick Ligase); NaOH/Brij-35 denaturation and neutralization; round 2 Adapter2 with Bst 2.0; "
         "150 µm filtration of fused BAGs; round 3 KAPA HiFi PCR with TG primer and Adapter3. The three rounds form a 17-nt BAG barcode."),
        ("Pool, AMPure-purify and tagment the amplicons", _scene("tagmented BAG amplicons", role("tag-bearing end"),
                                                       role("random Tn5 end")),
         "Two pooled 0.8× AMPure XP cleanups; ~1.5 ng amplicon is tagmented with Amplicon Tagment enzyme (55 °C, 5 min) and neutralized."),
        ("Add P5 and the N7 sample index by PCR", _scene("hybrid BAG-seq library", *tuple(final)),
         "P5-Adapter3 and a Nextera_N70x primer (10 cycles, Nextera PCR mix), then 0.8× AMPure XP. Paired-end 150-bp sequencing: the custom Read1-Adapter3 primer reads the tag-bearing end; the other read starts at a random tagmentation end."),
    ))

SPEC = SimpleNamespace(
    TITLE="Hybrid BAG-seq — joint single-nucleus DNA and RNA sequencing",
    NOTES="01_hybrid_bag_seq.html",
    SOURCE='Defining source: <a href="https://doi.org/10.1186/s13059-025-03790-5">Li et al., Genome Biology (2025)</a>; hybrid capture-primer sequences and split-pool chemistry are printed in the <a href="https://doi.org/10.17504/protocols.io.rm7vzqr92vx1/v1">linked author protocol</a>.',
    SUMMARY="Acrydite random-TG and poly(dT) primers copolymerized into one gel capture genomic DNA and RNA-derived cDNA from the same nucleus before three rounds of split-pool tagging.",
    CAVEAT=("The tag-bearing end is assembled from the author-protocol common arms; barcodes, varietal tags and the 2N stagger stay placeholders (a 3N round-3 family is also printed). "
            "The Read 1 order is derived from oligo polarity, not copied from a figure. The Tn5/N7/P7 end is a role because the Nextera adapter is not printed in the hybrid sources. "
            "The Exonuclease I/RNase H step is shown as performed, but its molecular targets are not stated by the protocol. N701 is a family representative. The older BAG capture primers differ from this hybrid version."),
    FINAL_LIBRARY=bag_final_library(),
    SEQ_PRIMERS=(sp.custom("Read 1", "Read1-Adapter3", CUSTOM_READ1_PRIMER,
                           "protocols.io rm7vzqr92vx1 v1 Primers and step 41"),),
    READ_LENGTHS={"Read 1": 150},
    ACQUIRED_ROLES=("Read 1",),
    SHOW_READ_CONTENT=False,
    READ_SPECS=(sp.read_spec("Read 1", 150, "BAG barcode 3, stagger, Adapter2, BAG barcode 2 + varietal tag 2, Adapter1, BAG barcode 1 + varietal tag 1, CATG, then template; article trims the first 80 bases after extracting barcode/varietal tags", "extract_then_include"),
                sp.read_spec("Read 2", 150, "random Tn5 end into template; Adapter1-derived GAGCGGACTCTGCG trimmed when present (primer not placed here)", "include"),
                sp.read_spec("Index 1 (i7)", None, "Nextera N70x sample index; acquisition length not stated in inspected sources", "undetermined"),
                sp.read_spec("Index 2 (i5)", None, "no i5 index in P5-Adapter3; acquisition not described", "undetermined")),
    FINAL_CAPTION="Tag-bearing end with source-printed common arms; Read 1 begins at BAG barcode round 3. The Tn5 end is a role.",
    SEQUENCING_INTRO="The custom Read1-Adapter3 primer (author protocol step 41) ends where BAG barcode round 3 begins. Read 2, from the random tagmentation end, and index reads are not placed in this model.",
    SEQUENCING_ENDING="NovaSeq 6000 paired-end 150-bp reads were used. Read 1 begins with an 80-base tag/adapter block containing the 17-nt BAG barcode and 12-nt varietal tag.",
    SEQUENCING_UNAVAILABLE=("Read 2 and the i7 index read start in the Nextera Tn5 end, whose adapter "
                            "sequence and read primers are not printed in the hybrid sources; they are not placed here."),
    workflow=bag_workflow,
    oligos=lambda: (f"5′-{BAG_DNA_CAPTURE_WRITTEN}-3′ — Acrydite random-TG DNA capture primer",
                             f"5′-{BAG_RNA_CAPTURE_WRITTEN}-3′ — Acrydite poly(dT) RNA capture primer",
                             f"5′-{HYBRID_TSO_WRITTEN}-3′ — hybrid template-switch oligo",
                             f"5′-{CUSTOM_READ1_PRIMER}-3′ — custom Read 1 primer (Read1-Adapter3; author protocol)",
                             f"5′-{P5_ADAPTER3}-3′ — final P5-Adapter3 PCR primer (author protocol)",
                             f"5′-{N701_PRIMER}-3′ — illustrative Nextera_N701 indexed PCR primer (author workbook)"),
)

def _bag_graph(spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
    nucleus = graph.input("nucleus", "single sorted nucleus containing genomic DNA and RNA",
                          MaterialKind.BIOLOGICAL,
                          state=_scene("single nucleus", role("genomic DNA"), role("RNA")),
                          primary=True)
    gel = graph.input("gel_condition", "hybrid-BAG polymerization and co-encapsulation condition",
                      MaterialKind.CONDITION)
    bag = graph.product("bag", "one nucleus co-encapsulated in a hybrid BAG",
                        MaterialKind.MOLECULAR,
                        state=_scene("hybrid BAG", role("polyacrylamide gel matrix"),
                                     role("Acrydite random-TG primer"),
                                     role("Acrydite poly(dT) primer"), role("single nucleus")))
    excluded = graph.product("excluded", "uncaptured nuclei or empty gels", MaterialKind.WASTE)
    graph.react("encapsulate", "Co-encapsulate and polymerize the hybrid BAG",
                StepKind.PARTITION,
                ((nucleus, InputRole.TEMPLATE), (gel, InputRole.CONDITION)),
                (bag, excluded), changes_molecule=False)
    graph.discard(excluded, "does not retain a one-nucleus hybrid gel particle")
    rt_primer = graph.input("rna_capture_primer", "Acrydite poly(dT) RNA capture primer",
                            MaterialKind.OLIGO, state=BAG_RNA_CAPTURE_WRITTEN)
    tso = graph.input("template_switch_oligo", "hybrid TSO carrying an NlaIII site",
                      MaterialKind.OLIGO, state=HYBRID_TSO_WRITTEN.replace("rG", "G"))
    rtase = graph.input("reverse_transcriptase", "Maxima H− reverse transcriptase (TSO added after 15 min; 42 °C 1 h)",
                        MaterialKind.CATALYST)
    rna_copy = graph.product("rna_copy", "gel-bound transcript cDNA",
                             MaterialKind.MOLECULAR, state=bag.state)
    graph.react("reverse_transcribe", "Reverse-transcribe gel-captured RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((bag, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (tso, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (rna_copy,))
    dna_primer = graph.input("dna_capture_primer", "Acrydite random-TG DNA capture primer",
                             MaterialKind.OLIGO, state=BAG_DNA_CAPTURE_WRITTEN)
    klenow = graph.input("klenow", "DNA polymerase I Large (Klenow) fragment added to the RT tube, 42 °C 30 min", MaterialKind.CATALYST)
    copies = graph.product("copies", "gel-bound genomic and transcript-derived copies",
                           MaterialKind.MOLECULAR, state=rna_copy.state)
    graph.react("copy_dna", "Klenow-extend genomic DNA in the same BAG",
                StepKind.EXTENSION,
                ((rna_copy, InputRole.TEMPLATE), (dna_primer, InputRole.STRUCTURAL),
                 (klenow, InputRole.CATALYST)), (copies,))
    exo = graph.input("exo_rnaseh", "Exonuclease I and RNase H (37 °C, 1.5 h)", MaterialKind.CATALYST)
    treated = graph.product("exo_treated", "Exo I/RNase H-treated gel-bound copies",
                            MaterialKind.MOLECULAR, state=copies.state)
    exo_waste = graph.product("exo_waste", "digestion products and released material removed by STOP-buffer washes (targets not stated)",
                              MaterialKind.WASTE)
    graph.react("exo_rnaseh_treat", "Treat BAGs with Exonuclease I and RNase H", StepKind.DIGESTION,
                ((copies, InputRole.TEMPLATE), (exo, InputRole.CATALYST)), (treated, exo_waste))
    graph.discard(exo_waste, "BAGs are centrifuged and washed with STOP-25/STOP-10; the protocol does not name the digested species")
    nla = graph.input("nlaiii", "NlaIII restriction endonuclease (two sequential reactions)", MaterialKind.CATALYST)
    cut = graph.product("cut", "NlaIII-cleaved gel-bound DNA and cDNA copies",
                        MaterialKind.MOLECULAR, state=copies.state)
    graph.react("digest", "Digest gel-bound copies with NlaIII", StepKind.DIGESTION,
                ((treated, InputRole.TEMPLATE), (nla, InputRole.CATALYST)), (cut,))
    current = cut
    reviewed = ["nucleus", "gel_condition", "encapsulate", "rna_capture_primer",
                "template_switch_oligo",
                "reverse_transcriptase", "reverse_transcribe", "dna_capture_primer",
                "klenow", "copy_dna", "exo_rnaseh", "exo_rnaseh_treat", "nlaiii", "digest"]
    for round_number in (1, 2):
        prefix = f"round{round_number}"
        if round_number == 2:
            naoh = graph.input("naoh_denaturation", "NaOH/Brij-35 denaturation (3 washes) then Tris/NaCl/EDTA/Tween neutralization",
                               MaterialKind.CONDITION)
            single = graph.product("denatured", "gel-anchored strands after alkaline denaturation", MaterialKind.MOLECULAR,
                                   state=current.state)
            released = graph.product("released_strands", "non-covalently bound strands released by NaOH (inferred fate)", MaterialKind.WASTE)
            graph.react("denature", "Denature with NaOH and neutralize before round 2", StepKind.STRAND_SEPARATION,
                        ((current, InputRole.TEMPLATE), (naoh, InputRole.CONDITION)), (single, released))
            graph.discard(released, "Supernatants are discarded after each NaOH wash; which strands leave is inferred from the covalent Acrydite anchor, not stated")
            reviewed.extend(("naoh_denaturation", "denature"))
            current = single
        adapter = graph.input(prefix + "_adapter", f"BAG barcode/varietal-tag adapter round {round_number}",
                              MaterialKind.ADAPTER)
        enzyme = graph.input(prefix + "_enzyme",
                             "Klenow exo− and Quick Ligase" if round_number == 1 else "Bst 2.0 DNA polymerase",
                             MaterialKind.CATALYST)
        tagged = graph.product(prefix + "_product", f"BAG-tagged molecules after round {round_number}",
                               MaterialKind.MOLECULAR, state=current.state)
        step_id = prefix + ("_ligation_extension" if round_number == 1 else "_linear_extension")
        graph.react(step_id,
                    "Add first BAG tag by ligation-extension" if round_number == 1 else
                    "Add second BAG tag by hybridization and Bst extension",
                    StepKind.LIGATION if round_number == 1 else StepKind.EXTENSION,
                    ((current, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                     (enzyme, InputRole.CATALYST)), (tagged,),
                    **({"geometry": AdapterLigationGeometry.unavailable(
                        "The linked author protocol prints the oligos, but their exact ligated molecular alignment is not yet resolved in this model.")}
                       if round_number == 1 else {}))
        reviewed.extend((prefix + "_adapter", prefix + "_enzyme", step_id))
        current = tagged
    strainer = graph.input("strainer", "150 µm PluriStrainer filtration", MaterialKind.SELECTION)
    filtered = graph.product("filtered", "unfused BAGs carrying two-round tags", MaterialKind.MOLECULAR,
                             state=current.state)
    fused = graph.product("fused_bags", "fused BAGs retained by the strainer", MaterialKind.WASTE)
    graph.react("filter_bags", "Filter out fused BAGs before round 3", StepKind.SELECTION,
                ((current, InputRole.TEMPLATE), (strainer, InputRole.SELECTION)), (filtered, fused))
    graph.discard(fused, "retained on the 150 µm strainer")
    reviewed.extend(("strainer", "filter_bags"))
    current = filtered
    round3_f = graph.input("round3_forward", "TG primer (TGTGTTGGGTGTGTTT*G*G; UP1)",
                           MaterialKind.OLIGO, state=ROUND3_TG_PRIMER)
    round3_r = graph.input("round3_reverse", "Adapter3-BC3-Adapter2-Z (2N member 1 shown)",
                           MaterialKind.OLIGO, state=ROUND3_ADAPTER3_EXAMPLE)
    polymerase = graph.input("round3_polymerase", "KAPA HiFi master mix, 10 cycles",
                             MaterialKind.CATALYST)
    amplified = graph.product("amplified", "three-round BAG-tagged amplicons",
                              MaterialKind.MOLECULAR, state=current.state)
    graph.react("round3_pcr", "Add the third BAG barcode by split PCR", StepKind.PCR,
                ((current, InputRole.TEMPLATE), (round3_f, InputRole.STRUCTURAL),
                 (round3_r, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(
                    spec.FINAL_LIBRARY, reason="This model has not assembled the source-printed barcode and primer families into an exact amplicon."))
    beads = graph.input("ampure_round3", "two pooled 0.8× AMPure XP cleanups", MaterialKind.SELECTION)
    purified = graph.product("purified", "pooled, purified three-round amplicons", MaterialKind.MOLECULAR,
                             state=amplified.state)
    unbound = graph.product("ampure_unbound", "primers and short products not bound at 0.8×", MaterialKind.WASTE)
    graph.react("cleanup_round3", "Pool by column and purify twice with 0.8× AMPure XP", StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (purified, unbound))
    graph.discard(unbound, "removed with the bead supernatant")
    transposome = graph.input("transposome", "Amplicon Tagment enzyme in Nextera TD buffer (55 °C, 5 min)", MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase (Nextera XT kit)", MaterialKind.CATALYST)
    tagmented = graph.product("tagmented", "tagmented BAG-tagged amplicons", MaterialKind.MOLECULAR,
                              state=amplified.state)
    waste = graph.product("tn5_waste", "fragments lacking the tag-bearing end (not amplified by P5-Adapter3; candidate)",
                          MaterialKind.WASTE)
    graph.react("tagment", "Tagment ~1.5 ng pooled amplicons",
                StepKind.TAGMENTATION,
                ((purified, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagmented, waste))
    graph.discard(waste, "not physically removed; lacks the Adapter3 site, so not exponentially amplified by the final PCR")
    p5_primer = graph.input("p5_adapter3", "P5-Adapter3 final PCR primer", MaterialKind.OLIGO,
                            state=P5_ADAPTER3)
    n7_primer = graph.input("n7_primer", "Nextera N70x indexed final PCR primer (N701 representative)", MaterialKind.OLIGO,
                            state=N701_PRIMER)
    final_polymerase = graph.input("final_polymerase", "Nextera PCR mix, 10 cycles", MaterialKind.CATALYST)
    pcr_product = graph.product("index_pcr_product", "indexed PCR product", MaterialKind.MOLECULAR,
                                state=spec.FINAL_LIBRARY)
    graph.react("index_pcr", "Add P5 and the N7 sample index by PCR", StepKind.PCR,
                ((tagmented, InputRole.TEMPLATE), (p5_primer, InputRole.STRUCTURAL),
                 (n7_primer, InputRole.STRUCTURAL), (final_polymerase, InputRole.CATALYST)),
                (pcr_product,), geometry=PCRGeometry.unavailable(
                    spec.FINAL_LIBRARY, reason="P5-Adapter3 binds the source-printed Adapter3 arm; the N70x site lies in the Nextera Tn5 adapter, which the hybrid sources do not print."))
    final_beads = graph.input("ampure_final", "0.8× AMPure XP cleanup", MaterialKind.SELECTION)
    library = graph.product("library", "hybrid BAG-seq library", MaterialKind.LIBRARY,
                            state=spec.FINAL_LIBRARY)
    final_unbound = graph.product("final_unbound", "primers and short products not bound at 0.8×", MaterialKind.WASTE)
    graph.react("cleanup_final", "Purify the indexed library with 0.8× AMPure XP", StepKind.CLEANUP,
                ((pcr_product, InputRole.TEMPLATE), (final_beads, InputRole.SELECTION)), (library, final_unbound))
    graph.discard(final_unbound, "removed with the bead supernatant")
    graph.sequence(library, SequencingPlan(library=spec.FINAL_LIBRARY, primers=tuple(spec.SEQ_PRIMERS),
                                           acquired_roles=spec.ACQUIRED_ROLES,
                                           read_lengths=spec.READ_LENGTHS))
    reviewed.extend(("round3_forward", "round3_reverse", "round3_polymerase", "round3_pcr",
                     "ampure_round3", "cleanup_round3",
                     "transposome", "tn5", "tagment", "p5_adapter3", "n7_primer",
                     "final_polymerase", "index_pcr", "ampure_final", "cleanup_final", "library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    return _bag_graph(SPEC)

SPEC.protocol_graph = protocol_graph
