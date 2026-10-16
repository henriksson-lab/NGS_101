"""Local source and reaction model for hybrid-bag-seq."""
from chemdraw import Construct, Workflow, feature
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry
from protocol_model import role, role_library
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

def bag_final_library() -> Construct:
    return role_library("hybrid BAG-seq library", (
        role("Read 1 common sequence"),
        role("BAG barcode round 3", length_bp=6, molecular_feature=BAG_BARCODE),
        role("BAG barcode round 2", length_bp=6, molecular_feature=BAG_BARCODE),
        role("varietal tag round 2", length_bp=4, molecular_feature=BAG_VARIETAL),
        role("BAG barcode round 1", length_bp=5, molecular_feature=BAG_BARCODE),
        role("varietal tag round 1", length_bp=4, molecular_feature=BAG_VARIETAL),
        role("NlaIII-proximal genomic tag", length_bp=4, molecular_feature=BAG_VARIETAL),
        role("CATG NlaIII site", length_bp=4),
        role("DNA- or RNA-derived insert"),
    ))

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
        ("Digest the bound copies with NlaIII", cut,
         "Cleavage exposes the sticky end used by the first split-pool extension."),
        ("Add a ligation-extension tag, a Bst-extension tag, then a split PCR tag", tagged,
         "The three rounds form a 17-nt BAG barcode and a 12-nt varietal tag."),
        ("Pool and tagment the amplicons", _scene("tagmented BAG amplicons", role("tag-bearing end"),
                                                       role("random Tn5 end")),
         "The author protocol incubates pooled amplicons with Nextera transposase before final indexing PCR."),
        ("Add P5 and the N7 sample index by PCR", _scene("hybrid BAG-seq library", *tuple(final)),
         "P5-Adapter3 and a Nextera_N70x primer complete the library. Paired-end 150-bp sequencing starts one read at the tag-bearing end and the other at a random tagmentation end."),
    ))

SPEC = SimpleNamespace(
    TITLE="Hybrid BAG-seq — joint single-nucleus DNA and RNA sequencing",
    NOTES="01_hybrid_bag_seq.html",
    SOURCE='Defining source: <a href="https://doi.org/10.1186/s13059-025-03790-5">Li et al., Genome Biology (2025)</a>; hybrid capture-primer sequences and split-pool chemistry are printed in the <a href="https://doi.org/10.17504/protocols.io.rm7vzqr92vx1/v1">linked author protocol</a>.',
    SUMMARY="Acrydite random-TG and poly(dT) primers copolymerized into one gel capture genomic DNA and RNA-derived cDNA from the same nucleus before three rounds of split-pool tagging.",
    CAVEAT="The linked author protocol prints the hybrid capture and tagging oligos, but their complete final adapter construct has not been assembled base by base here; unresolved regions are shown by role. N701 is a family representative, not the assigned index for every sample. The older BAG capture primers differ from this hybrid version.",
    FINAL_LIBRARY=bag_final_library(),
    SEQ_PRIMERS=(),
    READ_LENGTHS=None,
    FINAL_CAPTION="Role-level paired-end library; Read 1 begins with the split-pool tags and Read 2 begins at a random tagmentation end.",
    SEQUENCING_INTRO="",
    SEQUENCING_ENDING="NovaSeq 6000 paired-end 150-bp reads were used. Read 1 begins with an 80-base tag/adapter block containing the 17-nt BAG barcode and 12-nt varietal tag.",
    SEQUENCING_UNAVAILABLE=("The author protocol prints the custom Read 1 primer (Read1-Adapter3), "
                            "but this role-level model has not yet placed its binding site or assembled "
                            "the source-printed final primer and adapter families into a complete duplex."),
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
    rtase = graph.input("reverse_transcriptase", "reverse transcriptase",
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
    klenow = graph.input("klenow", "Klenow DNA polymerase", MaterialKind.CATALYST)
    copies = graph.product("copies", "gel-bound genomic and transcript-derived copies",
                           MaterialKind.MOLECULAR, state=rna_copy.state)
    graph.react("copy_dna", "Klenow-extend genomic DNA in the same BAG",
                StepKind.EXTENSION,
                ((rna_copy, InputRole.TEMPLATE), (dna_primer, InputRole.STRUCTURAL),
                 (klenow, InputRole.CATALYST)), (copies,))
    nla = graph.input("nlaiii", "NlaIII restriction endonuclease", MaterialKind.CATALYST)
    cut = graph.product("cut", "NlaIII-cleaved gel-bound DNA and cDNA copies",
                        MaterialKind.MOLECULAR, state=copies.state)
    graph.react("digest", "Digest gel-bound copies with NlaIII", StepKind.DIGESTION,
                ((copies, InputRole.TEMPLATE), (nla, InputRole.CATALYST)), (cut,))
    current = cut
    reviewed = ["nucleus", "gel_condition", "encapsulate", "rna_capture_primer",
                "template_switch_oligo",
                "reverse_transcriptase", "reverse_transcribe", "dna_capture_primer",
                "klenow", "copy_dna", "nlaiii", "digest"]
    for round_number in (1, 2):
        prefix = f"round{round_number}"
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
    round3_f = graph.input("round3_forward", "third-round split-PCR forward primer",
                           MaterialKind.OLIGO, state=ROUND3_TG_PRIMER)
    round3_r = graph.input("round3_reverse", "third-round split-PCR reverse primer",
                           MaterialKind.OLIGO, state=ROUND3_ADAPTER3_EXAMPLE)
    polymerase = graph.input("round3_polymerase", "third-round PCR polymerase",
                             MaterialKind.CATALYST)
    amplified = graph.product("amplified", "three-round BAG-tagged amplicons",
                              MaterialKind.MOLECULAR, state=current.state)
    graph.react("round3_pcr", "Add the third BAG barcode by split PCR", StepKind.PCR,
                ((current, InputRole.TEMPLATE), (round3_f, InputRole.STRUCTURAL),
                 (round3_r, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(
                    spec.FINAL_LIBRARY, reason="This model has not assembled the source-printed barcode and primer families into an exact amplicon."))
    transposome = graph.input("transposome", "library transposome", MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagmented = graph.product("tagmented", "tagmented BAG-tagged amplicons", MaterialKind.MOLECULAR,
                              state=amplified.state)
    waste = graph.product("tn5_waste", "same-end or non-library tagmentation products",
                          MaterialKind.WASTE)
    graph.react("tagment", "Pool and tagment the BAG-tagged amplicons",
                StepKind.TAGMENTATION,
                ((amplified, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagmented, waste))
    graph.discard(waste, "does not retain the tag-bearing end in a sequenceable library")
    p5_primer = graph.input("p5_adapter3", "P5-Adapter3 final PCR primer", MaterialKind.OLIGO,
                            state=P5_ADAPTER3)
    n7_primer = graph.input("n7_primer", "Nextera N70x indexed final PCR primer (N701 representative)", MaterialKind.OLIGO,
                            state=N701_PRIMER)
    final_polymerase = graph.input("final_polymerase", "Nextera PCR mix", MaterialKind.CATALYST)
    library = graph.product("library", "hybrid BAG-seq library", MaterialKind.LIBRARY,
                            state=spec.FINAL_LIBRARY)
    graph.react("index_pcr", "Add P5 and the N7 sample index by PCR", StepKind.PCR,
                ((tagmented, InputRole.TEMPLATE), (p5_primer, InputRole.STRUCTURAL),
                 (n7_primer, InputRole.STRUCTURAL), (final_polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    spec.FINAL_LIBRARY, reason="The final primer strings are printed in the author protocol, but this role-level model does not yet place both binding sites."))
    graph.sequence(library, SequencingPlan(limitation=spec.SEQUENCING_ENDING))
    reviewed.extend(("round3_forward", "round3_reverse", "round3_polymerase", "round3_pcr",
                     "transposome", "tn5", "tagment", "p5_adapter3", "n7_primer",
                     "final_polymerase", "index_pcr", "library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    return _bag_graph(SPEC)

SPEC.protocol_graph = protocol_graph
