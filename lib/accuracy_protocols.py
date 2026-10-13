"""Source-backed models for amplification and consensus-accuracy protocols.

The helpers keep circular consensus, endogenous endpoint families, and explicit UMI
families distinct.  They are not interchangeable merely because all reduce errors.
"""
from __future__ import annotations

from types import SimpleNamespace

import illumina as il
import nextera as nx
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from batch_ngs import nextera_library, seg, truseq_library
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature, revcomp
from circular import CircularizationGeometry
from endprep import (ThreePrimeOverhangAdapter, ligation_product_scene,
                     react_dA_tail_then_ligate,
                     repair_and_dA_tail)
from protocol_model import closed_circle, molecule, role
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)


# Source-transcribed facts used by the two protocol models below.  Keeping these
# separate from drawing geometry lets the focused protocol suites check the paper,
# while construction functions continue to own molecular invariants.
TRUEPRIME_PRIMER_LENGTH_NT = (7, 9)
TRUEPRIME_AMPLIFICATION_C = 30
TRUEPRIME_INACTIVATION_C = 65
BOTSEQS_TYPICAL_DILUTION_FOLD = 100_000
BOTSEQS_INSERT_BP = 350
BOTSEQS_READ_BP = 100
TRUSEQ_Y_END = ThreePrimeOverhangAdapter(
    "TruSeq Y adapter", il.STEM + "T", il.STEM_COMPLEMENT)
TRUSEQ_DA_LIGATION = reported_adapter_end(
    "TruSeq Y adapter", kind="cohesive", polarity="3-prime", terminal_sequence="T",
    statement="The documented TruSeq Y adapter presents a single 3-prime dT to a dA-tailed insert.")


def _workflow(initial: Scene, steps) -> Workflow:
    out = Workflow(MolecularState("Starting material", tuple(initial.rows())))
    for action, state, note in steps:
        out.react(action, state.rows(), name=action, note=note)
    return out


def _trueprime_workflow() -> Workflow:
    template = molecule("denatured genomic DNA", (role("single-stranded genomic template"),))
    primed = molecule("de novo primed template", (
        role("single-stranded genomic template"), role("7–9-nt DNA primer")))
    product = molecule("branched multiple-displacement product", (
        role("long genomic copies"), role("displaced single strands"), role("new DNA primers")))
    input_dna = molecule("genomic DNA", (role("genomic DNA"),))
    return _workflow(input_dna, (
        ("Lyse the cell, denature DNA and neutralize", template,
         "The single-cell workflow lyses on ice before neutralization; purified DNA is alkaline-denatured."),
        ("TthPrimPol synthesizes short DNA primers directly from dNTPs", primed,
         "The principal primer products are 7–9 nt; no synthetic random primer is added."),
        ("phi29 DNA polymerase extends the primers and displaces downstream strands", product,
         "Displacement exposes new single-stranded template for another priming round."),
        ("Repeat priming and strand-displacement synthesis at 30 °C", product,
         "The paper used 3 h for WGA and 3 or 6 h for single cells, then 10 min at 65 °C."),
    ))


def _trueprime_graph(final: Construct, authority: str) -> ProtocolGraph:
    graph = ProtocolGraph("TruePrime — primase-driven whole-genome amplification",
                          authority=authority)
    dna = graph.input("genomic_dna", "genomic DNA", MaterialKind.BIOLOGICAL,
                      state=molecule("genomic DNA", (role("genomic DNA"),)), primary=True)
    denaturation = graph.input("denaturation_condition", "alkaline denaturation and neutralization",
                               MaterialKind.CONDITION)
    single = graph.product("single_stranded", "denatured single-stranded genomic DNA",
                           MaterialKind.MOLECULAR,
                           state=molecule("denatured genomic DNA", (role("single-stranded genomic template"),)))
    graph.react("denature", "Denature genomic DNA and neutralize",
                StepKind.THERMAL_TREATMENT,
                ((dna, InputRole.TEMPLATE), (denaturation, InputRole.CONDITION)), (single,))
    primpol = graph.input("tth_primpol", "TthPrimPol", MaterialKind.CATALYST)
    primed = graph.product("primed", "templates bearing de-novo 7–9-nt DNA primers",
                           MaterialKind.MOLECULAR, state=single.state)
    graph.react("prime", "Synthesize short DNA primers directly from dNTPs",
                StepKind.EXTENSION,
                ((single, InputRole.TEMPLATE), (primpol, InputRole.CATALYST)), (primed,))
    phi29 = graph.input("phi29", "phi29 DNA polymerase", MaterialKind.CATALYST)
    amplified = graph.product("amplified", "branched multiple-displacement WGA product",
                              MaterialKind.MOLECULAR, state=final)
    graph.react("amplify", "Extend primers with repeated strand displacement",
                StepKind.EXTENSION,
                ((primed, InputRole.TEMPLATE), (phi29, InputRole.CATALYST)), (amplified,))
    fragmentation = graph.input("fragmentation_condition", "downstream fragmentation condition",
                                MaterialKind.CONDITION, inferred=True,
                                source=authority + " (downstream NEBNext implementation unspecified)")
    fragments = graph.product("fragments", "fragmented TruePrime-amplified DNA",
                              MaterialKind.MOLECULAR, state=amplified.state, inferred=True,
                              source=authority + " (downstream NEBNext implementation unspecified)")
    graph.react("fragment", "Fragment the amplified DNA for library preparation",
                StepKind.FRAGMENTATION,
                ((amplified, InputRole.TEMPLATE), (fragmentation, InputRole.CONDITION)),
                (fragments,), source=authority + " (downstream NEBNext kit unspecified)")
    adapter = graph.input("nebnext_adapter", "unspecified NEBNext sequencing adapter",
                          MaterialKind.ADAPTER, inferred=True,
                          source=authority + " (kit generation and oligos not identified)")
    ligase = graph.input("adapter_ligase", "NEBNext adapter ligase",
                         MaterialKind.CATALYST, inferred=True,
                         source=authority + " (kit generation not identified)")
    library = graph.product("sequencing_library", "unspecified NEBNext library",
                            MaterialKind.LIBRARY, inferred=True,
                            source=authority + " (adapter architecture unavailable)")
    graph.react("library_prep", "Prepare the unspecified NEBNext sequencing library",
                StepKind.LIGATION,
                ((fragments, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (library,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The paper does not identify the NEBNext kit generation or adapter-end geometry."),
                source=authority + " (downstream kit named only)")
    graph.sequence(library, SequencingPlan(
        acquired_roles=("Read 1", "Read 2"),
        limitation="The paper reports paired-end HiSeq 2500 sequencing but not the NEBNext kit generation, adapter oligos or primer sites."))
    graph.declare_reviewed_scope([
        "genomic_dna", "denaturation_condition", "denature", "tth_primpol", "prime",
        "phi29", "amplify", "fragmentation_condition", "fragment", "nebnext_adapter",
        "adapter_ligase", "library_prep", "sequencing_library"], authority)
    return graph.require_publishable()


def _botseqs_workflow(library: Construct) -> Workflow:
    dna = molecule("sheared genomic DNA", (role("genomic fragments with native endpoints"),))
    repaired = Construct([role("genomic insert")], name="end-repaired genomic insert")
    prepared = repair_and_dA_tail(repaired)
    ligated = ligation_product_scene(library, "genomic insert",
                                     label="TruSeq PCR-free molecules",
                                     adapter_group="Illumina Y adapter")
    ligated_plain = ligation_product_scene(
        library, "genomic insert", label="bottlenecked adapter-ligated library",
        adapter_group="Illumina Y adapter", highlight_components=False)
    families = molecule("redundantly sequenced bottleneck", (
        role("Watson read family"), role("same two genomic endpoints"), role("Crick read family")))
    w = Workflow(MolecularState("Sheared genomic DNA", tuple(dna.rows())))
    w.react("End-repair the sheared fragments",
            Scene.duplex(list(repaired), unpaired=("genomic insert'",)).rows(),
            name="end-repaired DNA")
    react_dA_tail_then_ligate(
        w, prepared, TRUSEQ_Y_END, ligated,
        dA_action="Add one 3-prime dA to each end",
        ligation_action="Ligate Illumina Y adapters",
        ligation_note="Adapter orientation distinguishes Watson and Crick descendants.")
    w.react("Dilute to about one haploid-genome equivalent", ligated_plain.rows(),
            name="bottlenecked adapter-ligated library",
            note="For 500 ng input the typical dilution was 100,000-fold; it is imposed before PCR.")
    w.react("PCR-amplify and paired-end sequence redundantly", families.rows(),
            name="Watson and Crick read families",
            note="The two genomic endpoints are the endogenous family key; no synthetic UMI is used.")
    return w


def _botseqs_final_rows(library: Construct):
    """Draw the actual adapter-ligation boundaries on the canonical TruSeq molecule."""
    sc = ligation_product_scene(library, "genomic insert", label=library.name,
                                adapter_group="Illumina Y adapter",
                                highlight_components=False)
    sc.labels("top")
    return sc.rows()


def _cirseq_workflow(library: Construct) -> Workflow:
    rna = molecule("purified viral RNA", (role("viral RNA"),))
    fragment = molecule("size-selected RNA", (role("85–100-nt RNA fragment", length_bp=92),))
    circ = closed_circle("circular RNA template", (role("RNA fragment"),),
                         closure_label="RNA intramolecular ligation")
    tandem = molecule("rolling-circle cDNA", (
        role("template copy 1"), role("template copy 2"), role("template copy 3")))
    ds_insert = Construct([role("tandem-repeat cDNA")],
                          name="double-stranded tandem-repeat cDNA")
    ds = Scene.duplex(list(ds_insert), label=ds_insert.name,
                      unpaired=("tandem-repeat cDNA'",))
    prepared = repair_and_dA_tail(ds_insert)
    ligated = ligation_product_scene(library, "tandem-repeat cDNA",
                                     label="TruSeq-ligated CirSeq library",
                                     adapter_group="TruSeq indexed adapter")
    w = Workflow(MolecularState("Starting material", tuple(rna.rows())))
    w.react("Fragment with Zn2+ and size-select", fragment.rows(), name="Fragment and size-select",
            note="Fragments are kept short enough for at least three copies in one read.")
    w.react("5′-phosphorylate and circularize RNA", circ.rows(), name="Circularize RNA",
            note="T4 RNA ligase 1 closes the native 3′→5′ junction.")
    w.react("Random-prime rolling-circle reverse transcription", tandem.rows(),
            name="Rolling-circle reverse transcription",
            note="Physically linked repeats derive independently from the same RNA template.")
    w.react("Synthesize the second strand and blunt the cDNA", ds.rows(),
            name="double-stranded tandem-repeat cDNA")
    react_dA_tail_then_ligate(
        w, prepared, TRUSEQ_Y_END, ligated,
        dA_action="Add one 3-prime dA to each cDNA end",
        ligation_action="Ligate TruSeq indexed adapters",
        ligation_note="Size selection retains tandem inserts long enough for within-read consensus.")
    return w


def _duplex_repair_workflow() -> Workflow:
    damaged = molecule("damaged or nicked DNA duplex", (role("damaged duplex with ragged ends"),))
    excised = molecule("lesion-excised duplex", (role("short gaps"), role("trimmed overhangs")))
    blunt = molecule("blunt repaired duplex", (role("restricted fill-in"), role("remaining nicks")))
    sealed = molecule("nick-sealed duplex", (role("original duplex interior"),))
    ata = molecule("dA-tailed duplex", (role("original duplex interior"), role("3′ dA", length_bp=1)))
    adapted_construct = Construct((
        role("custom duplex UMI adapter A", tag="umi",
             molecular_feature=feature("duplex_umi_a", "umi", "random")),
        role("repaired insert"),
        role("custom duplex UMI adapter B", tag="umi",
             molecular_feature=feature("duplex_umi_b", "umi", "random"))),
        name="duplex-UMI library")
    adapted = ligation_product_scene(
        adapted_construct, "repaired insert", adapter_group="duplex UMI adapter A",
        right_adapter_group="duplex UMI adapter B")
    return _workflow(damaged, (
        ("Excise damaged bases and remove single-stranded overhangs", excised,
         "EndoIV, Fpg, UDG, T4 PDG and EndoVIII act with ExoVII."),
        ("Blunt ends and restrict fill-in", blunt,
         "T4 PNK and T4 DNA polymerase fill short gaps and residual short 5′ overhangs."),
        ("Seal nicks with HiFi Taq DNA ligase", sealed,
         "Nick sealing precedes dA-tailing so polymerase does not start from open nicks."),
        ("Restricted dA-tail with dATP as the only nucleotide", ata,
         "Klenow fragment exo− and Taq add the ligation-compatible 3′ dA."),
        ("Ligate custom dual-index duplex-UMI adapters and PCR-amplify", adapted,
         "Complementary UMI families recover both original strands."),
    ))


TIME_B_PREFIX = "TGGGTGGAGGGTGG"
TIME_BARCODE = "AAAGA"                 # TS301, one published representative
TIME_TN5_REVERSE = "/5Phos/CTGTCTCTTATACACATCT"
TIME_METHYLATED_A = "T/iMe-dC/GT/iMe-dC/GG/iMe-dC/AG/iMe-dC/GT/iMe-dC/AGATGTGTATAAGAGA/iMe-dC/AG"
TIME_B_TS301 = TIME_B_PREFIX + TIME_BARCODE + nx.ME
TIME_BLOCK_B = "CTGTCTCTTATACACATCTHHHHHCCACCCTCCACCCA"
TIME_BLOCK_A = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
TIME_N501 = "AATGATACGGCGACCACCGAGATCTACACTAGATCGCTCGTCGGCAGCGTC"
TIME_N701 = "CAAGCAGAAGACGGCATACGAGATTCGCCTTAAAGCAGTGGTATCAACGCAGATCTGGGTGGAGGGTGG"
TIME_I5 = "TAGATCGC"
TIME_I7 = "TCGCCTTA"
TIME_CUSTOM_R2 = "AAGCAGTGGTATCAACGCAGATC" + TIME_B_PREFIX
TIME_CUSTOM_I1 = revcomp(TIME_CUSTOM_R2)


def _time_library():
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("i5 (N501)", TIME_I5, "cbc",
            feature=feature("pool_index_i5", "sample_index", "whitelist",
                            whitelist="TIME-seq Supplementary Table 5 indexed PCR primers")),
        seg("S5", nx.S5, "s5"), seg("left mosaic end", nx.ME, "me"),
        role("bisulfite-converted captured insert"),
        seg("right mosaic end reverse complement", nx.ME_RC, "me"),
        seg("Tn5 sample barcode reverse complement", revcomp(TIME_BARCODE), "cbc",
            feature=feature("tn5_sample_barcode", "inline_barcode", "whitelist",
                            whitelist="TIME-seq Supplementary Table 5 TS301–TS364")),
        seg("B-adaptor prefix reverse complement", revcomp(TIME_B_PREFIX)),
        seg("custom sequencing handle reverse complement",
            revcomp("AAGCAGTGGTATCAACGCAGATC"), "r2"),
        seg("i7 reverse complement (N701)", revcomp(TIME_I7), "cbc",
            feature=feature("pool_index_i7", "sample_index", "whitelist",
                            whitelist="TIME-seq Supplementary Table 5 indexed PCR primers")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="TIME-seq final library")
    primers = (sp.NEXTERA["R1"],
               sp.custom("Index 1 (i7)", "TIME-seq custom Index 1", TIME_CUSTOM_I1,
                         "Griffin et al. Supplementary Table 5"),
               sp.NEXTERA["I2"],
               sp.custom("Read 2", "TIME-seq custom Read 2", TIME_CUSTOM_R2,
                         "Griffin et al. Supplementary Table 5"))
    errors = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if errors:
        raise ValueError("TIME-seq final library: " + "; ".join(errors))
    return lib, primers


def _time_workflow() -> Workflow:
    dna = molecule("sample genomic DNA", (role("genomic DNA"),))
    tagged = molecule("sample-indexed tagmentation product", (
        role("methylated A adapter"), role("genomic fragment"),
        role("5-nt Tn5 sample barcode", tag="cbc",
             molecular_feature=feature("tn5_sample_barcode", "inline_barcode", "whitelist",
                                       whitelist="TIME-seq Supplementary Table 5 TS301–TS364")),
        role("B adapter")),
        junctions=(("methylated A adapter", "genomic fragment", "Tn5 strand transfer"),
                   ("genomic fragment", "5-nt Tn5 sample barcode", "Tn5 strand transfer")))
    pooled = molecule("pooled indexed DNA", (role("many Tn5-barcoded samples"),))
    captured = molecule("hybrid-captured pool", (role("biotin-RNA bait"), role("target loci")))
    converted = molecule("bisulfite-converted captured DNA", (role("converted target loci"),))
    return _workflow(dna, (
        ("Tagment each sample with methylated A and indexed B transposomes", tagged,
         "The 5-nt B-adapter barcode identifies the sample before pooling."),
        ("Pool samples and perform methylated end repair", pooled,
         "5-methyl-dCTP replaces dCTP so adapter cytosines survive bisulfite treatment."),
        ("Hybridize biotinylated RNA baits and enrich target loci", captured,
         "Capture occurs before bisulfite conversion."),
        ("Bisulfite-convert the captured DNA", converted,
         "Unmethylated cytosines become uracils while protected methylcytosines remain."),
        ("PCR with indexed P5 and P7 primers", molecule("final TIME-seq library", (
            role("P5/i5 end"), role("converted captured insert"), role("Tn5 sample barcode"),
            role("P7/i7 end"))),
         "Pool indices and the internal Tn5 barcode support two-level demultiplexing."),
    ))


def _botseqs_graph(library: Construct, primers) -> ProtocolGraph:
    source = "Hoang et al., PNAS (2016); doi:10.1073/pnas.1607794113"
    graph = ProtocolGraph("BotSeqS", authority=source)
    genomic = graph.input(
        "genomic_dna", "genomic DNA", MaterialKind.BIOLOGICAL,
        state=molecule("genomic DNA", (role("genomic DNA"),)), primary=True)
    shearing = graph.input("shearing_condition", "mechanical DNA-shearing condition",
                           MaterialKind.CONDITION)
    fragments = graph.product(
        "sheared_fragments", "sheared genomic DNA with native endpoints",
        MaterialKind.MOLECULAR,
        state=molecule("sheared genomic DNA", (role("genomic fragments with native endpoints"),)))
    graph.react("shear", "Shear genomic DNA", StepKind.FRAGMENTATION,
                ((genomic, InputRole.TEMPLATE), (shearing, InputRole.CONDITION)),
                (fragments,))
    repair = graph.input("end_repair_mix", "end-repair enzyme mix", MaterialKind.CATALYST)
    repaired_state = Construct([role("genomic insert")], name="end-repaired genomic insert")
    repaired = graph.product("repaired_fragments", "end-repaired genomic fragments",
                             MaterialKind.MOLECULAR, state=repaired_state)
    graph.react("end_repair", "End-repair the sheared fragments", StepKind.END_REPAIR,
                ((fragments, InputRole.TEMPLATE), (repair, InputRole.CATALYST)),
                (repaired,))
    tailing = graph.input("a_tailing_enzyme", "dA-tailing enzyme", MaterialKind.CATALYST)
    tailed_state = repair_and_dA_tail(repaired_state)
    tailed = graph.product("tailed_fragments", "dA-tailed genomic fragments",
                           MaterialKind.MOLECULAR, state=tailed_state)
    graph.react("a_tail", "Add one 3-prime dA to each end", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)),
                (tailed,))
    adapters = graph.input("truseq_y_adapters", "TruSeq Y adapters", MaterialKind.ADAPTER)
    ligase = graph.input("adapter_ligase", "DNA adapter ligase", MaterialKind.CATALYST)
    adapted = graph.product("adapted_pool", "TruSeq adapter-ligated genomic library pool",
                            MaterialKind.MOLECULAR, state=library)
    graph.react("adapter_ligation", "Ligate Illumina Y adapters", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,), geometry=TRUSEQ_DA_LIGATION)

    dilution = graph.input("dilution_bottleneck", "approximately one-genome-equivalent dilution",
                           MaterialKind.CONDITION)
    bottleneck = graph.product("bottleneck", "bottlenecked adapter-ligated molecules",
                               MaterialKind.MOLECULAR, state=library)
    excluded = graph.product("excluded_molecules", "adapter-ligated molecules excluded by dilution",
                             MaterialKind.WASTE)
    graph.react("bottleneck_sampling", "Dilute to about one haploid-genome equivalent",
                StepKind.SAMPLING,
                ((adapted, InputRole.TEMPLATE), (dilution, InputRole.CONDITION)),
                (bottleneck, excluded), changes_molecule=False)
    graph.discard(excluded, "falls outside the deliberately sampled molecular bottleneck")

    forward = graph.input("pcr_forward", "TruSeq library PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input("pcr_reverse", "TruSeq indexed library PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    final = graph.product("sequencing_library", "redundantly amplified BotSeqS library",
                          MaterialKind.LIBRARY, state=library)
    graph.react("library_pcr", "PCR-amplify the bottlenecked library", StepKind.PCR,
                ((bottleneck, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (final,), geometry=PCRGeometry.unavailable(
                    library,
                    reason="The paper uses ordinary TruSeq amplification; the chosen index oligo varies by library."))
    graph.sequence(final, SequencingPlan(
        library, tuple(primers), tuple(p.role for p in primers),
        {"Read 1": BOTSEQS_READ_BP, "Read 2": BOTSEQS_READ_BP}))
    graph.declare_reviewed_scope((
        "genomic_dna", "shearing_condition", "shear", "end_repair_mix", "end_repair",
        "a_tailing_enzyme", "a_tail", "truseq_y_adapters", "adapter_ligase",
        "adapter_ligation", "dilution_bottleneck", "bottleneck_sampling",
        "pcr_forward", "pcr_reverse", "pcr_polymerase", "library_pcr",
        "sequencing_library",
    ), source)
    return graph.require_publishable()


def _cirseq_graph(library: Construct, primers) -> ProtocolGraph:
    source = "Acevedo and Andino, Nature Protocols (2014); doi:10.1038/nprot.2014.118"
    graph = ProtocolGraph("CirSeq", authority=source)
    viral_rna = graph.input(
        "viral_rna", "purified viral RNA", MaterialKind.BIOLOGICAL,
        state=molecule("purified viral RNA", (role("viral RNA"),)), primary=True)
    zinc = graph.input("zinc_fragmentation", "Zn2+ fragmentation condition",
                       MaterialKind.CONDITION)
    fragments = graph.product(
        "rna_fragments", "fragmented viral RNA", MaterialKind.MOLECULAR,
        state=molecule("fragmented viral RNA", (role("short viral RNA fragments"),)))
    graph.react("fragment", "Fragment viral RNA with Zn2+", StepKind.FRAGMENTATION,
                ((viral_rna, InputRole.TEMPLATE), (zinc, InputRole.CONDITION)),
                (fragments,))
    size_rule = graph.input("size_selection", "85–100-nt size-selection rule",
                            MaterialKind.SELECTION)
    selected = graph.product(
        "selected_rna", "85–100-nt viral RNA fragments", MaterialKind.MOLECULAR,
        state=molecule("size-selected RNA", (role("85–100-nt RNA fragment", length_bp=92),)))
    rejected = graph.product("rejected_sizes", "RNA outside the 85–100-nt window",
                             MaterialKind.WASTE)
    graph.react("select_size", "Select 85–100-nt RNA fragments", StepKind.SELECTION,
                ((fragments, InputRole.TEMPLATE), (size_rule, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "does not fit at least three linked copies within the long read")

    kinase = graph.input("polynucleotide_kinase", "5-prime polynucleotide kinase",
                         MaterialKind.CATALYST)
    phosphorylated = graph.product("phosphorylated_rna", "5-prime-phosphorylated RNA fragments",
                                   MaterialKind.MOLECULAR, state=selected.state)
    graph.react("phosphorylate", "5-prime-phosphorylate selected RNA",
                StepKind.PHOSPHORYLATION,
                ((selected, InputRole.TEMPLATE), (kinase, InputRole.CATALYST)),
                (phosphorylated,))
    rna_ligase = graph.input("rna_ligase", "T4 RNA ligase 1", MaterialKind.CATALYST)
    circle_state = closed_circle("circular RNA template", (role("RNA fragment"),),
                                 closure_label="RNA intramolecular ligation")
    circles = graph.product("rna_circles", "circularized viral RNA fragments",
                            MaterialKind.MOLECULAR, state=circle_state)
    graph.react("circularize", "Circularize RNA fragments", StepKind.CIRCULARIZATION,
                ((phosphorylated, InputRole.TEMPLATE), (rna_ligase, InputRole.CATALYST)),
                (circles,), geometry=CircularizationGeometry.known(
                    "single-strand intramolecular ligation",
                    statement="The protocol circularizes phosphorylated RNA fragments with T4 RNA ligase 1."))

    random_rt = graph.input("random_rt_primer", "random reverse-transcription primer",
                            MaterialKind.OLIGO)
    reverse_transcriptase = graph.input("reverse_transcriptase", "rolling-circle reverse transcriptase",
                                       MaterialKind.CATALYST)
    tandem_state = molecule("rolling-circle cDNA", (
        role("template copy 1"), role("template copy 2"), role("template copy 3")))
    tandem = graph.product("tandem_cdna", "rolling-circle tandem-repeat cDNA",
                           MaterialKind.MOLECULAR, state=tandem_state)
    graph.react("rolling_circle_rt", "Random-prime rolling-circle reverse transcription",
                StepKind.REVERSE_TRANSCRIPTION,
                ((circles, InputRole.TEMPLATE), (random_rt, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (tandem,))
    second_primer = graph.input("second_strand_primer", "random second-strand primer",
                                MaterialKind.OLIGO)
    dna_polymerase = graph.input("second_strand_polymerase", "second-strand DNA polymerase",
                                 MaterialKind.CATALYST)
    ds_construct = Construct([role("tandem-repeat cDNA")],
                             name="double-stranded tandem-repeat cDNA")
    double_stranded = graph.product("double_stranded_cdna", ds_construct.name,
                                    MaterialKind.MOLECULAR, state=ds_construct)
    graph.react("second_strand", "Synthesize the second cDNA strand", StepKind.EXTENSION,
                ((tandem, InputRole.TEMPLATE), (second_primer, InputRole.STRUCTURAL),
                 (dna_polymerase, InputRole.CATALYST)), (double_stranded,))
    repair = graph.input("end_repair_mix", "blunt-end repair enzyme mix", MaterialKind.CATALYST)
    repaired = graph.product("repaired_cdna", "blunt tandem-repeat cDNA",
                             MaterialKind.MOLECULAR, state=ds_construct)
    graph.react("end_repair", "Blunt the double-stranded cDNA", StepKind.END_REPAIR,
                ((double_stranded, InputRole.TEMPLATE), (repair, InputRole.CATALYST)),
                (repaired,))
    tailing = graph.input("a_tailing_enzyme", "dA-tailing enzyme", MaterialKind.CATALYST)
    tailed = graph.product("tailed_cdna", "dA-tailed tandem-repeat cDNA",
                           MaterialKind.MOLECULAR,
                           state=repair_and_dA_tail(ds_construct))
    graph.react("a_tail", "Add one 3-prime dA to each cDNA end", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapters = graph.input("truseq_adapters", "TruSeq indexed adapters", MaterialKind.ADAPTER)
    ligase = graph.input("adapter_ligase", "DNA adapter ligase", MaterialKind.CATALYST)
    final = graph.product("sequencing_library", "TruSeq-ligated CirSeq library",
                          MaterialKind.LIBRARY, state=library)
    graph.react("adapter_ligation", "Ligate TruSeq indexed adapters", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (final,), geometry=TRUSEQ_DA_LIGATION)
    graph.sequence(final, SequencingPlan(
        library, tuple(primers), tuple(p.role for p in primers),
        {"Read 1": 300, "Index 1 (i7)": 8, "Read 2": 300},
        {"Read 1": "include", "Read 2": "include"}))
    graph.declare_reviewed_scope((
        "viral_rna", "zinc_fragmentation", "fragment", "size_selection", "select_size",
        "polynucleotide_kinase", "phosphorylate", "rna_ligase", "circularize",
        "random_rt_primer", "reverse_transcriptase", "rolling_circle_rt",
        "second_strand_primer", "second_strand_polymerase", "second_strand",
        "end_repair_mix", "end_repair", "a_tailing_enzyme", "a_tail",
        "truseq_adapters", "adapter_ligase", "adapter_ligation", "sequencing_library",
    ), source)
    return graph.require_publishable()


def _duplex_repair_graph(final: Construct) -> ProtocolGraph:
    source = "Xiong et al., Nucleic Acids Research (2022); doi:10.1093/nar/gkab855"
    graph = ProtocolGraph("Duplex-Repair", authority=source)
    damaged = graph.input(
        "damaged_dna", "damaged or nicked DNA duplex",
        MaterialKind.BIOLOGICAL,
        state=molecule("damaged or nicked DNA duplex", (
            role("damaged duplex with ragged ends"),)), primary=True)
    lesion_enzymes = graph.input(
        "lesion_enzymes", "EndoIV/Fpg/UDG/T4-PDG/EndoVIII/ExoVII repair cocktail",
        MaterialKind.CATALYST)
    excised = graph.product(
        "excised_dna", "lesion-excised duplex with trimmed overhangs",
        MaterialKind.MOLECULAR,
        state=molecule("lesion-excised duplex", (role("short gaps"), role("trimmed overhangs"))))
    graph.react(
        "lesion_excision", "Excise damaged bases and remove single-stranded overhangs",
        StepKind.CHEMICAL_CONVERSION,
        ((damaged, InputRole.TEMPLATE), (lesion_enzymes, InputRole.CATALYST)),
        (excised,))
    blunting = graph.input("blunting_enzymes", "T4 PNK and T4 DNA polymerase",
                           MaterialKind.CATALYST)
    blunt = graph.product(
        "blunt_dna", "blunt repaired duplex with restricted fill-in",
        MaterialKind.MOLECULAR,
        state=molecule("blunt repaired duplex", (
            role("restricted fill-in"), role("remaining nicks"))))
    graph.react("blunt_ends", "Blunt ends and restrict fill-in", StepKind.END_REPAIR,
                ((excised, InputRole.TEMPLATE), (blunting, InputRole.CATALYST)),
                (blunt,))
    nick_ligase = graph.input("nick_ligase", "HiFi Taq DNA ligase",
                              MaterialKind.CATALYST)
    sealed = graph.product(
        "sealed_dna", "nick-sealed repaired duplex", MaterialKind.MOLECULAR,
        state=molecule("nick-sealed duplex", (role("original duplex interior"),)))
    graph.react("seal_nicks", "Seal nicks before tailing", StepKind.NICK_SEALING,
                ((blunt, InputRole.TEMPLATE), (nick_ligase, InputRole.CATALYST)),
                (sealed,))
    tailing = graph.input(
        "restricted_tailing_enzyme", "Klenow fragment exo-minus and Taq polymerase",
        MaterialKind.CATALYST)
    datp = graph.input("datp_only", "dATP as the only tailing nucleotide",
                       MaterialKind.MOLECULAR)
    tailed = graph.product(
        "tailed_dna", "restricted dA-tailed repaired duplex",
        MaterialKind.MOLECULAR,
        state=molecule("dA-tailed duplex", (
            role("original duplex interior"), role("3-prime dA", length_bp=1))))
    graph.react("restricted_a_tail", "Restricted dA-tail with dATP only",
                StepKind.TAILING,
                ((sealed, InputRole.TEMPLATE), (datp, InputRole.STRUCTURAL),
                 (tailing, InputRole.CATALYST)), (tailed,))
    adapters = graph.input("duplex_umi_adapters", "custom dual-index duplex-UMI adapters",
                           MaterialKind.ADAPTER)
    adapter_ligase = graph.input("adapter_ligase", "DNA adapter ligase",
                                 MaterialKind.CATALYST)
    adapted = graph.product("adapted_dna", "duplex-UMI adapter-ligated DNA",
                            MaterialKind.MOLECULAR, state=final)
    graph.react("adapter_ligation", "Ligate custom duplex-UMI adapters",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (adapter_ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The paper identifies custom duplex-UMI adapters but does not publish their ligating-end sequences."))
    forward = graph.input("pcr_forward", "custom duplex-library PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input("pcr_reverse", "custom duplex-library PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("sequencing_library", "Duplex-Repair sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("library_pcr", "PCR-amplify the duplex-UMI library", StepKind.PCR,
                ((adapted, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The accessible source does not disclose the custom duplex-library primer sequences."))
    graph.sequence(library, SequencingPlan(
        limitation="Paired-end HiSeq 2500 or HiSeq X sequencing is reported, but custom adapter and sequencing-primer bases are unavailable."))
    graph.declare_reviewed_scope((
        "damaged_dna", "lesion_enzymes", "lesion_excision", "blunting_enzymes",
        "blunt_ends", "nick_ligase", "seal_nicks", "restricted_tailing_enzyme",
        "datp_only", "restricted_a_tail", "duplex_umi_adapters", "adapter_ligase",
        "adapter_ligation", "pcr_forward", "pcr_reverse", "pcr_polymerase",
        "library_pcr", "sequencing_library",
    ), source)
    return graph.require_publishable()


def _time_graph(library: Construct, primers) -> ProtocolGraph:
    source = ("Griffin et al., Nature Aging (2024), Supplementary Information and "
              "Supplementary Table 5; doi:10.1038/s43587-023-00555-2")
    graph = ProtocolGraph("TIME-seq", authority=source)
    batch = graph.input(
        "sample_batch", "batch of genomic-DNA samples",
        MaterialKind.BIOLOGICAL,
        state=molecule("sample genomic DNA batch", (role("independent genomic-DNA samples"),)),
        primary=True)
    sample_a = graph.product("sample_a", "representative genomic-DNA sample A",
                             MaterialKind.MOLECULAR, state=batch.state)
    sample_b = graph.product("sample_b", "representative genomic-DNA sample B",
                             MaterialKind.MOLECULAR, state=batch.state)
    graph.react("separate_samples", "Keep samples separate for Tn5 indexing",
                StepKind.SPLIT, ((batch, InputRole.TEMPLATE),), (sample_a, sample_b),
                changes_molecule=False)
    reviewed = ["sample_batch", "separate_samples"]

    def tagment(prefix: str, sample, barcode_label: str):
        adapters = graph.input(
            prefix + "_adapters",
            "methylated A and indexed B transposome adapters (" + barcode_label + ")",
            MaterialKind.ADAPTER)
        tn5 = graph.input(prefix + "_tn5", "Tn5 transposase for " + barcode_label,
                          MaterialKind.CATALYST)
        tagged = graph.product(
            prefix + "_tagged", barcode_label + " tagmented genomic DNA",
            MaterialKind.MOLECULAR,
            state=molecule("sample-indexed tagmentation product", (
                role("methylated A adapter"), role("genomic fragment"),
                role("5-nt Tn5 sample barcode", tag="cbc",
                     molecular_feature=feature(
                         "tn5_sample_barcode", "inline_barcode", "whitelist",
                         whitelist="TIME-seq Supplementary Table 5 TS301–TS364")),
                role("B adapter"))))
        losses = graph.product(prefix + "_losses", barcode_label + " non-library Tn5 products",
                               MaterialKind.WASTE)
        graph.react(
            prefix + "_tagment", "Tagment " + barcode_label + " with its indexed B transposome",
            StepKind.TAGMENTATION,
            ((sample, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
             (tn5, InputRole.CATALYST)), (tagged, losses))
        graph.discard(losses, "does not retain an amplifiable A/B tagmentation configuration")
        reviewed.extend((prefix + "_adapters", prefix + "_tn5", prefix + "_tagment"))
        return tagged

    tagged_a = tagment("a", sample_a, "sample A")
    tagged_b = tagment("b", sample_b, "sample B")
    pooled = graph.product(
        "pooled_samples", "pooled pre-indexed TIME-seq DNA",
        MaterialKind.MOLECULAR,
        state=molecule("pooled indexed DNA", (role("many Tn5-barcoded samples"),)))
    graph.react("pool_samples", "Pool samples after Tn5 barcode installation",
                StepKind.POOL,
                ((tagged_a, InputRole.TEMPLATE), (tagged_b, InputRole.TEMPLATE)),
                (pooled,), changes_molecule=False)
    reviewed.append("pool_samples")

    repair = graph.input("methylated_repair_enzymes", "methylated end-repair enzyme mix",
                         MaterialKind.CATALYST)
    methylated_dctp = graph.input("methylated_dctp", "5-methyl-dCTP repair nucleotide",
                                  MaterialKind.MOLECULAR)
    repaired = graph.product("repaired_pool", "methylated-end-repaired indexed DNA pool",
                             MaterialKind.MOLECULAR, state=pooled.state)
    graph.react("methylated_end_repair", "Perform methylated end repair",
                StepKind.END_REPAIR,
                ((pooled, InputRole.TEMPLATE), (methylated_dctp, InputRole.STRUCTURAL),
                 (repair, InputRole.CATALYST)), (repaired,))
    baits = graph.input("biotin_rna_baits", "biotinylated target-locus RNA baits",
                        MaterialKind.OLIGO)
    hybrids = graph.product(
        "bait_hybrids", "DNA hybridized to biotinylated target baits",
        MaterialKind.MOLECULAR,
        state=molecule("target-bait hybrids", (role("biotin-RNA bait"), role("target loci"))))
    graph.react("bait_hybridization", "Hybridize biotinylated RNA baits",
                StepKind.HYBRIDIZATION,
                ((repaired, InputRole.TEMPLATE), (baits, InputRole.STRUCTURAL)),
                (hybrids,))
    capture = graph.input("capture_beads", "streptavidin target-capture beads",
                          MaterialKind.SELECTION)
    captured = graph.product(
        "captured_pool", "hybrid-captured target DNA", MaterialKind.MOLECULAR,
        state=molecule("hybrid-captured pool", (role("biotin-RNA bait"), role("target loci"))))
    uncaptured = graph.product("uncaptured_pool", "uncaptured off-target DNA",
                               MaterialKind.WASTE)
    graph.react("target_capture", "Enrich bait-bound target loci", StepKind.AFFINITY_CAPTURE,
                ((hybrids, InputRole.TEMPLATE), (capture, InputRole.SELECTION)),
                (captured, uncaptured))
    graph.discard(uncaptured, "does not hybridize to the target-locus bait panel")
    bisulfite = graph.input("bisulfite", "bisulfite conversion condition",
                            MaterialKind.CONDITION)
    converted = graph.product(
        "converted_pool", "bisulfite-converted captured DNA", MaterialKind.MOLECULAR,
        state=molecule("bisulfite-converted captured DNA", (role("converted target loci"),)))
    graph.react("bisulfite_conversion", "Bisulfite-convert captured DNA",
                StepKind.CHEMICAL_CONVERSION,
                ((captured, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)),
                (converted,))
    n501 = graph.input("n501_primer", "published indexed PCR primer N501",
                       MaterialKind.OLIGO)
    n701 = graph.input("n701_primer", "published indexed PCR primer N701",
                       MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "bisulfite-library PCR polymerase",
                             MaterialKind.CATALYST)
    final = graph.product("sequencing_library", "final TIME-seq library",
                          MaterialKind.LIBRARY, state=library)
    geometry = PCRGeometry.exact(
        library, PCRPrimer("N501", TIME_N501), PCRPrimer("N701", TIME_N701),
        statement="Griffin et al. Supplementary Table 5 publishes both indexed PCR primers.")
    graph.react("indexed_pcr", "PCR with indexed P5 and P7 primers", StepKind.PCR,
                ((converted, InputRole.TEMPLATE), (n501, InputRole.STRUCTURAL),
                 (n701, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (final,), geometry=geometry)
    graph.sequence(final, SequencingPlan(
        library, tuple(primers), tuple(p.role for p in primers),
        {"Read 1": 145, "Index 1 (i7)": 8, "Index 2 (i5)": 8, "Read 2": 5},
        {"Read 1": "include", "Index 1 (i7)": "demultiplex",
         "Index 2 (i5)": "demultiplex", "Read 2": "extract"}))
    reviewed.extend((
        "methylated_repair_enzymes", "methylated_dctp", "methylated_end_repair",
        "biotin_rna_baits", "bait_hybridization", "capture_beads", "target_capture",
        "bisulfite", "bisulfite_conversion", "n501_primer", "n701_primer",
        "pcr_polymerase", "indexed_pcr", "sequencing_library",
    ))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def get(key: str):
    if key == "trueprime":
        final = Construct([role("TruePrime-amplified genomic DNA")], name="TruePrime WGA product")
        spec = SimpleNamespace(
            TITLE="TruePrime — primase-driven whole-genome amplification",
            NOTES="01_trueprime.html",
            SOURCE='<a href="https://doi.org/10.1038/ncomms13296">Picher et al., <i>Nature Communications</i> (2016)</a>.',
            SUMMARY="TthPrimPol makes 7–9-nt DNA primers directly on denatured genomic DNA; phi29 DNA polymerase extends them by multiple-displacement synthesis.",
            CAVEAT="TruePrime ends at amplified DNA. The paper used an unspecified NEBNext library preparation downstream, so no adapter sequence or sequencing-primer landing site is asserted here.",
            oligos=lambda: (), workflow=_trueprime_workflow, FINAL_LIBRARY=final,
            FINAL_CAPTION="The named protocol product is amplified genomic DNA, before sequencing adapters.",
            SEQ_PRIMERS=(), SEQUENCING_INTRO="", READ_LENGTHS=None,
            SEQUENCING_ENDING="The paper fragmented this product, prepared a NEBNext library and performed paired-end HiSeq 2500 sequencing.",
            SEQUENCING_UNAVAILABLE="The publication does not identify the NEBNext kit generation or publish its adapter and sequencing-primer sequences.")
        spec.protocol_graph = lambda: _trueprime_graph(final, spec.SOURCE)
        return spec
    if key == "botseqs":
        # The 2016 paper establishes standard TruSeq Y adapters but does not establish a
        # dual-index configuration.  Model the period-appropriate single-index molecule.
        lib, primers = truseq_library([role("genomic insert")],
                                      "BotSeqS TruSeq PCR-free library",
                                      dual_index=False)
        spec = SimpleNamespace(
            TITLE="BotSeqS — bottleneck sequencing system", NOTES="01_botseqs.html",
            SOURCE='<a href="https://doi.org/10.1073/pnas.1607794113">Hoang et al., <i>PNAS</i> (2016)</a>.',
            SUMMARY="A TruSeq PCR-free library is diluted after adapter ligation and before PCR; native fragment endpoints group redundant Watson and Crick read families without a synthetic UMI.",
            CAVEAT="", oligos=lambda: (), workflow=lambda: _botseqs_workflow(lib),
            FINAL_LIBRARY=lib, SEQ_PRIMERS=primers,
            READ_LENGTHS={"Read 1": BOTSEQS_READ_BP, "Read 2": BOTSEQS_READ_BP},
            final_rows=lambda: _botseqs_final_rows(lib),
            FINAL_CAPTION="Ordinary TruSeq paired-end molecule; error correction comes from the post-ligation bottleneck and shared genomic endpoints, not an added barcode.",
            SEQUENCING_INTRO="Paired reads establish both endpoints and P5/P7 orientation for Watson- and Crick-family consensus.")
        spec.protocol_graph = lambda: _botseqs_graph(lib, primers)
        return spec
    if key == "cirseq":
        lib, primers = truseq_library([role("tandem-repeat cDNA")], "CirSeq library", dual_index=False)
        spec = SimpleNamespace(
            TITLE="CirSeq — circular resequencing of RNA virus populations", NOTES="01_cirseq.html",
            SOURCE='<a href="https://doi.org/10.1038/nprot.2014.118">Acevedo and Andino, <i>Nature Protocols</i> (2014)</a>.',
            SUMMARY="Short viral RNA fragments are circularized and copied repeatedly by rolling-circle reverse transcription, placing at least three linked template copies inside one Illumina read.",
            CAVEAT="", oligos=lambda: (), workflow=lambda: _cirseq_workflow(lib),
            FINAL_LIBRARY=lib, SEQ_PRIMERS=primers, READ_LENGTHS={"Read 1": 300, "Index 1 (i7)": 8, "Read 2": 300},
            final_rows=lambda: ligation_product_scene(
                lib, "tandem-repeat cDNA", label=lib.name,
                adapter_group="TruSeq indexed adapter",
                highlight_components=False).rows(),
            FINAL_CAPTION="TruSeq adapters flank a double-stranded cDNA insert containing tandem copies of one circular RNA fragment.",
            SEQUENCING_INTRO="The protocol used single-end 300-cycle MiSeq reads; the binding view also shows the adapter's documented index and opposite-end primer sites.")
        spec.protocol_graph = lambda: _cirseq_graph(lib, primers)
        return spec
    if key == "duplex-repair":
        final = Construct([
            role("custom left duplex-UMI adapter", tag="umi",
                 molecular_feature=feature("duplex_umi_a", "umi", "random")),
            role("repaired insert"),
            role("custom right duplex-UMI adapter", tag="umi",
                 molecular_feature=feature("duplex_umi_b", "umi", "random")),
        ], name="Duplex-Repair sequencing library")
        spec = SimpleNamespace(
            TITLE="Duplex-Repair — damage-aware end repair for duplex sequencing", NOTES="01_duplex_repair.html",
            SOURCE='<a href="https://doi.org/10.1093/nar/gkab855">Xiong et al., <i>Nucleic Acids Research</i> (2022)</a>.',
            SUMMARY="Lesion excision, controlled blunting, nick sealing and dATP-only tailing are separated to minimize copying of one-strand damage into an apparent duplex mutation.",
            CAVEAT="The paper identifies custom dual-index duplex-UMI adapters but does not publish their sequences; adapter geometry and primer landing positions therefore remain unspecified.",
            oligos=lambda: (), workflow=_duplex_repair_workflow,
            FINAL_LIBRARY=final, FINAL_CAPTION="Published roles of the custom duplex-UMI library; exact adapter bases are unavailable.",
            SEQ_PRIMERS=(), READ_LENGTHS=None, SEQUENCING_INTRO="",
            SEQUENCING_ENDING="The authors sequenced custom duplex-UMI libraries on HiSeq 2500 or HiSeq X in paired-end mode.",
            SEQUENCING_UNAVAILABLE="The custom dual-index duplex-UMI adapter and primer oligos are not given in the accessible paper or supplement.")
        spec.protocol_graph = lambda: _duplex_repair_graph(final)
        return spec
    if key == "time-seq":
        lib, primers = _time_library()
        spec = SimpleNamespace(
            TITLE="TIME-seq — tagmentation-based indexing for methylation sequencing", NOTES="01_time_seq.html",
            SOURCE='<a href="https://doi.org/10.1038/s43587-023-00555-2">Griffin et al., <i>Nature Aging</i> (2024)</a>, Supplementary Information and Supplementary Table 5.',
            SUMMARY="Bisulfite-resistant Tn5 adapters install a 5-nt sample barcode before pooling; methylated end repair, hybrid capture, bisulfite conversion and indexed PCR produce the targeted methylation library.",
            CAVEAT="The displayed B-adapter uses published representative TS301; the paper publishes 64 alternative 5-nt sample indexes.",
            oligos=lambda: (
                f"<p><b>Tn5 reverse</b>: 5′-{TIME_TN5_REVERSE}-3′</p>",
                f"<p><b>Methylated adaptor A</b>: 5′-{TIME_METHYLATED_A}-3′</p>",
                f"<p><b>Indexed adaptor B, TS301</b>: 5′-{TIME_B_TS301}-3′</p>",
                f"<p><b>Adaptor B blocking oligo</b>: 5′-{TIME_BLOCK_B}-3′</p>",
                f"<p><b>Adaptor A blocking oligo</b>: 5′-{TIME_BLOCK_A}-3′</p>",
                f"<p><b>Indexed PCR primer N501</b>: 5′-{TIME_N501}-3′</p>",
                f"<p><b>Indexed PCR primer N701</b>: 5′-{TIME_N701}-3′</p>",
                f"<p><b>Custom Read 2 sequencing primer</b>: 5′-{TIME_CUSTOM_R2}-3′</p>",
                f"<p><b>Custom Index 1 sequencing primer</b>: 5′-{TIME_CUSTOM_I1}-3′</p>",
            ), workflow=_time_workflow,
            FINAL_LIBRARY=lib, SEQ_PRIMERS=primers,
            READ_LENGTHS={"Read 1": 145, "Index 1 (i7)": 8, "Index 2 (i5)": 8, "Read 2": 5},
            SHOW_READ_CONTENT=True,
            FINAL_CAPTION="The internal five-base Tn5 barcode identifies the pre-pooling sample; i5 and i7 identify the PCR pool. TS301 is drawn as one published representative.",
            SEQUENCING_INTRO="Read 2 and Index 1 use the two published custom sequencing primers; Read 1 and Index 2 use the matching Nextera sites. The economical run described in the Supplementary Information uses 145–153 Read 1 cycles, 8 i7 cycles, optional 8 i5 cycles and 5 Read 2 cycles; this view uses the 145-cycle Read 1 endpoint.")
        spec.protocol_graph = lambda: _time_graph(lib, primers)
        return spec
    raise KeyError(key)
