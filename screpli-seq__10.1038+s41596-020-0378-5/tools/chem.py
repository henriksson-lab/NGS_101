
from conformation_replication_specs import _three_c
from pcr import PCRPrimer
from chemdraw import Segment
from conformation_replication_specs import AdapterLigationGeometry, Construct, Definition, InputRole, MaterialKind, NascentPulse, PCRGeometry, ProtocolGraph, SequencingPlan, StepKind, da_tailed, molecule, role, sp
DEFINITION = Definition(
        "scRepli-seq — single-cell DNA replication sequencing", "10.1038/s41596-020-0378-5",
        "Miura et al., Nature Protocols (2020)",
        "A single S-phase cell is whole-genome amplified and sequenced; twofold copy-number differences distinguish replicated from unreplicated domains.",
        (role("SeqPlex-amplified genomic fragment (Read 1 may begin in a SEQXE handle)"),),
        "Each sample index identifies one single-cell library; libraries are built individually and the index is carried by the ligated NEXTflex Illumina index adaptor (NEXTflex DNA Barcodes - 48 in Procedure Step 58; dual-index NEXTflex sets are allowed). The article sequences pooled libraries (~40 cells per lane) on a HiSeq 1500 with 80-bp single reads (Library preparation and sequencing; Step 79 recommends >80 bp), and the pipeline maps Read 1 only (bwa aln/samse, Step 88) after trimming the Illumina adapter and any 5-prime SEQXE handle. The index read length is not stated; a 6-nt window follows the single-index trimming pattern printed in Step 82. Read 1 coverage, normalized to G1 controls, reports replicated and unreplicated domains. GEO GSE138634 also lists 84- and 127-nt and paired-end runs for individual deposited samples.", True,
        "Wet-lab route from the defining article (Miura et al. 2020, Procedure Steps 17-79): single cells sorted into Single Cell Lysis & Fragmentation Buffer with proteinase K (55 °C 1 h, then 99 °C 4 min, which also fragments the gDNA); SeqPlex Enhanced DNA Amplification Kit (SEQXE-50XN; DOP-PCR-based) pre-amplification then 24 cycles of amplification; 1.7x AMPure XP; Eco57I primer-removal digestion of the universal SEQXE adaptor (37 °C 1 h, 65 °C 20 min); 1.8x AMPure XP removing the ~40-bp digested adaptors; KAPA LTP (KK8232) end repair, 1.7x AMPure, on-bead A-tailing, 1.8x PEG/NaCl, NEXTflex adaptor ligation (15 min 20 °C), 1.0x PEG/NaCl; library PCR with the TPC1/TPC2 primer mix and KAPA HiFi HotStart ReadyMix for 5-7 cycles chosen by real-time PCR against KK2702 standards (published datasets 5 or 6 cycles); 0.9x AMPure XP to remove ~100-bp adaptor dimers. "
        "Not printed in the article: the SeqPlex primer and SEQXE adaptor sequences (the SEQXE pattern on the page is predicted from reads in Step 82, not a kit oligo), the NEXTflex adaptor sequences and their end geometry, and the index read length. The TruSeq-type adaptor ends drawn here are a model assumption (inferred styling); the TPC1/TPC2 primers are printed and their landing on those assumed P5/P7 ends is computed. The microscope/micromanipulator collection route and the GEO-listed Ion Proton samples are not drawn. The later DOP-PCR barcoded protocol (Bartlett et al., Curr. Protoc. 2022) is a different implementation and is not substituted.")


def _replication_graph(key: str, definition, final: Construct,
                       primers: tuple[sp.SeqPrimer, ...]) -> ProtocolGraph:
    """Typed replication-timing graphs retaining every comparison library."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    reviewed: list[str] = []
    plan = (SequencingPlan(final, tuple(primers), ("Read 1", "Index 1 (i7)"),
                           read_lengths={"Read 1": 80})
            if key == "screpli-seq" else
            SequencingPlan(final, tuple(primers), tuple(p.role for p in primers)))
    adapter_geometry = AdapterLigationGeometry.unavailable(
        "Miura 2020 Step 58 ligates a NEXTflex Illumina index adaptor (NEXTflex DNA Barcodes - 48) to the dA-tailed WGA DNA with KAPA LTP DNA ligase, but prints no adaptor sequence or end structure.")
    tpc1 = PCRPrimer("TPC1", "AATGATACGGCGACCACCGAG")
    tpc2 = PCRPrimer("TPC2", "CAAGCAGAAGACGGCATACGAG")

    def finish_short_read_branch(prefix: str, template, label: str, *, wga: bool) -> None:
        current = template
        if wga:
            lysis = graph.input(prefix + "_lysis", label + " Single Cell Lysis & Fragmentation Buffer with proteinase K; 55 °C 1 h, 99 °C 4 min (Steps 17-19)",
                                MaterialKind.CONDITION)
            lysed = graph.product(prefix + "_lysed", label + " heat-fragmented genomic DNA",
                                  MaterialKind.MOLECULAR,
                                  state=molecule(label + " fragmented gDNA", ("genome fragment",)))
            graph.react(prefix + "_lyse", "Lyse " + label + " and heat-fragment its gDNA",
                        StepKind.FRAGMENTATION,
                        ((current, InputRole.TEMPLATE), (lysis, InputRole.CONDITION)), (lysed,))
            current = lysed
            reviewed.extend((prefix + "_lysis", prefix + "_lyse"))
            dop_forward = graph.input(prefix + "_dop_forward", label + " SeqPlex Library Preparation Buffer/Enzyme priming (DOP-PCR-based; sequence not printed)",
                                      MaterialKind.OLIGO)
            dop_reverse = graph.input(prefix + "_dop_reverse", label + " SeqPlex amplification priming, second end (sequence not printed)",
                                      MaterialKind.OLIGO)
            wga_polymerase = graph.input(prefix + "_wga_polymerase", label + " DNA Polymerase for SeqPlex; 24 cycles 94 °C 15 s / 70 °C 5 min (Steps 22-25)",
                                         MaterialKind.CATALYST)
            amplified = graph.product(prefix + "_wga", label + " SeqPlex WGA product",
                                      MaterialKind.MOLECULAR,
                                      state=molecule(label + " WGA product", (
                                          "left SEQXE amplification handle",
                                          "genome copy",
                                          "right SEQXE amplification handle")))
            graph.react(prefix + "_wga_step", "Pre-amplify and whole-genome amplify " + label + " with SeqPlex",
                        StepKind.PCR,
                        ((current, InputRole.TEMPLATE),
                         (dop_forward, InputRole.STRUCTURAL),
                         (dop_reverse, InputRole.STRUCTURAL),
                         (wga_polymerase, InputRole.CATALYST)), (amplified,),
                        geometry=PCRGeometry.unavailable(
                            final, reason="SeqPlex primers are proprietary and not printed in the inspected sources; WGA products are not one fixed amplicon."))
            current = amplified
            reviewed.extend((prefix + "_dop_forward", prefix + "_dop_reverse",
                             prefix + "_wga_polymerase", prefix + "_wga_step"))
        else:
            extraction = graph.input(prefix + "_dna_extraction", label + " DNA purification condition",
                                     MaterialKind.CONDITION)
            purified = graph.product(prefix + "_purified_dna", label + " purified nascent DNA",
                                     MaterialKind.MOLECULAR, state=current.state)
            graph.react(prefix + "_extract", "Purify " + label + " DNA",
                        StepKind.EXTRACTION,
                        ((current, InputRole.TEMPLATE),
                         (extraction, InputRole.CONDITION)), (purified,))
            current = purified
            reviewed.extend((prefix + "_dna_extraction", prefix + "_extract"))

        if wga:
            beads = graph.input(prefix + "_ampure", label + " 1.7x AMPure XP purification (Steps 26-32)",
                                MaterialKind.SELECTION)
            cleaned = graph.product(prefix + "_wga_cleaned", label + " purified WGA DNA",
                                    MaterialKind.MOLECULAR, state=current.state)
            bead_waste = graph.product(prefix + "_ampure_waste", label + " AMPure supernatant and washes",
                                       MaterialKind.WASTE)
            graph.react(prefix + "_wga_cleanup", "Purify " + label + " WGA DNA",
                        StepKind.CLEANUP,
                        ((current, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                        (cleaned, bead_waste))
            graph.discard(bead_waste, "supernatant and washes removed from the beads")
            sonication = graph.input(prefix + "_sonication",
                                     label + " SeqPlex Primer Removal Enzyme Eco57I; 37 °C 1 h, 65 °C 20 min (Steps 34-35)",
                                     MaterialKind.CATALYST)
            digested = graph.product(prefix + "_digested", label + " Eco57I-digested WGA DNA",
                                     MaterialKind.MOLECULAR, state=molecule(label + " adapter-digested WGA fragment",
                                         ("possible residual SEQXE handle", "genome copy",
                                          "possible residual SEQXE handle")))
            cut_adaptors = graph.product(prefix + "_cut_adaptors", label + " digested SEQXE adaptor fragments (~40 bp)",
                                         MaterialKind.MOLECULAR, state=molecule("digested SEQXE adaptor", ("SEQXE adaptor fragment",)))
            graph.react(prefix + "_fragment", "Digest the universal SEQXE adaptor with Eco57I",
                        StepKind.DIGESTION,
                        ((cleaned, InputRole.TEMPLATE),
                         (sonication, InputRole.CATALYST)), (digested, cut_adaptors))
            removal = graph.input(prefix + "_adaptor_removal", label + " 1.8x AMPure XP (Steps 36-42; repeated if ~40-bp adaptors remain, Step 43)",
                                  MaterialKind.SELECTION)
            fragmented = graph.product(prefix + "_fragments", label + " adaptor-removed WGA DNA (100-500 bp smear, ~150 bp peak)",
                                       MaterialKind.MOLECULAR, state=digested.state)
            removed = graph.product(prefix + "_adaptor_waste", label + " supernatant with ~40-bp digested SEQXE adaptors",
                                    MaterialKind.WASTE)
            graph.react(prefix + "_adaptor_cleanup", "Remove digested SEQXE adaptors by AMPure size selection",
                        StepKind.CLEANUP,
                        ((digested, InputRole.TEMPLATE), (cut_adaptors, InputRole.TEMPLATE),
                         (removal, InputRole.SELECTION)), (fragmented, removed))
            graph.discard(removed, "short digested adaptors stay in the 1.8x AMPure supernatant")
            reviewed.extend((prefix + "_ampure", prefix + "_wga_cleanup", prefix + "_adaptor_removal",
                             prefix + "_adaptor_cleanup"))
        else:
            sonication = graph.input(prefix + "_sonication", label + " sonication condition",
                                     MaterialKind.CONDITION)
            fragmented = graph.product(prefix + "_fragments", label + " DNA fragments",
                                       MaterialKind.MOLECULAR, state=current.state)
            graph.react(prefix + "_fragment", "Fragment " + label + " DNA",
                        StepKind.FRAGMENTATION,
                        ((current, InputRole.TEMPLATE),
                         (sonication, InputRole.CONDITION)), (fragmented,))
        repair = graph.input(prefix + "_end_repair", label + (" KAPA LTP End Repair Enzyme Mix, 20 °C 30 min, then 1.7x AMPure (Steps 44-50)" if wga else " end-repair enzymes"),
                             MaterialKind.CATALYST)
        repaired = graph.product(prefix + "_repaired", label + " end-repaired DNA",
                                 MaterialKind.MOLECULAR, state=fragmented.state)
        graph.react(prefix + "_repair_step", "End-repair " + label + " DNA",
                    StepKind.END_REPAIR,
                    ((fragmented, InputRole.TEMPLATE),
                     (repair, InputRole.CATALYST)), (repaired,))
        tailing = graph.input(prefix + "_a_tailing", label + (" KAPA LTP A-Tailing Enzyme on beads, 30 °C 30 min, then 1.8x PEG/NaCl (Steps 51-57)" if wga else " dA-tailing enzyme"),
                              MaterialKind.CATALYST)
        tailed = graph.product(prefix + "_tailed", label + " dA-tailed DNA",
                               MaterialKind.MOLECULAR,
                               state=da_tailed(label + " dA-tailed DNA", label + " fragment"))
        graph.react(prefix + "_tail_step", "Add one 3-prime dA to each end",
                    StepKind.TAILING,
                    ((repaired, InputRole.TEMPLATE),
                     (tailing, InputRole.CATALYST)), (tailed,))
        adapters = graph.input(prefix + "_adapters", label + (" NEXTflex Illumina index adaptor, 5 µM (NEXTflex DNA Barcodes - 48; sequences not printed; Step 58)" if wga else " Illumina adapters"),
                               MaterialKind.ADAPTER, inferred=True)
        ligase = graph.input(prefix + "_adapter_ligase", label + (" KAPA LTP DNA Ligase, 20 °C 15 min, then 1.0x PEG/NaCl (Steps 58-65)" if wga else " DNA ligase"),
                             MaterialKind.CATALYST)
        adapted = graph.product(prefix + "_adapted", label + " adapter-ligated DNA",
                                MaterialKind.MOLECULAR, state=final, inferred=True)
        graph.react(prefix + "_ligation", "Ligate Illumina adapters to " + label,
                    StepKind.LIGATION,
                    ((tailed, InputRole.TEMPLATE),
                     (adapters, InputRole.STRUCTURAL),
                     (ligase, InputRole.CATALYST)), (adapted,),
                    geometry=adapter_geometry)
        reviewed.extend((
            prefix + "_sonication", prefix + "_fragment", prefix + "_end_repair",
            prefix + "_repair_step", prefix + "_a_tailing", prefix + "_tail_step",
            prefix + "_adapters", prefix + "_adapter_ligase", prefix + "_ligation",
        ))

        if not wga:
            denaturing = graph.input(prefix + "_denaturing", label + " DNA denaturation condition",
                                     MaterialKind.CONDITION)
            labelled_strand = graph.product(prefix + "_brdu_strand", label + " BrdU-labelled strand",
                                            MaterialKind.MOLECULAR, state=adapted.state)
            partner = graph.product(prefix + "_unlabelled_partner", label + " unlabelled partner strand",
                                    MaterialKind.WASTE)
            graph.react(prefix + "_denature", "Denature " + label + " adapter-ligated DNA",
                        StepKind.STRAND_SEPARATION,
                        ((adapted, InputRole.TEMPLATE),
                         (denaturing, InputRole.CONDITION)), (labelled_strand, partner))
            graph.discard(partner, "the BrdU-containing nascent strand is the selected analyte")
            antibody = graph.input(prefix + "_brdu_antibody", label + " anti-BrdU capture reagent",
                                   MaterialKind.SELECTION)
            selected = graph.product(prefix + "_selected", label + " BrdU-enriched library template",
                                     MaterialKind.MOLECULAR, state=adapted.state)
            unselected = graph.product(prefix + "_unselected", label + " unselected DNA",
                                       MaterialKind.WASTE)
            graph.react(prefix + "_brdu_capture", "Immunoprecipitate BrdU-labelled " + label + " DNA",
                        StepKind.AFFINITY_CAPTURE,
                        ((labelled_strand, InputRole.TEMPLATE),
                         (antibody, InputRole.SELECTION)), (selected, unselected))
            graph.discard(unselected, "not retained by anti-BrdU immunoprecipitation")
            current = selected
            reviewed.extend((
                prefix + "_denaturing", prefix + "_denature", prefix + "_brdu_antibody",
                prefix + "_brdu_capture",
            ))
        else:
            current = adapted

        forward = graph.input(prefix + "_index_forward", label + (" TPC1 primer 5'-AATGATACGGCGACCACCGA*G-3' (* phosphorothioate)" if wga else " indexed PCR forward primer"),
                              MaterialKind.OLIGO,
                              **({"state": Construct([Segment("TPC1", tpc1.sequence, "p5")], name="TPC1")} if wga else {}))
        reverse = graph.input(prefix + "_index_reverse", label + (" TPC2 primer 5'-CAAGCAGAAGACGGCATACGA*G-3' (* phosphorothioate)" if wga else " indexed PCR reverse primer"),
                              MaterialKind.OLIGO,
                              **({"state": Construct([Segment("TPC2", tpc2.sequence, "p7")], name="TPC2")} if wga else {}))
        pcr_polymerase = graph.input(prefix + "_pcr_polymerase", label + (" KAPA HiFi HotStart ReadyMix; 5-7 cycles set by real-time PCR against KK2702 standards (Steps 66-71)" if wga else " PCR polymerase"),
                                     MaterialKind.CATALYST)
        library = graph.product(prefix + "_library", label + " sequencing library",
                                MaterialKind.LIBRARY, state=final)
        dimers = None
        if wga:
            amplified_library = graph.product(prefix + "_pcr_product", label + " amplified library with adaptor dimers",
                                              MaterialKind.MOLECULAR, state=final)
        graph.react(prefix + "_index_pcr", ("Amplify the adaptor-indexed " if wga else "PCR-index the ") + label + " library",
                    StepKind.PCR,
                    ((current, InputRole.TEMPLATE),
                     (forward, InputRole.STRUCTURAL),
                     (reverse, InputRole.STRUCTURAL),
                     (pcr_polymerase, InputRole.CATALYST)), (amplified_library if wga else library,),
                    geometry=(PCRGeometry.exact(
                        final, tpc1, tpc2,
                        statement="TPC1/TPC2 are printed in Miura 2020 (Reagents); their sites are computed on the assumed TruSeq-type P5/P7 adaptor ends, since the NEXTflex adaptor sequence is not printed.")
                        if wga else PCRGeometry.unavailable(
                        final, reason="Sample index choice varies while the final TruSeq endpoint is canonical.")))
        if wga:
            post_pcr = graph.input(prefix + "_post_pcr_ampure", label + " 0.9x AMPure XP (Steps 72-77; repeated if ~100-bp adaptor dimers remain, Step 78)",
                                   MaterialKind.SELECTION)
            dimers = graph.product(prefix + "_dimer_waste", label + " supernatant with adaptor dimers",
                                   MaterialKind.WASTE)
            graph.react(prefix + "_library_cleanup", "Remove adaptor dimers from the " + label + " library",
                        StepKind.CLEANUP,
                        ((amplified_library, InputRole.TEMPLATE), (post_pcr, InputRole.SELECTION)),
                        (library, dimers))
            graph.discard(dimers, "adaptor dimers stay in the 0.9x AMPure supernatant")
            reviewed.extend((prefix + "_post_pcr_ampure", prefix + "_library_cleanup"))
        graph.sequence(library, plan)
        reviewed.extend((
            prefix + "_index_forward", prefix + "_index_reverse",
            prefix + "_pcr_polymerase", prefix + "_index_pcr", prefix + "_library",
        ))

    if key == "el-repli-seq":
        cells = graph.input("cycling_cells", "asynchronously cycling cells",
                            MaterialKind.BIOLOGICAL,
                            state=molecule("asynchronously cycling cells",
                                           ("replicating chromosomes",)), primary=True)
        brdu = graph.input("brdu", "BrdU pulse-label nucleotide",
                           MaterialKind.MOLECULAR)
        replication_polymerase = graph.input("replication_polymerase",
                                             "active cellular DNA polymerases",
                                             MaterialKind.CATALYST)
        labelled = graph.product("brdu_labelled_cells", "BrdU pulse-labelled S-phase cells",
                                 MaterialKind.MOLECULAR,
                                 state=NascentPulse("BrdU", ("early S phase", "late S phase")).labelled_scene())
        graph.react("pulse_label", "Pulse-label nascent DNA with BrdU",
                    StepKind.EXTENSION,
                    ((cells, InputRole.TEMPLATE), (brdu, InputRole.STRUCTURAL),
                     (replication_polymerase, InputRole.CATALYST)), (labelled,))
        gates = graph.input("facs_gates", "early/late S-phase FACS gates",
                            MaterialKind.SELECTION)
        early = graph.product("early_cells", "early-S BrdU-labelled cells",
                              MaterialKind.MOLECULAR, state=labelled.state)
        late = graph.product("late_cells", "late-S BrdU-labelled cells",
                             MaterialKind.MOLECULAR, state=labelled.state)
        graph.react("sort_fractions", "Sort early- and late-S fractions",
                    StepKind.SELECTION,
                    ((labelled, InputRole.TEMPLATE),
                     (gates, InputRole.SELECTION)), (early, late))
        reviewed.extend(("cycling_cells", "brdu", "replication_polymerase",
                         "pulse_label", "facs_gates", "sort_fractions"))
        finish_short_read_branch("early", early, "early-S", wga=False)
        finish_short_read_branch("late", late, "late-S", wga=False)
    else:
        cells = graph.input("sorted_cells", "S-phase cell plus parallel G1 controls",
                            MaterialKind.BIOLOGICAL,
                            state=molecule("flow-sorted cells", ("S-phase cell", "G1 controls")),
                            primary=True)
        s_cell = graph.product("s_phase_cell", "one S-phase cell",
                               MaterialKind.MOLECULAR, state=cells.state)
        g1_cell = graph.product("g1_control", "parallel G1 control cell",
                                MaterialKind.MOLECULAR, state=cells.state)
        graph.react("separate_cells", "Keep S-phase and G1 cells as separate libraries",
                    StepKind.SPLIT, ((cells, InputRole.TEMPLATE),), (s_cell, g1_cell))
        reviewed.extend(("sorted_cells", "separate_cells"))
        finish_short_read_branch("s_phase", s_cell, "S-phase cell", wga=True)
        finish_short_read_branch("g1", g1_cell, "G1 control", wga=True)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = SPEC.FINAL_LIBRARY
    primers = SPEC.SEQ_PRIMERS
    key = 'screpli-seq'
    return _replication_graph(key, d, final, primers)


from types import SimpleNamespace
from conformation_replication_specs import chain, _three_c
from chromosome_capture import ThreeCLibrary
from replication_labeling import CopyNumberReplication

def workflow():
    dosage = CopyNumberReplication("mid-S-phase", "SeqPlex whole-genome amplification")
    amplified = molecule("single-cell WGA product",
                         ("left SEQXE amplification handle", "cell genome copy",
                          "right SEQXE amplification handle"))
    adapted = molecule("single-cell Illumina library",
                       ("P5 / sample index", "WGA genomic fragment", "sample index / P7"),
                       junctions=(("P5 / sample index", "WGA genomic fragment", "adapter ligation"),
                                  ("WGA genomic fragment", "sample index / P7", "adapter ligation")))
    return chain(molecule("one flow-sorted S-phase cell", ("partly replicated genome",)), (
        ("Use copy number as the replication mark", dosage.states(),
         "Replicated domains have twice the DNA dosage of unreplicated domains in the same cell."),
        ("Lyse with proteinase K, heat-fragment, then whole-genome amplify with SeqPlex (SEQXE-50XN)", amplified,
         "Steps 17-25: 55 °C 1 h and 99 °C 4 min in lysis and fragmentation buffer, SeqPlex pre-amplification, then 24 amplification cycles (DOP-PCR-based kit). G1 cells processed in parallel provide the amplification and mappability reference."),
        ("AMPure-purify, digest the SEQXE adaptor with Eco57I, then remove the ~40-bp adaptors with 1.8x AMPure",
         molecule("adapter-removed WGA fragment", ("possible residual SEQXE handle", "WGA genomic fragment",
                                                   "possible residual SEQXE handle")),
         "Steps 26-43: 100-500 bp smear with a ~150 bp peak. Some reads still begin in a SEQXE-like sequence (Step 82), so removal is drawn as possibly incomplete."),
        ("End-repair and add one 3-prime dA to each end (KAPA LTP)",
         da_tailed("dA-tailed single-cell WGA fragment", "WGA genomic fragment"),
         "Steps 44-57: KAPA LTP Library Preparation Kit KK8232, slightly modified to reduce cost (ref. 43), with bead cleanups."),
        ("Ligate a NEXTflex Illumina index adaptor", adapted,
         "Step 58 names NEXTflex DNA Barcodes - 48 but prints no adaptor sequence; the TruSeq-style end is a model assumption."),
        ("Amplify 5-7 cycles with TPC1/TPC2, remove adaptor dimers, sequence 80-bp single reads on HiSeq 1500", adapted,
         "Steps 66-79. The index assigns reads to one individually built cell library; it is not an inline cell barcode."),
    ))

d = DEFINITION
if d.standard_truseq:
    from batch_ngs import truseq_library
    # NEXTflex adaptor sequence is not printed (Step 58): TruSeq-type ends are inferred;
    # the 6-nt i7 window follows the single-index trimming pattern printed in Step 82.
    FINAL_LIBRARY, SEQ_PRIMERS = truseq_library(list(d.final_parts), d.title + " final library",
                                                dual_index=False, inferred_adapters=True, i7_nt=6)
else:
    FINAL_LIBRARY = Construct([role("Illumina-compatible left end"), *d.final_parts, role("Illumina-compatible right end")], name=d.title + " final library")
    SEQ_PRIMERS = ()
OLIGOS = (
    "<b>TPC1 primer</b> — 5′-AATGATACGGCGACCACCGA*G-3′ (* phosphorothioate; HPLC grade; Miura 2020 Reagents, ref. 64)",
    "<b>TPC2 primer</b> — 5′-CAAGCAGAAGACGGCATACGA*G-3′ (* phosphorothioate; HPLC grade; Miura 2020 Reagents, ref. 64)",
)
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES="01_screpli_seq.html",
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=FINAL_LIBRARY,
    SEQ_PRIMERS=SEQ_PRIMERS, FINAL_CAPTION=d.sequencing,
    SEQUENCING_INTRO=d.sequencing, READ_LENGTHS={"Read 1": 80}, ACQUIRED_ROLES=("Read 1", "Index 1 (i7)"),
    SEQUENCING_ENDING=d.sequencing,
    SEQUENCING_UNAVAILABLE=("The defining protocol uses locus-specific or custom primer sequence that cannot be represented by one universal binding sequence."
                            if not SEQ_PRIMERS else ""),
    workflow=workflow, oligos=lambda: OLIGOS,
)
SPEC.protocol_graph = protocol_graph
