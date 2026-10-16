
from conformation_replication_specs import _three_c
from conformation_replication_specs import AdapterLigationGeometry, Construct, Definition, InputRole, MaterialKind, NascentPulse, PCRGeometry, ProtocolGraph, SequencingPlan, StepKind, da_tailed, molecule, role, sp
DEFINITION = Definition(
        "E/L Repli-seq — early/late replication timing", "10.1038/nprot.2017.148",
        "Marchal et al., Nature Protocols (2018)",
        "BrdU-labelled nascent DNA is separated into early- and late-S populations, adapter-ligated, immunoprecipitated and sequenced as two quantitative fractions.",
        (role("BrdU-enriched genomic fragment"),),
        "Ordinary Illumina genomic reads measure coverage in separate early- and late-S libraries; replication timing is their normalized log2 ratio.", True,
        "The published Supplementary Data 1 prints the NEBNext adaptor with 5′ phosphate, internal deoxyU and a terminal phosphorothioate bond, plus the i5/i7 index-primer families. Their source-printed sequences are shown below; exact annealed adaptor-end geometry is still under review.")


def _replication_graph(key: str, definition, final: Construct,
                       primers: tuple[sp.SeqPrimer, ...]) -> ProtocolGraph:
    """Typed replication-timing graphs retaining every comparison library."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    reviewed: list[str] = []
    plan = SequencingPlan(final, tuple(primers), tuple(p.role for p in primers))
    adapter_geometry = AdapterLigationGeometry.unavailable(
        "Supplementary Data 1 prints the NEBNext adaptor oligo and index-primer families, but the annealed ligation-end geometry has not yet been checked from the original source.")

    def finish_short_read_branch(prefix: str, template, label: str, *, wga: bool) -> None:
        current = template
        if wga:
            dop_forward = graph.input(prefix + "_dop_forward", label + " DOP-PCR forward primer",
                                      MaterialKind.OLIGO)
            dop_reverse = graph.input(prefix + "_dop_reverse", label + " DOP-PCR reverse primer",
                                      MaterialKind.OLIGO)
            wga_polymerase = graph.input(prefix + "_wga_polymerase", label + " WGA polymerase",
                                         MaterialKind.CATALYST)
            amplified = graph.product(prefix + "_wga", label + " SeqPlex WGA product",
                                      MaterialKind.MOLECULAR,
                                      state=molecule(label + " WGA product", (
                                          "left SEQXE amplification handle",
                                          "genome copy",
                                          "right SEQXE amplification handle")))
            graph.react(prefix + "_wga_step", "Whole-genome amplify " + label,
                        StepKind.PCR,
                        ((current, InputRole.TEMPLATE),
                         (dop_forward, InputRole.STRUCTURAL),
                         (dop_reverse, InputRole.STRUCTURAL),
                         (wga_polymerase, InputRole.CATALYST)), (amplified,),
                        geometry=PCRGeometry.unavailable(
                            final, reason="SeqPlex uses a degenerate DOP-PCR family rather than one fixed amplicon."))
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

        sonication = graph.input(prefix + "_sonication", label + " sonication condition",
                                 MaterialKind.CONDITION)
        fragmented = graph.product(prefix + "_fragments", label + " DNA fragments",
                                   MaterialKind.MOLECULAR, state=current.state)
        graph.react(prefix + "_fragment", "Fragment " + label + " DNA",
                    StepKind.FRAGMENTATION,
                    ((current, InputRole.TEMPLATE),
                     (sonication, InputRole.CONDITION)), (fragmented,))
        repair = graph.input(prefix + "_end_repair", label + " end-repair enzymes",
                             MaterialKind.CATALYST)
        repaired = graph.product(prefix + "_repaired", label + " end-repaired DNA",
                                 MaterialKind.MOLECULAR, state=fragmented.state)
        graph.react(prefix + "_repair_step", "End-repair " + label + " DNA",
                    StepKind.END_REPAIR,
                    ((fragmented, InputRole.TEMPLATE),
                     (repair, InputRole.CATALYST)), (repaired,))
        tailing = graph.input(prefix + "_a_tailing", label + " dA-tailing enzyme",
                              MaterialKind.CATALYST)
        tailed = graph.product(prefix + "_tailed", label + " dA-tailed DNA",
                               MaterialKind.MOLECULAR,
                               state=da_tailed(label + " dA-tailed DNA", label + " fragment"))
        graph.react(prefix + "_tail_step", "Add one 3-prime dA to each end",
                    StepKind.TAILING,
                    ((repaired, InputRole.TEMPLATE),
                     (tailing, InputRole.CATALYST)), (tailed,))
        adapters = graph.input(prefix + "_adapters", label + " NEBNext adaptor printed in Supplementary Data 1",
                               MaterialKind.ADAPTER,
                               source="Marchal et al. Supplementary Data 1; doi:10.1038/nprot.2017.148")
        ligase = graph.input(prefix + "_adapter_ligase", label + " DNA ligase",
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

        forward = graph.input(prefix + "_index_forward", label + " indexed PCR forward primer",
                              MaterialKind.OLIGO)
        reverse = graph.input(prefix + "_index_reverse", label + " indexed PCR reverse primer",
                              MaterialKind.OLIGO)
        pcr_polymerase = graph.input(prefix + "_pcr_polymerase", label + " PCR polymerase",
                                     MaterialKind.CATALYST)
        library = graph.product(prefix + "_library", label + " sequencing library",
                                MaterialKind.LIBRARY, state=final)
        graph.react(prefix + "_index_pcr", "PCR-index the " + label + " library",
                    StepKind.PCR,
                    ((current, InputRole.TEMPLATE),
                     (forward, InputRole.STRUCTURAL),
                     (reverse, InputRole.STRUCTURAL),
                     (pcr_polymerase, InputRole.CATALYST)), (library,),
                    geometry=PCRGeometry.unavailable(
                        final, reason="Sample index choice varies while the final TruSeq endpoint is canonical."))
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
    key = 'el-repli-seq'
    return _replication_graph(key, d, final, primers)


from types import SimpleNamespace
from conformation_replication_specs import chain, _three_c
from chromosome_capture import ThreeCLibrary
from replication_labeling import CopyNumberReplication

def workflow():
    pulse = NascentPulse("BrdU", ("early S phase", "late S phase"))
    adapted = molecule("adapter-ligated nascent DNA",
                       ("left Illumina adapter", "BrdU-containing nascent DNA", "right Illumina adapter"),
                       junctions=(("left Illumina adapter", "BrdU-containing nascent DNA", "adapter ligation"),
                                  ("BrdU-containing nascent DNA", "right Illumina adapter", "adapter ligation")))
    return chain(molecule("asynchronously cycling cells", ("replicating chromosomes",)), (
        ("Pulse-label nascent DNA with BrdU", pulse.labelled_scene(),
         "Only DNA synthesized during the two-hour pulse contains BrdU."),
        ("Fix and sort early- and late-S fractions by DNA content", pulse.sorted_scene(),
         "Early and late fractions remain separate libraries."),
        ("Purify, sonicate, end-repair and ligate universal Illumina adapters", adapted,
         "Adapters are installed before immunoprecipitation."),
        ("Denature, immunoprecipitate BrdU DNA, then index and PCR-amplify", adapted,
         "Sequencing coverage is compared as an early/late ratio."),
    ))

d = DEFINITION
if d.standard_truseq:
    from batch_ngs import truseq_library
    FINAL_LIBRARY, SEQ_PRIMERS = truseq_library(list(d.final_parts), d.title + " final library", dual_index=True)
else:
    FINAL_LIBRARY = Construct([role("Illumina-compatible left end"), *d.final_parts, role("Illumina-compatible right end")], name=d.title + " final library")
    SEQ_PRIMERS = ()
OLIGOS = (
    "<b>NEBNext Adaptor for Illumina (Supplementary Data 1)</b> — 5′-/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3′ (*: terminal phosphorothioate linkage)",
    "<b>BrdU-IP adaptor qPCR primers (Supplementary Methods)</b> — F 5′-ACACTCTTTCCCTACACGACGC-3′; R 5′-GACTGGAGTTCAGACGTGTGC-3′",
    "<b>Re-amplification primers (Supplementary Data 2)</b> — TS-Oligo 1 5′-AATGATACGGCGACCACCGAGA-3′; TS-Oligo 2 5′-CAAGCAGAAGACGGCATACGAG-3′",
)
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES="01_el_repli_seq.html",
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=FINAL_LIBRARY,
    SEQ_PRIMERS=SEQ_PRIMERS, FINAL_CAPTION=d.sequencing,
    SEQUENCING_INTRO=d.sequencing, READ_LENGTHS=None,
    SEQUENCING_ENDING=d.sequencing,
    SEQUENCING_UNAVAILABLE=("The defining protocol uses locus-specific or custom primer sequence that cannot be represented by one universal binding sequence."
                            if not SEQ_PRIMERS else ""),
    workflow=workflow, oligos=lambda: OLIGOS,
)
SPEC.protocol_graph = protocol_graph
