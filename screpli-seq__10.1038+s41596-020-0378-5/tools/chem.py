
from conformation_replication_specs import _three_c
from conformation_replication_specs import AdapterLigationGeometry, Construct, Definition, InputRole, MaterialKind, NascentPulse, PCRGeometry, ProtocolGraph, SequencingPlan, StepKind, da_tailed, molecule, role, sp
DEFINITION = Definition(
        "scRepli-seq — single-cell DNA replication sequencing", "10.1038/s41596-020-0378-5",
        "Miura et al., Nature Protocols (2020)",
        "A single S-phase cell is whole-genome amplified and sequenced; twofold copy-number differences distinguish replicated from unreplicated domains.",
        (role("SeqPlex-amplified genomic fragment"),),
        "Each PCR index identifies one single-cell library. Genomic Read 1/Read 2 coverage, normalized to G1 controls, reports replicated and unreplicated domains.", True,
        "The published Supplementary Software 1/2 prints a sample-index read-trimming pattern and a SEQXE motif predicted from example FASTQ reads. The latter is an analysis observation, not a disclosed synthesized SeqPlex reagent sequence; its exact molecular ends remain unresolved here. The later open DOP-PCR protocol uses a different whole-genome-amplification primer and is not substituted for this SeqPlex implementation.")


def _replication_graph(key: str, definition, final: Construct,
                       primers: tuple[sp.SeqPrimer, ...]) -> ProtocolGraph:
    """Typed replication-timing graphs retaining every comparison library."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    reviewed: list[str] = []
    plan = SequencingPlan(final, tuple(primers), tuple(p.role for p in primers))
    adapter_geometry = AdapterLigationGeometry.unavailable(
        "The protocol identifies the supported Illumina library preparation but does not print one immutable adapter-end reagent.")

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
        adapters = graph.input(prefix + "_adapters", label + " Illumina adapters",
                               MaterialKind.ADAPTER, inferred=True)
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
    key = 'screpli-seq'
    return _replication_graph(key, d, final, primers)


from types import SimpleNamespace
from conformation_replication_specs import chain, _three_c
from chromosome_capture import ThreeCLibrary
from replication_labeling import CopyNumberReplication

def workflow():
    dosage = CopyNumberReplication("mid-S-phase", "SeqPlex DOP-PCR whole-genome amplification")
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
        ("Whole-genome amplify each cell by SeqPlex DOP-PCR", amplified,
         "G1 cells processed in parallel provide the amplification and mappability reference."),
        ("Fragment WGA DNA, end-repair and add one 3-prime dA to each end",
         da_tailed("dA-tailed single-cell WGA fragment", "WGA genomic fragment"),
         "The amplified genomic fragment is shown before adapter ligation."),
        ("Ligate Illumina adapters", adapted,
         "The protocol permits KAPA LTP or NEBNext Ultra II library construction."),
        ("PCR-index each single-cell library and sequence", adapted,
         "The index assigns reads to one cell; it is not an inline cell barcode."),
    ))

d = DEFINITION
if d.standard_truseq:
    from batch_ngs import truseq_library
    FINAL_LIBRARY, SEQ_PRIMERS = truseq_library(list(d.final_parts), d.title + " final library", dual_index=False)
else:
    FINAL_LIBRARY = Construct([role("Illumina-compatible left end"), *d.final_parts, role("Illumina-compatible right end")], name=d.title + " final library")
    SEQ_PRIMERS = ()
OLIGOS = ()
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES="01_screpli_seq.html",
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
