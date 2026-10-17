
from conformation_replication_specs import _three_c
import illumina as il
from chemdraw import Segment
from adapter_ligation import reported_adapter_end
from conformation_replication_specs import AdapterLigationGeometry, Construct, Definition, InputRole, MaterialKind, NascentPulse, PCRGeometry, ProtocolGraph, SequencingPlan, StepKind, da_tailed, molecule, role, sp
DEFINITION = Definition(
        "E/L Repli-seq — early/late replication timing", "10.1038/nprot.2017.148",
        "Marchal et al., Nature Protocols (2018)",
        "BrdU-labelled nascent DNA is separated into early- and late-S populations, adapter-ligated, immunoprecipitated and sequenced as two quantitative fractions.",
        (role("BrdU-enriched genomic fragment"),),
        "Separate early- and late-S libraries are pooled with unique i5/i7 pairs and sequenced on a HiSeq; Step 81 says 50 bp single-end reads are generally sufficient (longer or paired-end reads are optional, e.g. for SNP phasing). Replication timing is the normalized early/late log2 coverage ratio.", True,
        "Supplementary Data 1 prints the NEBNext adaptor with 5′ phosphate, internal deoxyU and a terminal phosphorothioate bond, plus the i5/i7 index-primer families. Its hairpin fold (12-bp stem, 3′-T overhang) and USER loop opening are computed from that printed oligo, and the single 3′ dA from the one End Prep reaction is inferred from the kit design; the article prints neither. Intermediate adaptor states are drawn as labelled regions, not base-resolved duplexes.")


# Step 81: "Generally, 50bp single-end reads is sufficient"; the dual-index design
# (unique i5/i7 combination per library, Step 57) needs both index reads.  Index-read
# cycle counts are not printed (the Supplementary Data 1 indexes are 8 nt).
ACQUIRED_ROLES = ("Read 1", "Index 1 (i7)", "Index 2 (i5)")
READ_LENGTHS = {"Read 1": 50}
ARTICLE = "Marchal et al. 2018 Nature Protocols PROCEDURE Steps 1-81; doi:10.1038/nprot.2017.148"
SUPP = "Marchal et al. 2018 publisher Supplementary Data 1/2 and Supplementary Methods 1; doi:10.1038/nprot.2017.148"
# Supplementary Data 1 prints the NEBNext Adaptor for Illumina (#E7337A) as
# 5'-/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3'.
SUPP_DATA1_ADAPTOR = "GATCGGAAGAGCACACGTCTGAACTCCAGTCUACACTCTTTCCCTACACGACGCTCTTCCGATCT"
if SUPP_DATA1_ADAPTOR != il.NEBNEXT_HAIRPIN:
    raise ValueError("Supplementary Data 1 adaptor transcription differs from the NEBNext hairpin model")
_HP = il.NEBNEXT_HAIRPIN_MODEL
_DU = il.NEBNEXT_HAIRPIN_DU_POS


def adaptor_oligo() -> Construct:
    """The Supplementary Data 1 adaptor, split at its computed stem and dU positions."""
    seq, n = SUPP_DATA1_ADAPTOR, _HP.stem_len
    end = len(seq) - _HP.overhang_3_len
    return Construct([
        Segment("5-prime stem", seq[:n], "r2"),
        Segment("Read-2-side loop", seq[n:_DU], "r2"),
        Segment("deoxyU", seq[_DU:_DU + 1], "w1"),
        Segment("Read-1-side loop", seq[_DU + 1:end - n], "r1"),
        Segment("3-prime stem", seq[end - n:end], "r1"),
        Segment("3-prime dT (phosphorothioate-linked)", seq[end:], None),
    ], name="NEBNext Adaptor for Illumina (Supplementary Data 1)")


def _replication_graph(key: str, definition, final: Construct,
                       primers: tuple[sp.SeqPrimer, ...]) -> ProtocolGraph:
    """E/L Repli-seq graph following PROCEDURE Steps 1-81 for each S-phase fraction."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    reviewed: list[str] = []
    plan = SequencingPlan(final, tuple(primers), ACQUIRED_ROLES, read_lengths=READ_LENGTHS)

    def finish_fraction(prefix: str, template, label: str) -> None:
        lysis = graph.input(prefix + "_lysis", label + " SDS-PK buffer with 0.2 mg/ml proteinase K, 56 °C 2 h",
                            MaterialKind.CONDITION, source=ARTICLE)
        extraction = graph.input(prefix + "_dna_extraction", label + " Zymo Quick-DNA Microprep column",
                                 MaterialKind.SELECTION, source=ARTICLE)
        purified = graph.product(prefix + "_purified_dna", label + " purified genomic DNA",
                                 MaterialKind.MOLECULAR,
                                 state=molecule(label + " genomic DNA", ("BrdU-containing nascent strand", "unlabelled strands")))
        extract_waste = graph.product(prefix + "_extraction_waste", label + " column flow-through and washes",
                                      MaterialKind.WASTE, source=ARTICLE)
        graph.react(prefix + "_extract", "Lyse and column-purify " + label + " DNA (Steps 19-24)",
                    StepKind.EXTRACTION,
                    ((template, InputRole.TEMPLATE), (lysis, InputRole.CONDITION),
                     (extraction, InputRole.SELECTION)), (purified, extract_waste), source=ARTICLE)
        graph.discard(extract_waste, "column flow-through and washes are removed")
        sonication = graph.input(prefix + "_sonication", label + " Covaris 175 W, 10% duty, 200 cycles/burst, 120 s (~200 bp)",
                                 MaterialKind.CONDITION, source=ARTICLE)
        fragmented = graph.product(prefix + "_fragments", label + " sheared DNA",
                                   MaterialKind.MOLECULAR, state=purified.state)
        graph.react(prefix + "_fragment", "Sonicate " + label + " DNA (Steps 25-28)",
                    StepKind.FRAGMENTATION,
                    ((purified, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)),
                    (fragmented,), source=ARTICLE)
        end_prep = graph.input(prefix + "_end_prep", label + " NEBNext End Prep Enzyme Mix, 20 °C 30 min then 65 °C 30 min",
                               MaterialKind.CATALYST, source=ARTICLE)
        tailed = graph.product(prefix + "_tailed", label + " end-prepped dA-tailed DNA",
                               MaterialKind.MOLECULAR,
                               state=da_tailed(label + " dA-tailed DNA", label + " fragment"), inferred=True)
        graph.react(prefix + "_end_prep_step", "End-repair and dA-tail in one End Prep reaction (Steps 30-32)",
                    StepKind.TAILING,
                    ((fragmented, InputRole.TEMPLATE), (end_prep, InputRole.CATALYST)), (tailed,),
                    source=ARTICLE,
                    reason="The article runs one NEBNext Ultra End Prep reaction; the single 3-prime dA "
                           "is inferred from the kit's dA-tailing/TA-ligation design, not printed in the article.")
        adapters = graph.input(prefix + "_adapters", label + " NEBNext Adaptor for Illumina, 0.45 µM (1:10 below 100 ng DNA)",
                               MaterialKind.ADAPTER, state=adaptor_oligo(), source=SUPP)
        ligase = graph.input(prefix + "_adapter_ligase", label + " Blunt/TA Ligase Master Mix plus Ligation Enhancer, 20 °C 15 min",
                             MaterialKind.CATALYST, source=ARTICLE)
        closed = graph.product(prefix + "_adapted", label + " hairpin-adaptor-ligated DNA",
                               MaterialKind.MOLECULAR, inferred=True,
                               state=molecule(label + " hairpin-ligated DNA",
                                              ("closed NEBNext hairpin (dU loop)", "nascent DNA fragment",
                                               "closed NEBNext hairpin (dU loop)"),
                                              junctions=(("closed NEBNext hairpin (dU loop)", "nascent DNA fragment", "T/A ligation"),)))
        graph.react(prefix + "_ligation", "Ligate the hairpin adaptor to " + label + " DNA (Steps 33-34)",
                    StepKind.LIGATION,
                    ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                     (ligase, InputRole.CATALYST)), (closed,), source=ARTICLE,
                    geometry=reported_adapter_end(
                        "NEBNext Adaptor for Illumina", kind="cohesive", polarity="3-prime",
                        terminal_sequence=SUPP_DATA1_ADAPTOR[-1], body_sequence=SUPP_DATA1_ADAPTOR[:-1],
                        substrate_scenes=(da_tailed(label + " dA-tailed DNA", label + " fragment"),),
                        statement="Supplementary Data 1 prints the full 5-prime-phosphorylated adaptor with internal deoxyU "
                                  "and a phosphorothioate-linked 3-prime T; its 12-bp terminal self-complementary stem is "
                                  "computed, so the ligating end is a single 3-prime T overhang (fold derived, not printed)."))
        user = graph.input(prefix + "_user", label + " USER enzyme 3 µl, 37 °C 15 min",
                           MaterialKind.CATALYST, source=ARTICLE)
        opened = graph.product(prefix + "_opened", label + " USER-opened adaptor-ligated DNA",
                               MaterialKind.MOLECULAR, inferred=True,
                               state=molecule(label + " opened adaptor library",
                                              ("Read 2-side / Read 1-side opened adaptor arms", "nascent DNA fragment",
                                               "Read 1-side / Read 2-side opened adaptor arms")))
        user_waste = graph.product(prefix + "_user_waste", label + " excised uracil and loop-cleavage residues",
                                   MaterialKind.WASTE)
        graph.react(prefix + "_user_open", "Open the adaptor loops at deoxyU with USER (Steps 35-36)",
                    StepKind.DIGESTION,
                    ((closed, InputRole.TEMPLATE), (user, InputRole.CATALYST)), (opened, user_waste),
                    source=ARTICLE)
        clean1 = graph.input(prefix + "_ligation_column", label + " DNA Clean & Concentrator-5",
                             MaterialKind.SELECTION, source=ARTICLE)
        cleaned = graph.product(prefix + "_cleaned", label + " column-purified adaptor-ligated DNA",
                                MaterialKind.MOLECULAR, state=opened.state)
        clean1_waste = graph.product(prefix + "_ligation_column_waste", label + " ligation reagents, free adaptor and USER residues",
                                     MaterialKind.WASTE)
        graph.react(prefix + "_ligation_cleanup", "Column-purify " + label + " ligation (Step 37)",
                    StepKind.CLEANUP,
                    ((opened, InputRole.TEMPLATE), (clean1, InputRole.SELECTION)), (cleaned, clean1_waste),
                    source=ARTICLE, changes_molecule=False)
        graph.discard(user_waste, "removed with the Step 37 column flow-through")
        graph.discard(clean1_waste, "column flow-through is removed")

        denaturing = graph.input(prefix + "_denaturing", label + " 95 °C 5 min, ice 2 min",
                                 MaterialKind.CONDITION, source=ARTICLE)
        single = graph.product(prefix + "_brdu_strand", label + " BrdU-containing single strand",
                               MaterialKind.MOLECULAR, state=cleaned.state)
        partner = graph.product(prefix + "_unlabelled_strand", label + " unlabelled single strand (kept in the IP input)",
                                MaterialKind.MOLECULAR, state=cleaned.state)
        graph.react(prefix + "_denature", "Heat-denature " + label + " adaptor-ligated DNA (Step 41)",
                    StepKind.STRAND_SEPARATION,
                    ((cleaned, InputRole.TEMPLATE), (denaturing, InputRole.CONDITION)), (single, partner),
                    source=ARTICLE)
        antibody = graph.input(prefix + "_brdu_antibody", label + " mouse anti-BrdU (BD 555627) then rabbit anti-mouse IgG; 16,000 g pellet, IP-buffer wash",
                               MaterialKind.SELECTION, source=ARTICLE)
        selected = graph.product(prefix + "_selected", label + " BrdU-immunoprecipitated strands",
                                 MaterialKind.MOLECULAR, state=cleaned.state)
        unselected = graph.product(prefix + "_unselected", label + " supernatant with unlabelled strands",
                                   MaterialKind.WASTE)
        graph.react(prefix + "_brdu_capture", "Immunoprecipitate BrdU-containing " + label + " strands (Steps 42-51)",
                    StepKind.AFFINITY_CAPTURE,
                    ((single, InputRole.TEMPLATE), (partner, InputRole.TEMPLATE),
                     (antibody, InputRole.SELECTION)), (selected, unselected),
                    source=ARTICLE)
        graph.discard(unselected, "supernatant removed after centrifugation and wash (Steps 48-51)")
        digest = graph.input(prefix + "_ip_digest", label + " proteinase K digestion overnight 37 °C then 56 °C 60 min",
                             MaterialKind.CONDITION, source=ARTICLE)
        clean2 = graph.input(prefix + "_ip_column", label + " DNA Clean & Concentrator-5 (16 µl)",
                             MaterialKind.SELECTION, source=ARTICLE)
        ip_dna = graph.product(prefix + "_ip_dna", label + " purified BrdU-IP library template",
                               MaterialKind.MOLECULAR, state=cleaned.state)
        ip_waste = graph.product(prefix + "_ip_column_waste", label + " digested antibody and column flow-through",
                                 MaterialKind.WASTE)
        graph.react(prefix + "_ip_cleanup", "Release and column-purify " + label + " IP DNA (Steps 52-55)",
                    StepKind.CLEANUP,
                    ((selected, InputRole.TEMPLATE), (digest, InputRole.CONDITION),
                     (clean2, InputRole.SELECTION)), (ip_dna, ip_waste),
                    source=ARTICLE, changes_molecule=False)
        graph.discard(ip_waste, "column flow-through is removed")

        forward = graph.input(prefix + "_index_forward", label + " NEBNext i5 primer (i501-i508, Supplementary Data 1)",
                              MaterialKind.OLIGO, source=SUPP)
        reverse = graph.input(prefix + "_index_reverse", label + " NEBNext i7 primer (i701-i712, Supplementary Data 1)",
                              MaterialKind.OLIGO, source=SUPP)
        pcr_polymerase = graph.input(prefix + "_pcr_polymerase", label + " NEBNext Q5 Hot Start HiFi Master Mix; 14 cycles 98 °C 10 s / 65 °C 75 s",
                                     MaterialKind.CATALYST, source=ARTICLE)
        amplified = graph.product(prefix + "_amplified", label + " indexed PCR product",
                                  MaterialKind.MOLECULAR, state=final)
        graph.react(prefix + "_index_pcr", "Dual-index PCR of the " + label + " IP DNA (Steps 57-58)",
                    StepKind.PCR,
                    ((ip_dna, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                     (reverse, InputRole.STRUCTURAL), (pcr_polymerase, InputRole.CATALYST)), (amplified,),
                    source=ARTICLE,
                    geometry=PCRGeometry.unavailable(
                        final, reason="Each library gets a unique i5/i7 combination; the article's 14 cycles "
                                      "are a suggestion adjusted to IP yield, with Box 1 re-amplification as fallback."))
        beads = graph.input(prefix + "_ampure", label + " AMPure XP 90 µl per 100 µl (0.9x), three 80% ethanol washes",
                            MaterialKind.SELECTION, source=ARTICLE)
        library = graph.product(prefix + "_library", label + " sequencing library",
                                MaterialKind.LIBRARY, state=final)
        bead_waste = graph.product(prefix + "_ampure_supernatant", label + " bead supernatant with primers and primer dimers",
                                   MaterialKind.WASTE)
        graph.react(prefix + "_pcr_cleanup", "AMPure-purify the " + label + " library (Steps 60-72)",
                    StepKind.CLEANUP,
                    ((amplified, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (library, bead_waste),
                    source=ARTICLE, changes_molecule=False)
        graph.discard(bead_waste, "supernatant removed from magnetically held beads")
        graph.sequence(library, plan)
        reviewed.extend((
            prefix + "_lysis", prefix + "_dna_extraction", prefix + "_extract",
            prefix + "_sonication", prefix + "_fragment", prefix + "_end_prep",
            prefix + "_end_prep_step", prefix + "_adapters", prefix + "_adapter_ligase",
            prefix + "_ligation", prefix + "_user", prefix + "_user_open",
            prefix + "_ligation_column", prefix + "_ligation_cleanup",
            prefix + "_denaturing", prefix + "_denature", prefix + "_brdu_antibody",
            prefix + "_brdu_capture", prefix + "_ip_digest", prefix + "_ip_column",
            prefix + "_ip_cleanup", prefix + "_index_forward", prefix + "_index_reverse",
            prefix + "_pcr_polymerase", prefix + "_index_pcr", prefix + "_ampure",
            prefix + "_pcr_cleanup", prefix + "_library",
        ))

    cells = graph.input("cycling_cells", "asynchronously cycling cells",
                        MaterialKind.BIOLOGICAL,
                        state=molecule("asynchronously cycling cells",
                                       ("replicating chromosomes",)), primary=True)
    brdu = graph.input("brdu", "BrdU, 100 µM for 2 h at 37 °C", MaterialKind.MOLECULAR, source=ARTICLE)
    replication_polymerase = graph.input("replication_polymerase",
                                         "active cellular DNA polymerases",
                                         MaterialKind.CATALYST)
    labelled = graph.product("brdu_labelled_cells", "BrdU pulse-labelled S-phase cells",
                             MaterialKind.MOLECULAR,
                             state=NascentPulse("BrdU", ("early S phase", "late S phase")).labelled_scene())
    graph.react("pulse_label", "Pulse-label nascent DNA with BrdU (Steps 1-2)",
                StepKind.EXTENSION,
                ((cells, InputRole.TEMPLATE), (brdu, InputRole.STRUCTURAL),
                 (replication_polymerase, InputRole.CATALYST)), (labelled,), source=ARTICLE)
    gates = graph.input("facs_gates", "75% ethanol fixation, PI/RNase A staining, early/late S-phase FACS gates",
                        MaterialKind.SELECTION, source=ARTICLE)
    early = graph.product("early_cells", "early-S BrdU-labelled cells",
                          MaterialKind.MOLECULAR, state=labelled.state)
    late = graph.product("late_cells", "late-S BrdU-labelled cells",
                         MaterialKind.MOLECULAR, state=labelled.state)
    graph.react("sort_fractions", "Fix, PI-stain and sort early- and late-S fractions (Steps 3-18)",
                StepKind.SELECTION,
                ((labelled, InputRole.TEMPLATE),
                 (gates, InputRole.SELECTION)), (early, late), source=ARTICLE)
    reviewed.extend(("cycling_cells", "brdu", "replication_polymerase",
                     "pulse_label", "facs_gates", "sort_fractions"))
    finish_fraction("early", early, "early-S")
    finish_fraction("late", late, "late-S")
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
    SEQUENCING_INTRO=d.sequencing,
    SEQUENCING_ENDING=d.sequencing,
    SEQUENCING_UNAVAILABLE=("The defining protocol uses locus-specific or custom primer sequence that cannot be represented by one universal binding sequence."
                            if not SEQ_PRIMERS else ""),
    workflow=workflow, oligos=lambda: OLIGOS,
)
SPEC.protocol_graph = protocol_graph
