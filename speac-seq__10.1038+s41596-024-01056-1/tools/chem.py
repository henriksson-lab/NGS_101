"""Protocol-owned source, molecular workflow and reaction graph."""
from chemdraw import Construct, Scene, Segment, feature
from pcr import PCRGeometry
from protocol_graph import (AnnotationSlot, Disposition, Highlight, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, SourceClaim, SourceDiscrepancy,
                            StepKind, Terminal)
from protocol_model import chain, molecule, role, role_library
from extended_protocols import Definition, SAMPLE, GUIDE

DEFINITION = Definition('SPEAC-seq', '10.1038/s41596-024-01056-1', 'Faust Akl et al., Nature Protocols (2024)', 'CRISPR-perturbed and reporter cells are paired and cultured in droplets; functional fluorescence sorting precedes guide-cassette PCR and sequencing.', 'Guide counts are compared between functionally sorted and control droplet fractions; this is a pooled screen, not a transcriptome library.', '01_speac_seq.html')
PAYLOAD = (GUIDE, role("guide-flanking amplicon"))

# Verbatim from the defining protocol's "Oligonucleotides" list (Nature Protocols
# author manuscript PMC11805652, extracted lines 153-175; primers synthesized at IDT).
P5_STAGGERED = (
    ("0 nt stagger", "ATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTGATGTCCACGAGGTCTCT"),
    ("1 nt stagger", "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCGATGTCCACGAGGTCTCT"),
    ("2 nt stagger", "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTGCGATGTCCACGAGGTCTCT"),
    ("3 nt stagger", "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTAGCGATGTCCACGAGGTCTCT"),
    ("4 nt stagger", "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCAACGATGTCCACGAGGTCTCT"),
    ("5 nt stagger", "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTTGCACCGATGTCCACGAGGTCTCT"),
    ("6 nt stagger", "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTACGCAACGATGTCCACGAGGTCTCT"),
    ("7 nt stagger", "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTGAAGACCCGATGTCCACGAGGTCTCT"),
    ("8 nt stagger", "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTGAAGACCCTTGTGGAAAGGACGAAACACCG"),
)
P7_INDEX_PRIMER = "CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCCAATTCCCACTCCTTTCAAGACCT"
# Analysis anchor printed in the same protocol (line 703): Cutadapt -g ACACCG ... GTTTTAG.
GRNA_TRIM_ANCHORS = ("ACACCG", "GTTTTAG")

TITLE = DEFINITION.title
NOTES = DEFINITION.notes
SOURCE = f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.citation}</a>.'
SUMMARY = DEFINITION.summary
CAVEAT = ("The lentiCRISPRv2 vector sequence between the two primer sites is not printed in "
          "the defining protocol, so the amplicon interior and the final library are drawn as "
          "roles; the two ordering primers themselves are verbatim.")
FINAL_LIBRARY = role_library(TITLE + " final library roles", PAYLOAD, index=SAMPLE)
SEQ_PRIMERS = ()
FINAL_CAPTION = DEFINITION.endpoint
SEQUENCING_INTRO = DEFINITION.endpoint
READ_LENGTHS = None
SEQUENCING_ENDING = DEFINITION.endpoint
SEQUENCING_UNAVAILABLE = ("The P5 staggered primer mix carries the TruSeq Read 1 sequence and the P7 index "
                          "primer carries an eight-base index, but the defining protocol states only a depth "
                          "target (at least 100 reads per gRNA) — no instrument, read lengths or cycle plan — "
                          "so no read geometry is asserted on the role-level library.")


def workflow():
    return chain(molecule("pooled CRISPR-perturbed effector cells", (GUIDE, role("reporter cells"))), (
        ("Pair one perturbed cell with one reporter cell in a picolitre droplet",
         molecule("co-culture droplet", (GUIDE, role("perturbed cell"), role("fluorescent reporter cell"))),
         "Different dyes gate one-to-one pairing before the interaction develops."),
        ("Sort droplets by the functional fluorescent response",
         molecule("phenotype-enriched droplet fraction", (GUIDE, role("reporter state"))),
         "The screen's phenotype is a sorted population, not a single-cell transcriptome."),
        ("Break droplets with 20% PFO, recover genomic DNA and PCR-amplify the integrated guide cassette",
         molecule("SPEAC-seq library", (GUIDE, SAMPLE, role("guide-flanking amplicon"))),
         "Guide counts in response-positive and control fractions identify interaction regulators."),
    ))


def oligos():
    rows = ["<b>P5 staggered primer mix for Illumina sequencing</b> (IDT; one species per stagger, pooled at 100 µM)"]
    rows += [f"{name} — 5′-{seq}-3′" for name, seq in P5_STAGGERED]
    rows.append("<b>P7 index primer for lentiCRISPRv2</b> (unique per reaction; N<sub>8</sub> is the sample index)")
    rows.append(f"5′-{P7_INDEX_PRIMER}-3′")
    return tuple(rows)


def _p5_mix_state():
    """The ordering mix is nine separate species, each transcribed whole."""
    sc = Scene()
    for name, seq in P5_STAGGERED:
        # The "0 nt" row is printed with internal spaces in the author manuscript, so
        # reading it as one continuous 75-mer is this model's transcription, not a
        # verbatim string; it is therefore drawn and checked as inferred.
        spaced = name.startswith("0 ")
        sc.strand(name, [Segment(name, seq, "r1", inferred=spaced)], label=name)
        sc.labels(name)
    return sc


def _p7_state():
    return Construct((
        Segment("P7 and i7 landing", "CAAGCAGAAGACGGCATACGAGAT"),
        Segment("eight-base sample index, unassigned per reaction", "N" * 8, "cbc",
                length_bp=8,
                feature=feature("speac_sample_index", "sample_index", "unknown")),
        Segment("TruSeq Read 2 and vector-annealing 3′ end",
                "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCCAATTCCCACTCCTTTCAAGACCT"),
    ), name="P7 index primer")


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input(
        "cell_pools", "CRISPR-perturbed effector cells and fluorescent reporter cells",
        MaterialKind.BIOLOGICAL,
        state=molecule("pooled screen cells", (GUIDE, role("effector cells"),
                                                role("reporter cells"))), primary=True)
    droplets = graph.input("pairing_condition", "picolitre one-to-one droplet pairing condition",
                           MaterialKind.CONDITION)
    paired = graph.product("paired_droplets", "one-effector/one-reporter droplets",
                           MaterialKind.MOLECULAR,
                           state=molecule("co-culture droplet", (
                               GUIDE, role("perturbed cell"), role("fluorescent reporter cell"))))
    pairing_failures = graph.product("pairing_failures", "empty or incorrectly paired droplets",
                                     MaterialKind.WASTE)
    graph.react("pair_cells", "Pair one perturbed and one reporter cell per droplet",
                StepKind.PARTITION,
                ((cells, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (paired, pairing_failures))
    graph.discard(pairing_failures, "does not contain the required one-to-one cell pair")
    response_gate = graph.input("response_gate", "functional fluorescence sort gates",
                                MaterialKind.SELECTION)
    positive = graph.product("positive_droplets", "response-positive droplet fraction",
                             MaterialKind.MOLECULAR, state=paired.state)
    control = graph.product("control_droplets", "response-control droplet fraction",
                            MaterialKind.MOLECULAR, state=paired.state)
    graph.react("sort_response", "Sort droplets by functional response",
                StepKind.SELECTION,
                ((paired, InputRole.TEMPLATE), (response_gate, InputRole.SELECTION)),
                (positive, control))
    reviewed = ["cell_pools", "pairing_condition", "pair_cells", "response_gate",
                "sort_response"]

    def recover(prefix: str, template, label: str):
        pfo = graph.input(prefix + "_pfo", "20% (vol/vol) PFO in HFE-7500 droplet-breaking reagent",
                          MaterialKind.CONDITION)
        aqueous = graph.product(prefix + "_aqueous", label + " recovered aqueous cell phase",
                                MaterialKind.MOLECULAR,
                                state=molecule(label + " broken-droplet cells",
                                               (GUIDE, role("sorted cells"))))
        oil = graph.product(prefix + "_oil", label + " fluorinated oil phase",
                            MaterialKind.WASTE)
        graph.react(prefix + "_break", "Break the thawed emulsion with 20% PFO and take the aqueous layer",
                    StepKind.SOLUBILIZATION,
                    ((template, InputRole.TEMPLATE), (pfo, InputRole.CONDITION)),
                    (aqueous, oil))
        graph.discard(oil, "the fluorinated phase is separated from the aqueous cells by a 1,000g spin")
        extraction = graph.input(prefix + "_extraction",
                                 "Qiagen Blood and Tissue DNA isolation kit (two 200 µL elutions)",
                                 MaterialKind.CONDITION)
        dna = graph.product(prefix + "_dna", label + " integrated guide-cassette gDNA",
                            MaterialKind.MOLECULAR,
                            state=molecule(label + " guide DNA", (
                                GUIDE, role("guide-flanking genomic DNA"))))
        graph.react(prefix + "_recover", "Isolate genomic DNA from the recovered cells",
                    StepKind.EXTRACTION,
                    ((aqueous, InputRole.TEMPLATE), (extraction, InputRole.CONDITION)),
                    (dna,))
        beads = graph.input(prefix + "_precleanup", "AMPure XP beads at 2× beads-to-sample, 40 µL elution",
                            MaterialKind.CONDITION)
        concentrated = graph.product(prefix + "_dna_clean", label + " concentrated gDNA template",
                                     MaterialKind.MOLECULAR, state=dna.state)
        preclean_waste = graph.product(prefix + "_preclean_waste",
                                       label + " unbound gDNA cleanup supernatant", MaterialKind.WASTE)
        graph.react(prefix + "_concentrate", "Concentrate the eluted gDNA on AMPure XP beads",
                    StepKind.CLEANUP,
                    ((dna, InputRole.TEMPLATE), (beads, InputRole.CONDITION)),
                    (concentrated, preclean_waste))
        graph.discard(preclean_waste, "removed with the cleanup supernatant")
        reviewed.extend((prefix + "_pfo", prefix + "_break", prefix + "_oil",
                         prefix + "_extraction", prefix + "_recover",
                         prefix + "_precleanup", prefix + "_concentrate"))

        forward = graph.input(prefix + "_forward", "P5 staggered primer mix (0.5 µM final)",
                              MaterialKind.OLIGO, state=_p5_mix_state())
        reverse = graph.input(prefix + "_reverse",
                              "P7 index primer for lentiCRISPRv2, unique per reaction (0.5 µM final)",
                              MaterialKind.OLIGO, state=_p7_state())
        polymerase = graph.input(prefix + "_polymerase",
                                 "TaKaRa ExTaq polymerase, 75 mU per 100 µL, 28–35 cycles of 95/53/72 °C",
                                 MaterialKind.CATALYST)
        amplicon = graph.product(prefix + "_amplicon", label + " indexed guide amplicon",
                                 MaterialKind.LIBRARY, state=final)
        graph.react(prefix + "_pcr", "Amplify and index the integrated gRNA cassette", StepKind.PCR,
                    ((concentrated, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                     (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                    (amplicon,), geometry=PCRGeometry.unavailable(
                        final,
                        reason=("Both ordering primers are published verbatim, but the defining protocol "
                                "does not print the lentiCRISPRv2 sequence they anneal to, so the exact "
                                "amplicon ends cannot be located by computation.")))
        postbeads = graph.input(prefix + "_postcleanup",
                                "AMPure XP beads at 2× beads-to-sample, 40 µL elution",
                                MaterialKind.CONDITION)
        library = graph.product(prefix + "_library", label + " purified indexed library",
                                MaterialKind.LIBRARY, state=final)
        pcr_waste = graph.product(prefix + "_pcr_waste",
                                  label + " unincorporated primers and PCR buffer", MaterialKind.WASTE)
        graph.react(prefix + "_cleanup", "Purify the amplified library on AMPure XP beads",
                    StepKind.CLEANUP,
                    ((amplicon, InputRole.TEMPLATE), (postbeads, InputRole.CONDITION)),
                    (library, pcr_waste))
        graph.discard(pcr_waste, "removed with the AMPure supernatant before quantification")
        graph.sequence(library, SequencingPlan(
            limitation=(label + ": " + d.endpoint + " Libraries are checked on a Bioanalyzer high-sensitivity "
                        "chip (about 285 bp average for the lentiCRISPRv2 Brie library), quantified with a KAPA "
                        "library quantification kit, pooled at 4 nM each and sequenced to at least 100 reads "
                        "per gRNA; the protocol names no instrument or read length.")))
        reviewed.extend((prefix + "_forward", prefix + "_reverse", prefix + "_polymerase",
                         prefix + "_pcr", prefix + "_amplicon", prefix + "_postcleanup",
                         prefix + "_cleanup", prefix + "_library", prefix + "_pcr_waste"))

    recover("positive", positive, "response-positive")
    recover("control", control, "control")
    graph.annotate(Highlight(
        "positive_pcr", "One PCR installs both Illumina ends",
        "There is no separate adapter ligation or second indexing PCR: the P5 staggered mix adds P5 plus the "
        "TruSeq Read 1 sequence and a variable stagger, and the per-reaction P7 index primer adds the eight-base "
        "sample index plus P7, directly onto the integrated lentiCRISPRv2 cassette amplified from genomic DNA.",
        slot=AnnotationSlot.AFTER_WORKFLOW))
    graph.annotate(Highlight(
        "positive_forward", "The stagger is a phasing device, not an identifier",
        "The mix exists to offset the gRNA position between clusters on a low-diversity amplicon. The model "
        "attaches no identifier feature to it; only the P7 N<sub>8</sub> window carries a sample index.",
        slot=AnnotationSlot.AFTER_WORKFLOW))
    graph.annotate(SourceDiscrepancy(
        "positive_forward", "Printed stagger labels do not match the inserted bases",
        (SourceClaim("Nature Protocols author manuscript PMC11805652, Oligonucleotides list, extracted lines 155–171",
                     "Nine species are labelled 0–8 nt stagger, but relative to the shared "
                     "…CTTCCGATCT end the inserts are 0, 0, G, AG, CAA, TGCAC, ACGCAA, GAAGACCC bases for labels "
                     "0–7, the \"0 nt\" row is printed with internal spaces and lacks the first A of P5, and the \"8 nt\" row ends in "
                     "TTGTGGAAAGGACGAAACACCG instead of the CGATGTCCACGAGGTCTCT 3′ end shared by the others"),
         SourceClaim("Same protocol, gRNA amplification PCR table (extracted line 657)",
                     "The nine species enter the reaction only as one pooled \"P5 staggered primer mix\", "
                     "100 µM stock, 0.5 µM final")),
        "Transcribe all nine species verbatim and model the forward input as the pooled mix with an unassigned stagger.",
        "The primer block was read in the author-manuscript full text; the differences are visible between adjacent printed rows.",
        "Whether the labels, the missing leading A or the divergent 3′ end of the last species are typesetting "
        "losses in this author manuscript or real ordering sequences is unresolved; the publisher's typeset "
        "Nature Protocols version was not inspected, and no base is repaired here."))
    graph.annotate(SourceDiscrepancy(
        "positive_pcr", "P7 index primer stock concentration is internally inconsistent",
        (SourceClaim("Nature Protocols author manuscript PMC11805652, gRNA amplification PCR table, extracted line 659",
                     "\"P7 index primer (unique per reaction) — 5 mM — 0.5 µM — 10 µL\" per 100 µL reaction"),
         SourceClaim("Same table, extracted line 657 (the P5 row)",
                     "0.5 µL of a 100 µM stock gives the stated 0.5 µM final, so the table's own arithmetic "
                     "would need a 5 µM P7 stock for 10 µL to reach 0.5 µM")),
        "Model the P7 index primer at the stated 0.5 µM final working concentration and leave the stock unassigned.",
        "Both rows of the printed reaction table were read in the author-manuscript full text.",
        "Whether \"5 mM\" is a unit typo for 5 µM or a different intended stock is unresolved; the publisher's "
        "typeset version was not inspected and no value is silently repaired."))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
