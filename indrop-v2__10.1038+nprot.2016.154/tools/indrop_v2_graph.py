"""Zilionis et al. inDrop v2 reaction path and source-owned oligo inputs."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "indrop-v1__10.1016+j.cell.2015.04.044" / "tools"),
                str(HERE.parents[1] / "lib")]

import illumina as il
import indrop as I
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind,
                            Terminal)

SOURCE = "Zilionis et al. 2017, doi:10.1038/nprot.2016.154; Supplementary Table 1"
MAIN = ("Zilionis et al. 2017 Nature Protocols main text (publisher PDF), Reagent "
        "Setup and Procedure Steps 30-166, Fig. 3 and Table 1")
QC_SOURCE = ("Zilionis et al. 2017 Supplementary Method 2 (MOESM454), procedure step 9 "
             "and Anticipated result")
# Main text Table 1 lists instrument-specific cycle examples (Read 1 / index / Read 2);
# Step 166 requires Read 2 (cell barcode + UMI) of at least 50 bp.
TABLE1_CYCLES = (
    ("MiSeq, multiplexed", "MiSeq Reagent Kit v3 (150 cycles)", 89, 6, 55),
    ("NextSeq, multiplexed", "NextSeq 500/550 High-Output v2 (75 cycles)", 37, 6, 51),
    ("NextSeq, not multiplexed", "NextSeq 500/550 High-Output v2 (75 cycles)", 43, 0, 51),
    ("HiSeq 2500, multiplexed", "HiSeq SBS Kit v4 (2 x 50 cycles)", 61, 6, 51),
)
READ_LIMITATION = ("Acquired cycles are instrument-dependent: main text Table 1 gives four "
                   "examples and Step 166 requires a Read 2 of at least 50 bp; no single "
                   "configuration is modelled.")
PE2_N6 = I.PE2 + "N" * 6
PE2_PCR = il.P5 + "GGTC" + I.PE2[:18]
PE1_INDEX1 = il.P7 + "CGTGAT" + I.PE1[:16]


def protocol_graph() -> ProtocolGraph:
    """Trace v2 through bead synthesis, droplets, IVT, conversion and indexed PCR."""
    graph = ProtocolGraph(I.V2, authority=SOURCE)
    seed = graph.input("bead_seed", "acrydite photocleavable bead primer",
                       MaterialKind.OLIGO, state=I.SOURCE_ACRYDITE_BASES, source=SOURCE)
    plate = graph.input("barcode_plate_oligos", "two split-pool barcode plate oligo sets",
                        MaterialKind.OLIGO, source=SOURCE)
    bead_pol = graph.input("bead_polymerase", "Bst 2.0 DNA polymerase (60 °C; Steps 43-47)",
                           MaterialKind.CATALYST, source=MAIN)
    extended = graph.product("extended_beads", "bead primers after two split-pool extensions "
                             "(full-length plus round-1-only and unextended species)",
                             MaterialKind.MOLECULAR, source=QC_SOURCE)
    graph.react("bead_synthesis", "Build the bead barcode in two split-pool Bst 2.0 "
                "extensions, each stopped (STOP-25) and denatured with NaOH (Steps 30-63)",
                StepKind.EXTENSION,
                ((seed, InputRole.TEMPLATE), (plate, InputRole.STRUCTURAL),
                 (bead_pol, InputRole.CATALYST)), (extended,), source=MAIN)
    hinfi_bead = graph.input("bead_hinfi", "FastDigest HinfI (37 °C, 2 h; Step 66)",
                             MaterialKind.CATALYST, source=MAIN)
    digested = graph.product("hinfi_treated_beads", "HinfI-treated bead primers",
                             MaterialKind.MOLECULAR, source=MAIN)
    graph.react("bead_hinfi_digest", "HinfI digestion of barcoded beads (purpose not stated "
                "for this step)", StepKind.DIGESTION,
                ((extended, InputRole.TEMPLATE), (hinfi_bead, InputRole.CATALYST)),
                (digested,), source=MAIN)
    ba19 = graph.input("ba19_oligo", "BA19 oligo 5'-BAAAAAAAAAAAAAAAAAAA-3'",
                       MaterialKind.OLIGO, state="B" + "A" * 19, source=MAIN)
    protected = graph.product("ba19_protected_beads",
                              "bead primers with BA19 annealed to full-length T19V tails",
                              MaterialKind.MOLECULAR, source=MAIN)
    graph.react("ba19_hybridization", "Hybridize BA19 to bead primers (RT, 30 min; Steps 70-72)",
                StepKind.HYBRIDIZATION,
                ((digested, InputRole.TEMPLATE), (ba19, InputRole.STRUCTURAL)),
                (protected,), source=MAIN)
    exo = graph.input("exonuclease_i", "Exonuclease I (0.27 U/µl, RT, 2 h; Step 74)",
                      MaterialKind.CATALYST, source=MAIN)
    trimmed = graph.product("exo_treated_beads",
                            "BA19-protected full-length bead primers after ExoI",
                            MaterialKind.MOLECULAR, source=MAIN)
    incomplete = graph.product("incomplete_bead_primers",
                               "ExoI hydrolysis products of shorter bead primers "
                               "(~102/101 and 65 nt precursor bands)",
                               MaterialKind.MOLECULAR, source=QC_SOURCE)
    graph.react("bead_cleanup", "Exonuclease I removes primers missing poly-T tails; "
                "stopped with STOP-25 (Steps 74-75)",
                StepKind.DIGESTION,
                ((protected, InputRole.TEMPLATE), (exo, InputRole.CATALYST)),
                (trimmed, incomplete), source=MAIN)
    graph.mark_terminal(Terminal(incomplete.id, Disposition.UNRESOLVED,
        "Step 74 states ExoI is added 'to remove primers missing poly-T tails' and the "
        "Fig. 6 legend calls the lower bands synthesis intermediates hydrolysed by ExoI; "
        "Step 76 subsequently washes the beads, but this graph does not resolve "
        "the soluble hydrolysis products versus bead-attached residual fragments. "
        "Hydrolysis is not itself physical discard."))
    bead = graph.product("barcoded_bead", "finished split-pool barcode bead primer",
                         MaterialKind.OLIGO, state=I.bead_primer(), source=SOURCE)
    removed_ba19 = graph.product("removed_ba19", "BA19 oligo stripped from the beads",
                                 MaterialKind.WASTE, source=MAIN)
    alkali = graph.input("bead_denaturation", "NaOH denaturation solution, then "
                         "neutralization buffer (Steps 57-60, repeated at Step 76)",
                         MaterialKind.CONDITION, source=MAIN)
    graph.react("ba19_removal", "STOP-10 washes and alkaline denaturation remove BA19 "
                "(Steps 76-79)", StepKind.STRAND_SEPARATION,
                ((trimmed, InputRole.TEMPLATE), (alkali, InputRole.CONDITION)),
                (bead, removed_ba19), source=MAIN)
    graph.mark_terminal(Terminal(removed_ba19.id, Disposition.DISCARD,
        "Removed with the wash supernatants (Steps 76 and 53-60)."))
    graph.annotate(SourceDiscrepancy(
        "bead_cleanup", "Published bead-QC intermediate band size differs",
        (SourceClaim("Supplementary Method 2, procedure step 9",
                     "bands at 135, 102 and 65 nt before clean-up"),
         SourceClaim("Supplementary Method 2, Anticipated result",
                     "expected sizes 135, 101 and 65 nt before ExoI treatment")),
        "No single intermediate length is asserted; the bead primer keeps the published "
        "variable 8-11 nt barcode 1, so the round-1-only species varies in length.",
        "Both values are printed in the same publisher supplement (MOESM454). The main "
        "text (Fig. 6 legend) gives the full-length primer as 134-136 nt and no "
        "intermediate size.",
        "Which value the authors intended is unresolved."))

    cells = graph.input("cells", "single-cell suspension", MaterialKind.BIOLOGICAL,
                        primary=True)
    loading = graph.input("droplet_loading", "one-cell/one-bead droplet loading condition",
                          MaterialKind.CONDITION)
    droplets = graph.product("productive_droplets", "one-cell/one-bead droplets",
                             MaterialKind.MOLECULAR)
    excluded = graph.product("excluded_droplets", "empty or multiplet droplets",
                             MaterialKind.WASTE)
    graph.react("partition", "Co-encapsulate cells and beads and photorelease primers",
                StepKind.PARTITION,
                ((cells, InputRole.TEMPLATE), (bead, InputRole.SPIKE_IN),
                 (loading, InputRole.CONDITION)), (droplets, excluded))
    graph.mark_terminal(Terminal(excluded.id, Disposition.UNRESOLVED,
        "Physical fate not stated in the main text or inspected supplements; the published "
        "pipeline (MOESM460 README step 2) separates cell barcodes from low-abundance background "
        "computationally by a read-count threshold."))
    rt = graph.input("reverse_transcriptase", "SuperScript III reverse transcriptase "
                     "(in-drop, after 365 nm UV primer release; 50 °C 2 h, 70 °C 15 min; "
                     "Steps 127-128)", MaterialKind.CATALYST, source=MAIN)
    first = graph.product("barcoded_cdna", "bead-barcoded first-strand cDNA:mRNA hybrid "
                          "with single-stranded T7 promoter", MaterialKind.MOLECULAR,
                          source=MAIN)
    graph.react("reverse_transcription", "Reverse-transcribe from released T19V primers",
                StepKind.REVERSE_TRANSCRIPTION,
                ((droplets, InputRole.TEMPLATE), (bead, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (first,), source=MAIN)
    digest = graph.input("exo_hinfi_mix", "ExoI (1 U/µl) + FastDigest HinfI after PFO "
                         "emulsion breaking; 37 °C 30 min (Steps 131-137)",
                         MaterialKind.CATALYST, source=MAIN)
    cleaned = graph.product("cleaned_cdna", "cDNA:mRNA hybrid after digestion and 1.2x "
                            "AMPure (Step 138)", MaterialKind.MOLECULAR, source=MAIN)
    leftovers = graph.product("digested_primers", "unused ssDNA primers (ExoI) and "
                              "T7-promoter primer dimers/hairpins (HinfI)",
                              MaterialKind.WASTE, source=MAIN)
    spent_beads = graph.product("spent_hydrogel_beads", "hydrogel beads retained on the "
                                "0.45 µm filter", MaterialKind.WASTE, source=MAIN)
    graph.react("post_rt_digest", "Digest unused primers (ExoI) and primer dimers (HinfI); "
                "remove beads by filtration; 1.2x AMPure", StepKind.DIGESTION,
                ((first, InputRole.TEMPLATE), (digest, InputRole.CATALYST)),
                (cleaned, leftovers, spent_beads), source=MAIN)
    graph.mark_terminal(Terminal(leftovers.id, Disposition.UNRESOLVED,
        "Fig. 3 legend: the target is immune because it is a cDNA:mRNA hybrid with a "
        "single-stranded T7 promoter. The 1.2x AMPure clean-up follows, but the source does "
        "not state that the digestion products are removed by it."))
    graph.mark_terminal(Terminal(spent_beads.id, Disposition.DISCARD,
        "Step 136: the filter membrane containing the barcoded hydrogel beads is discarded."))
    second_pol = graph.input("second_strand_polymerase", "NEBNext mRNA Second Strand "
                             "Synthesis Module (16 °C 2.5 h, 65 °C 20 min; Steps 139-140)",
                             MaterialKind.CATALYST, source=MAIN)
    ds = graph.product("t7_templates", "double-stranded T7-bearing cDNA",
                       MaterialKind.MOLECULAR, source=MAIN)
    graph.react("second_strand", "Second-strand synthesis makes the T7 promoter "
                "double-stranded", StepKind.EXTENSION,
                ((cleaned, InputRole.TEMPLATE), (second_pol, InputRole.CATALYST)), (ds,),
                source=MAIN)
    t7 = graph.input("t7_polymerase", "HiScribe T7 High Yield kit (37 °C 13-15 h; "
                     "Steps 141-142)", MaterialKind.CATALYST, source=MAIN)
    arna = graph.product("arna", "linearly amplified antisense RNA after 1.3x AMPure "
                         "(Step 143)", MaterialKind.MOLECULAR, source=MAIN)
    graph.react("ivt", "Transcribe pooled cDNA from the T7 promoter",
                StepKind.TRANSCRIPTION,
                ((ds, InputRole.TEMPLATE), (t7, InputRole.CATALYST)), (arna,), source=MAIN)
    frag = graph.input("zinc_fragmentation", "Ambion RNA Fragmentation Reagents "
                       "(zinc ion; 70 °C for 3 min, then stop solution; Steps 146-150)",
                       MaterialKind.CATALYST, source=MAIN)
    frag_arna = graph.product("fragmented_arna", "fragmented aRNA after AMPure",
                              MaterialKind.MOLECULAR, source=MAIN)
    graph.react("arna_fragmentation", "Fragment aRNA by zinc-ion hydrolysis",
                StepKind.FRAGMENTATION,
                ((arna, InputRole.TEMPLATE), (frag, InputRole.CATALYST)), (frag_arna,),
                source=MAIN)

    pe2n6 = graph.input("pe2_n6", "PE2-N6 random primer", MaterialKind.OLIGO,
                        state=PE2_N6, source=SOURCE)
    rt2 = graph.input("primescript_rt", "PrimeScript reverse transcriptase (primer "
                      "annealing 70 °C 3 min; 30 °C 10 min, 42 °C 1 h, 70 °C 15 min; "
                      "Steps 151-154)", MaterialKind.CATALYST, source=MAIN)
    converted = graph.product("sense_cdna", "PE2-tagged sense cDNA after 1.2x AMPure "
                              "(Step 156)", MaterialKind.MOLECULAR, source=MAIN)
    graph.react("second_rt", "Random-prime fragmented aRNA with PE2-N6",
                StepKind.REVERSE_TRANSCRIPTION,
                ((frag_arna, InputRole.TEMPLATE), (pe2n6, InputRole.STRUCTURAL),
                 (rt2, InputRole.CATALYST)), (converted,), source=MAIN)
    forward = graph.input("pcr_forward", "PE2 primer (P5 side)",
                          MaterialKind.OLIGO, state=PE2_PCR, source=SOURCE)
    reverse = graph.input("pcr_reverse", "PE1 primer index 1 (P7 side)",
                          MaterialKind.OLIGO, state=PE1_INDEX1, source=SOURCE)
    pcr_pol = graph.input("pcr_polymerase", "Kapa HiFi HotStart (2 cycles at 55 °C "
                          "annealing, then qPCR-determined 65 °C cycles, typically 8-12; "
                          "Steps 157-161)", MaterialKind.CATALYST, source=MAIN)
    final = I.v2_library()
    library = graph.product("library", final.name + " (0.7x AMPure; 400-600 bp average, "
                            "Steps 163-164)", MaterialKind.LIBRARY, state=final, source=MAIN)
    graph.react("library_pcr", "Add P5, P7 and sample index by enrichment PCR",
                StepKind.PCR,
                ((converted, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (pcr_pol, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="The pre-PCR IVT-conversion intermediate is not represented base-for-base."),
                source=MAIN)
    graph.annotate(SourceDiscrepancy(
        "library_pcr", "Indexed PCR primer named PE1 and PE2 in the same passage",
        (SourceClaim("Main text Reagent Setup, 'PE1/PE2 primer mix'",
                     "'PE1 primer: 5'-CAAGCAGAAGACGGCATACGAGATXXXXXXCTCTTTCCCTACACGA' ... "
                     "'See Supplementary Table 2 for a list of 24 PE2 primer variants'"),
         SourceClaim("Supplementary Table 2 (MOESM455)",
                     "indexed rows are 'PE1 primer index 1-24'; 'PE2 primer' is the "
                     "unindexed P5 primer")),
        "The model indexes the PE1 (P7-side) primer, matching the printed bases and "
        "Supplementary Table 2.",
        "The printed sequences agree; only the name in one sentence differs.",
        "None for the bases; the sentence is read as a naming slip."))
    graph.sequence(library, SequencingPlan(
        library=final, primers=I.V2_PRIMERS,
        acquired_roles=tuple(p.role for p in I.V2_PRIMERS),
        read_lengths={},
        read_uses={"Read 1": "include", "Index 1 (i7)": "demultiplex",
                   "Read 2": "extract"},
        limitation=READ_LIMITATION))
    graph.declare_reviewed_scope((
        "bead_seed", "barcode_plate_oligos", "bead_polymerase", "bead_synthesis",
        "extended_beads", "bead_hinfi", "bead_hinfi_digest", "ba19_oligo",
        "ba19_hybridization", "exonuclease_i", "bead_cleanup", "bead_denaturation",
        "ba19_removal",
        "cells", "droplet_loading", "partition", "reverse_transcriptase",
        "reverse_transcription", "exo_hinfi_mix", "post_rt_digest",
        "second_strand_polymerase", "second_strand", "t7_polymerase", "ivt",
        "zinc_fragmentation", "arna_fragmentation", "pe2_n6", "primescript_rt",
        "second_rt", "pcr_forward", "pcr_reverse", "pcr_polymerase", "library_pcr",
        "library"), MAIN)
    return graph.require_publishable()
