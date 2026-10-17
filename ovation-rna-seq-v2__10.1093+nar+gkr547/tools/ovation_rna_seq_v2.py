"""Source-backed models for named isothermal-amplification sequencing protocols."""
from __future__ import annotations

from types import SimpleNamespace

import illumina as il
from batch_ngs import nextera_library, truseq_amplicon_library, truseq_library
from chemdraw import Construct, Segment, feature, oligo, revcomp
from isothermal import (lamp_product, lamp_topology, mda_product, proximity_extension_pair,
                        proximity_extension_product, spia_product)
from nanopore_library import (motor_library, motor_sequencing_plan, no_primer_ending,
                              role as nanopore_role)
from pcr import PCRGeometry, PCRPrimer
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import (Highlight, InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)
from protocol_model import molecule, role
import seqprimers as sp
from proximity_topology import ProximityExtensionTopology




OVATION_LIBRARY, OVATION_PRIMERS = truseq_library([
    role("fragment of Ribo-SPIA amplified cDNA")
], "Ovation-derived Illumina library", dual_index=False, inferred_adapters=True)


def _ovation_graph(source: str) -> ProtocolGraph:
    g = ProtocolGraph("Ovation RNA-Seq V2 / Ribo-SPIA", authority=source)
    rna = g.input("rna", "total RNA", MaterialKind.BIOLOGICAL,
                  state=molecule("total RNA", (role("RNA transcript"),)), primary=True)
    first_primer = g.input(
        "first_primer", "random-hexamer and poly-T chimeric first-strand primer mix",
        MaterialKind.OLIGO, inferred=True,
        source=source + " (Introduction and Methods name the primer classes; bases not printed)")
    rt = g.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    first = g.product("first_strand", "first-strand cDNA/RNA hybrid",
                      MaterialKind.MOLECULAR,
                      state=molecule("first-strand product", (
                          role("RNA template"), role("first-strand cDNA"),
                          role("chimeric primer handle"))))
    g.react("first_strand_rt", "Prime RNA with random hexamers and poly-T chimeric primer; reverse transcribe",
            StepKind.REVERSE_TRANSCRIPTION,
            ((rna, InputRole.TEMPLATE), (first_primer, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (first,))
    heat = g.input("template_heating", "heating step that partially degrades the RNA template",
                   MaterialKind.CONDITION)
    polymerase = g.input("dna_polymerase", "DNA polymerase (not named in the paper)",
                         MaterialKind.CATALYST)
    duplex = g.product("cdna_duplex", "double-stranded cDNA with a DNA/RNA heteroduplex end",
                       MaterialKind.MOLECULAR,
                       state=molecule("SPIA-competent cDNA", (
                           role("double-stranded transcript copy"),
                           role("DNA/RNA heteroduplex at the 5' end"))))
    g.react("second_strand", "Partially degrade the RNA template by heating; synthesize the second strand",
            StepKind.EXTENSION,
            ((first, InputRole.TEMPLATE), (heat, InputRole.CONDITION),
             (polymerase, InputRole.CATALYST)), (duplex,))
    rnaclean = g.input("rnaclean_xp", "Agencourt RNAClean XP beads", MaterialKind.SELECTION)
    purified = g.product("purified_duplex", "bead-purified double-stranded cDNA",
                         MaterialKind.MOLECULAR, state=duplex.state)
    bead_waste = g.product("rnaclean_supernatant", "RNAClean XP supernatant and washes",
                           MaterialKind.WASTE)
    g.react("purify_duplex", "Purify double-stranded cDNA with RNAClean XP beads",
            StepKind.CLEANUP,
            ((duplex, InputRole.TEMPLATE), (rnaclean, InputRole.SELECTION)),
            (purified, bead_waste))
    g.discard(bead_waste, "unbound reaction components removed by the bead purification")
    spia_primer = g.input("spia_primer", "SPIA primer", MaterialKind.OLIGO,
                          inferred=True, source=source + " (primer bases not printed)")
    rnase_h = g.input("rnase_h", "RNase H", MaterialKind.CATALYST)
    amplified = g.product("spia_cdna", "linearly amplified Ribo-SPIA cDNA population",
                          MaterialKind.MOLECULAR, state=spia_product())
    g.react("spia", "RNase H removes the heteroduplex RNA; SPIA primer extension displaces the forward strand",
            StepKind.EXTENSION,
            ((purified, InputRole.TEMPLATE), (spia_primer, InputRole.STRUCTURAL),
             (rnase_h, InputRole.CATALYST), (polymerase, InputRole.CATALYST)), (amplified,))
    hexamers = g.input("second_strand_hexamers", "random hexamers for linear second-strand amplification",
                       MaterialKind.OLIGO, inferred=True,
                       source=source + " (primer class named; bases and modifications not printed)")
    ds_spia = g.product("ds_spia_cdna", "double-stranded amplified cDNA",
                        MaterialKind.MOLECULAR,
                        state=molecule("double-stranded SPIA cDNA", (
                            role("SPIA cDNA strand"), role("random-hexamer-primed complementary strand"))))
    g.react("hexamer_second_strand", "Amplify the second-strand cDNA linearly with random hexamers",
            StepKind.EXTENSION,
            ((amplified, InputRole.TEMPLATE), (hexamers, InputRole.STRUCTURAL),
             (polymerase, InputRole.CATALYST)), (ds_spia,))
    shear = g.input("shear_condition", "Covaris S1 sonication: duty cycle 5, intensity 3, 200 cycles/burst, 180 s",
                    MaterialKind.CONDITION)
    fragmented = g.product("fragmented_cdna", "sheared amplified cDNA (100–360 bp)", MaterialKind.MOLECULAR,
                           state=molecule("sheared cDNA", (role("transcript-derived cDNA fragment"),)))
    g.react("fragment", "Shear ~0.5–1 µg double-stranded cDNA (sheared branch)",
            StepKind.FRAGMENTATION,
            ((ds_spia, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (fragmented,))
    gel = g.input("agarose_gel", "2.5% agarose gel, 250–350 bp excision for Illumina",
                  MaterialKind.SELECTION)
    sized = g.product("sized_cdna", "250–350 bp cDNA fragments", MaterialKind.MOLECULAR,
                      state=fragmented.state)
    gel_rest = g.product("gel_rejected", "cDNA outside the excised size range", MaterialKind.WASTE)
    g.react("size_select", "Gel size-select 250–350 bp fragments", StepKind.SELECTION,
            ((fragmented, InputRole.TEMPLATE), (gel, InputRole.SELECTION)), (sized, gel_rest))
    g.discard(gel_rest, "left in the gel outside the excised size range")
    end_repair = g.input("end_repair_enzymes", "NEB library-kit end-repair enzymes",
                         MaterialKind.CATALYST, source=source + " (kit and enzyme mix not named)")
    blunt = g.product("blunt_cdna", "blunt-ended cDNA fragments", MaterialKind.MOLECULAR,
                      state=sized.state)
    g.react("end_repair", "Blunt the fragments by end repair", StepKind.END_REPAIR,
            ((sized, InputRole.TEMPLATE), (end_repair, InputRole.CATALYST)), (blunt,))
    tailing = g.input("da_tailing_enzyme", "NEB library-kit dA-tailing enzyme",
                      MaterialKind.CATALYST, source=source + " (kit and enzyme not named)")
    tailed = g.product("da_tailed_cdna", "dA-tailed cDNA fragments", MaterialKind.MOLECULAR,
                       state=sized.state)
    g.react("da_tail", "dA-tail the fragments (Illumina libraries)", StepKind.TAILING,
            ((blunt, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapters = g.input("library_adapters", "platform-specific double-stranded bar-coded adapters (NEB kit)",
                       MaterialKind.ADAPTER, inferred=True,
                       source=source + " (exact kit and adapter oligos unresolved)")
    ligase = g.input("adapter_ligase", "NEB library-kit ligase", MaterialKind.CATALYST,
                     source=source + " (kit and ligase not named)")
    ligated = g.product("adapter_ligated", "adapter-ligated cDNA fragments", MaterialKind.MOLECULAR)
    g.react("ligate_adapters", "Ligate bar-coded adapters", StepKind.LIGATION,
            ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
             (ligase, InputRole.CATALYST)), (ligated,),
            geometry=AdapterLigationGeometry.unavailable(
                "The defining paper names platform-specific bar-coded NEB adapters after dA-tailing "
                "but does not identify the kit or print the adapter strands and ends."))
    purification = g.input("ligation_purification", "automated post-reaction purification (method not named)",
                           MaterialKind.SELECTION)
    cleaned = g.product("purified_ligation", "purified adapter-ligated fragments", MaterialKind.MOLECULAR)
    cleaned_waste = g.product("purification_waste", "unligated adapters and reaction components",
                              MaterialKind.WASTE)
    g.react("purify_ligation", "Purify the ligation (automated)", StepKind.CLEANUP,
            ((ligated, InputRole.TEMPLATE), (purification, InputRole.SELECTION)),
            (cleaned, cleaned_waste))
    g.discard(cleaned_waste, "removed by the automated purification")
    gel2 = g.input("library_gel", "manual gel purification (size range not stated)",
                   MaterialKind.SELECTION)
    gel_library = g.product("gel_purified_ligation", "gel-purified adapter-ligated fragments",
                            MaterialKind.MOLECULAR)
    gel2_rest = g.product("library_gel_rejected", "material outside the excised band", MaterialKind.WASTE)
    g.react("gel_purify", "Gel-purify the ligated library (manual)", StepKind.SELECTION,
            ((cleaned, InputRole.TEMPLATE), (gel2, InputRole.SELECTION)), (gel_library, gel2_rest))
    g.discard(gel2_rest, "left in the gel outside the excised band")
    forward = g.input("pcr_forward", "library PCR forward primer", MaterialKind.OLIGO, inferred=True,
                      source=source + " (primer not named or printed)")
    reverse = g.input("pcr_reverse", "library PCR reverse primer", MaterialKind.OLIGO, inferred=True,
                      source=source + " (primer not named or printed)")
    pcr_polymerase = g.input("pcr_polymerase", "library PCR polymerase (not named)", MaterialKind.CATALYST)
    library = g.product("library", "Illumina library from Ribo-SPIA cDNA",
                        MaterialKind.LIBRARY, state=OVATION_LIBRARY, inferred=True,
                        source=source + " (canonical shell inferred from compatibility)")
    g.react("library_pcr", "Amplify the library (15 cycles, manual)", StepKind.PCR,
            ((gel_library, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
             (reverse, InputRole.STRUCTURAL), (pcr_polymerase, InputRole.CATALYST)), (library,),
            geometry=PCRGeometry.unavailable(
                None, reason="The defining paper does not print the library PCR primers."))
    g.sequence(library, SequencingPlan(
        OVATION_LIBRARY, OVATION_PRIMERS,
        ("Read 1", "Index 1 (i7)", "Read 2"),
        {"Read 1": 50, "Read 2": 50},
        {"Read 1": "include", "Index 1 (i7)": "demultiplex", "Read 2": "include"}))
    g.annotate(SourceDiscrepancy(
        "second_strand", "The article describes second-strand synthesis in two ways",
        (SourceClaim("Introduction, Ovation paragraph (PMC3185437 JATS line 26)",
                     "the RNA template is partially degraded in a heating step and the second "
                     "strand is synthesized along the first-strand cDNA using DNA polymerase"),
         SourceClaim("Materials and methods, cDNA synthesis (PMC3185437 JATS line 56)",
                     "double-stranded DNA is generated by fragmentation of the mRNA template "
                     "strand using RNA-dependant DNA polymerase")),
        "The model follows the Introduction: heat-fragmented RNA template, then an unnamed "
        "DNA polymerase copies the first-strand cDNA. No second-strand primer is drawn.",
        "A second strand copied on a cDNA template is DNA-templated; the Methods enzyme "
        "wording does not describe a separate reagent or step.",
        "The commercial V2 user guide and Kurn et al. (2005), which would identify the "
        "polymerase and any priming of the second strand, remain uninspected."))
    g.annotate(Highlight(
        "fragment", "Unsheared and mixed branches are not drawn",
        "The study also prepared libraries from unsheared cDNA (100–550 bp) and an equal "
        "mix of sheared and unsheared cDNA (100–500 bp). Optional DNase treatment of the "
        "amplified material in the comparison series is not drawn either."))
    g.declare_reviewed_scope((
        "rna", "first_primer", "reverse_transcriptase", "first_strand_rt",
        "template_heating", "dna_polymerase", "second_strand", "rnaclean_xp",
        "purify_duplex", "spia_primer", "rnase_h", "spia", "second_strand_hexamers",
        "hexamer_second_strand", "shear_condition", "fragment", "agarose_gel",
        "size_select", "end_repair_enzymes", "end_repair", "da_tailing_enzyme", "da_tail",
        "library_adapters", "adapter_ligase", "ligate_adapters", "ligation_purification",
        "purify_ligation", "library_gel", "gel_purify", "pcr_forward", "pcr_reverse",
        "pcr_polymerase", "library_pcr", "library"), source)
    return g.require_publishable()



source = ('<a href="https://doi.org/10.1093/nar/gkr547">'
          'Tariq et al., <i>Nucleic Acids Research</i> (2011)</a>; the NuGEN '
          'Ovation RNA-Seq V2 user guide it follows remains uninspected.')
spec = SimpleNamespace(
    TITLE="Ovation RNA-Seq V2 — Ribo-SPIA whole-transcriptome amplification",
    NOTES="01_ovation-rna-seq-v2.html", SOURCE=source,
    SUMMARY=("Random hexamers and a poly-T chimeric primer initiate first-strand cDNA. "
             "After heat-mediated partial RNA degradation and second-strand synthesis, "
             "RNase H repeatedly exposes the SPIA priming site and a strand-displacing "
             "polymerase releases linear cDNA copies; random hexamers then copy the "
             "second strand before conventional Illumina library construction."),
    CAVEAT=("The article does not print the chimeric, SPIA or random primer bases, and "
            "the commercial V2 user guide and the cited SPIA predecessor remain "
            "uninspected. The library used unnamed New England Biolabs kits with "
            "platform-specific bar-coded adapters; the canonical TruSeq-compatible ends, "
            "the i7 index position and its length are inferred, not published."),
    oligos=lambda: (), FINAL_LIBRARY=OVATION_LIBRARY,
    FINAL_CAPTION=("INFERRED shell — sequencing library made from sheared Ribo-SPIA cDNA. "
                   "The transcript-derived insert is source-established; exact commercial "
                   "adapter/index oligos are not printed in the defining paper."),
    SEQ_PRIMERS=OVATION_PRIMERS,
    ACQUIRED_ROLES=("Read 1", "Index 1 (i7)", "Read 2"),
    READ_LENGTHS={"Read 1": 50, "Read 2": 50},
    READ_USES={"Read 1": "include", "Index 1 (i7)": "demultiplex", "Read 2": "include"},
    SEQUENCING_INTRO=("The study acquired paired-end 2 × 50-base reads on a Genome "
                      "Analyzer IIx. The libraries were bar-coded, but the paper does not "
                      "state the barcode position or index-read length; the index read "
                      "below follows the inferred shell and has no published length."))


def protocol_graph() -> ProtocolGraph:
    return _ovation_graph(source)

spec.protocol_graph = protocol_graph
globals().update(vars(spec))
