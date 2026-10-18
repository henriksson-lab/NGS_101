"""Source-backed models for named isothermal-amplification sequencing protocols."""
from __future__ import annotations

from types import SimpleNamespace

from chemdraw import Construct, Scene, Segment, complement_segments, feature
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import (Disposition, Terminal, Highlight, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind)
from protocol_model import molecule, role
import seqprimers as sp


def _duplex(name):
    parts = [role(name)]
    return Scene.duplex(parts, unpaired=(name + "'",), label=name)


def _population(*names):
    sc = Scene()
    for i, name in enumerate(names):
        sc.strand("component_" + str(i), [role(name)], label=name)
    return sc


def _tailed():
    body = role("cDNA fragment")
    bottom = complement_segments([body])[0]
    sc = Scene()
    sc.strand("top", [body, Segment("top 3-prime dA", "A")])
    sc.anneal("bottom", [bottom, Segment("bottom 3-prime dA", "A")],
              to="top", pair=(bottom.name, body.name), unpaired=(bottom.name,))
    sc.overhang("top", "top 3-prime dA", polarity="3-prime", owner="substrate")
    sc.overhang("bottom", "bottom 3-prime dA", polarity="3-prime", owner="substrate")
    return sc


OVATION_LIBRARY = Construct([
    role("native NEB barcoded adapter region: placement unassigned", inferred=True,
         molecular_feature=feature("ovation_sample_barcode", "sample_index", "unknown",
             note="Bar-coded adapter is source-positive; represented region does not assign physical side, barcode width or read.")),
    role("fragment of Ribo-SPIA amplified cDNA"),
    role("opposite native NEB adapter end: architecture unassigned", inferred=True)
], name="Ovation-derived Illumina library")
OVATION_PRIMERS = ()


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
                      state=_population("RNA template (hybrid geometry unassigned)",
                                        "chimeric-primer-derived first-strand cDNA"))
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
                       state=_population("first cDNA strand with primer-derived terminal RNA region",
                                         "opposite DNA strand (heteroduplex geometry unassigned)"))
    heat_products = g.product("heat_rna_fragments", "RNA products of partial heat degradation",
                              MaterialKind.MOLECULAR,
                              state=_population("partially degraded RNA products"))
    g.react("second_strand", "Partially degrade the RNA template by heating; synthesize the second strand",
            StepKind.EXTENSION,
            ((first, InputRole.TEMPLATE), (heat, InputRole.CONDITION),
             (polymerase, InputRole.CATALYST)), (duplex, heat_products))
    g.mark_terminal(Terminal(heat_products.id, Disposition.UNRESOLVED,
                            "Heat degradation is not physical separation; later bead partition unassigned."))
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
                          MaterialKind.MOLECULAR, state=_population("retained SPIA template population",
                                            "separate displaced single-strand cDNA copies"))
    hydrolysate = g.product("rnase_h_products", "RNA hydrolysis products from RNase H",
                           MaterialKind.MOLECULAR, state=_population("RNA hydrolysis products"))
    g.react("spia", "RNase H removes the heteroduplex RNA; SPIA primer extension displaces the forward strand",
            StepKind.EXTENSION,
            ((purified, InputRole.TEMPLATE), (spia_primer, InputRole.STRUCTURAL),
             (rnase_h, InputRole.CATALYST), (polymerase, InputRole.CATALYST)), (amplified, hydrolysate))
    g.mark_terminal(Terminal(hydrolysate.id, Disposition.UNRESOLVED,
                            "RNase H hydrolysis does not itself remove the products."))
    hexamers = g.input("second_strand_hexamers", "random hexamers for linear second-strand amplification",
                       MaterialKind.OLIGO, inferred=True,
                       source=source + " (primer class named; bases and modifications not printed)")
    ds_spia = g.product("ds_spia_cdna", "double-stranded amplified cDNA",
                        MaterialKind.MOLECULAR,
                        state=_duplex("random-hexamer-copied SPIA cDNA"))
    g.react("hexamer_second_strand", "Amplify the second-strand cDNA linearly with random hexamers",
            StepKind.EXTENSION,
            ((amplified, InputRole.TEMPLATE), (hexamers, InputRole.STRUCTURAL),
             (polymerase, InputRole.CATALYST)), (ds_spia,))
    shear = g.input("shear_condition", "Covaris S1 sonication: duty cycle 5, intensity 3, 200 cycles/burst, 180 s",
                    MaterialKind.CONDITION)
    fragmented = g.product("fragmented_cdna", "sheared amplified cDNA (100–360 bp)", MaterialKind.MOLECULAR,
                           state=_duplex("transcript-derived cDNA fragment"))
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
                       state=_tailed())
    g.react("da_tail", "dA-tail the fragments (Illumina libraries)", StepKind.TAILING,
            ((blunt, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapters = g.input("library_adapters", "platform-specific double-stranded bar-coded adapters (NEB kit)",
                       MaterialKind.ADAPTER, inferred=True,
                       source=source + " (exact kit and adapter oligos unresolved)")
    ligase = g.input("adapter_ligase", "NEB library-kit ligase", MaterialKind.CATALYST,
                     source=source + " (kit and ligase not named)")
    ligated = g.product("adapter_ligated", "adapter-ligated cDNA fragments", MaterialKind.MOLECULAR,
                        state=molecule("native ligation product; duplex/nicks unassigned",
                                       tuple(OVATION_LIBRARY)))
    g.react("ligate_adapters", "Ligate bar-coded adapters", StepKind.LIGATION,
            ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
             (ligase, InputRole.CATALYST)), (ligated,),
            geometry=AdapterLigationGeometry.unavailable(
                "The defining paper names platform-specific bar-coded NEB adapters after dA-tailing "
                "but does not identify the kit or print the adapter strands and ends."))
    purification = g.input("ligation_purification", "automated post-reaction purification (method not named)",
                           MaterialKind.SELECTION)
    cleaned = g.product("purified_ligation", "purified adapter-ligated fragments", MaterialKind.MOLECULAR,
                        state=ligated.state)
    cleaned_waste = g.product("purification_waste", "unligated adapters and reaction components",
                              MaterialKind.WASTE)
    g.react("purify_ligation", "Purify the ligation (automated)", StepKind.CLEANUP,
            ((ligated, InputRole.TEMPLATE), (purification, InputRole.SELECTION)),
            (cleaned, cleaned_waste))
    g.discard(cleaned_waste, "removed by the automated purification")
    gel2 = g.input("library_gel", "manual gel purification (size range not stated)",
                   MaterialKind.SELECTION)
    gel_library = g.product("gel_purified_ligation", "gel-purified adapter-ligated fragments",
                            MaterialKind.MOLECULAR, state=cleaned.state)
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
                        source=source + " (native adapter strands and index position unassigned)")
    g.react("library_pcr", "Amplify the library (15 cycles, manual)", StepKind.PCR,
            ((gel_library, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
             (reverse, InputRole.STRUCTURAL), (pcr_polymerase, InputRole.CATALYST)), (library,),
            geometry=PCRGeometry.unavailable(
                gel_library.state, reason="Native NEB PCR primer identities and strand geometry remain unassigned in the inspected source scope."))
    g.sequence(library, SequencingPlan(
        OVATION_LIBRARY, OVATION_PRIMERS,
        ("Read 1", "Read 2"),
        {"Read 1": 50, "Read 2": 50},
        {"Read 1": "include", "Read 2": "include"},
        limitation="Published paired 2×50; native primer sites and barcode acquisition unassigned."))
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
        "gel_purify", "Manual gel/PCR order is a selected schematic order",
        "Methods list manual gel purification and 15-cycle library amplification together. "
        "The graph places the gel first as a declared model ordering; their experimental "
        "relative order is not established by the inspected article. Barcodes are source-positive, "
        "but position, width and acquisition are unassigned; no default index is appended."))
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
            "platform-specific bar-coded adapters. Their native strand architecture, barcode "
            "position/width and PCR/sequence primers remain unassigned in this source scope. "
            "INFERRED role ends do not assert a canonical TruSeq shell."),
    oligos=lambda: (), FINAL_LIBRARY=OVATION_LIBRARY,
    FINAL_CAPTION=("INFERRED native end roles — sequencing library made from sheared Ribo-SPIA cDNA. "
                   "The transcript-derived insert is source-established; exact commercial "
                   "adapter/index oligos are not printed in the defining paper."),
    SEQ_PRIMERS=OVATION_PRIMERS,
    ACQUIRED_ROLES=("Read 1", "Read 2"),
    READ_LENGTHS={"Read 1": 50, "Read 2": 50},
    READ_USES={"Read 1": "include", "Read 2": "include"},
    SHOW_READ_CONTENT=False,
    READ_SPECS=(
        sp.read_spec("Read 1", 50, "cDNA fragment; exact entry site unassigned", "include"),
        sp.read_spec("Read 2", 50, "opposite cDNA end; exact entry site unassigned", "include"),
        sp.read_spec("Index 1 (i7)", "unassigned", "Barcode acquisition/position unassigned", "undetermined"),
        sp.read_spec("Index 2 (i5)", "unassigned", "Barcode acquisition/position unassigned", "undetermined")),
    SEQUENCING_INTRO=("Published results describe paired 2 × 50-base reads on Genome Analyzer IIx. "
                      "Bar-coded libraries are established; index acquisition and cycle counts "
                      "are unassigned here. This is not a finding that index reads were absent."))



def protocol_graph() -> ProtocolGraph:
    return _ovation_graph(source)

spec.SEQUENCING_ENDING = spec.SEQUENCING_INTRO
spec.protocol_graph = protocol_graph
globals().update(vars(spec))
