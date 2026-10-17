"""Source-owned LIGR-seq definition and reaction graph.

Order follows the Sharma et al. (2016) Experimental Procedures and Supplemental
Experimental Procedures: in vivo AMT/365 nm crosslinking, RNA extraction and DNase,
Ribo-Zero Gold depletion, limited S1 nuclease digestion, circRNA-ligase proximity
ligation (85 C denaturation, 60 C ligation), RNase R enrichment, 254 nm crosslink
reversal, then a SMARTer stranded library.  Ligation precedes RNase R.  The SMARTer
kit oligos are not printed in the inspected sources, so the library is a role
endpoint rather than an inferred TruSeq shell.
"""
from __future__ import annotations

from chemdraw import Construct, feature
from interaction_capture import ReversibleDuplexCrosslink
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import chain, molecule, role
from protocol_specs import Definition, get

RANDOM_BARCODE = feature("ligr.random_barcode", "inline_barcode", "random",
                         note="5-prime random barcode trimmed from single-end reads; "
                              "length and origin within the SMARTer kit not stated")
INDEX = feature("ligr.index", "sample_index", "unknown",
                note="barcode index carried by the kit reverse PCR primer; selected index not stated")

FINAL_LIBRARY = Construct([
    role("SMARTer forward-primer end (kit sequence not printed)"),
    role("5-prime random barcode", molecular_feature=RANDOM_BARCODE),
    role("RNA partner A-derived cDNA"),
    role("RNA partner B-derived cDNA"),
    role("SMARTer reverse-primer end (kit sequence not printed)"),
    role("barcode index (position within the reverse primer not stated)", molecular_feature=INDEX),
], name="LIGR-seq SMARTer stranded library (role endpoint)")

SEQUENCING = ("Single-end, strand-specific reads; the 5-prime random barcode is trimmed "
              "before Aligater chimera detection. Read length, platform primer sequences "
              "and the selected index are not given in the inspected sources. Partner order "
              "in a read depends on which RNA end was ligated, and circular ligation "
              "products can appear as concatemers of the two RNAs (Figure S1C).")

DEFINITION = Definition(
    "LIGR-seq — ligation of interacting RNA", "10.1016/j.molcel.2016.04.030",
    "Sharma et al., Molecular Cell (2016)",
    "In vivo AMT crosslinks hold RNA duplexes together while S1-trimmed ends are joined by circRNA ligase; RNase R then enriches crosslinked, ligated RNA before 254 nm reversal and SMARTer stranded library construction.",
    "role", ("5-prime random barcode", "RNA A-derived cDNA", "RNA B-derived cDNA"),
    SEQUENCING, "RNA",
    "Model correction (2026-10-17): earlier versions enriched with RNase R before ligation and drew a dual-indexed TruSeq shell with a dA ligation junction; the source performs ligation first, and the SMARTer kit is a template-switching library whose oligo sequences are not printed in the article or supplement. The library is therefore shown by role. The supplement names the SMARTer stranded kit generically (no catalogue number or revision); the 2013 Clontech components list and a third-party transcript of the 2014 manual (revision 050714) mark the SMARTer and N6 oligos proprietary. The study's kit generation, selected index and read length remain unverified."
)


def _workflow():
    x = ReversibleDuplexCrosslink("AMT", "RNase R digests uncrosslinked RNA")
    return chain(molecule("cellular RNA in 293T cells", (role("RNA molecules"),)), (
        ("Crosslink paired RNA in cells with AMT and 365 nm UV", x.duplex(),
         "20 µg/ml AMT on ice, 365 nm for 30 min in 5-min intervals; TRI-reagent extraction and TURBO DNase I follow."),
        ("Deplete rRNA, then limited S1 nuclease digestion", x.duplex(),
         "Ribo-Zero Gold on 4 µg total RNA; 400 ng digested with diluted S1 nuclease 30 min at room temperature, stopped by phenol–chloroform. S1 trims single-stranded RNA next to duplexes to create ligatable free ends."),
        ("Proximity-ligate free ends with circRNA ligase", x.chimera(),
         "50 ng RNA denatured 2 min at 85 °C, chilled, then ATP and circRNA ligase for 1 h at 60 °C. The elevated temperature limits in vitro rehybridization while the AMT crosslink keeps partners together."),
        ("Enrich with RNase R, then reverse crosslinks at 254 nm", molecule(
            "photoreversed RNA chimera", (role("RNA partner A fragment"), role("RNA partner B fragment")),
            junctions=(("RNA partner A fragment", "RNA partner B fragment", "RNA proximity ligation"),)),
         "RNase R (10 min, 37 °C) degrades linear RNA lacking crosslink protection; circular ligation products also resist it. After phenol–chloroform and ethanol, 254 nm for 10 min reverses AMT, which is removed by ethanol precipitation."),
        ("Build a SMARTer stranded library and PAGE-select >200 nt", molecule(
            "LIGR-seq library", tuple(FINAL_LIBRARY)),
         "First-strand primer adds a 5-prime universal sequence; the reverse transcriptase adds non-templated nucleotides that anneal a universal forward primer for second-strand synthesis. PCR (20 or 17 cycles) uses an indexed reverse primer; AMPure 1:2 cleanups; 6% TBE PAGE, products above 200 nt excised."),
    ))


SPEC = get("ligr-seq", definition=DEFINITION, final_override=FINAL_LIBRARY,
           workflow_fn=_workflow,
           sequencing_unavailable="The SMARTer kit primer sequences are not printed in the inspected article, supplement or kit components list, so read-primer landing is shown by role only.")


def protocol_graph():
    d = DEFINITION
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    topology = ReversibleDuplexCrosslink("AMT", "RNase R digests uncrosslinked RNA")
    reviewed = []

    def add(*ids):
        reviewed.extend(ids)

    rna = graph.input("rna", "293T cellular RNA containing paired regions",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("cellular RNA", (role("RNA molecules"),)), primary=True)
    amt = graph.input("crosslinker", "20 µg/ml AMT, 365 nm UV 30 min on ice",
                      MaterialKind.CONDITION)
    crosslinked = graph.product("crosslinked", "AMT-crosslinked RNA duplexes in cells",
                                MaterialKind.MOLECULAR, state=topology.duplex())
    graph.react("crosslink", "Photocrosslink paired RNA in vivo with AMT",
                StepKind.CROSS_LINKING,
                ((rna, InputRole.TEMPLATE), (amt, InputRole.CONDITION)), (crosslinked,))
    add("rna", "crosslinker", "crosslink")

    extraction = graph.input("tri_dnase", "TRI-reagent extraction and TURBO DNase I",
                             MaterialKind.SELECTION)
    total = graph.product("total_rna", "crosslinked total RNA", MaterialKind.MOLECULAR,
                          state=crosslinked.state)
    dna_waste = graph.product("dna_waste", "DNase-degraded DNA and non-aqueous fractions",
                              MaterialKind.WASTE)
    graph.react("extract", "Extract RNA and digest DNA", StepKind.EXTRACTION,
                ((crosslinked, InputRole.TEMPLATE), (extraction, InputRole.SELECTION)),
                (total, dna_waste))
    graph.discard(dna_waste, "TRI-reagent phase separation and DNase I digestion")
    add("tri_dnase", "extract", "dna_waste")

    ribozero = graph.input("ribozero", "Ribo-Zero Gold rRNA depletion (4 µg input)",
                           MaterialKind.SELECTION)
    depleted = graph.product("depleted", "rRNA-depleted crosslinked RNA",
                             MaterialKind.MOLECULAR, state=crosslinked.state)
    rrna = graph.product("rrna", "captured rRNA", MaterialKind.WASTE)
    graph.react("deplete_rrna", "Deplete ribosomal RNA", StepKind.SELECTION,
                ((total, InputRole.TEMPLATE), (ribozero, InputRole.SELECTION)),
                (depleted, rrna))
    graph.discard(rrna, "removed by Ribo-Zero Gold capture")
    add("ribozero", "deplete_rrna", "rrna")

    s1 = graph.input("s1_nuclease", "S1 nuclease (1:100), 30 min room temperature; phenol–chloroform stop",
                     MaterialKind.CATALYST)
    trimmed = graph.product("trimmed", "S1-trimmed crosslinked duplexes with free RNA ends",
                            MaterialKind.MOLECULAR, state=crosslinked.state)
    graph.react("s1_digest", "Limited S1 single-strand digestion of 400 ng RNA",
                StepKind.DIGESTION,
                ((depleted, InputRole.TEMPLATE), (s1, InputRole.CATALYST)), (trimmed,))
    add("s1_nuclease", "s1_digest")

    ligase = graph.input("proximity_ligase",
                         "circRNA ligase with ATP (85 °C denaturation, 1 h 60 °C)",
                         MaterialKind.CATALYST)
    chimera = graph.product("chimera", "crosslinked, proximity-ligated RNA chimera",
                            MaterialKind.MOLECULAR, state=topology.chimera())
    graph.react("proximity_ligation", "Ligate proximal free RNA ends of 50 ng RNA",
                StepKind.PROXIMITY_LIGATION,
                ((trimmed, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (chimera,))
    add("proximity_ligase", "proximity_ligation")

    rnase_r = graph.input("rnase_r", "RNase R, 10 min 37 °C; phenol–chloroform, ethanol",
                          MaterialKind.CATALYST)
    enriched = graph.product("enriched", "RNase R-resistant crosslinked or circular RNA",
                             MaterialKind.MOLECULAR, state=chimera.state)
    linear = graph.product("linear_digested", "RNase R-degraded uncrosslinked linear RNA",
                           MaterialKind.WASTE)
    graph.react("rnase_r_digest", "Digest uncrosslinked RNA with RNase R",
                StepKind.DIGESTION,
                ((chimera, InputRole.TEMPLATE), (rnase_r, InputRole.CATALYST)),
                (enriched, linear))
    graph.discard(linear, "nucleotides from RNase R degradation are removed by phenol–chloroform extraction and ethanol precipitation")
    add("rnase_r", "rnase_r_digest", "linear_digested")

    uv254 = graph.input("photoreversal", "254 nm UV 10 min on ice; ethanol precipitation removes AMT",
                        MaterialKind.CONDITION)
    reversed_chimera = graph.product("reversed_chimera", "photoreversed RNA chimera",
                                     MaterialKind.MOLECULAR, state=molecule(
                                         "photoreversed RNA chimera",
                                         (role("RNA partner A fragment"), role("RNA partner B fragment")),
                                         junctions=(("RNA partner A fragment", "RNA partner B fragment",
                                                     "RNA proximity ligation"),)))
    graph.react("photoreverse", "Reverse AMT crosslinks at 254 nm",
                StepKind.CHEMICAL_CONVERSION,
                ((enriched, InputRole.TEMPLATE), (uv254, InputRole.CONDITION)),
                (reversed_chimera,))
    add("photoreversal", "photoreverse")

    rt_primer = graph.input("rt_primer", "SMARTer first-strand primer adding a 5-prime universal sequence (sequence not printed)",
                            MaterialKind.OLIGO)
    rt = graph.input("rt_enzyme", "SMARTer reverse transcriptase that adds non-templated nucleotides",
                     MaterialKind.CATALYST)
    first = graph.product("first_strand", "first-strand cDNA with universal 5-prime sequence and added 3-prime nucleotides",
                          MaterialKind.MOLECULAR, state=molecule("first-strand cDNA", (
                              role("universal sequence"), role("chimera cDNA"),
                              role("non-templated nucleotides"))))
    graph.react("reverse_transcribe", "First-strand cDNA synthesis", StepKind.REVERSE_TRANSCRIPTION,
                ((reversed_chimera, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (first,))
    add("rt_primer", "rt_enzyme", "reverse_transcribe")

    forward = graph.input("universal_forward", "SMARTer universal forward primer annealing to the added nucleotides (sequence not printed)",
                          MaterialKind.OLIGO)
    double = graph.product("ds_cdna", "double-stranded SMARTer cDNA", MaterialKind.MOLECULAR,
                           state=molecule("ds cDNA", (role("universal forward sequence"),
                                                      role("chimera cDNA"), role("universal sequence"))))
    graph.react("second_strand", "Second-strand synthesis from the universal forward primer",
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)),
                (double,))
    add("universal_forward", "second_strand")

    ampure = graph.input("ampure_rt", "AMPure SPRI 1:2 sample:beads", MaterialKind.SELECTION)
    cleaned = graph.product("cleaned_cdna", "bead-recovered cDNA", MaterialKind.MOLECULAR,
                            state=double.state)
    sup = graph.product("rt_supernatant", "SPRI supernatant", MaterialKind.WASTE)
    graph.react("cleanup_rt", "SPRI cleanup of the RT reaction", StepKind.CLEANUP,
                ((double, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)), (cleaned, sup))
    graph.discard(sup, "unbound material removed with the bead supernatant")
    add("ampure_rt", "cleanup_rt", "rt_supernatant")

    pcr_f = graph.input("pcr_forward", "SMARTer forward PCR primer (sequence not printed)", MaterialKind.OLIGO)
    pcr_r = graph.input("pcr_reverse", "SMARTer indexed reverse PCR primer (index not stated)", MaterialKind.OLIGO)
    pol = graph.input("pcr_polymerase", "SMARTer kit PCR polymerase (20 cycles rep 1, 17 cycles rep 2)",
                      MaterialKind.CATALYST)
    amplified = graph.product("amplified", "indexed SMARTer amplicons", MaterialKind.MOLECULAR,
                              state=molecule("amplified library", tuple(FINAL_LIBRARY)))
    graph.react("library_pcr", "Amplify and index the library", StepKind.PCR,
                ((cleaned, InputRole.TEMPLATE), (pcr_f, InputRole.STRUCTURAL),
                 (pcr_r, InputRole.STRUCTURAL), (pol, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY, reason="SMARTer kit primer sequences are not printed in the inspected sources."))
    add("pcr_forward", "pcr_reverse", "pcr_polymerase", "library_pcr")

    gel = graph.input("page_select", "AMPure 1:2, then 6% TBE PAGE; excise >200 nt, elute overnight",
                      MaterialKind.SELECTION)
    library = graph.product("library", "LIGR-seq sequencing library", MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY)
    short = graph.product("short_products", "products below 200 nt left in the gel",
                          MaterialKind.WASTE)
    graph.react("size_select", "Bead cleanup and PAGE size selection", StepKind.SELECTION,
                ((amplified, InputRole.TEMPLATE), (gel, InputRole.SELECTION)), (library, short))
    graph.discard(short, "not excised from the PAGE gel")
    add("page_select", "size_select", "short_products", "library")

    graph.sequence(library, SequencingPlan(limitation=SEQUENCING))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


SPEC.protocol_graph = protocol_graph
