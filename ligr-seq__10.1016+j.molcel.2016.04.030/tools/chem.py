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

from types import SimpleNamespace
import seqprimers as sp
from chemdraw import Construct, Scene, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind, Disposition, Terminal, SourceClaim, SourceDiscrepancy)
from protocol_model import molecule, role

RANDOM_BARCODE = feature("ligr.random_barcode", "inline_barcode", "random",
                         note="5-prime random barcode trimmed from single-end reads; "
                              "length and origin within the SMARTer kit not stated")
INDEX = feature("ligr.index", "sample_index", "unknown",
                note="barcode index carried by the kit reverse PCR primer; selected index not stated")

FINAL_LIBRARY = Construct([
    role("SMARTer forward-primer end; unassigned", inferred=True),
    role("5-prime random barcode", molecular_feature=RANDOM_BARCODE),
    role("RNA partner A-derived cDNA"),
    role("RNA partner B-derived cDNA"),
    role("indexed SMARTer reverse-primer end; internal index position unassigned", molecular_feature=INDEX, inferred=True),
], name="LIGR-seq SMARTer stranded library (role endpoint)")

SEQUENCING = ("Single-end, strand-specific reads; the 5-prime random barcode is trimmed "
              "before Aligater chimera detection. Read length, platform primer sequences "
              "and the selected index remain unassigned in this bounded model. Partner order "
              "in a read depends on which RNA end was ligated, and circular ligation "
              "products can appear as concatemers of the two RNAs (Figure S1C).")

TITLE = "LIGR-seq — ligation of interacting RNA"
NOTES = "01_ligr_seq.html"
SOURCE = '<a href="https://doi.org/10.1016/j.molcel.2016.04.030">Sharma et al., Molecular Cell (2016)</a>, Experimental Procedures and Supplemental Experimental Procedures.'
SUMMARY = "AMT holds RNA partners through S1 trimming and proximity ligation; RNase R enrichment precedes photoreversal and SMARTer stranded library construction."
CAVEAT = "INFERRED ends — the native SMARTer kit version, molecular ends and primer species remain unassigned in this bounded model. RNA partner lengths and crosslinked pairing extent are unassigned. Linear and circular ligation classes may coexist; their proportions and RT products are not fixed. The 5′ random barcode is source-reported read content, not a proven UMI."
SEQ_PRIMERS = ()
READ_LENGTHS = None
ACQUIRED_ROLES = ("Read 1",)
READ_SPECS = (sp.read_spec("Read 1", "unassigned", "Source-reported 5′ random barcode followed by chimera-derived sequence; partner order and possible concatemers vary", "extract_then_include"),)
SEQUENCING_INTRO = SEQUENCING
SEQUENCING_ENDING = SEQUENCING
SEQUENCING_UNAVAILABLE = "Kit version, primer sites, read length and actual index acquisition remain unassigned; a compatible later kit does not establish them."
FINAL_CAPTION = "INFERRED ends — role-level SMARTer endpoint; no fixed partner lengths, index length or primer geometry is assigned."


def crosslinked_state():
    sc = Scene()
    for name in ("RNA partner A", "RNA partner B"):
        sc.strand(name, [role(name + "; sequence and paired extent unassigned")], label=name)
        sc.note(name, "AMT interstrand crosslink retains the partners; no base-resolved fold assigned")
    return sc


def chimera_state():
    sc = Scene()
    sc.strand("linear chimera", [role("RNA partner A fragment"), role("RNA partner B fragment")], label="representative ligated RNA class")
    sc.junction("linear chimera", "RNA partner A fragment", "RNA partner B fragment", "proximity-ligation junction; bases unassigned")
    sc.strand("circular class", [role("circular ligation products; native strand order/closure unassigned")], label="separate source-reported circular class; not a linear continuation")
    sc.note("circular class", "Figure S1A/C: circular products resist RNase R and can yield concatenated RNA-partner alignments")
    return sc


def protocol_graph():
    source = SOURCE
    graph = ProtocolGraph(TITLE, authority=source)
    rna = graph.input("rna", "293T cellular RNA containing paired regions",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("cellular RNA", (role("RNA molecules"),)), primary=True)
    amt = graph.input("crosslinker", "20 µg/ml AMT, 365 nm UV 30 min on ice",
                      MaterialKind.CONDITION)
    crosslinked = graph.product("crosslinked", "AMT-crosslinked RNA duplexes in cells",
                                MaterialKind.MOLECULAR, state=crosslinked_state())
    graph.react("crosslink", "Photocrosslink paired RNA in vivo with AMT",
                StepKind.CROSS_LINKING,
                ((rna, InputRole.TEMPLATE), (amt, InputRole.CONDITION)), (crosslinked,))

    extraction = graph.input("tri_dnase", "TRI-reagent extraction and TURBO DNase I",
                             MaterialKind.SELECTION)
    total = graph.product("total_rna", "crosslinked total RNA", MaterialKind.MOLECULAR,
                          state=crosslinked.state)
    dna_waste = graph.product("dna_waste", "DNase-derived DNA material; later recovery unassigned",
                              MaterialKind.MOLECULAR)
    graph.react("extract", "Extract RNA and digest DNA", StepKind.EXTRACTION,
                ((crosslinked, InputRole.TEMPLATE), (extraction, InputRole.SELECTION)),
                (total, dna_waste))
    graph.mark_terminal(Terminal(dna_waste.id, Disposition.UNRESOLVED, "DNA hydrolysis does not establish physical disposal; the exact partition through subsequent RNA recovery is unassigned."))

    ribozero = graph.input("ribozero", "Ribo-Zero Gold rRNA depletion (4 µg input)",
                           MaterialKind.SELECTION)
    depleted = graph.product("depleted", "rRNA-depleted crosslinked RNA",
                             MaterialKind.MOLECULAR, state=crosslinked.state)
    rrna = graph.product("rrna", "captured rRNA", MaterialKind.WASTE)
    graph.react("deplete_rrna", "Deplete ribosomal RNA", StepKind.SELECTION,
                ((total, InputRole.TEMPLATE), (ribozero, InputRole.SELECTION)),
                (depleted, rrna))
    graph.discard(rrna, "removed by Ribo-Zero Gold capture")

    s1 = graph.input("s1_nuclease", "S1 nuclease (1:100), 30 min room temperature; phenol–chloroform stop",
                     MaterialKind.CATALYST)
    trimmed = graph.product("trimmed", "S1-trimmed crosslinked duplexes with free RNA ends",
                            MaterialKind.MOLECULAR, state=crosslinked.state)
    s1_products = graph.product("s1_products", "RNA material hydrolysed during S1 trimming", MaterialKind.MOLECULAR)
    graph.react("s1_digest", "Limited S1 single-strand digestion of 400 ng RNA",
                StepKind.DIGESTION,
                ((depleted, InputRole.TEMPLATE), (s1, InputRole.CATALYST)), (trimmed, s1_products))
    graph.mark_terminal(Terminal(s1_products.id, Disposition.UNRESOLVED,
        "S1 hydrolysis does not establish physical disposal; partition through the source phenol–chloroform recovery is unassigned."))

    ligase = graph.input("proximity_ligase",
                         "circRNA ligase with ATP (85 °C denaturation, 1 h 60 °C)",
                         MaterialKind.CATALYST)
    chimera = graph.product("chimera", "crosslinked, proximity-ligated RNA chimera",
                            MaterialKind.MOLECULAR, state=chimera_state())
    unligated = graph.product("unligated", "RNA outside the tracked chimera class", MaterialKind.MOLECULAR)
    graph.react("proximity_ligation", "Ligate proximal free RNA ends of 50 ng RNA",
                StepKind.PROXIMITY_LIGATION,
                ((trimmed, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (chimera, unligated))
    graph.mark_terminal(Terminal(unligated.id, Disposition.UNRESOLVED,
        "These classes are accounted separately, not physically purified here. RNaseR resistance and recovery of nonchimeric RNA through later steps remain unassigned."))

    rnase_r = graph.input("rnase_r", "RNase R, 10 min 37 °C; phenol–chloroform, ethanol",
                          MaterialKind.CATALYST)
    enriched = graph.product("enriched", "RNase R-resistant crosslinked or circular RNA",
                             MaterialKind.MOLECULAR, state=chimera.state)
    linear = graph.product("linear_digested", "RNase R hydrolysis products from susceptible linear RNA",
                           MaterialKind.MOLECULAR)
    graph.react("rnase_r_digest", "Digest uncrosslinked RNA with RNase R",
                StepKind.DIGESTION,
                ((chimera, InputRole.TEMPLATE), (rnase_r, InputRole.CATALYST)),
                (enriched, linear))
    graph.mark_terminal(Terminal(linear.id, Disposition.UNRESOLVED, "RNase R hydrolyses RNA. Phenol–chloroform and ethanol recovery follow, but quantitative exclusion of these molecular products is not modelled."))

    uv254 = graph.input("photoreversal", "254 nm UV 10 min on ice; ethanol precipitation removes AMT",
                        MaterialKind.CONDITION)
    reversed_chimera = graph.product("reversed_chimera", "photoreversed RNA chimera",
                                     MaterialKind.MOLECULAR, state=chimera_state())
    graph.react("photoreverse", "Reverse AMT crosslinks at 254 nm",
                StepKind.CHEMICAL_CONVERSION,
                ((enriched, InputRole.TEMPLATE), (uv254, InputRole.CONDITION)),
                (reversed_chimera,))

    rt_primer = graph.input("rt_primer", "SMARTer first-strand primer adding a 5-prime universal sequence (sequence not printed)",
                            MaterialKind.OLIGO)
    rt = graph.input("rt_enzyme", "SMARTer reverse transcriptase that adds non-templated nucleotides",
                     MaterialKind.CATALYST)
    first = graph.product("first_strand", "first-strand cDNA with universal 5-prime sequence and added 3-prime nucleotides",
                          MaterialKind.MOLECULAR, state=molecule("first-strand cDNA", (
                              role("universal sequence"), role("chimera cDNA"),
                              role("non-templated nucleotides"))))
    residual_rna = graph.product("rt_rna", "RNA templates after RT; fate unassigned", MaterialKind.MOLECULAR, state=reversed_chimera.state)
    graph.react("reverse_transcribe", "First-strand cDNA synthesis", StepKind.REVERSE_TRANSCRIPTION,
                ((reversed_chimera, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (first, residual_rna))
    graph.mark_terminal(Terminal(residual_rna.id, Disposition.UNRESOLVED, "RT does not physically remove its RNA templates; their recovery through the kit cleanup remains unassigned."))

    forward = graph.input("universal_forward", "SMARTer universal forward primer annealing to the added nucleotides (sequence not printed)",
                          MaterialKind.OLIGO)
    ds_parts = (role("universal forward sequence"), role("chimera cDNA"), role("universal sequence"))
    double = graph.product("ds_cdna", "double-stranded SMARTer cDNA", MaterialKind.MOLECULAR,
                           state=Scene.duplex(ds_parts, unpaired=tuple(x.name + "'" for x in ds_parts)))
    graph.react("second_strand", "Second-strand synthesis from the universal forward primer",
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)),
                (double,))

    ampure = graph.input("ampure_rt", "AMPure SPRI 1:2 sample:beads", MaterialKind.SELECTION)
    cleaned = graph.product("cleaned_cdna", "bead-recovered cDNA", MaterialKind.MOLECULAR,
                            state=double.state)
    sup = graph.product("rt_supernatant", "SPRI supernatant", MaterialKind.WASTE)
    graph.react("cleanup_rt", "SPRI cleanup of the RT reaction", StepKind.CLEANUP,
                ((double, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)), (cleaned, sup))
    graph.discard(sup, "unbound material removed with the bead supernatant")

    pcr_f = graph.input("pcr_forward", "SMARTer forward PCR primer (sequence not printed)", MaterialKind.OLIGO)
    pcr_r = graph.input("pcr_reverse", "SMARTer indexed reverse PCR primer (index not stated)", MaterialKind.OLIGO)
    pol = graph.input("pcr_polymerase", "SMARTer kit PCR polymerase (20 cycles rep 1, 17 cycles rep 2)",
                      MaterialKind.CATALYST)
    amplified = graph.product("amplified", "indexed SMARTer amplicons", MaterialKind.MOLECULAR,
                              state=molecule("amplified library", tuple(FINAL_LIBRARY)))
    graph.react("library_pcr", "INFERRED ends — amplify and index the library", StepKind.PCR,
                ((cleaned, InputRole.TEMPLATE), (pcr_f, InputRole.STRUCTURAL),
                 (pcr_r, InputRole.STRUCTURAL), (pol, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(
                    cleaned.state, reason="Actual preceding cDNA state; native kit primer strands and first-cycle copying remain unassigned, not a completed library template."))

    gel = graph.input("page_select", "AMPure 1:2, then 6% TBE PAGE; excise >200 nt, elute overnight",
                      MaterialKind.SELECTION)
    library = graph.product("library", "LIGR-seq sequencing library", MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY)
    short = graph.product("short_products", "products below 200 nt left in the gel",
                          MaterialKind.WASTE)
    graph.react("size_select", "INFERRED ends — bead cleanup and PAGE size selection", StepKind.SELECTION,
                ((amplified, InputRole.TEMPLATE), (gel, InputRole.SELECTION)), (library, short))
    graph.discard(short, "not excised from the PAGE gel")

    graph.sequence(library, SequencingPlan(library=FINAL_LIBRARY, acquired_roles=ACQUIRED_ROLES, limitation=SEQUENCING))
    graph.annotate(SourceDiscrepancy("cleanup_rt", "SPRI size wording and final PAGE window",
        (SourceClaim("Supplemental Experimental Procedures, Library Preparation, lines846–850", "1:2 sample:beads to obtain fragments <100 bp"),
         SourceClaim("Same section, lines851–853", "products above 200 nt excised from 6% TBE PAGE")),
        "Use the printed 1:2 bead ratio without assigning a SPRI size cutoff; retain the explicit final PAGE threshold above200nt.",
        "Both statements appear in the same published library-preparation paragraph; the model keeps their different scopes visible.",
        "The retained-versus-excluded meaning of the <100bp wording is unresolved; it is not silently corrected."))
    graph.declare_reviewed_scope(("rna", "crosslinker", "crosslink", "tri_dnase", "extract", "ribozero", "deplete_rrna", "s1_nuclease", "s1_digest", "proximity_ligase", "proximity_ligation", "rnase_r", "rnase_r_digest", "photoreversal", "photoreverse", "rt_primer", "rt_enzyme", "reverse_transcribe", "universal_forward", "second_strand", "ampure_rt", "cleanup_rt", "pcr_forward", "pcr_reverse", "pcr_polymerase", "library_pcr", "page_select", "size_select", "library"), source)
    return graph.require_publishable()


SPEC = SimpleNamespace(TITLE=TITLE, NOTES=NOTES, SOURCE=SOURCE, SUMMARY=SUMMARY,
    CAVEAT=CAVEAT, FINAL_LIBRARY=FINAL_LIBRARY, SEQ_PRIMERS=SEQ_PRIMERS,
    READ_LENGTHS=READ_LENGTHS, ACQUIRED_ROLES=ACQUIRED_ROLES, READ_SPECS=READ_SPECS,
    SEQUENCING_INTRO=SEQUENCING_INTRO, SEQUENCING_ENDING=SEQUENCING_ENDING,
    SEQUENCING_UNAVAILABLE=SEQUENCING_UNAVAILABLE, FINAL_CAPTION=FINAL_CAPTION,
    oligos=lambda: (), protocol_graph=protocol_graph)
