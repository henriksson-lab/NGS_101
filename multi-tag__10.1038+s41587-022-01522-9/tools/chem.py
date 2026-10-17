"""Protocol-owned source, molecular workflow and reaction graph.

Source scope (round 2026-10-17-breadth-first-13): published Meers et al. 2022,
PMC10188359 JATS Methods ("Conjugates for MulTI-Tag", "Bulk MulTI-Tag protocol",
"Single-cell MulTI-Tag", "Sequencing and data pre-processing") and published
Supplementary Table 1 (MOESM2; byte-identical to the local preprint media-2.xlsx).
The graph draws the single-cell branch; the bulk branch shares the conjugate
tagmentation rounds and replaces the 96-well i7 transposome / ICELL8 steps with one
i7 transposome and tube PCR.
"""
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role, role_library
from extended_protocols import Definition, ident

DEFINITION = Definition('MulTI-Tag', '10.1038/s41587-022-01522-9', 'Meers et al., Nature Biotechnology (2022)', 'Sequential antibody–adapter-conjugate pA–Tn5 tagmentation inserts target-specific i5 adapter barcodes; a secondary-antibody i7 tagmentation in 96 barcoded wells and ICELL8 microwell PCR add combinatorial cell identity.', 'Each genomic fragment carries a target barcode and a three-part combinatorial cell identity (i7 transposome well, ICELL8 i5 and i7 PCR primers).', '01_multi_tag.html')

# Identifier roles (Supplementary Table 1: every barcode column is eight bases).
TN5_WELL = ident("i7 transposome-well barcode", "i7_transposome_well", "cell_barcode", length=8)
I5_PCR = ident("ICELL8 i5 primer barcode", "icell8_i5_primer", "cell_barcode", length=8)
I7_PCR = ident("ICELL8 i7 primer barcode", "icell8_i7_primer", "cell_barcode", length=8)
TARGET_BC = ident("target barcode", "target_barcode", "feature_barcode", length=8)
INSERT = role("target-proximal genomic insert")

# Supplementary Table 1 example rows. Each is one member of a printed family; the
# experiment chooses the row (target, well or nanowell), so none is universal.
TN5MEREV = "CTGTCTCTTATACACATCT"            # README: [phos]CTGTCTCTTATACACATCT
P5_I5_1 = "TCGTCGGCAGCGTCTCCACGCTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG"  # H3K27me3, TATAGCCT
P7_I7_13 = "GTCTCGTGGGCTCGGCTGTCCCTGTCCGTCTACATCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG"  # first single-cell row, GTCTACAT
SCI_AD1_1 = "AATGATACGGCGACCACCGAGATCTACACTAGATCGCTCGTCGGCAGCGTCTCCACGC"  # i5 primer, TAGATCGC
SCI_AD2_1 = "CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGCTGTCCCTGTCC"   # i7 primer, TCGCCTTA

TITLE = DEFINITION.title
NOTES = DEFINITION.notes
SOURCE = f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.citation}</a>.'
SUMMARY = DEFINITION.summary
CAVEAT = ('Published Supplementary Table 1 prints the target-coded P5_i5 and '
          'well-coded P7_i7 adapters (both annealed to 5′-phosphorylated Tn5MErev) '
          'and the 72 i5 and 72 i7 amplification primers; the oligos shown are '
          'example rows of those families. The diagram draws the single-cell branch '
          'at role level: the paper does not state which sequencing primers or '
          'index-read lengths acquired the target and cell barcodes, so no read '
          'landing is drawn. Paired-end 25 × 25 HiSeq reads are published.')
FINAL_LIBRARY = role_library(TITLE + " final library roles",
                             (TARGET_BC, INSERT, TN5_WELL, I7_PCR), index=I5_PCR)
SEQ_PRIMERS = ()
FINAL_CAPTION = ('Role order: P5 end – ICELL8 i5 primer barcode – P5_i5 target barcode – '
                 'mosaic end – insert – mosaic end – P7_i7 transposome-well barcode – '
                 'ICELL8 i7 primer barcode – P7 end. Barcode order derived from the '
                 'Supplementary Table 1 oligo layouts.')
SEQUENCING_INTRO = DEFINITION.endpoint
READ_LENGTHS = None
SEQUENCING_ENDING = ('Illumina HiSeq, paired-end 25 × 25 reads (Methods, "Sequencing and '
                     'data pre-processing"). Index-read lengths and the primers that read '
                     'the target and cell barcodes are not stated in the inspected sources.')
SEQUENCING_UNAVAILABLE = ("The adapters place the target barcode between "
                          "TCGTCGGCAGCGTCTCCACGC and the mosaic end, so a standard Nextera "
                          "Read 1 primer is not implied; the paper does not name the "
                          "sequencing primers used.")


def workflow():
    return chain(molecule("bead-bound nuclei with chromatin targets",
                          (role("target 1"), role("target 2"), role("target 3"))), (
        ("Assemble antibody–i5 adapter conjugate transposome; bind, wash, tagment (Mg2+), EDTA wash",
         molecule("first target tagged", (TARGET_BC, role("target-1 chromatin fragment"))),
         "The conjugated P5_i5 adapter barcode identifies the antibody target."),
        ("Repeat conjugate binding and tagmentation for each further target",
         molecule("multiplex i5-tagged chromatin", (TARGET_BC, role("target-proximal genomic DNA"))),
         "Order chosen empirically; H3K27me3 conjugate last in every published order."),
        ("Secondary antibody, then i7 transposome tagmentation in 96 barcoded wells",
         molecule("i5/i7-tagged chromatin", (TARGET_BC, INSERT, TN5_WELL)),
         "pA–Tn5 loaded only with P7_i7 adapters adds the reverse adapter and a well barcode."),
        ("Pool, DNase, strain, dispense into ICELL8 nanowells, SDS/Triton, PCR",
         molecule("MulTI-Tag library", (I5_PCR, TARGET_BC, INSERT, TN5_WELL, I7_PCR)),
         "i5 and i7 primer barcodes per nanowell complete combinatorial cell identity."),
    ))


def oligos():
    return ('Tn5MErev — 5′-/5Phos/CTGTCTCTTATACACATCT-3′',
            'P5_i5_1 (H3K27me3; 5′-aminated for conjugation, modification not printed in the table) — '
            '5′-' + P5_I5_1 + '-3′',
            'P7_i7_13 (single-cell well adapter, example row) — 5′-' + P7_I7_13 + '-3′',
            'sci_Ad1.1_TAGATCGC (i5 primer, example row) — 5′-' + SCI_AD1_1 + '-3′',
            'sci_Ad2.1_TAAGGCGA (i7 primer, example row) — 5′-' + SCI_AD2_1 + '-3′')


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    source = f"{d.citation}; doi:{d.doi}; PMC10188359 Methods and Supplementary Table 1"
    graph = ProtocolGraph(d.title, authority=source)
    reviewed: list[str] = []

    def inp(identifier, label, kind, **kw):
        reviewed.append(identifier)
        return graph.input(identifier, label, kind, **kw)

    def step(identifier, label, kind, inputs, outputs, **kw):
        reviewed.append(identifier)
        return graph.react(identifier, label, kind, inputs, outputs, **kw)

    def waste(identifier, label, reason):
        item = graph.product(identifier, label, MaterialKind.WASTE)
        return item, reason

    # Adapter preparation ("Conjugates for MulTI-Tag").
    i5_oligo = inp("i5_adapter", "P5_i5 target adapter, 5′-aminated (example row P5_i5_1, H3K27me3)",
                   MaterialKind.OLIGO, state=P5_I5_1)
    merev_i5 = inp("tn5merev_i5", "Tn5MErev, 5′-phosphorylated (for i5 adapters)",
                   MaterialKind.OLIGO, state=TN5MEREV)
    annealed_i5 = graph.product("annealed_i5", "Tn5MErev-annealed P5_i5 adapter (100 µM, 1× PBS)",
                                MaterialKind.MOLECULAR,
                                state=molecule("annealed i5 adapter", (TARGET_BC, role("mosaic-end duplex"))))
    step("anneal_i5", "Anneal 200 µM P5_i5 adapter to equimolar Tn5MErev in 1× PBS",
         StepKind.HYBRIDIZATION, ((i5_oligo, InputRole.TEMPLATE), (merev_i5, InputRole.STRUCTURAL)),
         (annealed_i5,))
    for_conj = graph.product("i5_for_conjugation", "annealed i5 adapter for antibody conjugation",
                             MaterialKind.MOLECULAR, state=annealed_i5.state)
    free_i5 = graph.product("i5_free", "unconjugated annealed i5 adapter of the same sequence (16 pmol per transposome)",
                            MaterialKind.MOLECULAR, state=annealed_i5.state)
    step("split_i5", "Reserve unconjugated annealed adapter", StepKind.SPLIT,
         ((annealed_i5, InputRole.TEMPLATE),), (for_conj, free_i5), changes_molecule=False)
    antibody = inp("primary_antibodies", "target-specific primary antibodies (PBS; no BSA, azide or glycerol)",
                   MaterialKind.SELECTION)
    conj_kit = inp("conjugation_kit", "Abcam Oligonucleotide Conjugation Kit; random amine conjugation, ~10:1 adapter excess",
                   MaterialKind.CONDITION)
    conjugates = graph.product("conjugates", "antibody–i5 adapter conjugates (one target barcode per antibody)",
                               MaterialKind.MOLECULAR,
                               state=molecule("antibody–adapter conjugate", (role("primary antibody"), TARGET_BC)))
    step("conjugate", "Conjugate 5′-amine annealed adapters to 100 µg primary antibody",
         StepKind.CHEMICAL_CONVERSION,
         ((for_conj, InputRole.TEMPLATE), (antibody, InputRole.SELECTION),
          (conj_kit, InputRole.CONDITION)), (conjugates,))
    first_conj = graph.product("first_conjugate", "first target conjugate", MaterialKind.MOLECULAR,
                               state=conjugates.state)
    later_conj = graph.product("later_conjugates", "remaining target conjugates", MaterialKind.MOLECULAR,
                               state=conjugates.state)
    step("split_conjugates", "Use one conjugate per sequential round", StepKind.SPLIT,
         ((conjugates, InputRole.TEMPLATE),), (first_conj, later_conj), changes_molecule=False)
    first_free = graph.product("first_free_i5", "free i5 adapter for the first round", MaterialKind.MOLECULAR,
                               state=annealed_i5.state)
    later_free = graph.product("later_free_i5", "free i5 adapter for later rounds", MaterialKind.MOLECULAR,
                               state=annealed_i5.state)
    step("split_free_i5", "Match free adapter to each conjugate", StepKind.SPLIT,
         ((free_i5, InputRole.TEMPLATE),), (first_free, later_free), changes_molecule=False)

    # Nuclei on WGA beads ("Single-cell MulTI-Tag").
    nuclei = graph.input("nuclei", "nuclei (1.8 million for single-cell)", MaterialKind.BIOLOGICAL,
                         state=molecule("nuclei", (role("target 1"), role("target 2"), role("target 3"))),
                         primary=True)
    reviewed.append("nuclei")
    beads = inp("wga_beads", "biotin-WGA on streptavidin T1 Dynabeads (bulk: ConA beads)", MaterialKind.SELECTION)
    on_beads = graph.product("bead_nuclei", "bead-bound nuclei", MaterialKind.BIOLOGICAL, state=nuclei.state)
    step("bind_beads", "Pre-bind nuclei to WGA beads", StepKind.BINDING,
         ((nuclei, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (on_beads,))

    def conjugate_round(prefix, template, conj, free, label):
        tn5 = inp(prefix + "_patn5", "pA–Tn5 (~5 µM, 1 µl per conjugate)", MaterialKind.CATALYST)
        mg = inp(prefix + "_mg", "Wash-300 + 2 mM EDTA binding; Wash-300 + 10 mM MgCl2, 1 h 37 °C", MaterialKind.CONDITION)
        tagged = graph.product(prefix + "_tagged", label + " i5-tagged chromatin", MaterialKind.MOLECULAR,
                               state=molecule(label + " tagged chromatin", (TARGET_BC, role("target-proximal genomic DNA"))))
        unbound, reason = waste(prefix + "_unbound", label + " unbound conjugate transposome (Wash-300 washes, then 5 mM EDTA wash)",
                                "removed in the three Wash-300 washes and the post-tagmentation EDTA wash")
        step(prefix + "_tagment",
             "Assemble " + label + " conjugate transposome (conjugate + free adapter + pA–Tn5, 30–60 min RT), "
             "bind 1 h RT or overnight 4 °C, wash ×3, tagment, EDTA wash",
             StepKind.TAGMENTATION,
             ((template, InputRole.TEMPLATE), (conj, InputRole.STRUCTURAL), (free, InputRole.STRUCTURAL),
              (tn5, InputRole.CATALYST), (mg, InputRole.CONDITION)), (tagged, unbound))
        graph.discard(unbound, reason)
        return tagged

    first = conjugate_round("first", on_beads, first_conj, first_free, "first-target")
    multiplex = conjugate_round("later", first, later_conj, later_free, "remaining-target")

    secondary = inp("secondary_antibodies", "species-matched secondary antibodies, 1 h RT in Wash Buffer",
                    MaterialKind.SELECTION)
    with_secondary = graph.product("secondary_bound", "secondary-antibody-bound tagged chromatin",
                                   MaterialKind.MOLECULAR, state=multiplex.state)
    sec_waste, sec_reason = waste("secondary_unbound", "unbound secondary antibody", "removed by two Wash Buffer washes")
    step("bind_secondary", "Bind secondary antibodies to the primary conjugates", StepKind.BINDING,
         ((multiplex, InputRole.TEMPLATE), (secondary, InputRole.SELECTION)), (with_secondary, sec_waste))
    graph.discard(sec_waste, sec_reason)

    i7_oligo = inp("i7_adapter", "P7_i7 well adapter (96 single-cell rows; example P7_i7_13)",
                   MaterialKind.OLIGO, state=P7_I7_13)
    merev_i7 = inp("tn5merev_i7", "Tn5MErev, 5′-phosphorylated (for i7 adapters)", MaterialKind.OLIGO, state=TN5MEREV)
    annealed_i7 = graph.product("annealed_i7", "Tn5MErev-annealed P7_i7 adapter", MaterialKind.MOLECULAR,
                                state=molecule("annealed i7 adapter", (role("mosaic-end duplex"), TN5_WELL)))
    step("anneal_i7", "Anneal P7_i7 adapter to Tn5MErev", StepKind.HYBRIDIZATION,
         ((i7_oligo, InputRole.TEMPLATE), (merev_i7, InputRole.STRUCTURAL)), (annealed_i7,))
    i7_tn5 = inp("i7_patn5", "pA–Tn5 loaded only with i7 adapter (1 µl + 32 pmol adapter per transposome)",
                 MaterialKind.CATALYST)
    i7_mg = inp("i7_mg", "nuclei split equally across 96 i7 transposome wells; 1 h RT in Wash-300, wash ×3, 10 mM MgCl2 1 h 37 °C",
                MaterialKind.CONDITION)
    dual = graph.product("dual_tagged", "i5/i7-tagged chromatin carrying an i7 well barcode", MaterialKind.MOLECULAR,
                         state=molecule("dual-tagged chromatin", (TARGET_BC, INSERT, TN5_WELL)))
    i7_waste, i7_reason = waste("i7_unbound", "unbound i7 transposome", "removed by three Wash-300 washes")
    step("i7_tagment", "Secondary CUT&Tag with i7 transposome in 96 barcoded wells, then re-aggregate",
         StepKind.TAGMENTATION,
         ((with_secondary, InputRole.TEMPLATE), (annealed_i7, InputRole.STRUCTURAL),
          (i7_tn5, InputRole.CATALYST), (i7_mg, InputRole.CONDITION)), (dual, i7_waste))
    graph.discard(i7_waste, i7_reason)

    dnase = inp("rq1_dnase", "RQ1 DNase, 10 min at 0 °C, then RQ1 Stop Buffer", MaterialKind.CATALYST)
    treated = graph.product("dnase_nuclei", "DNase-treated pooled nuclei", MaterialKind.MOLECULAR, state=dual.state)
    digest, digest_reason = waste("dnase_digest", "DNase-digested extranuclear DNA",
                                  "digested and washed away with 10 mM TAPS")
    step("dnase", "DNase-treat pooled nuclei on a cold block", StepKind.DIGESTION,
         ((dual, InputRole.TEMPLATE), (dnase, InputRole.CATALYST)), (treated, digest))
    graph.discard(digest, digest_reason)
    strainer = inp("strainers", "two 20-µm cell strainers, 300 g 10 min", MaterialKind.SELECTION)
    strained = graph.product("strained_nuclei", "strained nuclei", MaterialKind.MOLECULAR, state=dual.state)
    clumps, clump_reason = waste("strainer_retained", "material retained on the strainers", "not passed through the strainers")
    step("strain", "Force nuclei through 20-µm strainers", StepKind.SELECTION,
         ((treated, InputRole.TEMPLATE), (strainer, InputRole.SELECTION)), (strained, clumps))
    graph.discard(clumps, clump_reason)
    dispense = inp("icell8", "DAPI + ICELL8 Second Diluent; 857 nuclei/µl, 640 µl dispensed into an ICELL8 chip",
                   MaterialKind.CONDITION)
    in_wells = graph.product("nanowell_nuclei", "nuclei in ICELL8 nanowells", MaterialKind.MOLECULAR, state=dual.state)
    undispensed, und_reason = waste("not_dispensed", "nuclei not dispensed into nanowells",
                                    "not dispensed into the ICELL8 chip")
    step("dispense", "Dispense nuclei into an ICELL8 microfluidic chip", StepKind.PARTITION,
         ((strained, InputRole.TEMPLATE), (dispense, InputRole.CONDITION)), (in_wells, undispensed),
         changes_molecule=False)
    graph.discard(undispensed, und_reason)
    sds = inp("sds_triton", "SDS denaturation then Triton X-100 neutralization in microwells (bulk: TAPS-SDS 58 °C 1 h)",
              MaterialKind.CONDITION)
    released = graph.product("released_dna", "released tagmented DNA", MaterialKind.MOLECULAR, state=dual.state)
    step("denature", "Release tagmented DNA with SDS; neutralize with Triton X-100", StepKind.DEPROTEINIZATION,
         ((in_wells, InputRole.TEMPLATE), (sds, InputRole.CONDITION)), (released,))

    fwd = inp("i5_primer", "i5 primer, one of 72 dispensed per nanowell (example sci_Ad1.1_TAGATCGC)",
              MaterialKind.OLIGO, state=SCI_AD1_1)
    rev = inp("i7_primer", "i7 primer, one of 72 dispensed per nanowell (example sci_Ad2.1_TAAGGCGA)",
              MaterialKind.OLIGO, state=SCI_AD2_1)
    pol = inp("pcr_mix", "PCR master mix (bulk: 2× NEBNext, 58 °C 5 min / 72 °C 5 min gap fill then cycling; "
              "single-cell microwell conditions cited to ref. 62)", MaterialKind.CATALYST)
    amplified = graph.product("amplified", "nanowell-indexed amplified library", MaterialKind.MOLECULAR, state=final)
    step("pcr", "Gap-fill and amplify in nanowells with i5/i7 primers", StepKind.PCR,
         ((released, InputRole.TEMPLATE), (fwd, InputRole.STRUCTURAL), (rev, InputRole.STRUCTURAL),
          (pol, InputRole.CATALYST)), (amplified,),
         geometry=PCRGeometry.unavailable(final, reason=(
             "Primer 3′ ends equal the 5′ regions of the adapters (i5: TCGTCGGCAGCGTCTCCACGC; "
             "i7: GTCTCGTGGGCTCGGCTGTCCCTGTCC), but the extension product is drawn at role "
             "level and the single-cell cycling is cited to ref. 62, not restated.")))
    beads2 = inp("highprep", "two rounds of 1.3× HighPrep beads; elute 20 µl 10 mM Tris-HCl/1 mM EDTA",
                 MaterialKind.SELECTION)
    library = graph.product("library", "MulTI-Tag sequencing library", MaterialKind.LIBRARY, state=final)
    sup, sup_reason = waste("bead_supernatant", "bead-cleanup supernatant", "not retained by 1.3× HighPrep beads")
    step("cleanup", "Re-aggregate microwell contents and bead-purify", StepKind.CLEANUP,
         ((amplified, InputRole.TEMPLATE), (beads2, InputRole.SELECTION)), (library, sup))
    graph.discard(sup, sup_reason)
    graph.sequence(library, SequencingPlan(limitation=SEQUENCING_ENDING))
    reviewed.append("library")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
