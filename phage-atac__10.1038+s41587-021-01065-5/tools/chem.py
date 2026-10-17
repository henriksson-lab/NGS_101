"""Protocol-owned source, molecular workflow and reaction graph.

Source authority: Fiskin et al., Nat. Biotechnol. 40, 374 (2022), author manuscript
PMC10549956 Methods (PHAGE-ATAC cell staining, library preparation, PDT library and
sequencing passages) and the published Supplementary Information (MOESM1) Tables 1
and 4.  The bioRxiv preprint (10.1101/2020.10.01.322420) is superseded and is not
used for active claims.
"""
from chemdraw import Construct, Segment, feature
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role
from extended_protocols import Definition, CELL, TARGET

DEFINITION = Definition(
    'PHAGE-ATAC', '10.1038/s41587-021-01065-5', 'Fiskin et al., Nature Biotechnology (2022)',
    'Nanobody-displaying phagemids carry 10x-compatible PAC tags so one droplet barcode '
    'joins protein-derived tags to ATAC and mitochondrial DNA.',
    'Separate protein-derived-tag (PDT) and ATAC libraries share the 10x cell barcode.',
    '01_phage_atac.html')

# Published Supplementary Table 1 (MOESM1 text lines 75-81), 5'->3', no modifications.
EF147 = "AATGATACGGCGACCACCGAGA"
EF91 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGATACCGCGGTGTATTATTGC"
EF91_HANDLE = "GATACCGCGGTGTATTATTGC"   # 3' end; the FR3 handle introduced by EF79

# Sample index identity: 10x i7 sample-index primers; their sequences and length are
# not printed in the paper.
SAMPLE_I7 = role("i7 sample index", tag="cbc",
                 molecular_feature=feature("sample_index", "sample_index", "unknown"))

TITLE = DEFINITION.title
NOTES = DEFINITION.notes
SOURCE = (f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.citation}</a> '
          '(author manuscript PMC10549956 Methods; published Supplementary Tables 1 and 4).')
SUMMARY = DEFINITION.summary
CAVEAT = ('Published Supplementary Table 1 prints the PDT-specific primers EF147 and EF91, '
          'and Supplementary Table 4 prints the CDR3 barcode readouts (read by R3 or R1). '
          'The 10x Chromium Single Cell ATAC gel-bead oligo, Amp-Mix, SI-PCR Primer B and '
          'i7 sample-index primers are named but not printed; the ATAC library itself is '
          'prepared "according to the Chromium Single Cell ATAC Solution protocol", whose '
          'version-matched vendor guide has not been read, so those reactions are drawn '
          'at role level. The order of regions in the PDT library below is '
          'INFERRED from the PAC-tag (RD1) and EF79/EF91 handle positions; it is not drawn '
          'as a published sequence.')

PDT_LIBRARY = Construct([
    role("P5 end from SI-PCR Primer B / EF147", inferred=True),
    CELL,
    role("Read 1 (PAC-tag RD1 captured by bead oligo)", inferred=True),
    TARGET,
    role("FR3 handle + ME/Read 2 (EF91 tail)", inferred=True),
    SAMPLE_I7,
    role("P7 end from i7 sample-index primer", inferred=True),
], name="PHAGE-ATAC PDT library (role order inferred)")
FINAL_LIBRARY = PDT_LIBRARY
SEQ_PRIMERS = ()
FINAL_CAPTION = ('INFERRED — PDT library role order. Supplementary Table 4 states the CDR3 '
                 'barcode is read in R3 (starting with the EF91 handle) for most clones and in '
                 'R1 for hashtags PH-C/PH-D. The ATAC library follows the 10x protocol and is '
                 'not drawn.')
READ_LENGTHS = {"Read 1": 34, "Read 2": 34}
ACQUIRED_ROLES = ("Read 1", "Read 2")
SEQUENCING_INTRO = ('PDT and ATAC libraries were pooled and paired-end sequenced, 2 × 34 cycles, '
                    'on a NextSeq 550 (High Output); index-read cycles are not stated. Demultiplexed '
                    'with CellRanger-ATAC mkfastq into the 10x scATAC R1/R2/R3 convention.')
SEQUENCING_ENDING = DEFINITION.endpoint
SEQUENCING_UNAVAILABLE = ("2 × 34 paired-end cycles are published; the 10x primer and index "
                          "sequences, and the acquired index cycles, are not printed in the paper.")


def workflow():
    return chain(molecule("cell stained with nanobody-displaying M13 phage", (role("surface or intracellular epitope"), TARGET)), (
        ("Fix (1% formaldehyde), lyse (0.1% NP-40) and tagment in bulk",
         molecule("phage-bound tagmented cell", (TARGET, role("PAC-tag phagemid"), role("ATAC fragments"))),
         "A milder lysis buffer replaces the 10x scATAC lysis buffer, which lost pNb staining."),
        ("Co-encapsulate with 10x ATAC gel beads; linear amplification",
         molecule("droplet-barcoded products", (CELL, TARGET, role("phagemid tag"), role("ATAC fragment"))),
         "The bead oligo captures the RD1 in the PAC tag alongside ATAC fragments."),
        ("Split 5 µl for PDT and 35 µl for ATAC libraries",
         molecule("PHAGE-ATAC readouts", (CELL, TARGET, role("ATAC or phage insert"))),
         "The common cell barcode links protein abundance, accessibility and mitochondrial variants."),
    ))


def oligos():
    return ('EF147 PDT PCR — 5′-' + EF147 + '-3′',
            'EF91 PDT PCR — 5′-' + EF91 + '-3′')


def _oligo(name, seq):
    return Construct([Segment(name, seq)], name=name)


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    source = (f"{d.citation}; doi:{d.doi}; PMC10549956 Methods (library preparation, "
              "PDT library, sequencing); published Supplementary Tables 1 and 4")
    vendor = ("Chromium Single Cell ATAC Solution protocol (10x Genomics), named in the "
              "Methods; version-matched guide not read in this pass")
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input(
        "stained_cells", "cells stained with nanobody-displaying M13 phage, fixed and lysed",
        MaterialKind.BIOLOGICAL,
        state=molecule("phage-stained cell", (
            role("surface or intracellular epitope"), TARGET, role("PAC-tag phagemid"))),
        primary=True)
    transposome = graph.input("atac_transposome", "10x ATAC transposition mix",
                              MaterialKind.ADAPTER, source=vendor)
    tn5 = graph.input("tn5", "ATAC Tn5 transposase", MaterialKind.CATALYST, source=vendor)
    tagmented = graph.product(
        "tagmented_cells", "phage-bound cells with tagmented accessible chromatin",
        MaterialKind.MOLECULAR,
        state=molecule("phage-bound tagmented cell", (
            TARGET, role("PAC-tag phagemid"), role("ATAC fragments"))))
    losses = graph.product("tn5_losses", "non-library ATAC products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment accessible chromatin (10x ATAC protocol)",
                StepKind.TAGMENTATION,
                ((cells, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagmented, losses), source=vendor)
    graph.discard(losses, "does not retain an amplifiable accessible-chromatin fragment")
    droplets = graph.input("droplet_condition", "10x ATAC chip and Chromium controller",
                           MaterialKind.CONDITION)
    occupied = graph.product("occupied_droplets", "cells co-encapsulated with ATAC gel beads",
                             MaterialKind.MOLECULAR, state=tagmented.state)
    empty = graph.product("empty_droplets", "empty or collision droplets", MaterialKind.WASTE)
    graph.react("partition", "Co-encapsulate cells and ATAC gel beads",
                StepKind.PARTITION,
                ((tagmented, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (occupied, empty), changes_molecule=False)
    graph.discard(empty, "does not contain one interpretable cell/bead combination")
    bead_primer = graph.input("bead_primer", "10x ATAC cell-barcode bead oligo (sequence not printed)",
                              MaterialKind.OLIGO, source=vendor)
    polymerase = graph.input("linear_polymerase", "10x droplet linear-amplification mix",
                             MaterialKind.CATALYST, source=vendor)
    barcoded = graph.product(
        "barcoded_products", "cell-barcoded ATAC and phagemid products",
        MaterialKind.MOLECULAR,
        state=molecule("droplet-barcoded products", (
            CELL, TARGET, role("phagemid tag"), role("ATAC fragment"))))
    graph.react("linear_amplification", "Linear PCR applies the droplet barcode to ATAC fragments and PAC-tag RD1",
                StepKind.EXTENSION,
                ((occupied, InputRole.TEMPLATE), (bead_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (barcoded,))
    silane = graph.input("silane_beads", "MyOne silane beads", MaterialKind.SELECTION)
    eluate = graph.product("silane_eluate", "barcoded products in 40 µl elution buffer I",
                           MaterialKind.MOLECULAR, state=barcoded.state)
    silane_waste = graph.product("silane_waste", "broken-emulsion supernatant", MaterialKind.WASTE)
    graph.react("silane_cleanup", "Break emulsion; MyOne silane cleanup, elute 40 µl",
                StepKind.CLEANUP, ((barcoded, InputRole.TEMPLATE), (silane, InputRole.SELECTION)),
                (eluate, silane_waste))
    graph.discard(silane_waste, "removed from bead-bound barcoded products")
    atac = graph.product("atac_branch", "35 µl eluate for scATAC library",
                         MaterialKind.MOLECULAR, state=barcoded.state)
    protein = graph.product("protein_branch", "5 µl eluate for PDT library",
                            MaterialKind.MOLECULAR, state=barcoded.state)
    graph.react("split_readouts", "Split 5 µl (PDT) and 35 µl (ATAC)",
                StepKind.SPLIT, ((eluate, InputRole.TEMPLATE),), (atac, protein),
                changes_molecule=False)

    # ATAC branch: vendor protocol, sequences not printed.
    atac_fwd = graph.input("atac_forward", "10x SI-PCR Primer B", MaterialKind.OLIGO, source=vendor)
    atac_rev = graph.input("atac_reverse", "10x i7 sample-index primer", MaterialKind.OLIGO, source=vendor)
    atac_pol = graph.input("atac_polymerase", "10x Amp Mix", MaterialKind.CATALYST, source=vendor)
    atac_lib = graph.product("atac_library", "PHAGE-ATAC chromatin sequencing library",
                             MaterialKind.LIBRARY,
                             state=Construct([role("P5 end"), CELL, role("Tn5 ATAC fragment"),
                                              SAMPLE_I7, role("P7 end")], name="scATAC library roles"))
    graph.react("atac_pcr", "Sample-index PCR per the Chromium Single Cell ATAC protocol",
                StepKind.PCR,
                ((atac, InputRole.TEMPLATE), (atac_fwd, InputRole.STRUCTURAL),
                 (atac_rev, InputRole.STRUCTURAL), (atac_pol, InputRole.CATALYST)), (atac_lib,),
                source=vendor,
                geometry=PCRGeometry.unavailable(
                    atac_lib.state,
                    reason="Performed per the 10x Chromium Single Cell ATAC protocol; the paper "
                           "prints no primer sequences and the vendor guide is not read here."))
    graph.sequence(atac_lib, SequencingPlan(limitation=SEQUENCING_UNAVAILABLE,
                                            acquired_roles=ACQUIRED_ROLES,
                                            read_lengths=READ_LENGTHS))

    # PDT branch: published Methods "PDT library preparation".
    ef147 = graph.input("pdt_forward", "EF147 (Supplementary Table 1)", MaterialKind.OLIGO,
                        state=_oligo("EF147", EF147))
    ef91 = graph.input("pdt_reverse", "EF91 (Supplementary Table 1)", MaterialKind.OLIGO,
                       state=_oligo("EF91", EF91))
    kapa = graph.input("pdt_polymerase", "2× KAPA polymerase", MaterialKind.CATALYST)
    pdt_amp = graph.product("pdt_amplicon", "PDT-specific amplicon", MaterialKind.MOLECULAR,
                            state=molecule("PDT amplicon", (
                                role("EF147 / P5 end", inferred=True), CELL,
                                role("Read 1 (PAC-tag RD1)", inferred=True), TARGET,
                                role("FR3 handle + EF91 ME/Read 2 tail", inferred=True))))
    graph.react("pdt_pcr", "PDT PCR with EF147/EF91: 95 °C 3 min; 20 × (95 °C 20 s, 60 °C 30 s, 72 °C 20 s); 72 °C 5 min",
                StepKind.PCR,
                ((protein, InputRole.TEMPLATE), (ef147, InputRole.STRUCTURAL),
                 (ef91, InputRole.STRUCTURAL), (kapa, InputRole.CATALYST)), (pdt_amp,),
                geometry=PCRGeometry.unavailable(
                    pdt_amp.state,
                    reason="EF147 and EF91 are published, but the barcoded template depends on "
                           "the unprinted 10x bead oligo; EF91's 3' end matches the EF79 FR3 handle "
                           "(Supplementary Table 1) and the binding order is inferred."))
    spri1 = graph.input("pdt_spri_065", "65 µl SPRIselect (first side)", MaterialKind.SELECTION)
    pdt_sup = graph.product("pdt_supernatant", "160 µl supernatant retained", MaterialKind.MOLECULAR,
                            state=pdt_amp.state)
    pdt_large = graph.product("pdt_large_bound", "large fragments bound to first-side beads",
                              MaterialKind.WASTE)
    graph.react("pdt_spri_upper", "Add 65 µl SPRIselect; keep 160 µl supernatant",
                StepKind.CLEANUP, ((pdt_amp, InputRole.TEMPLATE), (spri1, InputRole.SELECTION)),
                (pdt_sup, pdt_large))
    graph.discard(pdt_large, "bound to the first SPRIselect addition; supernatant is kept")
    spri2 = graph.input("pdt_spri_192", "192 µl SPRIselect (second side)", MaterialKind.SELECTION)
    pdt_clean = graph.product("pdt_clean", "size-selected PDT library in 40 µl EB",
                              MaterialKind.MOLECULAR, state=pdt_amp.state)
    pdt_small = graph.product("pdt_small_unbound", "unbound small products and primers",
                              MaterialKind.WASTE)
    graph.react("pdt_spri_lower", "Add 192 µl SPRIselect; 2 × 80% ethanol; elute 40 µl EB",
                StepKind.CLEANUP, ((pdt_sup, InputRole.TEMPLATE), (spri2, InputRole.SELECTION)),
                (pdt_clean, pdt_small))
    graph.discard(pdt_small, "left in the supernatant of the second SPRIselect addition")
    si_b = graph.input("pdt_si_pcr_b", "10x SI-PCR Primer B (7.5 µl)", MaterialKind.OLIGO, source=vendor)
    i7 = graph.input("pdt_i7", "10x i7 sample-index primer (2.5 µl)", MaterialKind.OLIGO, source=vendor)
    amp_mix = graph.input("pdt_amp_mix", "10x Amp-Mix (50 µl)", MaterialKind.CATALYST, source=vendor)
    pdt_indexed = graph.product("pdt_indexed", "indexed PDT library", MaterialKind.MOLECULAR,
                                state=PDT_LIBRARY)
    graph.react("pdt_index_pcr", "15 ng PDT library; indexing PCR 98 °C 45 s; 6 × (98 °C 20 s, 67 °C 30 s, 72 °C 20 s); 72 °C 1 min",
                StepKind.PCR,
                ((pdt_clean, InputRole.TEMPLATE), (si_b, InputRole.STRUCTURAL),
                 (i7, InputRole.STRUCTURAL), (amp_mix, InputRole.CATALYST)), (pdt_indexed,),
                geometry=PCRGeometry.unavailable(
                    PDT_LIBRARY,
                    reason="SI-PCR Primer B and the i7 sample-index primers are named 10x "
                           "reagents; their sequences are not printed in the paper."))
    spri3 = graph.input("pdt_spri_120", "120 µl SPRIselect", MaterialKind.SELECTION)
    pdt_lib = graph.product("pdt_library", "PHAGE-ATAC PDT sequencing library", MaterialKind.LIBRARY,
                            state=PDT_LIBRARY)
    pdt_waste = graph.product("pdt_index_waste", "unbound primers and reaction components",
                              MaterialKind.WASTE)
    graph.react("pdt_final_spri", "120 µl SPRIselect; elute 40 µl EB; Qubit and Bioanalyzer QC",
                StepKind.CLEANUP, ((pdt_indexed, InputRole.TEMPLATE), (spri3, InputRole.SELECTION)),
                (pdt_lib, pdt_waste))
    graph.discard(pdt_waste, "left in the SPRIselect supernatant")
    graph.sequence(pdt_lib, SequencingPlan(limitation=SEQUENCING_UNAVAILABLE,
                                           acquired_roles=ACQUIRED_ROLES,
                                           read_lengths=READ_LENGTHS))
    graph.declare_reviewed_scope(list(graph.materials) + list(graph.steps), source)
    return graph.require_publishable()
