from types import SimpleNamespace
from chemdraw import feature
from protocol_model import chain, role_library

from chemdraw import Segment
from readout_protocols import Construct, Definition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, _authority, _unavailable, molecule, role
DEFINITION = Definition("SDR-seq — single-cell DNA–RNA sequencing",
        "10.1038/s41592-025-02805-0", "Lindenhofer et al., Nature Methods (2025)",
        "In-situ barcoded poly(dT) RT is followed by Tapestri droplet multiplex PCR; distinct R2N and R2 overhangs split targeted gDNA and RNA into matched sequencing libraries.",
        "The gDNA branch spans targeted variants; the RNA branch reports transcript sequence plus cell barcode, RT sample barcode and UMI.",
        "01_sdr_seq.html")


def _sdr_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = _authority(d)
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input(
        "fixed_cells", "fixed permeabilized cells containing genomic DNA and poly(A) RNA",
        MaterialKind.BIOLOGICAL,
        state=molecule("fixed permeabilized cell", (role("genomic DNA"), role("poly(A) RNA"))),
        primary=True)
    rt_primer = graph.input("barcoded_rt_primer",
                            "V2 RT_primer (Supp. Table 2): 23-nt CS handle + 8-nt sample BC + N8 UMI + T16VN",
                            MaterialKind.OLIGO)
    reverse_transcriptase = graph.input("reverse_transcriptase", "Maxima H Minus reverse transcriptase",
                                       MaterialKind.CATALYST)
    tagged = graph.product(
        "tagged_cells", "cells retaining gDNA and tagged cDNA", MaterialKind.MOLECULAR,
        state=molecule("cell retaining gDNA and tagged cDNA", (role("genomic DNA"),
            *RT_HANDLE)))
    graph.react(
        "in_situ_rt", "Reverse-transcribe RNA in situ with barcoded poly(dT)",
        StepKind.REVERSE_TRANSCRIPTION,
        ((cells, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
         (reverse_transcriptase, InputRole.CATALYST)), (tagged,))

    droplets = graph.input("droplet_condition", "Tapestri first-droplet encapsulation (Mission Bio v2 device)",
                           MaterialKind.CONDITION)
    occupied = graph.product("occupied_droplets", "cell-containing first droplets",
                             MaterialKind.MOLECULAR, state=tagged.state)
    empty = graph.product("empty_droplets", "empty or unusable droplets", MaterialKind.WASTE)
    graph.react("partition", "Encapsulate single cells in first droplets", StepKind.PARTITION,
                ((tagged, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (occupied, empty))
    graph.discard(empty, "does not contain a usable cell template")
    proteinase = graph.input("proteinase_k", "proteinase K", MaterialKind.CATALYST)
    lysed = graph.product("lysed_droplets", "lysed, proteinase K-treated cells in first droplets",
                          MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("lyse_digest", "Lyse cells and digest protein with proteinase K in first droplets",
                StepKind.DEPROTEINIZATION,
                ((occupied, InputRole.TEMPLATE), (proteinase, InputRole.CATALYST)), (lysed,))

    target_forward = graph.input("target_forward_pool",
                                 "gDNA forward pool: Droplet_gDNA_primer_fwd GTACTCGCAGTAGTC-gDNA_specific",
                                 MaterialKind.OLIGO)
    rna_forward = graph.input("rna_forward", "Droplet_amplification_fwd universal RNA forward primer on the RT handle (V2)",
                              MaterialKind.OLIGO,
                              state=(Segment("Droplet_amplification_fwd", RNA_FORWARD),))
    branch_reverse = graph.input(
        "branch_reverse_pool", "reverse pool added in the first droplet: R2N-gDNA_specific and R2-cDNA_specific primers",
        MaterialKind.OLIGO)
    bead_oligos = graph.input("barcode_bead_oligos", "Mission Bio barcoding-bead cell-BC oligos with CS overhang (bases not published)",
                              MaterialKind.OLIGO)
    pcr_polymerase = graph.input("droplet_pcr_polymerase", "droplet PCR polymerase",
                                 MaterialKind.CATALYST)
    dna_amplicons = graph.product(
        "dna_amplicons", "cell-barcoded targeted gDNA amplicons", MaterialKind.MOLECULAR,
        state=molecule("gDNA target amplicons", (*CELL_BARCODE,
            role("gDNA target"), role("R2N overhang"))))
    rna_amplicons = graph.product(
        "rna_amplicons", "cell-barcoded targeted cDNA amplicons", MaterialKind.MOLECULAR,
        state=molecule("cDNA target amplicons", (*CELL_BARCODE, role("GACACGTC (3′ part of the 23-nt RNA handle)"), *RT_HANDLE[1:],
            role("R2 overhang"))))
    graph.react(
        "droplet_multiplex_pcr", "Second droplet: add forward primers, PCR reagents and barcoding bead; multiplex PCR",
        StepKind.PCR,
        ((lysed, InputRole.TEMPLATE), (target_forward, InputRole.STRUCTURAL),
         (rna_forward, InputRole.STRUCTURAL),
         (branch_reverse, InputRole.STRUCTURAL), (bead_oligos, InputRole.STRUCTURAL),
         (pcr_polymerase, InputRole.CATALYST)), (dna_amplicons, rna_amplicons),
        geometry=_unavailable(
            final, "The panel-specific target primers and commercial bead oligo are not one fixed public amplicon."))
    breaker = graph.input("emulsion_breaking", "emulsion breaking and amplicon recovery (Mission Bio protocol)",
                          MaterialKind.CONDITION)
    dna_recovered = graph.product("dna_recovered", "recovered gDNA amplicons", MaterialKind.MOLECULAR,
                                  state=dna_amplicons.state)
    rna_recovered = graph.product("rna_recovered", "recovered cDNA amplicons", MaterialKind.MOLECULAR,
                                  state=rna_amplicons.state)
    emulsion_waste = graph.product("emulsion_waste", "oil phase and droplet carry-over", MaterialKind.WASTE)
    graph.react("break_emulsion", "Break emulsions and recover pooled amplicons", StepKind.CLEANUP,
                ((dna_amplicons, InputRole.TEMPLATE), (rna_amplicons, InputRole.TEMPLATE),
                 (breaker, InputRole.CONDITION)), (dna_recovered, rna_recovered, emulsion_waste))
    graph.discard(emulsion_waste, "removed when emulsions are broken; recovery details are in the Mission Bio protocol")

    index_polymerase = graph.input("index_polymerase", "library indexing PCR polymerase",
                                   MaterialKind.CATALYST)
    dna_forward = graph.input("dna_index_forward", "Nextera_R1N_S50x i5 primer (gDNA library)",
                              MaterialKind.OLIGO)
    dna_reverse = graph.input("dna_index_reverse", "Nextera_R2N_N70x i7 primer (gDNA library)",
                              MaterialKind.OLIGO)
    dna_library = graph.product("dna_library", "SDR-seq targeted DNA library",
                                MaterialKind.LIBRARY, state=final)
    graph.react(
        "dna_library_pcr", "Build the R2N-marked DNA library", StepKind.PCR,
        ((dna_recovered, InputRole.TEMPLATE), (dna_forward, InputRole.STRUCTURAL),
         (dna_reverse, InputRole.STRUCTURAL), (index_polymerase, InputRole.CATALYST)),
        (dna_library,), geometry=_unavailable(
            final, "The public role-level description does not establish a universal R2N indexing pair."))
    rna_index_forward = graph.input("rna_index_forward", "Nextera_R1N_S50x i5 primer (RNA library)",
                              MaterialKind.OLIGO)
    rna_reverse = graph.input("rna_index_reverse", "Illumina_R2_N70x i7 primer (RNA library)",
                              MaterialKind.OLIGO)
    rna_library = graph.product("rna_library", "SDR-seq targeted RNA library",
                                MaterialKind.LIBRARY, state=final)
    graph.react(
        "rna_library_pcr", "Build the R2-marked RNA library", StepKind.PCR,
        ((rna_recovered, InputRole.TEMPLATE), (rna_index_forward, InputRole.STRUCTURAL),
         (rna_reverse, InputRole.STRUCTURAL), (index_polymerase, InputRole.CATALYST)),
        (rna_library,), geometry=_unavailable(
            final, "The public role-level description does not establish a universal R2 indexing pair."))
    graph.sequence(dna_library, SequencingPlan(
        limitation="DNA branch: the genomic-target read is acquired together with its cell barcode; exact panel-dependent cycles are not universal."))
    graph.sequence(rna_library, SequencingPlan(
        limitation="RNA branch: transcript, cell barcode, RT sample barcode and UMI are acquired; exact panel-dependent cycles are not universal."))
    graph.declare_reviewed_scope((
        "fixed_cells", "barcoded_rt_primer", "reverse_transcriptase", "in_situ_rt",
        "droplet_condition", "partition", "proteinase_k", "lyse_digest", "rna_forward",
        "emulsion_breaking", "break_emulsion", "target_forward_pool", "branch_reverse_pool",
        "barcode_bead_oligos", "droplet_pcr_polymerase", "droplet_multiplex_pcr",
        "index_polymerase", "dna_index_forward", "dna_index_reverse", "dna_library_pcr",
        "rna_index_forward", "rna_index_reverse", "rna_library_pcr", "dna_library",
        "rna_library",
    ), source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    key = 'sdr-seq'
    return _sdr_graph(d, final)

# Methods "SDR-seq data analysis": cell BC1, variable-length 14-17-bp linker, cell BC2,
# 15-bp constant linker, then (RNA only) sample BC and UMI.  Whether the 15-bp constant
# linker is the 15-nt capture sequence GTACTCGCAGTAGTC is inferred, not stated.
CELL_BARCODE = (
    role("cell BC1", tag="cbc", molecular_feature=feature(
        "cell_barcode_1", "cell_barcode", "unknown", group="cell_barcode", part="1")),
    role("variable 14-17-bp linker"),
    role("cell BC2", tag="cbc", molecular_feature=feature(
        "cell_barcode_2", "cell_barcode", "unknown", group="cell_barcode", part="2")),
    role("15-bp constant linker (capture-sequence junction)", inferred=True),
)
# Supplementary Table 2 (V2): Droplet_amplification_fwd and the 5' handle of every RT_primer.
RNA_FORWARD = "GTACTCGCAGTAGTCGACACGTC"
UMI = role("UMI", tag="umi", length_bp=8,
           molecular_feature=feature("umi", "umi", "unknown"))
RT_SAMPLE_BARCODE = role("RT sample barcode", tag="cbc", length_bp=8,
                         molecular_feature=feature("rt_sample_barcode", "sample_index", "unknown"))
RT_HANDLE = (role("23-nt CS/amplification handle"), RT_SAMPLE_BARCODE, UMI,
             role("T16VN"), role("cDNA target"))

def _sdr():
    return chain(molecule("fixed permeabilized cell", (role("genomic DNA"),
                                                        role("poly(A) RNA"))), (
        ("Reverse-transcribe RNA in situ with barcoded poly(dT)",
         molecule("cell retaining gDNA and tagged cDNA", (role("genomic DNA"), *RT_HANDLE)),
         "The RT primer installs a UMI and sample barcode before droplet processing."),
        ("Multiplex-amplify targeted gDNA and cDNA in Tapestri droplets",
         molecule("target amplicons", (*CELL_BARCODE,
                                        role("gDNA or cDNA target"),
                                        role("R2N or R2 overhang"))),
         "Capture-sequence complementarity copies the bead cell barcode onto both molecular branches."),
        ("Use distinct reverse overhangs to build separate DNA and RNA libraries",
         molecule("matched SDR-seq libraries", (role("gDNA: Nextera-R2N library"),
                                                 *CELL_BARCODE,
                                                 role("RNA: TruSeq-R2 library"),
                                                 RT_SAMPLE_BARCODE, UMI)),
         "The DNA read spans the variant target; the RNA read retains cell, sample and UMI identity."),
    ))

CAVEAT = ("The page models the V2 design (separate gDNA R2N and RNA R2 libraries); V1 (proof-of-concept) "
          "used R2N on both branches and built libraries with the Mission Bio user guide. "
          "Published discrepancy: Supplementary Table 3 lists the V1 proof-of-concept RT primers as "
          "RT_primer_1_8 to RT_primer_4_8 (the V2 names in Supplementary Table 2, with the GACACGTC spacer), "
          "whereas Supplementary Table 1 lists V1 RT_primer_1 to RT_primer_4 without that spacer. "
          "The model does not depend on the V1 primer; which sequence the POP run used remains unresolved. "
          "The bead oligo bases, whether the 15-bp constant linker of the read structure is the 15-nt capture sequence, and read cycles are not given in the inspected sources.")
d = DEFINITION
FINAL_LIBRARY = role_library("SDR-seq matched library roles", (*CELL_BARCODE, RT_SAMPLE_BARCODE, UMI, role("targeted gDNA or cDNA insert")))
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES=d.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary, CAVEAT=CAVEAT, FINAL_LIBRARY=FINAL_LIBRARY, SEQ_PRIMERS=(),
    FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint, READ_LENGTHS=None,
    SEQUENCING_ENDING=d.endpoint,
    SEQUENCING_UNAVAILABLE=("Supplementary Tables 1 and 2 publish experiment-specific RT and library primers. The V1 proof-of-concept used a common R2N arm; later V2 experiments separate genomic and RNA libraries. The Tapestri bead oligo panel still requires vendor-version verification."),
    workflow=_sdr, oligos=lambda: (
        "<p>Published Supplementary Table 1, V1 RT_primer_1 (5′→3′): GTACTCGCAGTAGTC<cbc>TCGCCTTA</cbc><umi>NNNNNNNN</umi>TTTTTTTTTTTTTTTTVN</p>",
        "<p>Published Supplementary Table 2, V2 RT_primer_1_8 (5′→3′): GTACTCGCAGTAGTCGACACGTC<cbc>TCGCCTTA</cbc><umi>NNNNNNNN</umi>TTTTTTTTTTTTTTTTVN; Droplet_amplification_fwd: " + RNA_FORWARD + "</p>",
        "<p>Published Supplementary Table 2, V2 example genomic i7 Nextera_R2N_N701 (5′→3′): CAAGCAGAAGACGGCATACGAGAT<cbc>TCGCCTTA</cbc>GTCTCGTGGGCTCGG</p>",
        "<p>Published Supplementary Table 2, V2 example RNA i7 Illumina_R2_N701 (5′→3′): CAAGCAGAAGACGGCATACGAGAT<cbc>TCGCCTTA</cbc>GTGACTGGAGTTCAGACGTGT</p>",
    ),
)
SPEC.protocol_graph = protocol_graph
globals().update(vars(SPEC))
