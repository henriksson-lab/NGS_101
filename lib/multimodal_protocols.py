"""Shared molecular models for single-cell multimodal identifier transfer.

The central invariant is identifier provenance: every feature says whether it identifies
a cell, molecule, sample/pool, antibody or genomic target.  Branching workflows keep the
same cell identity while allowing physically distinct sequencing libraries.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

from chemdraw import feature
from circular import CircularizationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import chain, molecule, role, role_library


def ident(name, identifier, kind, *, length=None, tag="cbc"):
    return role(name, tag=tag, length_bp=length,
                molecular_feature=feature(identifier, kind, "unknown"))


CELL = ident("cell barcode", "cell_barcode", "cell_barcode")
UMI = ident("UMI", "umi", "umi", tag="umi")
SAMPLE = ident("sample index", "sample_index", "sample_index", length=8)
POOL = ident("staining-pool barcode", "pool_barcode", "sample_index", length=10)
ANTIBODY = ident("antibody barcode", "antibody_barcode", "feature_barcode")
DNA_TYPE = ident("DNA-branch barcode", "dna_branch", "inline_barcode", length=6)
RNA_TYPE = ident("RNA-branch barcode", "rna_branch", "inline_barcode", length=6)
DNA_UMI = ident("DNA UMI", "dna_umi", "umi", length=8, tag="umi")
RNA_UMI = ident("RNA UMI", "rna_umi", "umi", length=8, tag="umi")
GUIDE = ident("12-nt internal guide barcode", "guide_ibar", "guide_barcode", length=12)

SCONE_ONE_TN5 = "GTCTCGTGGGCTCGGTCATGNNNNNNNNAGATGTGTATAAGAGACAG"
SCONE_P5_TN5 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
SCONE_MOSAIC = "CTGTCTCTTATACACATCT"
SCONE_RT_POLYT = "GTCTCGTGGGCTCGGAACTTNNNNNNNNTTTTTTTTTTTTTTTTTTTTVN"
SCONE_RT_GGG = "GTCTCGTGGGCTCGGAACTTNNNNNNNNGGGHN"
SCONE_RT_TTT = "GTCTCGTGGGCTCGGAACTTNNNNNNNNTTTVN"
SCONE_3P_ADAPTER = SCONE_RT_GGG
SCONE_AMP_PRIMER = "GATGTGTGGAGGTCTCGTGGGCTCGG"
SCONE_ILLU_R1 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
SCONE_ILLU_R2 = "CCGAGCCCACGAGACCTGTCTCTTATACACATCT"
SCONE_ILLU_I7 = "AGATGTGTATAAGAGACAGGTCTCGTGGGCTCGG"

PHOSPHO_BOA = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNVTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
PHOSPHO_BOB = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTTGCTAGGACCGGCCTTAAAGC"
PHOSPHO_P5 = "AATGATACGGCGACCACCGA"
PHOSPHO_P7_TSA = "CAAGCAGAAGACGGCATACGAGATXXXXXXXXGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
PHOSPHO_P7_TSB = "CAAGCAGAAGACGGCATACGAGATXXXXXXXXGTGACTGGAGTTCAGACGTGTGC"
PHOSPHO_HTO_ADDITIVE = "CCTTGGCACCCGAGAATTCC"
PHOSPHO_ADT_ADDITIVE = "GTGACTGGAGTTCAGACGTGTGCTC"
PHOSPHO_MULTIOME_P5 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC"


def _scone_oligos():
    """Exact oligos printed in scONE-seq Supplementary Table S2.

    Index windows remain symbolic because the table defines primer families, while the
    three custom Illumina sequencing primers are invariant concrete sequences.
    """
    return (
        "<b>One-Tn5 adaptor</b> — 5′-[biotin]" + SCONE_ONE_TN5 + "-3′",
        "<b>P5-Tn5 adaptor</b> — 5′-" + SCONE_P5_TN5 + "-3′",
        "<b>Mosaic</b> — 5′-[phos]" + SCONE_MOSAIC + "-3′",
        "<b>RT-polyT</b> — 5′-[biotin]" + SCONE_RT_POLYT + "-3′",
        "<b>RT-GGG</b> — 5′-[biotin]" + SCONE_RT_GGG + "-3′",
        "<b>RT-TTT</b> — 5′-[biotin]" + SCONE_RT_TTT + "-3′",
        "<b>3′ adaptor</b> — 5′-" + SCONE_3P_ADAPTER + "-3′",
        "<b>Amp-primer</b> — 5′-" + SCONE_AMP_PRIMER + "-3′",
        "<b>P7 index-primer family</b> — 5′-CAAGCAGAAGACGGCATACGAGAT[i7]AGATGTGTATAAGAGACAGGTCTCGTGGGCTCGG-3′",
        "<b>P5 index-primer family</b> — 5′-AATGATACGGCGACCACCGAGATCTACAC[i5]TCGTCGGCAGCGTC-3′",
        "<b>ONE-Illu Read 1</b> — 5′-" + SCONE_ILLU_R1 + "-3′",
        "<b>ONE-Illu Read 2</b> — 5′-" + SCONE_ILLU_R2 + "-3′",
        "<b>ONE-Illu I7 Index Primer</b> — 5′-" + SCONE_ILLU_I7 + "-3′",
    )


def _phospho_oligos():
    """Published Phospho-seq tag, bridge and amplification oligos.

    The bridge sequences are inherited explicitly from ASAP-seq, the method cited for
    this molecular conversion. Antibody-specific 15-base barcodes and sample-specific
    i7 choices remain variable windows rather than invented representative bases.
    """
    return (
        "<b>TotalSeq-A-style hashing tag</b> — "
        "5′-/5AmMC12/CCTTGGCACCCGAGAATTCCA[15-nt antibody barcode]B"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA*A*A-3′",
        "<b>TotalSeq-B-style antibody tag</b> — "
        "5′-/5AmMC12/GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT[10-nt UMI]"
        "[15-nt antibody barcode][9-nt UMI]GCTTTAAGGCCGGTCCTAGC*A*A-3′",
        "<b>30-nt blocking oligo</b> — 5′-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/3ddC/-3′",
        "<b>Bridge oligo A (BOA)</b> — 5′-" + PHOSPHO_BOA + "/3InvdT/-3′",
        "<b>Bridge oligo B (BOB)</b> — 5′-" + PHOSPHO_BOB + "/3InvdT/-3′",
        "<b>P5 tag-library primer</b> — 5′-" + PHOSPHO_P5 + "-3′",
        "<b>TotalSeq-A P7 index-primer family</b> — 5′-" + PHOSPHO_P7_TSA + "-3′",
        "<b>TotalSeq-B P7 index-primer family</b> — 5′-" + PHOSPHO_P7_TSB + "-3′",
        "<b>Multiome HTO additive</b> — 5′-" + PHOSPHO_HTO_ADDITIVE + "*C*C-3′",
        "<b>Multiome ADT additive</b> — 5′-GTGACTGGAGTTCAGACGTGTGC*T*C-3′",
        "<b>Multiome SI-PCR P5</b> — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGC*T*C-3′",
    )


def _sccc():
    return chain(molecule("cell carrying self-reporting transposons", (
        role("TF-directed piggyBac transposase"), role("integrated SRT"))), (
        ("Transcribe through the SRT terminal repeat into flanking genomic DNA",
         molecule("self-reporting RNA", (role("reporter transcript"),
                                          role("transposon terminal repeat"),
                                          role("genomic insertion flank"), role("poly(A)"))),
         "Removing the transposon's polyadenylation signal makes the insertion junction part of the RNA 3′ UTR."),
        ("Run droplet RT and split the barcoded first-strand product",
         molecule("SRT cDNA", (CELL, UMI, role("insertion-flank cDNA"),
                                role("terminal-repeat cDNA"))),
         "One aliquot becomes the expression library; the other is amplified specifically for SRT transcripts."),
        ("Circularize the SRT amplicon, remove linear DNA and shear",
         molecule("junction-proximal fragment", (CELL, UMI,
                                                  role("genomic insertion flank"),
                                                  role("terminal repeat"))),
         "Circularization moves the 3′ cell barcode and UMI next to the transposon–genome junction before short-read library construction."),
        ("Capture the junction-bearing branch and add custom Illumina adapters",
         molecule("scCC library", (role("custom P5 end"), CELL, UMI,
                                    role("insertion flank"), role("terminal repeat"),
                                    role("custom P7 end"))),
         "The shared droplet barcode joins insertion sites to the separately sequenced transcriptome."),
    ))


def _cotech():
    return chain(molecule("antibody-bound permeabilized cells", (
        role("target chromatin epitope"), role("poly(A) RNA"))), (
        ("Bind uniquely barcoded protein-A–Tn5 and tagment target-proximal DNA",
         molecule("tagged chromatin branch", (role("T5 well barcode"),
                                               role("epitope-proximal DNA"),
                                               role("T7 well barcode"))),
         "The antibody localizes pA–Tn5; the first well identifies the chromatin branch."),
        ("Reverse-transcribe mRNA with a matching well-barcoded poly(dT) primer",
         molecule("paired intracellular branches", (role("barcoded chromatin fragments"),
                                                      role("RT well barcode"),
                                                      role("transcript cDNA"))),
         "Chromatin and RNA receive a shared first-round well identity before pooling."),
        ("Pool, redistribute and preamplify both branches with a second well index",
         molecule("two-index cell identity", (ident("first-round well barcode", "cotech_round1", "cell_barcode"),
                                               ident("second-round well barcode", "cotech_round2", "cell_barcode"),
                                               role("DNA and cDNA products"))),
         "The ordered pair of wells identifies a cell after combinatorial indexing."),
        ("Split the product into chromatin and transcriptome libraries",
         molecule("matched CoTECH libraries", (role("TruSeq-like chromatin library"),
                                                role("shared cell identity"),
                                                role("Nextera-like RNA library"))),
         "The two physical libraries remain linked by the same combinatorial barcode."),
    ))


def _gotcha():
    return chain(molecule("single-cell ATAC reaction", (role("tagmented accessible DNA"),
                                                         role("target genomic locus"))), (
        ("Spike a locus-specific primer into droplet barcoding extension",
         molecule("barcoded ATAC and target molecules", (CELL,
             role("ordinary accessible-chromatin fragments"), role("target-locus amplicon"))),
         "The target primer copies the same droplet cell barcode used by the scATAC library onto the genomic locus."),
        ("Split off an aliquot and perform biotinylated hemi-nested target PCR",
         molecule("enriched genotype amplicon", (CELL, role("variant-spanning genomic insert"),
                                                  role("biotinylated nested end"))),
         "Nested amplification enriches the selected locus while ordinary ATAC fragments are only linearly carried through."),
        ("Streptavidin-capture and index on-bead",
         molecule("GoT-ChA library", (role("P5 / Read 1 end"), CELL,
                                       role("variant-spanning insert"), SAMPLE,
                                       role("Read 2 / P7 end"))),
         "The genotype library can be pooled with the parent scATAC library and joined computationally by cell barcode."),
    ))


def _scone():
    return chain(molecule("lysed single cell or nucleus", (role("genomic DNA"),
                                                            role("total RNA"))), (
        ("Tagment genomic DNA with a custom molecularly barcoded adapter",
         molecule("tagged DNA", (DNA_TYPE, DNA_UMI, role("tagmented genomic fragment"))),
         "The six-base branch barcode distinguishes DNA; the eight-base UMI identifies an original tagged molecule."),
        ("Internally prime total RNA with a separately barcoded RT primer",
         molecule("first-strand cDNA", (RNA_TYPE, RNA_UMI, role("full-transcript cDNA"))),
         "Internal priming captures polyadenylated and non-polyadenylated RNA."),
        ("Poly(C)-tail cDNA and add the opposite amplification adapter",
         molecule("co-amplifiable DNA and cDNA", (role("common amplification end"),
                                                   role("DNA or RNA insert"),
                                                   role("branch barcode and UMI"))),
         "DNA-specific and RNA-specific identifiers permit both products to share one-tube amplification."),
        ("Amplify together and complete the indexed sequencing library",
         molecule("scONE-seq library", (role("P5 end"), SAMPLE,
                                         role("branch barcode"), role("UMI"),
                                         role("DNA or cDNA insert"), role("P7 end"))),
         "Reads are separated by the molecular branch barcode after sequencing."),
    ))


def _scito():
    return chain(molecule("universally oligo-conjugated antibody", (
        role("antibody"), role("20-nt antibody handle"))), (
        ("Hybridize a pool-specific splint carrying antibody and pool identity",
         molecule("SCITO antibody tag", (role("antibody handle complement"),
                                          role("TruSeq Read 2 handle"), ANTIBODY,
                                          POOL, role("10x feature-capture sequence"))),
         "One antibody conjugation can be reused; the splint supplies both target and staining-pool identity."),
        ("Stain split cell pools, combine them and overload droplets",
         molecule("tagged cell", (ANTIBODY, POOL, role("surface-protein counts"))),
         "Pool identity makes many otherwise colliding cells distinguishable."),
        ("Capture tags on gel beads and copy droplet barcode plus UMI",
         molecule("SCITO-seq library", (role("P5 end"), CELL, UMI,
                                         ANTIBODY, POOL, role("Read 2 / P7 end"))),
         "Antibody-plus-pool barcode and droplet barcode form the combinatorial cell identity."),
    ))


def _phospho():
    return chain(molecule("fixed permeabilized hashed cell", (
        role("accessible chromatin"), role("intracellular phosphoproteins"))), (
        ("Stain intracellular targets with DNA-tagged antibodies",
         molecule("tagged cell", (ANTIBODY, role("phosphoprotein epitope"),
                                   ident("hash tag", "sample_hash", "sample_index"))),
         "Click-conjugated antibodies, including TotalSeq-A and TotalSeq-B formats, carry distinct feature barcodes."),
        ("Tagment chromatin and add bridge oligos before 10x ATAC GEM barcoding",
         molecule("capturable branches", (role("tagmented accessible DNA"),
                                           role("bridge-bound antibody tag"))),
         "Bridge oligos make the antibody tags compatible with the gel-bead capture sequence."),
        ("Separate ATAC DNA from the tag-containing supernatant after SPRI cleanup",
         molecule("matched Phospho-seq libraries", (CELL, role("ATAC insert"), UMI,
                                                     ANTIBODY, role("phosphoprotein tag"))),
         "ATAC fragments remain on beads; hashing and antibody tags are amplified from the saved supernatant."),
    ))


def _cropseq_multi():
    return chain(molecule("CROPseq-multi 3-prime LTR cassette", (
        role("U6 promoter, reverse orientation"), role("tRNA"), role("sgRNA 1"),
        role("tRNA"), role("sgRNA 2"), GUIDE)), (
        ("Duplicate the engineered 3′ LTR during lentiviral integration",
         molecule("integrated perturbation cassette", (role("5-prime LTR copy"),
             role("selection-gene transcript"), role("3-prime LTR copy"))),
         "Keeping the Pol III cassette on the minus strand prevents tRNA processing from cleaving lentiviral genomic RNA or the Pol II selection transcript."),
        ("Transcribe and tRNA-process the multiplexed Pol III guide array",
         molecule("two mature perturbation guides", (role("sgRNA 1 with iBAR"),
                                                       role("sgRNA 2 with iBAR"))),
         "RNases P and Z release separate guides; each 12-base iBAR lies only 19 bp from its linked spacer."),
        ("Capture perturbation identity by mRNA, direct guide capture or T7-IVT",
         molecule("CROPseq-multi readout", (CELL, UMI, GUIDE,
                                             role("linked guide spacer"))),
         "The same iBAR is transcribed antisense in polyadenylated mRNA and within the Pol III sgRNA; CSMv2 additionally includes a T7 promoter."),
    ))


def _chair():
    return chain(molecule("crosslinked nucleus", (role("chromatin contacts"),
                                                   role("accessible chromatin"),
                                                   role("poly(A) RNA"))), (
        ("Digest chromatin in situ with a restriction enzyme",
         molecule("restricted crosslinked chromatin", (role("restriction fragment A"),
                                                        role("crosslinked restriction fragment B"))),
         "Restriction cleavage exposes ends while crosslinks preserve spatial neighbours."),
        ("Proximity-ligate crosslinked restriction ends",
         molecule("3D contact junction", (role("restriction fragment A"),
                                           role("restriction fragment B")),
                  junctions=(("restriction fragment A", "restriction fragment B", "proximity ligation"),)),
         "The ligation junction records a chromatin contact."),
        ("Tagment accessible chromatin while retaining contact products",
         molecule("ChAIR DNA precursor", (role("Tn5 end"), role("accessible anchor"),
                                           role("contact partner"), role("Tn5 end"))),
         "The DNA branch yields ordinary accessibility fragments and paired-end contact tags anchored at accessible chromatin."),
        ("Bar-code DNA and RNA together in 10x Multiome droplets",
         molecule("three matched readouts", (CELL, UMI, role("ChAIR-RNA"),
                                              role("ChAIR-ATAC"), role("ChAIR-PET"))),
         "One droplet cell barcode links gene expression, accessibility and proximity-ligation junctions."),
    ))


@dataclass(frozen=True)
class Definition:
    title: str; doi: str; citation: str; summary: str; endpoint: str; notes: str


DEFS = {
    "sccc": Definition("scCC — single-cell calling cards", "10.1016/j.cell.2020.06.037",
        "Moudgil et al., Cell (2020)", "Self-reporting transposons encode their genomic insertion flank in RNA; circularization brings droplet identifiers next to that junction so TF-binding history can be joined to gene expression.",
        "The calling-card read resolves the insertion flank and retains the same cell barcode and UMI as the parent 10x cDNA.", "01_sccc.html"),
    "cotech": Definition("CoTECH — chromatin occupancy and transcriptome", "10.1038/s41592-021-01129-z",
        "Xiong et al., Nature Methods (2021)", "Barcoded antibody-targeted pA–Tn5 and barcoded poly(dT) RT mark chromatin and RNA in the same first well; redistribution supplies a second combinatorial cell index.",
        "Separate chromatin and RNA libraries share the ordered pair of well barcodes.", "01_cotech.html"),
    "got-cha": Definition("GoT-ChA — genotyping of targeted loci with chromatin accessibility", "10.1101/2022.05.11.491515",
        "Namikawa et al., bioRxiv (2022)", "A locus-specific primer participates in single-cell ATAC droplet barcoding, after which hemi-nested biotin capture enriches the genomic variant without losing the ATAC cell barcode.",
        "The genotype read spans the chosen variant; the cell barcode links it to the parent accessibility profile.", "01_got_cha.html"),
    "scone-seq": Definition("scONE-seq — one-tube single-cell DNA and RNA sequencing", "10.1126/sciadv.abp8901",
        "Yu et al., Science Advances (2023)", "Custom Tn5 adapters label genomic DNA while internally primed RT labels total RNA with a distinct branch barcode and UMI; both branches are amplified together.",
        "A six-base branch barcode separates DNA from RNA and an eight-base UMI identifies molecules within each branch.", "01_scone_seq.html"),
    "scito-seq": Definition("SCITO-seq — single-cell combinatorial indexed cytometry sequencing", "10.1038/s41592-021-01222-3",
        "Hwang et al., Nature Methods (2021)", "Hybridized splints add antibody and staining-pool barcodes before droplet feature capture adds a second cell barcode and UMI.",
        "Protein identity and staining-pool identity are read with the droplet barcode and UMI.", "01_scito_seq.html"),
    "phospho-seq": Definition("Phospho-seq — intracellular protein dynamics with chromatin accessibility", "10.1101/2023.03.27.534442",
        "Chung et al., bioRxiv (2023; published 2025)", "DNA-tagged intracellular antibodies are bridged into a 10x scATAC capture reaction; tag and ATAC branches are separated after GEM barcoding while retaining cell identity.",
        "The ATAC library and antibody/hash tag libraries share a droplet cell barcode; antibody barcodes identify total or phosphorylated targets.", "01_phospho_seq.html"),
    "cropseq-multi": Definition("CROPseq-multi — multiplexed perturbation cassette", "10.1101/2024.03.17.585235",
        "Walton et al., bioRxiv (2024)", "A reverse-oriented tRNA–sgRNA array in the lentiviral 3′ LTR expresses multiple guides while short internal barcodes support mRNA, direct-guide and optical readouts.",
        "A 12-base internal guide barcode remains linked to its spacer and can be read from poly(A) mRNA, Pol III guide RNA or the CSMv2 T7-IVT product.", "01_cropseq_multi.html"),
    "chair": Definition("ChAIR — chromatin accessibility, interaction and RNA", "10.1038/s41592-025-02658-7",
        "Chai et al., Nature Methods (2025)", "In-situ restriction digestion and proximity ligation record chromatin contacts, Tn5 marks accessible anchors, and 10x Multiome barcodes the DNA and RNA branches together.",
        "ChAIR-RNA, ChAIR-ATAC and ChAIR-PET share a droplet cell barcode; paired DNA ends distinguish contact junctions from accessibility fragments.", "01_chair.html"),
}


def _finish_pcr(graph, prefix, template, final, endpoint, reviewed, *, label):
    forward = graph.input(prefix + "_forward", label + " forward PCR primer",
                          MaterialKind.OLIGO)
    reverse = graph.input(prefix + "_reverse", label + " reverse PCR primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input(prefix + "_polymerase", label + " PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr", "Amplify and index the " + label + " library",
                StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="The source defines a primer family or kit-level endpoint rather than one universal amplicon."))
    graph.sequence(library, SequencingPlan(limitation=endpoint))
    reviewed.extend((prefix + "_forward", prefix + "_reverse",
                     prefix + "_polymerase", prefix + "_pcr",
                     prefix + "_library"))
    return library


def _finish_rt_pcr(graph, prefix, template, final, endpoint, reviewed, *, label):
    primer = graph.input(prefix + "_rt_primer", label + " reverse-transcription primer",
                         MaterialKind.OLIGO)
    enzyme = graph.input(prefix + "_rt_enzyme", label + " reverse transcriptase",
                         MaterialKind.CATALYST)
    cdna = graph.product(prefix + "_cdna", label + " cDNA", MaterialKind.MOLECULAR,
                         state=molecule(label + " cDNA", (role("cDNA insert"),)))
    graph.react(prefix + "_rt", "Reverse-transcribe the " + label + " branch",
                StepKind.REVERSE_TRANSCRIPTION,
                ((template, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
                 (enzyme, InputRole.CATALYST)), (cdna,))
    reviewed.extend((prefix + "_rt_primer", prefix + "_rt_enzyme", prefix + "_rt"))
    return _finish_pcr(graph, prefix, cdna, final, endpoint, reviewed, label=label)


def _sccc_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cell = graph.input("cell", "cell carrying an integrated self-reporting transposon",
                       MaterialKind.BIOLOGICAL, state=_sccc().initial, primary=True)
    rna_polymerase = graph.input("rna_polymerase", "cellular RNA polymerase",
                                 MaterialKind.CATALYST)
    transcript = graph.product("srt_transcript", "self-reporting transcript containing the insertion flank",
                               MaterialKind.MOLECULAR,
                               state=molecule("SRT transcript", (role("reporter transcript"),
                                   role("terminal repeat"), role("genomic insertion flank"), role("poly(A)"))))
    graph.react("transcribe_srt", "Transcribe through the SRT into flanking genomic DNA",
                StepKind.TRANSCRIPTION,
                ((cell, InputRole.TEMPLATE), (rna_polymerase, InputRole.CATALYST)),
                (transcript,))
    bead_primer = graph.input("bead_primer", "droplet oligo(dT) cell-barcode/UMI primer",
                              MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product("barcoded_cdna", "cell-barcoded SRT cDNA",
                         MaterialKind.MOLECULAR,
                         state=molecule("barcoded SRT cDNA", (CELL, UMI,
                             role("insertion-flank cDNA"), role("terminal-repeat cDNA"))))
    graph.react("reverse_transcribe", "Reverse-transcribe and barcode the SRT transcript",
                StepKind.REVERSE_TRANSCRIPTION,
                ((transcript, InputRole.TEMPLATE), (bead_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    expression = graph.product("expression_branch", "expression-library aliquot",
                               MaterialKind.MOLECULAR, state=cdna.state)
    junction = graph.product("junction_branch", "SRT-junction aliquot",
                             MaterialKind.MOLECULAR, state=cdna.state)
    graph.react("split_branches", "Split expression and SRT-junction branches",
                StepKind.SPLIT, ((cdna, InputRole.TEMPLATE),), (expression, junction),
                changes_molecule=False)
    reviewed = ["cell", "rna_polymerase", "transcribe_srt", "bead_primer",
                "reverse_transcriptase", "reverse_transcribe", "split_branches"]
    expression_final = role_library("scCC expression-library roles",
                                    (CELL, UMI, role("transcript cDNA")), index=SAMPLE)
    _finish_pcr(graph, "expression", expression, expression_final,
                "Expression branch: droplet transcriptome reads linked by cell barcode.",
                reviewed, label="scCC expression")
    circle_ligase = graph.input("circle_ligase", "SRT amplicon circularization ligase",
                                MaterialKind.CATALYST)
    circle = graph.product("circle", "circular SRT amplicon", MaterialKind.MOLECULAR,
                           state=junction.state)
    graph.react("circularize", "Circularize the SRT amplicon", StepKind.CIRCULARIZATION,
                ((junction, InputRole.TEMPLATE), (circle_ligase, InputRole.CATALYST)), (circle,),
                geometry=CircularizationGeometry.known(
                    "double-strand intramolecular ligation",
                    statement="The SRT branch circularizes its amplicon before exonuclease selection."))
    exonuclease = graph.input("exonuclease", "linear-DNA exonuclease",
                              MaterialKind.CATALYST)
    cleaned = graph.product("clean_circle", "exonuclease-resistant SRT circles",
                            MaterialKind.MOLECULAR, state=circle.state)
    linear_waste = graph.product("linear_waste", "digested linear DNA", MaterialKind.WASTE)
    graph.react("remove_linear", "Remove residual linear DNA", StepKind.DIGESTION,
                ((circle, InputRole.TEMPLATE), (exonuclease, InputRole.CATALYST)),
                (cleaned, linear_waste))
    graph.discard(linear_waste, "does not contain retained circular SRT junction molecules")
    shear = graph.input("shear_condition", "DNA shearing condition", MaterialKind.CONDITION)
    fragments = graph.product("junction_fragments", "junction-proximal circle fragments",
                              MaterialKind.MOLECULAR, state=cleaned.state)
    graph.react("shear", "Shear the circularized product", StepKind.FRAGMENTATION,
                ((cleaned, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (fragments,))
    capture = graph.input("capture_probe", "junction capture reagent", MaterialKind.SELECTION)
    captured = graph.product("captured_junction", "captured insertion-junction fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    uncaptured = graph.product("uncaptured", "uncaptured circle fragments", MaterialKind.WASTE)
    graph.react("capture_junction", "Capture the junction-bearing branch",
                StepKind.AFFINITY_CAPTURE,
                ((fragments, InputRole.TEMPLATE), (capture, InputRole.SELECTION)),
                (captured, uncaptured))
    graph.discard(uncaptured, "does not carry the selected SRT junction")
    reviewed.extend(("circle_ligase", "circularize", "exonuclease", "remove_linear",
                     "shear_condition", "shear", "capture_probe", "capture_junction"))
    _finish_pcr(graph, "junction", captured, final, d.endpoint, reviewed,
                label="scCC insertion junction")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _cotech_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input("cells", "antibody-bound permeabilized cells containing poly(A) RNA",
                        MaterialKind.BIOLOGICAL, state=_cotech().initial, primary=True)
    transposome = graph.input("barcoded_patn5", "first-well-barcoded pA–Tn5",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "protein-A–Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged_cells", "cells with target-proximal barcoded chromatin",
                           MaterialKind.MOLECULAR, state=cells.state)
    tn5_waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment antibody-targeted chromatin",
                StepKind.TAGMENTATION,
                ((cells, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tn5_waste))
    graph.discard(tn5_waste, "does not retain a target-proximal barcoded fragment")
    rt_primer = graph.input("rt_primer", "matching first-well-barcoded poly(dT) primer",
                            MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    paired = graph.product("paired_branches", "cells containing first-well-barcoded chromatin and cDNA",
                           MaterialKind.MOLECULAR,
                           state=molecule("paired intracellular branches", (
                               role("barcoded chromatin fragments"),
                               role("RT well barcode"), role("transcript cDNA"))))
    graph.react("reverse_transcribe", "Reverse-transcribe RNA with the matching well barcode",
                StepKind.REVERSE_TRANSCRIPTION,
                ((tagged, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (paired,))
    second_well = graph.input("second_well", "second-round combinatorial-index partition",
                              MaterialKind.CONDITION)
    indexed = graph.product("indexed_products", "two-well-indexed chromatin and cDNA products",
                            MaterialKind.MOLECULAR,
                            state=molecule("two-index cell identity", (
                                ident("first-round well barcode", "cotech_round1", "cell_barcode"),
                                ident("second-round well barcode", "cotech_round2", "cell_barcode"),
                                role("DNA and cDNA products"))))
    collisions = graph.product("collisions", "ambiguous second-well collisions", MaterialKind.WASTE)
    graph.react("redistribute", "Pool and redistribute products into second-index wells",
                StepKind.PARTITION,
                ((paired, InputRole.TEMPLATE), (second_well, InputRole.CONDITION)),
                (indexed, collisions), changes_molecule=False)
    graph.discard(collisions, "does not have an unambiguous ordered pair of cell indexes")
    chromatin = graph.product("chromatin_branch", "CoTECH chromatin branch",
                              MaterialKind.MOLECULAR, state=indexed.state)
    rna = graph.product("rna_branch", "CoTECH transcriptome branch",
                        MaterialKind.MOLECULAR, state=indexed.state)
    graph.react("split_libraries", "Split chromatin and transcriptome libraries",
                StepKind.SPLIT, ((indexed, InputRole.TEMPLATE),), (chromatin, rna),
                changes_molecule=False)
    reviewed = ["cells", "barcoded_patn5", "tn5", "tagment", "rt_primer",
                "reverse_transcriptase", "reverse_transcribe", "second_well",
                "redistribute", "split_libraries"]
    chrom_final = role_library("CoTECH chromatin-library roles",
                               (CELL, role("target-proximal chromatin insert")), index=SAMPLE)
    rna_final = role_library("CoTECH transcriptome-library roles",
                             (CELL, role("transcript cDNA")), index=SAMPLE)
    _finish_pcr(graph, "chromatin", chromatin, chrom_final,
                "Chromatin branch: " + d.endpoint, reviewed, label="CoTECH chromatin")
    _finish_pcr(graph, "rna", rna, rna_final,
                "RNA branch: " + d.endpoint, reviewed, label="CoTECH RNA")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _gotcha_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input("cells", "single-cell ATAC reaction containing a target genomic locus",
                        MaterialKind.BIOLOGICAL, state=_gotcha().initial, primary=True)
    bead_primer = graph.input("bead_primer", "droplet cell-barcode extension primer",
                              MaterialKind.OLIGO)
    locus_primer = graph.input("locus_primer", "locus-specific spike-in primer",
                               MaterialKind.OLIGO)
    polymerase = graph.input("extension_polymerase", "droplet extension polymerase",
                             MaterialKind.CATALYST)
    barcoded = graph.product("barcoded_products", "cell-barcoded ATAC and target-locus products",
                             MaterialKind.MOLECULAR,
                             state=molecule("barcoded ATAC and target molecules", (
                                 CELL, role("accessible-chromatin fragments"),
                                 role("target-locus amplicon"))))
    graph.react("barcode_extension", "Copy the droplet barcode onto ATAC and target molecules",
                StepKind.EXTENSION,
                ((cells, InputRole.TEMPLATE), (bead_primer, InputRole.STRUCTURAL),
                 (locus_primer, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (barcoded,))
    atac = graph.product("atac_branch", "parent scATAC branch", MaterialKind.MOLECULAR,
                         state=barcoded.state)
    target = graph.product("target_branch", "target-enrichment aliquot", MaterialKind.MOLECULAR,
                           state=barcoded.state)
    graph.react("split_branches", "Split parent ATAC and target-genotype branches",
                StepKind.SPLIT, ((barcoded, InputRole.TEMPLATE),), (atac, target),
                changes_molecule=False)
    reviewed = ["cells", "bead_primer", "locus_primer", "extension_polymerase",
                "barcode_extension", "split_branches"]
    atac_final = role_library("GoT-ChA parent ATAC roles",
                              (CELL, role("accessible-chromatin insert")), index=SAMPLE)
    _finish_pcr(graph, "atac", atac, atac_final,
                "Parent scATAC branch linked by the same cell barcode.", reviewed,
                label="GoT-ChA parent ATAC")
    nested_forward = graph.input("nested_forward", "biotinylated hemi-nested target primer",
                                 MaterialKind.OLIGO)
    nested_reverse = graph.input("nested_reverse", "target PCR reverse primer",
                                 MaterialKind.OLIGO)
    nested_polymerase = graph.input("nested_polymerase", "target PCR polymerase",
                                    MaterialKind.CATALYST)
    enriched = graph.product("enriched_target", "biotinylated variant-spanning genotype amplicon",
                             MaterialKind.MOLECULAR,
                             state=molecule("enriched genotype amplicon", (
                                 CELL, role("variant-spanning genomic insert"),
                                 role("biotinylated nested end"))))
    graph.react("nested_pcr", "Perform biotinylated hemi-nested target PCR", StepKind.PCR,
                ((target, InputRole.TEMPLATE), (nested_forward, InputRole.STRUCTURAL),
                 (nested_reverse, InputRole.STRUCTURAL),
                 (nested_polymerase, InputRole.CATALYST)), (enriched,),
                geometry=PCRGeometry.unavailable(
                    enriched.state, reason="The locus-specific panel determines the target-primer geometry."))
    beads = graph.input("streptavidin", "streptavidin capture beads", MaterialKind.SELECTION)
    captured = graph.product("captured_target", "streptavidin-captured genotype amplicon",
                             MaterialKind.MOLECULAR, state=enriched.state)
    uncaptured = graph.product("uncaptured_target", "uncaptured target-PCR material",
                               MaterialKind.WASTE)
    graph.react("capture", "Capture biotinylated target amplicons on streptavidin",
                StepKind.AFFINITY_CAPTURE,
                ((enriched, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (captured, uncaptured))
    graph.discard(uncaptured, "does not carry the selected biotinylated target amplicon")
    reviewed.extend(("nested_forward", "nested_reverse", "nested_polymerase",
                     "nested_pcr", "streptavidin", "capture"))
    _finish_pcr(graph, "genotype", captured, final, d.endpoint, reviewed,
                label="GoT-ChA genotype")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _scone_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cell = graph.input("cell", "lysed single cell or nucleus containing genomic DNA and total RNA",
                       MaterialKind.BIOLOGICAL, state=_scone().initial, primary=True)
    dna = graph.product("dna_branch", "genomic-DNA branch", MaterialKind.MOLECULAR,
                        state=cell.state)
    rna = graph.product("rna_branch", "total-RNA branch", MaterialKind.MOLECULAR,
                        state=cell.state)
    graph.react("separate_modalities", "Expose genomic DNA and total RNA as parallel substrates",
                StepKind.SPLIT, ((cell, InputRole.TEMPLATE),), (dna, rna),
                changes_molecule=False)
    transposome = graph.input("dna_transposome", "DNA-branch-barcode/UMI-loaded transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged_dna = graph.product("tagged_dna", "branch-barcoded molecularly tagged genomic DNA",
                               MaterialKind.MOLECULAR,
                               state=molecule("tagged DNA", (DNA_TYPE, DNA_UMI,
                                   role("tagmented genomic fragment"))))
    dna_waste = graph.product("dna_waste", "non-library tagmentation products",
                              MaterialKind.WASTE)
    graph.react("tagment_dna", "Tagment genomic DNA with the DNA-branch adapter",
                StepKind.TAGMENTATION,
                ((dna, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged_dna, dna_waste))
    graph.discard(dna_waste, "does not retain a molecularly tagged genomic insert")
    rt_primer = graph.input("rna_rt_primer", "RNA-branch-barcode/UMI internal RT primer",
                            MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product("cdna", "internally primed branch-barcoded cDNA",
                         MaterialKind.MOLECULAR,
                         state=molecule("first-strand cDNA", (RNA_TYPE, RNA_UMI,
                             role("full-transcript cDNA"))))
    graph.react("reverse_transcribe", "Internally prime total RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    tailing_enzyme = graph.input("tailing_enzyme", "poly(C)-tailing enzyme",
                                 MaterialKind.CATALYST)
    tailed = graph.product("tailed_cdna", "poly(C)-tailed cDNA with the opposite amplification end",
                           MaterialKind.MOLECULAR, state=cdna.state)
    graph.react("tail_cdna", "Poly(C)-tail the cDNA", StepKind.TAILING,
                ((cdna, InputRole.TEMPLATE), (tailing_enzyme, InputRole.CATALYST)), (tailed,))
    pooled = graph.product("pooled_templates", "co-amplifiable tagged DNA and cDNA",
                           MaterialKind.MOLECULAR,
                           state=molecule("co-amplifiable DNA and cDNA", (
                               role("common amplification end"), role("DNA or RNA insert"),
                               role("branch barcode and UMI"))))
    graph.react("pool_branches", "Pool DNA and RNA products for one-tube amplification",
                StepKind.POOL,
                ((tagged_dna, InputRole.TEMPLATE), (tailed, InputRole.TEMPLATE)), (pooled,))
    reviewed = ["cell", "separate_modalities", "dna_transposome", "tn5",
                "tagment_dna", "rna_rt_primer", "reverse_transcriptase",
                "reverse_transcribe", "tailing_enzyme", "tail_cdna", "pool_branches"]
    _finish_pcr(graph, "combined", pooled, final, d.endpoint, reviewed,
                label="scONE-seq combined DNA/RNA")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _scito_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input("cell_pools", "cells divided into staining pools",
                        MaterialKind.BIOLOGICAL,
                        state=molecule("cell pools", (role("surface proteins"),)), primary=True)
    antibodies = graph.input("antibodies", "universally oligo-conjugated antibodies",
                             MaterialKind.MOLECULAR)
    splints = graph.input("splints", "antibody- and pool-identifying splint oligos",
                          MaterialKind.OLIGO)
    tags = graph.product("antibody_tags", "SCITO antibody tags carrying target and pool identity",
                         MaterialKind.SELECTION,
                         state=molecule("SCITO antibody tag", (ANTIBODY, POOL,
                             role("10x feature-capture sequence"))))
    graph.react("hybridize_splints", "Hybridize pool-specific splints to universal antibody handles",
                StepKind.HYBRIDIZATION,
                ((antibodies, InputRole.TEMPLATE), (splints, InputRole.STRUCTURAL)), (tags,))
    stained = graph.product("stained_cells", "SCITO-tagged cells", MaterialKind.MOLECULAR,
                            state=molecule("tagged cell", (ANTIBODY, POOL,
                                role("surface-protein counts"))))
    unbound = graph.product("unbound_tags", "unbound antibody tags", MaterialKind.WASTE)
    graph.react("stain", "Stain cells with target- and pool-barcoded antibodies",
                StepKind.BINDING,
                ((cells, InputRole.TEMPLATE), (tags, InputRole.SELECTION)),
                (stained, unbound))
    graph.discard(unbound, "does not report a cell-associated protein measurement")
    droplets = graph.input("droplets", "10x feature-barcoding droplet condition",
                           MaterialKind.CONDITION)
    occupied = graph.product("occupied_droplets", "overloaded droplets containing pool-resolved cells",
                             MaterialKind.MOLECULAR, state=stained.state)
    empty = graph.product("empty_droplets", "empty or unresolved droplets", MaterialKind.WASTE)
    graph.react("partition", "Combine stained pools and partition cells into droplets",
                StepKind.PARTITION,
                ((stained, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (occupied, empty), changes_molecule=False)
    graph.discard(empty, "does not contain a resolvable pool/cell combination")
    bead_primer = graph.input("bead_primer", "10x cell-barcode/UMI feature-capture primer",
                              MaterialKind.OLIGO)
    polymerase = graph.input("extension_polymerase", "feature-tag extension polymerase",
                             MaterialKind.CATALYST)
    captured = graph.product("captured_tags", "droplet-barcoded SCITO feature tags",
                             MaterialKind.MOLECULAR,
                             state=molecule("captured SCITO tags", (CELL, UMI,
                                 ANTIBODY, POOL)))
    graph.react("capture_tags", "Copy droplet barcode and UMI onto antibody tags",
                StepKind.EXTENSION,
                ((occupied, InputRole.TEMPLATE), (bead_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (captured,))
    reviewed = ["cell_pools", "antibodies", "splints", "hybridize_splints", "stain",
                "droplets", "partition", "bead_primer", "extension_polymerase",
                "capture_tags"]
    _finish_pcr(graph, "feature", captured, final, d.endpoint, reviewed,
                label="SCITO feature-tag")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _phospho_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input("cells", "fixed permeabilized hashed cells with accessible chromatin",
                        MaterialKind.BIOLOGICAL, state=_phospho().initial, primary=True)
    antibodies = graph.input("antibodies", "DNA-tagged intracellular antibodies and hash tags",
                             MaterialKind.SELECTION)
    stained = graph.product("stained_cells", "antibody-tagged hashed cells",
                            MaterialKind.MOLECULAR,
                            state=molecule("tagged cell", (ANTIBODY,
                                ident("hash tag", "sample_hash", "sample_index"),
                                role("accessible chromatin"))))
    unbound = graph.product("unbound_antibodies", "unbound antibody tags", MaterialKind.WASTE)
    graph.react("stain", "Stain intracellular phosphoprotein targets",
                StepKind.BINDING,
                ((cells, InputRole.TEMPLATE), (antibodies, InputRole.SELECTION)),
                (stained, unbound))
    graph.discard(unbound, "does not report a cell-associated intracellular target")
    transposome = graph.input("transposome", "10x scATAC loaded transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagmented = graph.product("tagmented_cells", "tagged cells with tagmented accessible chromatin",
                              MaterialKind.MOLECULAR, state=stained.state)
    tn5_waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment accessible chromatin", StepKind.TAGMENTATION,
                ((stained, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagmented, tn5_waste))
    graph.discard(tn5_waste, "does not retain an amplifiable accessibility fragment")
    bridge = graph.input("bridge_oligos", "antibody-tag-to-10x-capture bridge oligos",
                         MaterialKind.OLIGO)
    bridged = graph.product("bridged_cells", "cells with capturable ATAC and antibody-tag branches",
                            MaterialKind.MOLECULAR, state=_phospho().reactions[1].after)
    graph.react("bridge_tags", "Hybridize antibody tags to 10x-compatible bridge oligos",
                StepKind.HYBRIDIZATION,
                ((tagmented, InputRole.TEMPLATE), (bridge, InputRole.STRUCTURAL)), (bridged,))
    droplets = graph.input("droplets", "10x ATAC GEM barcoding condition", MaterialKind.CONDITION)
    barcoded = graph.product("barcoded_products", "cell-barcoded ATAC and antibody-tag products",
                             MaterialKind.MOLECULAR,
                             state=molecule("barcoded Phospho-seq products", (
                                 CELL, role("ATAC insert"), UMI, ANTIBODY,
                                 role("phosphoprotein tag"))))
    failed = graph.product("failed_partitions", "empty or ambiguous GEM products", MaterialKind.WASTE)
    graph.react("partition", "Partition and barcode both molecular branches in GEMs",
                StepKind.PARTITION,
                ((bridged, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (barcoded, failed), changes_molecule=False)
    graph.discard(failed, "does not retain an interpretable cell barcode")
    atac = graph.product("atac_branch", "ATAC bead-bound branch", MaterialKind.MOLECULAR,
                         state=barcoded.state)
    tags = graph.product("tag_branch", "antibody/hash-tag supernatant branch",
                         MaterialKind.MOLECULAR, state=barcoded.state)
    graph.react("split_branches", "Separate ATAC DNA from antibody/hash-tag supernatant",
                StepKind.SPLIT, ((barcoded, InputRole.TEMPLATE),), (atac, tags),
                changes_molecule=False)
    reviewed = ["cells", "antibodies", "stain", "transposome", "tn5", "tagment",
                "bridge_oligos", "bridge_tags", "droplets", "partition", "split_branches"]
    atac_final = role_library("Phospho-seq ATAC-library roles",
                              (CELL, role("accessible-chromatin insert")), index=SAMPLE)
    tag_final = role_library("Phospho-seq tag-library roles",
                             (CELL, UMI, ANTIBODY, role("hash or phosphoprotein tag")), index=SAMPLE)
    _finish_pcr(graph, "atac", atac, atac_final, "ATAC branch: " + d.endpoint,
                reviewed, label="Phospho-seq ATAC")
    _finish_pcr(graph, "tags", tags, tag_final, "Antibody/hash branch: " + d.endpoint,
                reviewed, label="Phospho-seq tag")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _cropseq_multi_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cassette = graph.input("cassette", "CROPseq-multi lentiviral 3-prime-LTR cassette",
                           MaterialKind.MOLECULAR, state=_cropseq_multi().initial, primary=True)
    cells = graph.input("target_cells", "target cells", MaterialKind.BIOLOGICAL)
    integrase = graph.input("integrase", "lentiviral integration machinery",
                            MaterialKind.CATALYST)
    integrated = graph.product("integrated_cassette", "integrated duplicated-LTR perturbation cassette",
                               MaterialKind.MOLECULAR,
                               state=molecule("integrated perturbation cassette", (
                                   role("5-prime LTR copy"), role("selection-gene transcript"),
                                   role("3-prime LTR copy"))))
    graph.react("integrate", "Integrate and duplicate the engineered 3-prime LTR",
                StepKind.INTEGRATION,
                ((cassette, InputRole.TEMPLATE), (cells, InputRole.ANALYTE),
                 (integrase, InputRole.CATALYST)), (integrated,))
    rna_polymerase = graph.input("rna_polymerase", "Pol III transcription machinery",
                                 MaterialKind.CATALYST)
    guides = graph.product("guides", "tRNA-processed guides carrying linked iBARs",
                           MaterialKind.MOLECULAR,
                           state=molecule("processed perturbation guides", (
                               role("sgRNA 1 with iBAR"), role("sgRNA 2 with iBAR"))))
    graph.react("transcribe_guides", "Transcribe and tRNA-process the multiplexed guide array",
                StepKind.TRANSCRIPTION,
                ((integrated, InputRole.TEMPLATE), (rna_polymerase, InputRole.CATALYST)),
                (guides,))
    mrna = graph.product("mrna_branch", "poly(A)-mRNA perturbation-identity branch",
                         MaterialKind.MOLECULAR, state=guides.state)
    direct = graph.product("direct_branch", "direct-guide-capture branch",
                           MaterialKind.MOLECULAR, state=guides.state)
    t7 = graph.product("t7_branch", "CSMv2 T7-IVT branch", MaterialKind.MOLECULAR,
                       state=guides.state)
    graph.react("split_readouts", "Route perturbation identity to the three supported readouts",
                StepKind.SPLIT, ((guides, InputRole.TEMPLATE),), (mrna, direct, t7),
                changes_molecule=False)
    reviewed = ["cassette", "target_cells", "integrase", "integrate", "rna_polymerase",
                "transcribe_guides", "split_readouts"]
    mrna_final = role_library("CROPseq-multi mRNA-readout roles",
                              (CELL, UMI, GUIDE, role("linked guide spacer")), index=SAMPLE)
    direct_final = role_library("CROPseq-multi direct-guide roles",
                                (CELL, UMI, GUIDE, role("linked guide spacer")), index=SAMPLE)
    _finish_rt_pcr(graph, "mrna", mrna, mrna_final, "mRNA branch: " + d.endpoint,
                   reviewed, label="CROPseq-multi mRNA")
    _finish_rt_pcr(graph, "direct", direct, direct_final,
                   "direct-guide branch: " + d.endpoint, reviewed,
                   label="CROPseq-multi direct guide")
    t7_polymerase = graph.input("t7_rna_polymerase", "T7 RNA polymerase", MaterialKind.CATALYST)
    t7_rna = graph.product("t7_rna", "T7-amplified CSMv2 guide RNA", MaterialKind.MOLECULAR,
                           state=t7.state)
    graph.react("t7_transcription", "Amplify the CSMv2 guide-barcode product by T7 IVT",
                StepKind.TRANSCRIPTION,
                ((t7, InputRole.TEMPLATE), (t7_polymerase, InputRole.CATALYST)), (t7_rna,))
    reviewed.extend(("t7_rna_polymerase", "t7_transcription"))
    _finish_rt_pcr(graph, "t7", t7_rna, final, "T7-IVT branch: " + d.endpoint,
                   reviewed, label="CROPseq-multi T7-IVT")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _chair_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    nucleus = graph.input("nucleus", "crosslinked nucleus containing chromatin and poly(A) RNA",
                          MaterialKind.BIOLOGICAL, state=_chair().initial, primary=True)
    restriction = graph.input("restriction_enzyme", "in-situ restriction endonuclease",
                              MaterialKind.CATALYST)
    digested = graph.product("digested", "restricted crosslinked chromatin",
                             MaterialKind.MOLECULAR,
                             state=molecule("restricted chromatin", (
                                 role("restriction fragment A"),
                                 role("crosslinked restriction fragment B"), role("poly(A) RNA"))))
    graph.react("digest", "Digest chromatin in situ", StepKind.DIGESTION,
                ((nucleus, InputRole.TEMPLATE), (restriction, InputRole.CATALYST)), (digested,))
    ligase = graph.input("proximity_ligase", "proximity DNA ligase", MaterialKind.CATALYST)
    contacts = graph.product("contacts", "chromatin-contact junctions retained with RNA",
                             MaterialKind.MOLECULAR,
                             state=molecule("3D contact junction", (
                                 role("restriction fragment A"), role("restriction fragment B"),
                                 role("poly(A) RNA")),
                                 junctions=(("restriction fragment A", "restriction fragment B",
                                             "proximity ligation"),)))
    graph.react("proximity_ligation", "Proximity-ligate crosslinked restriction ends",
                StepKind.PROXIMITY_LIGATION,
                ((digested, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (contacts,))
    transposome = graph.input("transposome", "ATAC loaded transposome", MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "accessible-anchor and contact products retained with RNA",
                           MaterialKind.MOLECULAR, state=_chair().reactions[2].after)
    tn5_waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment accessible chromatin while retaining contact products",
                StepKind.TAGMENTATION,
                ((contacts, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tn5_waste))
    graph.discard(tn5_waste, "does not retain an amplifiable accessible or contact fragment")
    droplets = graph.input("droplets", "10x Multiome droplet-barcoding condition",
                           MaterialKind.CONDITION)
    barcoded = graph.product("barcoded", "cell-barcoded RNA, ATAC and contact products",
                             MaterialKind.MOLECULAR,
                             state=molecule("three matched readouts", (CELL, UMI,
                                 role("ChAIR-RNA"), role("ChAIR-ATAC"), role("ChAIR-PET"))))
    failed = graph.product("failed_partitions", "empty or ambiguous droplets", MaterialKind.WASTE)
    graph.react("partition", "Barcode DNA and RNA together in Multiome droplets",
                StepKind.PARTITION,
                ((tagged, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (barcoded, failed), changes_molecule=False)
    graph.discard(failed, "does not retain an interpretable cell barcode")
    rna = graph.product("rna_branch", "ChAIR-RNA branch", MaterialKind.MOLECULAR,
                        state=barcoded.state)
    atac = graph.product("atac_branch", "ChAIR-ATAC branch", MaterialKind.MOLECULAR,
                         state=barcoded.state)
    pet = graph.product("pet_branch", "ChAIR-PET contact branch", MaterialKind.MOLECULAR,
                        state=barcoded.state)
    graph.react("split_readouts", "Separate RNA, accessibility and contact readouts",
                StepKind.SPLIT, ((barcoded, InputRole.TEMPLATE),), (rna, atac, pet),
                changes_molecule=False)
    reviewed = ["nucleus", "restriction_enzyme", "digest", "proximity_ligase",
                "proximity_ligation", "transposome", "tn5", "tagment", "droplets",
                "partition", "split_readouts"]
    rna_final = role_library("ChAIR-RNA roles", (CELL, UMI, role("transcript cDNA")), index=SAMPLE)
    atac_final = role_library("ChAIR-ATAC roles", (CELL, role("accessible-chromatin insert")), index=SAMPLE)
    pet_final = role_library("ChAIR-PET roles", (CELL, role("contact-junction insert")), index=SAMPLE)
    _finish_rt_pcr(graph, "rna", rna, rna_final, "ChAIR-RNA: " + d.endpoint,
                   reviewed, label="ChAIR RNA")
    _finish_pcr(graph, "atac", atac, atac_final, "ChAIR-ATAC: " + d.endpoint,
                reviewed, label="ChAIR ATAC")
    _finish_pcr(graph, "pet", pet, pet_final, "ChAIR-PET: " + d.endpoint,
                reviewed, label="ChAIR contact")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def get(key: str):
    d = DEFS[key]
    workflows = {"sccc": _sccc, "cotech": _cotech, "got-cha": _gotcha,
                 "scone-seq": _scone, "scito-seq": _scito,
                 "phospho-seq": _phospho, "cropseq-multi": _cropseq_multi,
                 "chair": _chair}
    payloads = {
        "sccc": (CELL, UMI, role("transposon–genome junction")),
        "cotech": (CELL, role("chromatin or cDNA insert")),
        "got-cha": (CELL, role("variant-spanning genomic insert")),
        "scone-seq": (role("DNA/RNA branch barcode"), role("branch UMI"), role("insert")),
        "scito-seq": (CELL, UMI, ANTIBODY, POOL),
        "phospho-seq": (CELL, UMI, ANTIBODY, role("ATAC or protein insert")),
        "cropseq-multi": (CELL, UMI, GUIDE, role("linked guide spacer")),
        "chair": (CELL, UMI, role("RNA, ATAC or contact insert")),
    }
    final = role_library(d.title + " final library roles", payloads[key], index=SAMPLE)
    spec = SimpleNamespace(
        TITLE=d.title, NOTES=d.notes,
        SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
        SUMMARY=d.summary, CAVEAT="", FINAL_LIBRARY=final, SEQ_PRIMERS=(),
        FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint, READ_LENGTHS=None,
        SEQUENCING_ENDING=d.endpoint,
        SEQUENCING_UNAVAILABLE="Target panels, pool barcodes and kit revisions vary; the defining source establishes the read roles, while this page does not invent a universal primer sequence.",
        workflow=workflows[key], oligos=lambda: (),
    )
    if key == "scone-seq":
        spec.CAVEAT = ("Supplementary Table S2 publishes the custom molecular adapters, "
                       "amplification oligos and custom Illumina Read 1, Read 2 and i7 "
                       "sequencing primers. Its ONE-Illu i5 row instead directs users to "
                       "the sequencing-equipment instructions; that one primer is not "
                       "reconstructed here.")
        spec.SEQUENCING_UNAVAILABLE = ("Custom Read 1, Read 2 and i7 primer sequences are "
                                       "published and listed above. The source leaves the "
                                       "ONE-Illu i5 primer to the instrument manufacturer, "
                                       "so a complete four-primer landing model remains mixed.")
        spec.oligos = _scone_oligos
    elif key == "scito-seq":
        spec.CAVEAT = ("The defining paper fixes the splint architecture and lengths, "
                       "including a 20-nt antibody handle, TruSeq Read 2 region, 10-nt "
                       "compound antibody/pool barcode and 22-nt 10x feature-capture "
                       "region. It does not print the handle or barcode bases, and the "
                       "completed library follows the commercial 10x 3′ v3 Feature "
                       "Barcoding workflow.")
        spec.SEQUENCING_UNAVAILABLE = ("The original article, author manuscript, "
                                       "Supplementary Tables 1–4 and authors' public "
                                       "analysis repository do not disclose a universal "
                                       "SCITO handle/barcode set. Read geometry also "
                                       "inherits the commercial 10x 3′ v3 endpoint.")
    elif key == "phospho-seq":
        spec.CAVEAT = ("The Phospho-seq paper prints the complete home-conjugated tag "
                       "architectures, tag-library PCR primers and run lengths. It names "
                       "BOA and BOB but inherits their exact sequences from the cited "
                       "ASAP-seq bridge method; those published sequences are listed "
                       "above. Antibody barcodes and i7 bases are intentionally variable.")
        spec.SEQUENCING_UNAVAILABLE = ("The antibody-tag branch has published PCR "
                                       "oligos and a 34-base Read 1, 8-base i7, 16-base "
                                       "i5 and 34-base Read 2 run. Its cell-barcode end "
                                       "is supplied by the versioned commercial 10x "
                                       "scATAC gel-bead chemistry, so this shared page "
                                       "does not fabricate that proprietary oligo.")
        spec.oligos = _phospho_oligos
    factories = {
        "sccc": _sccc_graph,
        "cotech": _cotech_graph,
        "got-cha": _gotcha_graph,
        "scone-seq": _scone_graph,
        "scito-seq": _scito_graph,
        "phospho-seq": _phospho_graph,
        "cropseq-multi": _cropseq_multi_graph,
        "chair": _chair_graph,
    }
    if key in factories:
        spec.protocol_graph = lambda: factories[key](d, final)
    return spec
