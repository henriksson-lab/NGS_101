"""Shared molecular models for single-cell multimodal identifier transfer.

The central invariant is identifier provenance: every feature says whether it identifies
a cell, molecule, sample/pool, antibody or genomic target.  Branching workflows keep the
same cell identity while allowing physically distinct sequencing libraries.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

from chemdraw import feature
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
    return SimpleNamespace(
        TITLE=d.title, NOTES=d.notes,
        SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
        SUMMARY=d.summary, CAVEAT="", FINAL_LIBRARY=final, SEQ_PRIMERS=(),
        FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint, READ_LENGTHS=None,
        SEQUENCING_ENDING=d.endpoint,
        SEQUENCING_UNAVAILABLE="Target panels, pool barcodes and kit revisions vary; the defining source establishes the read roles, while this page does not invent a universal primer sequence.",
        workflow=workflows[key], oligos=lambda: (),
    )
