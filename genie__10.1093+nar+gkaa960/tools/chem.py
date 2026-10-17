"""GenIE targeted cDNA readout (Cooper et al., NAR 2020, doi:10.1093/nar/gkaa960)."""
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "GenIE — endogenous variant editing with targeted cDNA sequencing"
NOTES = "01_genie.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1093/nar/gkaa960">Cooper et al., <i>Nucleic Acids Research</i> (2020)</a>.'
SUMMARY = ("Arrayed CRISPR editing creates endogenous alleles; a gene-specific RT-PCR "
           "amplicon from RNA and a matched genomic-DNA amplicon from the same cells are "
           "barcoded and sequenced, and their allele ratios are compared.")

# Supplementary Table 3 prints one common tail per primer orientation, "added to the 5'
# end of the forward and reverse primers respectively".  The forward tail is the TruSeq
# Read 1 primer; the reverse tail is the Illumina paired-end / Genome-Analyzer-era arm of
# the cited 96-plex barcoding method (ref 8, Kozarewa & Turner 2011), NOT the current
# TruSeq Read 2 arm.
FORWARD_TAIL = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
REVERSE_TAIL = il.PE_READ2_ARM
assert FORWARD_TAIL == il.TRUSEQ_READ1 and REVERSE_TAIL != il.TRUSEQ_READ2

# One complete locus from Supplementary Table 3 keeps the drawn molecule concrete; every
# genomic part is an experiment-specific parameter.
TARGET = "MUL1_rs6700034"
RT_PRIMER = "CTTGAGCTCTTGGGAGACC"
FORWARD_LOCUS = "GAGGATGAACTGGCACAGC"
REVERSE_LOCUS = "CTGGGCGACTGATTTCGAG"
# Table 3 prints tails and locus primers in separate cells; the ordered oligo is their
# documented concatenation, so the graph states keep the verbatim locus halves.
FORWARD_OLIGO = FORWARD_TAIL + FORWARD_LOCUS
REVERSE_OLIGO = REVERSE_TAIL + REVERSE_LOCUS

CAVEAT = ("The adaptor tails, the gene-specific RT primer and the amplicon primers above "
          "are verbatim from Supplementary Table 3 for one locus; every locus has its own "
          "set. The paper delegates the barcoding PCR to ref 8 and prints neither its "
          "primers nor which amplicon end receives P5, so only the outer flow-cell and "
          "index shell is inferred.")

_UNPRINTED = ("the second PCR adds the MiSeq barcodes of ref 8 (Kozarewa & Turner 2011); "
              "the paper prints neither those primers nor which amplicon end receives P5")


def _library(insert_name: str, name: str) -> Construct:
    """Indexed amplicon: exact primer-installed arms, inferred terminal shell."""
    return Construct([
        seg("P5", il.P5, "p5", inferred=True),
        seg("i5", "J" * 8, "cbc", placeholder=True, inferred=True,
            feature=feature("sample_index_i5", "sample_index", "unknown")),
        seg("forward adaptor tail (TruSeq Read 1 arm)", FORWARD_TAIL, "r1"),
        seg(f"{TARGET} forward primer", FORWARD_LOCUS),
        seg(insert_name, "X" * 56, placeholder=True),
        seg(f"{TARGET} reverse primer complement", revcomp(REVERSE_LOCUS)),
        seg("reverse adaptor tail complement (Illumina PE arm)", revcomp(REVERSE_TAIL), "r2"),
        seg("i7 reverse complement", "I" * 8, "cbc", placeholder=True, inferred=True,
            feature=feature("sample_index_i7", "sample_index", "unknown")),
        seg("P7 reverse complement", il.P7_RC, "p7", inferred=True),
    ], name=name)


FINAL_LIBRARY = _library("edited-locus cDNA interval", "GenIE targeted cDNA library")
GDNA_LIBRARY = _library("edited-locus genomic interval", "GenIE matched gDNA library")
SEQ_PRIMERS = (
    sp.custom("Read 1", "forward adaptor tail primer (TruSeq Read 1 arm)", FORWARD_TAIL,
              "Cooper et al. 2020 Supplementary Table 3 'Forward primer adaptor sequence'"),
    sp.custom("Read 2", "reverse adaptor tail primer (Illumina PE arm)", REVERSE_TAIL,
              "Cooper et al. 2020 Supplementary Table 3 'Reverse primer adaptor sequence'",
              note="Printed as the common reverse-primer tail; the matching sequencing "
                   "primer itself is delegated to ref 8."),
)
FINAL_CAPTION = ("Representative " + TARGET + " cDNA amplicon. Both internal boundaries and "
                 "both adaptor tails are source-exact; the flow-cell and index shell is "
                 "inferred because " + _UNPRINTED + ".")
EXTRA_MOLECULAR_STATES = ((
    "Matched genomic-DNA library",
    GDNA_LIBRARY,
    "The same primer pair amplifies the locus from MagAttract-purified gDNA of the same "
    "cell pellets; cDNA and gDNA allele ratios are compared. For SDF4 only, Supplementary "
    "Table 3 gives a separate spliced cDNA forward primer."),)
READ_LENGTHS = {"Read 1": 150, "Read 2": 150}
SEQUENCING_INTRO = ("Amplicons are kept below 295 bp so that a 150-base paired-end MiSeq "
                    "run reads across the edited site from both ends.")


def oligos():
    return (
        f"<b>Forward primer adaptor tail</b> — 5′-{FORWARD_TAIL}-3′",
        f"<b>Reverse primer adaptor tail</b> — 5′-{REVERSE_TAIL}-3′",
        f"<b>{TARGET} gene-specific RT primer</b> — 5′-{RT_PRIMER}-3′",
        f"<b>{TARGET} forward oligo</b> — 5′-{FORWARD_OLIGO}-3′",
        f"<b>{TARGET} reverse oligo</b> — 5′-{REVERSE_OLIGO}-3′",
    )


def protocol_graph():
    source="Defining source: https://doi.org/10.1093/nar/gkaa960"
    inferred_source=source+" (barcoding PCR delegated to ref 8)"
    graph=ProtocolGraph("GenIE",authority=source)
    cells=graph.input("cells","arrayed hiPSC population",MaterialKind.BIOLOGICAL,primary=True)
    guide=graph.input("guide","full-length synthetic modified sgRNA (Synthego)",MaterialKind.OLIGO)
    donor=graph.input("hdr_donor","100-nt ssODN HDR template (IDT Ultramer)",MaterialKind.OLIGO)
    editor=graph.input("editor","purified SpCas9 or VRQR SpCas9 RNP",MaterialKind.CATALYST)
    edited=graph.product("edited_cells","nucleofected endogenous allelic-series population",
                         MaterialKind.BIOLOGICAL)
    graph.react("edit","Nucleofect Cas9 RNP with the ssODN donor",StepKind.GENOME_EDITING,
                ((cells,InputRole.TEMPLATE),(guide,InputRole.STRUCTURAL),
                 (donor,InputRole.STRUCTURAL),(editor,InputRole.CATALYST)),(edited,))
    # The same frozen cell pellets feed two parallel branches.
    rna_extraction=graph.input("rna_extraction",
                               "Direct-zol/TRI-Reagent extraction with in-column and TURBO DNase",
                               MaterialKind.CONDITION)
    rna=graph.product("rna","DNase-treated total RNA including nascent transcripts",
                      MaterialKind.MOLECULAR)
    graph.react("extract_rna","Extract RNA with a Trizol-based kit",StepKind.EXTRACTION,
                ((edited,InputRole.TEMPLATE),(rna_extraction,InputRole.CONDITION)),(rna,))
    gdna_extraction=graph.input("gdna_extraction","MagAttract HMW kit with proteinase K lysis",
                                MaterialKind.CONDITION)
    gdna=graph.product("gdna","genomic DNA from the same edited pellets",MaterialKind.MOLECULAR)
    graph.react("extract_gdna","Purify genomic DNA from matched pellets",StepKind.EXTRACTION,
                ((edited,InputRole.TEMPLATE),(gdna_extraction,InputRole.CONDITION)),(gdna,))
    rt_primer=graph.input("rt_primer",f"gene-specific RT primer ({TARGET}, 0.1 µM)",
                          MaterialKind.OLIGO,state=RT_PRIMER)
    rt=graph.input("reverse_transcriptase","SuperScript IV reverse transcriptase",
                   MaterialKind.CATALYST)
    cdna=graph.product("cdna","locus-containing first-strand cDNA",MaterialKind.MOLECULAR)
    graph.react("reverse_transcription",
                "Reverse-transcribe with a gene-specific primer 30-50 bp outside the amplicon",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(cdna,))
    fwd=graph.input("forward_primer",f"tailed {TARGET} forward amplicon primer",
                    MaterialKind.OLIGO,state=FORWARD_LOCUS)
    rev=graph.input("reverse_primer",f"tailed {TARGET} reverse amplicon primer",
                    MaterialKind.OLIGO,state=REVERSE_LOCUS)
    polymerase=graph.input("polymerase","PowerUp SYBR Green master mix",MaterialKind.CATALYST)
    cdna_amplicon=graph.product("cdna_amplicon","tailed cDNA amplicon across the edited site",
                                MaterialKind.MOLECULAR)
    graph.react("cdna_pcr","Amplify the edited locus from cDNA (PCR1, 30 cycles)",StepKind.PCR,
                ((cdna,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),
                 (rev,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(cdna_amplicon,),
                geometry=PCRGeometry.unavailable(FINAL_LIBRARY,
                    reason="The amplified genomic interval is locus-specific; the drawn molecule keeps one tabulated primer pair and a placeholder interval."))
    gdna_amplicon=graph.product("gdna_amplicon","tailed genomic amplicon across the edited site",
                                MaterialKind.MOLECULAR)
    graph.react("gdna_pcr","Amplify the same locus from genomic DNA (PCR1, 30 cycles)",
                StepKind.PCR,
                ((gdna,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),
                 (rev,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(gdna_amplicon,),
                geometry=PCRGeometry.unavailable(GDNA_LIBRARY,
                    reason="Same tailed primer pair as the cDNA branch, except for SDF4 where the supplement gives a separate spliced cDNA forward primer."))
    i5=graph.input("index_forward","MiSeq barcoding forward primer",MaterialKind.OLIGO,
                   inferred=True,source=inferred_source)
    i7=graph.input("index_reverse","MiSeq barcoding reverse primer",MaterialKind.OLIGO,
                   inferred=True,source=inferred_source)
    library=graph.product("cdna_library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY,inferred=True,source=inferred_source)
    graph.react("cdna_barcoding_pcr","Add MiSeq barcodes to the cDNA amplicon (PCR2)",
                StepKind.PCR,
                ((cdna_amplicon,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
                 (i7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.unavailable(FINAL_LIBRARY,
                    reason="The barcoding primers are delegated to ref 8, so only the terminal shell is inferred."))
    gdna_library=graph.product("gdna_library",GDNA_LIBRARY.name,MaterialKind.LIBRARY,
                               state=GDNA_LIBRARY,inferred=True,source=inferred_source)
    graph.react("gdna_barcoding_pcr","Add MiSeq barcodes to the genomic amplicon (PCR2)",
                StepKind.PCR,
                ((gdna_amplicon,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
                 (i7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(gdna_library,),
                geometry=PCRGeometry.unavailable(GDNA_LIBRARY,
                    reason="The barcoding primers are delegated to ref 8, so only the terminal shell is inferred."))
    plan=dict(primers=SEQ_PRIMERS,acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
              read_lengths=READ_LENGTHS,read_uses={"Read 1":"include","Read 2":"include"},
              limitation=("The terminal MiSeq index shell is inferred: "+_UNPRINTED+
                          ". Index cycles are therefore not represented as acquired reads."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,**plan))
    graph.sequence(gdna_library,SequencingPlan(library=GDNA_LIBRARY,**plan))
    graph.declare_reviewed_scope((
        "cells","guide","hdr_donor","editor","edit","rna_extraction","extract_rna",
        "gdna_extraction","extract_gdna","rt_primer","reverse_transcriptase",
        "reverse_transcription","forward_primer","reverse_primer","polymerase",
        "cdna_pcr","gdna_pcr","index_forward","index_reverse","cdna_barcoding_pcr",
        "gdna_barcoding_pcr","cdna_library","gdna_library"),source)
    return graph.require_publishable()
