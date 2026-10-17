"""Original paired-end eCLIP molecular workflow and exact ENCODE SOP adapters."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from adapter_ligation import reported_adapter_end
from batch_ngs import seg
from chemdraw import Construct, Scene, complement_segments, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)

TITLE = "eCLIP — enhanced crosslinking and immunoprecipitation"
NOTES = "01_eclip.html"
SOURCE = 'Defining paper and protocol: <a href="https://doi.org/10.1038/nmeth.3810">Van Nostrand et al. (2016)</a>; <a href="https://www.encodeproject.org/documents/842f7424-5396-424a-a1a3-3f18707c3222/@@download/attachment/eCLIP_SOP_v1.P_110915.pdf">ENCODE eCLIP SOP v1.P</a>.'
SUMMARY = "UV-crosslinked RBP–RNA complexes are immunoprecipitated and RNase-trimmed. A 3′ RNA adapter precedes reverse transcription; a second adapter is ligated to cDNA so both read-through and crosslink-truncated cDNAs can amplify."
CAVEAT = "The SOP adds two different barcoded RNA adapters to each sample (for example the colour-balanced X1A + X1B pair) and lists several D5/D7 index choices. The panels draw the X1A, D501 and D701 examples; positions, widths and primer family do not vary between those choices. Read lengths follow the defining paper (2 × 55 for the N10 cDNA adapter); acquired index-read cycles are not stated there."
I5_FEATURE=feature("sample_index_i5","sample_index","fixed")
I7_FEATURE=feature("sample_index_i7","sample_index","fixed")
INLINE_FEATURE=feature("inline_rna_barcode","inline_barcode","fixed")
RANDOMER_5=feature("molecular_barcode_rna_adapter","umi","random",group="molecular_barcode",part="RNA adapter")
RANDOMER_10=feature("molecular_barcode_cdna_adapter","umi","random",group="molecular_barcode",part="cDNA adapter")

RIL19 = "AGATCGGAAGAGCGTCGTG"
X1A_TAIL = "AGATCGGAAGAGCGTCGTGTAG"
RNA_X1A_BARCODE = "ATATAGG"
RNA_X1B_BARCODE = "AATAGCA"   # SOP: RNA_X1A and RNA_X1B are the colour-balanced pair
AR17 = "ACACGACGCTCTTCCGA"
RAND103TR3_FIXED = "AGATCGGAAGAGCACACGTCTG"
D501_I5 = "TATAGCCT"
D701_I7_OLIGO = "CGAGTAAT"
# Supplementary Protocol 1 (supplement-2 text lines 357-365)
PCR_F_D501 = il.P5 + D501_I5 + il.TRUSEQ_READ1
PCR_R_D701 = "CAAGCAGAAGACGGCATACGAGAT" + D701_I7_OLIGO + "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC"
RAND103TR3 = "N" * 10 + RAND103TR3_FIXED

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"),
    seg("D501 i5", D501_I5, "cbc", feature=I5_FEATURE),
    seg("Read 1 arm", il.TRUSEQ_READ1, "r1"),
    seg("X1A five-base randomer", "U" * 5, "umi", placeholder=True, feature=RANDOMER_5),
    seg("X1A inline RNA barcode read", revcomp(RNA_X1A_BARCODE), "cbc", feature=INLINE_FEATURE),
    seg("RBP-bound RNA cDNA", "X" * 34, placeholder=True),
    seg("rand103Tr3 ten-base randomer", "V" * 10, "umi", placeholder=True, feature=RANDOMER_10),
    seg("Read 2 arm first base", "A", "r2"),
    seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
    seg("D701 i7 reverse complement", "ATTACTCG", "cbc", feature=I7_FEATURE),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="paired-end dual-index eCLIP library")
SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"], sp.TRUSEQ["R2"])
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS):
    raise ValueError("eCLIP final library: " + "; ".join(problems))

FINAL_CAPTION = "Original paired-end eCLIP amplicon using the SOP's X1A/D501/D701 examples. Read 1 starts with the five-base randomer and reverse complement of the seven-base inline RNA barcode; Read 2 starts with the ten-base cDNA-adapter randomer."
READ_LENGTHS = {"Read 1": 55, "Read 2": 55}
SEQUENCING_INTRO = "Original eCLIP was sequenced on HiSeq 2500 or 4000 as paired-end 50-nt (N5 cDNA adapter) or 55-nt (N10 adapter, drawn here) reads, with D5/i5 and D7/i7 index reads whose cycle counts the paper does not state. Read 1 begins with five random bases then the reverse complement of the RNA-adapter barcode; Read 2 begins at rand103Tr3's ten-base randomer and approaches the RNA 5′ end/crosslink-stop boundary."


def crosslinked_rna_rows():
    rna = [
        seg("RNA 5-prime flank", "X" * 14, placeholder=True),
        seg("crosslinked nucleotide", "X", "umi", placeholder=True),
        seg("RNA 3-prime flank", "X" * 14, placeholder=True),
    ]
    sc = Scene()
    sc.strand("RBP-bound RNA", rna)
    sc.mark("RBP-bound RNA", "crosslinked nucleotide", "UV-crosslinked RBP")
    sc.note("RBP-bound RNA", "partial RNase + antibody-bead enrichment")
    return sc.rows()


def rna_adapter_pair():
    """SOP RNA_X1A + RNA_X1B: /5Phos/ barcode, rN5, 22-base tail, /3SpC3/ (RNA)."""
    sc = Scene()
    for name, barcode in (("RNA_X1A", RNA_X1A_BARCODE), ("RNA_X1B", RNA_X1B_BARCODE)):
        sc.strand(name, [
            seg(f"{name} inline barcode", barcode.replace("T", "U"), "cbc",
                feature=INLINE_FEATURE),
            seg(f"{name} five-base randomer", "N" * 5, "umi", placeholder=True,
                feature=RANDOMER_5),
            seg(f"{name} 3-prime tail", X1A_TAIL.replace("T", "U"), "r1"),
        ], mod5="p", mod3="SpC3")
    return sc


def rna_adapter_ligation_rows():
    """On-bead ligation of the 5'-phosphorylated, 3'-SpC3-blocked X1A RNA adapter."""
    parts = [
        seg("RBP-bound RNA fragment", "X" * 24, placeholder=True),
        seg("X1A inline barcode", RNA_X1A_BARCODE.replace("T", "U"), "cbc",
            feature=INLINE_FEATURE),
        seg("X1A five-base randomer", "N" * 5, "umi", placeholder=True, feature=RANDOMER_5),
        seg("X1A 3-prime tail", X1A_TAIL.replace("T", "U"), "r1"),
    ]
    sc = Scene()
    sc.strand("adapter-ligated RNA", parts, mod3="SpC3")
    sc.junction("adapter-ligated RNA", "RBP-bound RNA fragment", "X1A inline barcode",
                "T4 RNA ligase 1 junction")
    sc.note("adapter-ligated RNA", "fragment 3'-OH after FastAP + T4 PNK; adapter 5'-phosphate")
    return sc.rows()


def reverse_transcription_rows():
    rna = [
        seg("RNA 5-prime flank", "X" * 14, placeholder=True),
        seg("peptide-remnant crosslink", "X", "umi", placeholder=True),
        seg("RNA 3-prime fragment", "X" * 20, placeholder=True),
        seg("3-prime RNA adapter", "X" * 18, "r1", placeholder=True),
    ]
    sc = Scene()
    sc.strand("RNA", rna)
    truncated = complement_segments(rna[2:])
    read_through = complement_segments(rna)
    sc.anneal("crosslink-truncated cDNA", truncated, to="RNA",
              pair=("3-prime RNA adapter'", "3-prime RNA adapter"))
    sc.arrow("crosslink-truncated cDNA", "RT stops at peptide remnant")
    sc.anneal("read-through cDNA", read_through, to="RNA",
              pair=("3-prime RNA adapter'", "3-prime RNA adapter"))
    sc.arrow("read-through cDNA", "read-through RT")
    sc.mark("RNA", "peptide-remnant crosslink", "crosslink site")
    return sc.rows()


def cdna_adapter_ligation_rows():
    """rand103Tr3's 5'-phosphate joins the cDNA 3' end; its 3' end is SpC3-blocked."""
    def product(name, cdna_name):
        parts = [
            seg("AR17-primed RNA-adapter complement", "X" * 18, "r1", placeholder=True),
            seg(cdna_name, "X" * 28, placeholder=True),
            seg("rand103Tr3 ten-base randomer", "N" * 10, "umi", placeholder=True,
                feature=RANDOMER_10),
            seg("rand103Tr3 fixed region", RAND103TR3_FIXED, "r2"),
        ]
        sc = Scene()
        sc.strand(name, parts, mod3="SpC3")
        sc.junction(name, cdna_name, "rand103Tr3 ten-base randomer", "cDNA 3' ligation")
        return sc.rows()

    return [
        *product("crosslink-truncated product", "cDNA ending at crosslink"),
        *product("read-through product", "cDNA to RNA 5' end"),
    ]


def sections():
    return [
        ("Crosslink, trim and immunoprecipitate", crosslinked_rna_rows(), "UV fixes direct contacts; limited RNase and RBP immunoprecipitation define the captured RNA fragment."),
        ("Ligate the 3′ RNA adapter on beads", rna_adapter_ligation_rows(), "After FastAP/T4 PNK end repair, the barcoded RNA adapter is ligated on beads before gel and membrane size selection."),
        ("Reverse transcription can stop at the crosslink", reverse_transcription_rows(), "Both truncated and read-through cDNAs are retained."),
        ("Ligate rand103Tr3 to cDNA", cdna_adapter_ligation_rows(), "The cDNA 3′ ligation adds a ten-base randomer and bypasses iCLIP circularization, making both cDNA classes amplifiable."),
    ]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    complexes=graph.input("rbp_rna","RBP-bound cellular RNA",MaterialKind.BIOLOGICAL,
                          primary=True)
    uv=graph.input("uv","UV cross-linking condition",MaterialKind.CONDITION)
    crosslinked=graph.product("crosslinked","UV-crosslinked RBP-RNA complexes",
                              MaterialKind.MOLECULAR,state=tuple(crosslinked_rna_rows()))
    graph.react("crosslink","Covalently cross-link direct RBP-RNA contacts",
                StepKind.CROSS_LINKING,
                ((complexes,InputRole.TEMPLATE),(uv,InputRole.CONDITION)),(crosslinked,))
    rnase=graph.input("rnase","limited RNase I (Ambion)",MaterialKind.CATALYST)
    trimmed=graph.product("trimmed","RNase-trimmed crosslinked RBP-RNA fragments",
                          MaterialKind.MOLECULAR)
    digestion_waste=graph.product("digestion_waste","released RNA fragments",
                                   MaterialKind.WASTE)
    graph.react("trim","Partially digest unprotected RNA",StepKind.DIGESTION,
                ((crosslinked,InputRole.TEMPLATE),(rnase,InputRole.CATALYST)),
                (trimmed,digestion_waste))
    graph.discard(digestion_waste,"RNA outside the protected RBP footprint is digested")
    antibody=graph.input("antibody_beads","RBP-specific antibody beads",
                         MaterialKind.SELECTION)
    enriched=graph.product("enriched","immunoprecipitated RBP-RNA fragments",
                           MaterialKind.MOLECULAR)
    unbound=graph.product("unbound","unbound RNA-protein material",MaterialKind.WASTE)
    graph.react("immunoprecipitate","Immunoprecipitate the crosslinked RBP",
                StepKind.AFFINITY_CAPTURE,
                ((trimmed,InputRole.TEMPLATE),(antibody,InputRole.SELECTION)),
                (enriched,unbound))
    graph.discard(unbound,"material not retained by the RBP antibody is washed away")
    fastap=graph.input("fastap","FastAP alkaline phosphatase",MaterialKind.CATALYST)
    pnk=graph.input("t4_pnk","T4 PNK (pH 6.5 buffer)",MaterialKind.CATALYST)
    repaired=graph.product("repaired","on-bead RNA fragments with 3'-OH ends",
                           MaterialKind.MOLECULAR)
    graph.react("dephosphorylate","Dephosphorylate RNA fragment ends on beads",
                StepKind.DEPHOSPHORYLATION,
                ((enriched,InputRole.TEMPLATE),(fastap,InputRole.CATALYST),
                 (pnk,InputRole.CATALYST)),(repaired,))
    rna_adapter=graph.input("rna_adapter","X1A (+X1B pair) /5Phos/ /3SpC3/ barcoded RNA adapters",
                            MaterialKind.ADAPTER,state=rna_adapter_pair())
    rna_ligase=graph.input("rna_ligase","T4 RNA ligase 1 high concentration, PEG 8000, DMSO",
                           MaterialKind.CATALYST)
    adapter_rna=graph.product("adapter_rna","3-prime-adapter-ligated RBP-bound RNA",
                              MaterialKind.MOLECULAR,state=tuple(rna_adapter_ligation_rows()))
    free_adapter=graph.product("free_rna_adapter","unincorporated RNA adapter",
                               MaterialKind.WASTE)
    rna_geometry=reported_adapter_end(
        "X1A RNA adapter",kind="single-strand",polarity="5-prime",
        terminal_sequence=RNA_X1A_BARCODE+"N"*5+X1A_TAIL,
        statement="The ENCODE SOP prints the X1A barcode, five random bases and its 22-base tail at the RNA 3-prime boundary; RiL19 is a distinct 19-base adapter.")
    graph.react("ligate_rna_adapter","Ligate the barcoded RNA adapter on beads",
                StepKind.LIGATION,
                ((repaired,InputRole.TEMPLATE),(rna_adapter,InputRole.STRUCTURAL),
                 (rna_ligase,InputRole.CATALYST)),(adapter_rna,free_adapter),geometry=rna_geometry)
    graph.discard(free_adapter,"on-bead ligation lets the washes remove unincorporated adapter (paper Methods)")
    membrane=graph.input("membrane","SDS-PAGE, nitrocellulose transfer and band-to-+75 kDa excision",
                         MaterialKind.SELECTION)
    proteinase=graph.input("proteinase_k","Proteinase K",MaterialKind.CATALYST)
    released=graph.product("released_rna","size-selected adapter-ligated RNA with peptide remnant",
                           MaterialKind.MOLECULAR)
    off_window=graph.product("off_window","complexes outside the excised membrane region",
                             MaterialKind.WASTE)
    graph.react("membrane_select","Select the RBP band to 75 kDa above it and release RNA",
                StepKind.SELECTION,
                ((adapter_rna,InputRole.TEMPLATE),(membrane,InputRole.SELECTION),
                 (proteinase,InputRole.CATALYST)),(released,off_window))
    graph.discard(off_window,"membrane outside the excised lane region is not processed")
    rt_primer=graph.input("rt_primer","AR17 reverse-transcription primer",
                          MaterialKind.OLIGO,state=AR17)
    rt=graph.input("reverse_transcriptase","AffinityScript reverse transcriptase, 55 °C",MaterialKind.CATALYST)
    truncated=graph.product("truncated_cdna","crosslink-truncated cDNA",
                            MaterialKind.MOLECULAR)
    readthrough=graph.product("readthrough_cdna","crosslink-read-through cDNA",
                              MaterialKind.MOLECULAR)
    graph.react("reverse_transcribe","Reverse-transcribe through or up to the peptide remnant",
                StepKind.REVERSE_TRANSCRIPTION,
                ((released,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(truncated,readthrough))
    cdna_pool=graph.product("cdna_pool","pooled truncated and read-through cDNA",
                            MaterialKind.MOLECULAR)
    graph.react("pool_cdna","Retain both informative cDNA classes",StepKind.POOL,
                ((truncated,InputRole.TEMPLATE),(readthrough,InputRole.TEMPLATE)),
                (cdna_pool,),changes_molecule=False)
    exosap=graph.input("exosap","ExoSAP-IT, then NaOH 70 °C and HCl neutralisation",
                       MaterialKind.CATALYST)
    clean_cdna=graph.product("clean_cdna","single-stranded cDNA after primer and RNA removal",
                             MaterialKind.MOLECULAR)
    rt_waste=graph.product("rt_waste","excess AR17 primer and hydrolysed RNA template",
                           MaterialKind.WASTE)
    graph.react("cleanup_cdna","Remove excess primer, hydrolyse RNA, Silane cleanup",
                StepKind.CLEANUP,
                ((cdna_pool,InputRole.TEMPLATE),(exosap,InputRole.CATALYST)),
                (clean_cdna,rt_waste))
    graph.discard(rt_waste,"ExoSAP digests excess oligonucleotide; NaOH hydrolyses RNA; MyONE Silane cleanup")
    cdna_adapter=graph.input("cdna_adapter","5-prime-phosphorylated rand103Tr3 N10 adapter",
                             MaterialKind.ADAPTER,state=RAND103TR3)
    cdna_ligase=graph.input("cdna_ligase","T4 RNA ligase 1 high concentration, PEG 8000, DMSO, overnight",
                            MaterialKind.CATALYST)
    dual=graph.product("dual_adapter","dual-adapter eCLIP cDNA",MaterialKind.MOLECULAR,
                       state=tuple(cdna_adapter_ligation_rows()))
    cdna_geometry=reported_adapter_end(
        "rand103Tr3",kind="single-strand",polarity="5-prime",
        terminal_sequence="N"*10+RAND103TR3_FIXED,
        statement="The SOP prints the rand103Tr3 N10 adapter and specifies its ligation to single-stranded cDNA.")
    graph.react("ligate_cdna_adapter","Ligate rand103Tr3 to both cDNA classes",
                StepKind.LIGATION,
                ((clean_cdna,InputRole.TEMPLATE),(cdna_adapter,InputRole.STRUCTURAL),
                 (cdna_ligase,InputRole.CATALYST)),(dual,),geometry=cdna_geometry)
    d5=graph.input("d5_primer","D501 indexed PCR primer",MaterialKind.OLIGO,state=PCR_F_D501)
    d7=graph.input("d7_primer","D701 indexed PCR primer",MaterialKind.OLIGO,state=PCR_R_D701)
    polymerase=graph.input("pcr_polymerase","Q5 PCR master mix (cycles from qPCR Ct)",MaterialKind.CATALYST)
    amplified=graph.product("amplified","amplified eCLIP library",MaterialKind.MOLECULAR)
    graph.react("library_pcr","Amplify and add D5/D7 indexes",StepKind.PCR,
                ((dual,InputRole.TEMPLATE),(d5,InputRole.STRUCTURAL),
                 (d7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (amplified,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The exact final sites are modelled, but the compact module does not reconstruct the pre-PCR duplex orientation."))
    gel=graph.input("agarose_gel","agarose gel size selection",MaterialKind.SELECTION)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY)
    gel_waste=graph.product("gel_waste","products outside the excised size range",
                            MaterialKind.WASTE)
    graph.react("size_select","Size-select the amplified library on agarose",
                StepKind.SELECTION,
                ((amplified,InputRole.TEMPLATE),(gel,InputRole.SELECTION)),
                (library,gel_waste),changes_molecule=False)
    graph.discard(gel_waste,"gel region not excised; the excised range is not stated in the inspected Methods")
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS),
                                          read_lengths=READ_LENGTHS))
    graph.annotate(SourceDiscrepancy("membrane_select","Published RNA-size equivalent of the 75 kDa membrane window",
        (SourceClaim("Van Nostrand 2016 Online Methods (PMC4887338 JATS text line 123)","75 kDa (~150 nt of RNA) above the protein size"),
         SourceClaim("Van Nostrand 2016 Figure 1a legend (JATS text line 363)","75 kDa (~220 nt of RNA) above the protein size")),
        "Draw the window as the RBP band to 75 kDa above it, as both passages and the SOP state; no nucleotide length is assigned.",
        "Both values are printed in the defining article; Supplementary Protocol 1 states only the 75 kDa cut.",
        "Which nucleotide equivalent the authors intended is unresolved."))
    graph.annotate(SourceDiscrepancy("ligate_cdna_adapter","cDNA adapter name and step heading in Supplementary Protocol 1",
        (SourceClaim("Supplementary Protocol 1 oligo list (supplement-2 text lines 352–353)","rand103Tr3 /5Phos/N10AGATCGGAAGAGCACACGTCTG/3SpC3/"),
         SourceClaim("Supplementary Protocol 1 Day 3 step (text lines 1126–1128)","\u201c5\u2019 linker ligate cDNA\u201d; \u201cAdd 0.8 \u00b5l rand3Tr3 adapter\u201d")),
        "Use the listed rand103Tr3 oligo, ligated by its 5'-phosphate to the cDNA 3' end as the paper Methods state.",
        "The Methods say the second adapter is ligated to the cDNA fragment 3\u2032 end; rand3Tr3 is the only adapter named in the step and no separate rand3Tr3 sequence is listed.",
        "Identity of rand3Tr3 with rand103Tr3 is inferred from context, not stated."))
    graph.declare_reviewed_scope((
        "rbp_rna","uv","crosslink","rnase","trim","antibody_beads",
        "immunoprecipitate","fastap","t4_pnk","dephosphorylate","rna_adapter",
        "rna_ligase","ligate_rna_adapter","membrane","proteinase_k","membrane_select",
        "rt_primer","reverse_transcriptase","reverse_transcribe","pool_cdna",
        "exosap","cleanup_cdna",
        "cdna_adapter","cdna_ligase","ligate_cdna_adapter","d5_primer",
        "d7_primer","pcr_polymerase","library_pcr","agarose_gel","size_select",
        "library"),SOURCE)
    return graph.require_publishable()
