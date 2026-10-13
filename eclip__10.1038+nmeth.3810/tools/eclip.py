"""Original paired-end eCLIP molecular workflow and exact ENCODE SOP adapters."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from adapter_ligation import reported_adapter_end
from batch_ngs import seg
from chemdraw import Construct, Scene, complement_segments, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from rna_special import adapter_ligation_scene

TITLE = "eCLIP — enhanced crosslinking and immunoprecipitation"
NOTES = "01_eclip.html"
SOURCE = 'Defining paper and protocol: <a href="https://doi.org/10.1038/nmeth.3810">Van Nostrand et al. (2016)</a>; <a href="https://www.encodeproject.org/documents/842f7424-5396-424a-a1a3-3f18707c3222/@@download/attachment/eCLIP_SOP_v1.P_110915.pdf">ENCODE eCLIP SOP v1.P</a>.'
SUMMARY = "UV-crosslinked RBP–RNA complexes are immunoprecipitated and RNase-trimmed. A 3′ RNA adapter precedes reverse transcription; a second adapter is ligated to cDNA so both read-through and crosslink-truncated cDNAs can amplify."
CAVEAT = "The SOP provides several color-balanced RNA barcode and D5/D7 index choices. The final panel uses the published X1A, D501 and D701 examples; their positions, widths and primer family do not vary."
I5_FEATURE=feature("sample_index_i5","sample_index","fixed")
I7_FEATURE=feature("sample_index_i7","sample_index","fixed")
INLINE_FEATURE=feature("inline_rna_barcode","inline_barcode","fixed")
RANDOMER_5=feature("molecular_barcode_rna_adapter","umi","random",group="molecular_barcode",part="RNA adapter")
RANDOMER_10=feature("molecular_barcode_cdna_adapter","umi","random",group="molecular_barcode",part="cDNA adapter")

RIL19 = "AGATCGGAAGAGCGTCGTGTG"
RNA_X1A_BARCODE = "ATATAGG"
AR17 = "ACACGACGCTCTTCCGA"
RAND103TR3_FIXED = "AGATCGGAAGAGCACACGTCTG"
D501_I5 = "TATAGCCT"
D701_I7_OLIGO = "CGAGTAAT"

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"),
    seg("D501 i5", D501_I5, "cbc", feature=I5_FEATURE),
    seg("Read 1 arm", il.TRUSEQ_READ1, "r1"),
    seg("X1A inline RNA barcode", RNA_X1A_BARCODE, "cbc", feature=INLINE_FEATURE),
    seg("X1A five-base randomer", "U" * 5, "umi", placeholder=True, feature=RANDOMER_5),
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

FINAL_CAPTION = "Original paired-end eCLIP amplicon using the SOP's X1A/D501/D701 examples. Read 1 starts with the seven-base inline barcode and five-base randomer; Read 2 starts with the ten-base cDNA-adapter randomer."
SEQUENCING_INTRO = "Original eCLIP uses paired-end reads and both D5/i5 and D7/i7 index reads. Read 1 begins at the RNA-adapter barcode; Read 2 begins at rand103Tr3's ten-base randomer and approaches the RNA 5′ end/crosslink-stop boundary."


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


def reverse_transcription_rows():
    rna = [
        seg("RNA 5-prime flank", "X" * 14, placeholder=True),
        seg("peptide-remnant crosslink", "X", "umi", placeholder=True),
        seg("RNA 3-prime fragment", "X" * 20, placeholder=True),
        seg("3-prime RNA adapter", "X" * 18, "r2", placeholder=True),
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
    def product(name, cdna_name, *, read_through=False):
        parts = [
            seg("rand103Tr3 ten-base randomer", "N" * 10, "umi", placeholder=True,
                feature=RANDOMER_10),
            seg("rand103Tr3 fixed region", RAND103TR3_FIXED, "r2"),
            seg(cdna_name, "X" * 28, placeholder=True),
        ]
        if read_through:
            parts.append(seg("RNA-adapter complement", "X" * 18, "r2", placeholder=True))
        sc = Scene()
        sc.strand(name, parts, mod5="p")
        sc.junction(name, "rand103Tr3 fixed region", cdna_name, "cDNA ligation")
        return sc.rows()

    return [
        *product("crosslink-truncated product", "cDNA ending at crosslink"),
        *product("read-through product", "full-length cDNA", read_through=True),
    ]


def sections():
    return [
        ("Crosslink, trim and immunoprecipitate", crosslinked_rna_rows(), "UV fixes direct contacts; limited RNase and RBP immunoprecipitation define the captured RNA fragment."),
        ("Ligate the 3′ RNA adapter on beads", adapter_ligation_scene(fragment_name="RBP-bound RNA fragment").rows(), "The indexed RNA adapter is ligated before gel purification and reverse transcription."),
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
    rnase=graph.input("rnase","limited RNase",MaterialKind.CATALYST)
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
    rna_adapter=graph.input("rna_adapter","X1A barcoded randomer RNA adapter",
                            MaterialKind.ADAPTER)
    rna_ligase=graph.input("rna_ligase","RNA 3-prime adapter ligase",
                           MaterialKind.CATALYST)
    adapter_rna=graph.product("adapter_rna","3-prime-adapter-ligated RBP-bound RNA",
                              MaterialKind.MOLECULAR)
    rna_geometry=reported_adapter_end(
        "X1A RNA adapter",kind="single-strand",polarity="5-prime",
        terminal_sequence=RNA_X1A_BARCODE+"N"*5+RIL19,
        statement="The ENCODE SOP prints the X1A barcode/randomer and RIL19 adapter sequence used at the RNA 3-prime boundary.")
    graph.react("ligate_rna_adapter","Ligate the barcoded RNA adapter on beads",
                StepKind.LIGATION,
                ((enriched,InputRole.TEMPLATE),(rna_adapter,InputRole.STRUCTURAL),
                 (rna_ligase,InputRole.CATALYST)),(adapter_rna,),geometry=rna_geometry)
    rt_primer=graph.input("rt_primer","AR17 reverse-transcription primer",
                          MaterialKind.OLIGO)
    rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
    truncated=graph.product("truncated_cdna","crosslink-truncated cDNA",
                            MaterialKind.MOLECULAR)
    readthrough=graph.product("readthrough_cdna","crosslink-read-through cDNA",
                              MaterialKind.MOLECULAR)
    graph.react("reverse_transcribe","Reverse-transcribe through or up to the peptide remnant",
                StepKind.REVERSE_TRANSCRIPTION,
                ((adapter_rna,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(truncated,readthrough))
    cdna_pool=graph.product("cdna_pool","pooled truncated and read-through cDNA",
                            MaterialKind.MOLECULAR)
    graph.react("pool_cdna","Retain both informative cDNA classes",StepKind.POOL,
                ((truncated,InputRole.TEMPLATE),(readthrough,InputRole.TEMPLATE)),
                (cdna_pool,),changes_molecule=False)
    cdna_adapter=graph.input("cdna_adapter","5-prime-phosphorylated rand103Tr3 N10 adapter",
                             MaterialKind.ADAPTER)
    cdna_ligase=graph.input("cdna_ligase","single-stranded cDNA ligase",
                            MaterialKind.CATALYST)
    dual=graph.product("dual_adapter","dual-adapter eCLIP cDNA",MaterialKind.MOLECULAR,
                       state=tuple(cdna_adapter_ligation_rows()))
    cdna_geometry=reported_adapter_end(
        "rand103Tr3",kind="single-strand",polarity="5-prime",
        terminal_sequence="N"*10+RAND103TR3_FIXED,
        statement="The SOP prints the rand103Tr3 N10 adapter and specifies its ligation to single-stranded cDNA.")
    graph.react("ligate_cdna_adapter","Ligate rand103Tr3 to both cDNA classes",
                StepKind.LIGATION,
                ((cdna_pool,InputRole.TEMPLATE),(cdna_adapter,InputRole.STRUCTURAL),
                 (cdna_ligase,InputRole.CATALYST)),(dual,),geometry=cdna_geometry)
    d5=graph.input("d5_primer","D501 indexed PCR primer",MaterialKind.OLIGO)
    d7=graph.input("d7_primer","D701 indexed PCR primer",MaterialKind.OLIGO)
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY)
    graph.react("library_pcr","Amplify and add D5/D7 indexes",StepKind.PCR,
                ((dual,InputRole.TEMPLATE),(d5,InputRole.STRUCTURAL),
                 (d7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The exact final sites are modelled, but the compact module does not reconstruct the pre-PCR duplex orientation."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "rbp_rna","uv","crosslink","rnase","trim","antibody_beads",
        "immunoprecipitate","rna_adapter","rna_ligase","ligate_rna_adapter",
        "rt_primer","reverse_transcriptase","reverse_transcribe","pool_cdna",
        "cdna_adapter","cdna_ligase","ligate_cdna_adapter","d5_primer",
        "d7_primer","pcr_polymerase","library_pcr","library"),SOURCE)
    return graph.require_publishable()
