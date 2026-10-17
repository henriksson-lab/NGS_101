"""Molecular model for CH-ATAC-seq from Supplementary Table 1."""
from __future__ import annotations

import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct, Segment, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


MGI_P5 = "GAACGACATGGCTACGATCCGACTT"
MGI_P7 = "TGTGAGCCAAGGAGTTGTTGTCTTC"
HY_BARCODE = "TTCTCGCATG"       # representative CH-ATAC Barcoded_HY_oligo_1
TN5_BARCODE = "TCCTACCAGT"      # representative Tn5_barcode_primer_1
I7_OLIGO_INDEX = "TAGGTCCGAT"   # representative MGI_P7_index_1


def indexed_tn5_transfer() -> list[Segment]:
    """Tn5_barcode_primer_1 as Supplementary Table 1 prints it: 52 bases.

    The table writes it `A*AGCAGTGGTATCAACGCAGAGT...` with the Modification column
    reading `*:Phosphorthioate`. The `*` marks the phosphorothioate **bond** between
    the first and second base, in the same notation the table uses elsewhere; the
    base string therefore begins `AAGCAGTGG...` and is exactly `rt.SMART_HANDLE`.
    (Confirmed inside the same table: `Barcoded_HY_oligo_*` opens with
    `ACTCTGCGTTGATACCACTGCTT`, the 23-nt reverse complement of that handle.)
    """
    return [seg("SMART handle", rt.SMART_HANDLE, "tso"),
            seg("Tn5 barcode", TN5_BARCODE, "cbc"),
            seg("mosaic end", nx.ME, "me")]


def primer_c() -> list[Segment]:
    return [seg("s7", nx.S7, "s7"), seg("mosaic end", nx.ME, "me")]


def hy_head() -> list[Segment]:
    return [seg("s5", nx.S5, "s5"), seg("mosaic end", nx.ME, "me")]


def barcoded_hy() -> list[Segment]:
    return [seg("SMART-handle reverse complement", revcomp(rt.SMART_HANDLE), "tso"),
            seg("HY barcode", HY_BARCODE, "cbc"),
            seg("HY-head reverse complement", revcomp(nx.ADAPTOR_S5), "s5")]


SEQ_PRIMERS = (sp.NEXTERA["R1"], sp.NEXTERA["I1"], sp.NEXTERA["R2"])


def hy_tn5_junction() -> str:
    """The HY/Tn5 junction, computed as the unique overlap of the two printed oligos.

    `revcomp(Barcoded_HY_oligo)` reads `s5 · ME · HY-barcode' · SMART handle` and the
    indexed Tn5 transfer strand reads `SMART handle · Tn5 barcode · ME`. Their maximal
    overlap is the shared handle, so the order of the P5-side segments in the filled-in
    library follows from the oligo table alone and is not an inference. Raises if the
    overlap ever stops being unique or stops being the handle.
    """
    hy_top = revcomp("".join(s.top for s in barcoded_hy()))
    transfer = "".join(s.top for s in indexed_tn5_transfer())
    overlaps = [n for n in range(1, min(len(hy_top), len(transfer)) + 1)
                if hy_top[-n:] == transfer[:n]]
    if overlaps != [len(rt.SMART_HANDLE)]:
        raise ValueError("CH-ATAC HY/Tn5 oligos no longer share exactly one overlap")
    return rt.SMART_HANDLE


def final_library(insert_nt: int = 32) -> Construct:
    lib = Construct([
        seg("MGI P5", MGI_P5, "p5"), seg("s5", nx.S5, "s5"),
        seg("mosaic end", nx.ME, "me"),
        seg("HY barcode, read orientation", revcomp(HY_BARCODE), "cbc",
            feature=feature("cell_hy", "cell_barcode", "combinatorial", group="cell_id", part="HY")),
        seg("SMART handle shared by the HY and Tn5 oligos", hy_tn5_junction(), "tso"),
        seg("Tn5 barcode", TN5_BARCODE, "cbc",
            feature=feature("cell_tn5", "cell_barcode", "combinatorial", group="cell_id", part="Tn5")),
        seg("mosaic end after barcodes", nx.ME, "me"),
        seg("accessible genomic DNA", "X" * insert_nt, placeholder=True),
        seg("opposite mosaic end", nx.ME_RC, "me"),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
        seg("i7, read orientation", revcomp(I7_OLIGO_INDEX), "cbc",
            feature=feature("sample_i7", "sample_index", "whitelist", whitelist="published MGI P7 index set")),
        seg("MGI P7 reverse complement", revcomp(MGI_P7), "p7"),
    ], name="CH-ATAC sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS,
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems:
        raise ValueError("invalid CH-ATAC library: " + "; ".join(problems))
    return lib


def primer_landings(lib: Construct | None = None):
    lib = lib or final_library()
    return [(p, sp.locate(lib, p)) for p in SEQ_PRIMERS]


def read1_layout() -> list[tuple[str, str]]:
    return [("1&ndash;10", "HY barcode"), ("11&ndash;33", "SMART handle — dark"),
            ("34&ndash;43", "Tn5 barcode"), ("44&ndash;62", "mosaic end — dark"),
            ("63&ndash;", "accessible genomic DNA")]


METHODS = ("STAR Methods, 'Library construction of CH-ATAC-seq' "
           "(doi:10.1016/j.devcel.2024.01.015, p. e4)")

# STAR Methods, 'MGI library preparation and sequencing' (p. e5) prints the read 1 and
# read 2 cycle counts and the dark windows, but neither the index-read length nor the
# bases of the "custom TM (Tn5 modified) sequencing primers".
READ_RECIPE = ("The defining STAR Methods print the DNBSEQ-T7 recipe: 100 cycles of read 1 "
               "with a dark reaction from cycles 11-33 and 44-62, and 100 cycles of read 2, "
               "after circularization with a VAHTS Circularization Kit for MGI and DNB "
               "preparation. Those dark windows match the modelled read 1 (HY barcode 1-10, "
               "shared handle 11-33, Tn5 barcode 34-43, mosaic end 44-62). The custom TM "
               "sequencing primer bases are not printed: the Nextera primers drawn here are "
               "placed by computation on the modelled library, and the 10-cycle index read "
               "follows from the 10-nt i7 rather than from a printed cycle count.")

WASTE_REASON = ("not physically separated in the model: same-ended products are computed "
                "non-library fates (handle/handle ends acquire s5 at both ends and are "
                "expected to be PCR-suppressed; Primer C/Primer C ends lack the P5-side "
                "site); their carry-over is not quantified in any inspected source")


def protocol_graph() -> ProtocolGraph:
    source="Defining source: https://doi.org/10.1016/j.devcel.2024.01.015"
    graph=ProtocolGraph("CH-ATAC-seq",authority=source)
    nuclei=graph.input("nuclei","isolated nuclei",MaterialKind.BIOLOGICAL,primary=True)
    t_a=graph.product("tn5_well_a","representative Tn5-index well A nuclei",MaterialKind.BIOLOGICAL)
    t_b=graph.product("tn5_well_b","representative Tn5-index well B nuclei",MaterialKind.BIOLOGICAL)
    graph.react("tn5_split","Split nuclei across 384 indexed-Tn5 wells",StepKind.SPLIT,
                ((nuclei,InputRole.TEMPLATE),),(t_a,t_b))
    transposomes=graph.input("indexed_transposomes","indexed-handle and Primer-C transposomes",
                             MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged_a=graph.product("tagged_a","Tn5-indexed nuclei from well A",MaterialKind.MOLECULAR)
    tagged_b=graph.product("tagged_b","Tn5-indexed nuclei from well B",MaterialKind.MOLECULAR)
    waste_a=graph.product("tagment_waste_a","non-library products from well A",MaterialKind.WASTE)
    waste_b=graph.product("tagment_waste_b","non-library products from well B",MaterialKind.WASTE)
    graph.react("tagment_a","Tagment with indexed-handle and Primer-C transposomes in well A",
                StepKind.TAGMENTATION,
                ((t_a,InputRole.TEMPLATE),(transposomes,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged_a,waste_a))
    graph.react("tagment_b","Tagment with indexed-handle and Primer-C transposomes in well B",
                StepKind.TAGMENTATION,
                ((t_b,InputRole.TEMPLATE),(transposomes,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged_b,waste_b))
    graph.discard(waste_a,WASTE_REASON)
    graph.discard(waste_b,WASTE_REASON)
    pool1=graph.product("tn5_pool","pooled Tn5-indexed nuclei",MaterialKind.MOLECULAR)
    graph.react("tn5_pooling","Pool the indexed-Tn5 wells",StepKind.POOL,
                ((tagged_a,InputRole.TEMPLATE),(tagged_b,InputRole.TEMPLATE)),(pool1,))
    h_a=graph.product("hy_well_a","representative HY-index well A nuclei",MaterialKind.MOLECULAR)
    h_b=graph.product("hy_well_b","representative HY-index well B nuclei",MaterialKind.MOLECULAR)
    graph.react("hy_split","Redistribute nuclei across 768 HY-barcode wells",StepKind.SPLIT,
                ((pool1,InputRole.TEMPLATE),),(h_a,h_b))
    hy=graph.input("hy_duplex","pre-annealed barcoded HY capture duplex",MaterialKind.ADAPTER)
    captured_a=graph.product("hy_captured_a","HY-barcode-captured fragments from well A",
                             MaterialKind.MOLECULAR)
    captured_b=graph.product("hy_captured_b","HY-barcode-captured fragments from well B",
                             MaterialKind.MOLECULAR)
    graph.react("hybridize_a","Hybridize HY duplex to the SMART handle in well A",
                StepKind.HYBRIDIZATION,
                ((h_a,InputRole.TEMPLATE),(hy,InputRole.STRUCTURAL)),(captured_a,))
    graph.react("hybridize_b","Hybridize HY duplex to the SMART handle in well B",
                StepKind.HYBRIDIZATION,
                ((h_b,InputRole.TEMPLATE),(hy,InputRole.STRUCTURAL)),(captured_b,))
    # Supplementary Table 1, common primers: Block_tail_primer_oligo, "Block free
    # annealing hybridization primers". Its bases equal the SMART handle, so it can
    # occupy the free HY duplexes' handle-complementary overhang. STAR Methods
    # "Library construction of CH-ATAC-seq" (mmc3 p. e4): hybridization 37 C 90 min,
    # then 0.5 ul of 100 uM block tail primer oligo per well, 30 min at 37 C, then pool.
    block=graph.input("block_tail_oligo","Block_tail_primer_oligo (SMART-handle blocker)",
                      MaterialKind.OLIGO,state=(seg("SMART handle",rt.SMART_HANDLE,"tso"),))
    blocked_a=graph.product("hy_blocked_a","HY-captured nuclei from well A with free HY overhangs blocked",
                            MaterialKind.MOLECULAR)
    blocked_b=graph.product("hy_blocked_b","HY-captured nuclei from well B with free HY overhangs blocked",
                            MaterialKind.MOLECULAR)
    graph.react("block_a","Block free HY duplexes with Block_tail_primer_oligo in well A",
                StepKind.HYBRIDIZATION,
                ((captured_a,InputRole.TEMPLATE),(block,InputRole.STRUCTURAL)),(blocked_a,))
    graph.react("block_b","Block free HY duplexes with Block_tail_primer_oligo in well B",
                StepKind.HYBRIDIZATION,
                ((captured_b,InputRole.TEMPLATE),(block,InputRole.STRUCTURAL)),(blocked_b,))
    pool2=graph.product("hy_pool","pooled Tn5/HY dual-barcoded nuclei",MaterialKind.MOLECULAR)
    graph.react("hy_pooling","Pool the HY-index wells",StepKind.POOL,
                ((blocked_a,InputRole.TEMPLATE),(blocked_b,InputRole.TEMPLATE)),(pool2,))
    # STAR Methods (mmc3 p. e4): the pooled, washed nuclei get a T4 PNK + ATP mix,
    # 37 C 30 min, before they are re-split for the final indexed PCR.
    pnk=graph.input("t4_pnk","T4 polynucleotide kinase with ATP",MaterialKind.CATALYST)
    phosphorylated=graph.product("pnk_pool","pooled dual-barcoded nuclei after T4 PNK",
                                 MaterialKind.MOLECULAR)
    graph.react("pnk","Phosphorylate 5' ends with T4 PNK in the pool",StepKind.PHOSPHORYLATION,
                ((pool2,InputRole.TEMPLATE),(pnk,InputRole.CATALYST)),(phosphorylated,),
                source=METHODS)
    pcr_a=graph.product("pcr_well_a","representative final-index PCR well A input",MaterialKind.MOLECULAR)
    pcr_b=graph.product("pcr_well_b","representative final-index PCR well B input",MaterialKind.MOLECULAR)
    graph.react("pcr_split","Divide dual-barcoded material across final indexed-PCR wells",
                StepKind.SPLIT,((phosphorylated,InputRole.TEMPLATE),),(pcr_a,pcr_b))
    # STAR Methods (mmc3 p. e4): in each final well, a gap-fill mix of Ampligase DNA
    # ligase with its buffer, dNTP and T4 DNA polymerase, 37 C 60 min. Placement and
    # reagents are printed; which nicks the ligase seals is derived from the oligo
    # layout (see the note), so the product carries no exact state.
    t4_pol=graph.input("t4_dna_polymerase","T4 DNA polymerase with dNTP",MaterialKind.CATALYST)
    ampligase=graph.input("ampligase","Ampligase thermostable DNA ligase",MaterialKind.CATALYST)
    filled_a=graph.product("hy_filled_a","gap-filled and sealed fragments in PCR well A",
                           MaterialKind.MOLECULAR)
    filled_b=graph.product("hy_filled_b","gap-filled and sealed fragments in PCR well B",
                           MaterialKind.MOLECULAR)
    for suffix,template,filled in (("a",pcr_a,filled_a),("b",pcr_b,filled_b)):
        graph.react("hy_fill_"+suffix,
                    "Gap-fill with T4 DNA polymerase and seal with Ampligase in PCR well "+suffix.upper(),
                    StepKind.EXTENSION,
                    ((template,InputRole.TEMPLATE),(hy,InputRole.STRUCTURAL),
                     (t4_pol,InputRole.CATALYST),(ampligase,InputRole.CATALYST)),(filled,),
                    source=METHODS)
    # STAR Methods (mmc3 p. e4): proteinase K/SDS nuclei lysis 55 C 15 min, Tween
    # quench, then PCR mix in the same well.
    lysis=graph.input("nuclei_lysis","proteinase K/SDS nuclei lysis buffer",MaterialKind.CATALYST)
    released={}
    for suffix,filled in (("a",filled_a),("b",filled_b)):
        released[suffix]=graph.product("released_"+suffix,"released fragments in PCR well "+suffix.upper(),
                                       MaterialKind.MOLECULAR)
        graph.react("lyse_"+suffix,"Release fragments with proteinase K/SDS in PCR well "+suffix.upper(),
                    StepKind.DEPROTEINIZATION,
                    ((filled,InputRole.TEMPLATE),(lysis,InputRole.CATALYST)),(released[suffix],),
                    source=METHODS,changes_molecule=False)
    # Supplementary Table 1, "PCR P5 primer" (5'-phosphorylated) and MGI_P7_index_1.
    p5=graph.input("mgi_p5_primer","MGI P5-side PCR primer",MaterialKind.OLIGO,
                   state=(seg("MGI P5",MGI_P5,"p5"),seg("s5",nx.S5,"s5")))
    p7=graph.input("indexed_mgi_p7_primer","indexed MGI P7-side PCR primer",
                   MaterialKind.OLIGO,
                   state=(seg("MGI P7",MGI_P7,"p7"),seg("i7",I7_OLIGO_INDEX,"cbc"),
                          seg("s7",nx.S7,"s7")))
    kapa=graph.input("kapa_hifi","KAPA HiFi HotStart ReadyMix",MaterialKind.CATALYST)
    final=final_library()
    for suffix,template in (("a",released["a"]),("b",released["b"])):
        library=graph.product("library_"+suffix,final.name+" sublibrary "+suffix.upper(),
                              MaterialKind.LIBRARY,state=final)
        graph.react("library_pcr_"+suffix,"Add indexed MGI ends in a 12-cycle PCR in well "+suffix.upper(),
                    StepKind.PCR,
                    ((template,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
                     (p7,InputRole.STRUCTURAL),(kapa,InputRole.CATALYST)),(library,),
                    source=METHODS,
                    geometry=PCRGeometry.exact(final,
                        PCRPrimer("MGI P5 primer",MGI_P5+nx.S5,anneal_nt=len(nx.S5)),
                        PCRPrimer("indexed MGI P7 primer",MGI_P7+I7_OLIGO_INDEX+nx.S7,
                                  anneal_nt=len(nx.S7)),
                        statement="Supplementary Table 1 prints the PCR P5 primer and the 96 MGI_P7 indexed primers; each lands on its Nextera arm with its 3' s5 or s7 segment."))
        graph.sequence(library,SequencingPlan(
            library=final,primers=SEQ_PRIMERS,
            acquired_roles=("Read 1","Index 1 (i7)","Read 2"),
            read_lengths={"Read 1":100,"Index 1 (i7)":10,"Read 2":100},
            read_uses={"Read 1":"extract_then_include","Index 1 (i7)":"demultiplex",
                       "Read 2":"include"},
            limitation=READ_RECIPE))
    graph.declare_reviewed_scope((
        "nuclei","tn5_split","indexed_transposomes","tn5","tagment_a","tagment_b",
        "tn5_pooling","hy_split","hy_duplex","hybridize_a","hybridize_b","block_tail_oligo",
        "block_a","block_b","hy_pooling","t4_pnk","pnk","pcr_split","t4_dna_polymerase","ampligase",
        "hy_fill_a","hy_fill_b","nuclei_lysis","lyse_a","lyse_b","kapa_hifi","mgi_p5_primer",
        "indexed_mgi_p7_primer","library_pcr_a","library_pcr_b","library_a","library_b"),source)
    return graph.require_publishable()
