"""Molecular construct model for txci-ATAC-seq."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def kit_boundary(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    """10x-kit sequence known here only through the secondary reconstruction."""
    return seg(name, top, tag, inferred=True, **kw)


TN5_BARCODE = "GAACCGCG"       # representative A1 well, Supplementary Table S5
I7_OLIGO_INDEX = "TCGCCTTA"    # P7.S701, Supplementary Table S7
SHORT_SBS = il.TRUSEQ_READ2[-18:]


def tn5_a() -> list[Segment]:
    return [seg("s5", nx.S5, "s5"), seg("mosaic end", nx.ME, "me")]


def tn5_b() -> list[Segment]:
    return [seg("partial TruSeq Read 2", SHORT_SBS, "r2"),
            seg("Tn5 barcode", TN5_BARCODE, "cbc"),
            seg("mosaic end", nx.ME, "me")]


def tn5_bottom() -> list[Segment]:
    return [seg("mosaic-end reverse complement", nx.ME_RC, "me")]


def transposome_scene(which: str = "A") -> Scene:
    transfer = tn5_a() if which == "A" else tn5_b()
    sc = Scene()
    sc.strand("transfer", transfer, label=f"Tn5ME-{which}")
    sc.anneal("bottom", tn5_bottom(), to="transfer",
              pair=("mosaic-end reverse complement", "mosaic end"),
              label="Tn5MErev", mod5="p")
    return sc


def bead_oligo() -> list[Segment]:
    return [kit_boundary("P5", il.P5, "p5"),
            kit_boundary("GEM barcode", "N" * 16, "cbc", placeholder=True,
                         feature=feature("cell_gem", "cell_barcode", "whitelist", whitelist="10x ATAC whitelist", group="cell_id", part="GEM")),
            kit_boundary("s5", nx.S5, "s5")]


def p7_primer() -> list[Segment]:
    return [seg("P7", il.P7, "p7"), seg("i7 oligo index", I7_OLIGO_INDEX, "cbc"),
            seg("TruSeq Read 2", il.TRUSEQ_READ2, "r2")]


SEQ_PRIMERS = (sp.NEXTERA["R1"], sp.TRUSEQ["I1"],
               sp.NEXTERA["I2"], sp.TRUSEQ["R2"])


def final_library(insert_nt: int = 32) -> Construct:
    lib = Construct([
        *bead_oligo(), seg("mosaic end", nx.ME, "me"),
        seg("accessible genomic DNA", "X" * insert_nt, placeholder=True),
        seg("opposite mosaic end", nx.ME_RC, "me"),
        seg("Tn5 barcode, read orientation", revcomp(TN5_BARCODE), "cbc",
            feature=feature("cell_tn5", "cell_barcode", "combinatorial", group="cell_id", part="Tn5")),
        seg("TruSeq Read 2 reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
        seg("i7, read orientation", revcomp(I7_OLIGO_INDEX), "cbc",
            feature=feature("sample_i7", "sample_index", "whitelist", whitelist="published P7 index set")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="txci-ATAC sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid final txci-ATAC construct: " + "; ".join(problems))
    return lib


def primer_landings(lib: Construct | None = None):
    lib = lib or final_library()
    out = []
    for primer in SEQ_PRIMERS:
        hit = sp.locate(lib, primer)
        if hit is None:
            raise ValueError(f"{primer.name} does not land on {lib.name}")
        out.append((primer, hit))
    return out


def read2_layout() -> list[tuple[str, str]]:
    return [("1&ndash;8", "Tn5 barcode"), ("9&ndash;27", "mosaic end"),
            ("28&ndash;", "accessible genomic DNA")]


def protocol_graph() -> ProtocolGraph:
    source="Defining source: https://doi.org/10.1186/s13059-023-03150-1"
    kit_source=source+" (stock 10x gel-bead oligo follows the marked secondary reconstruction)"
    graph=ProtocolGraph("txci-ATAC-seq",authority=source)
    nuclei=graph.input("nuclei","isolated nuclei",MaterialKind.BIOLOGICAL,primary=True)
    well_a=graph.product("tn5_well_a","representative indexed-Tn5 well A nuclei",MaterialKind.BIOLOGICAL)
    well_b=graph.product("tn5_well_b","representative indexed-Tn5 well B nuclei",MaterialKind.BIOLOGICAL)
    graph.react("tn5_split","Split nuclei across 96 indexed-Tn5 wells",StepKind.SPLIT,
                ((nuclei,InputRole.TEMPLATE),),(well_a,well_b))
    transposomes=graph.input("indexed_transposomes","Tn5ME-A/Tn5ME-B indexed transposome set",
                             MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged_a=graph.product("tagged_a","pre-indexed nuclei from well A",MaterialKind.MOLECULAR)
    tagged_b=graph.product("tagged_b","pre-indexed nuclei from well B",MaterialKind.MOLECULAR)
    waste_a=graph.product("tagment_waste_a","non-library products from well A",MaterialKind.WASTE)
    waste_b=graph.product("tagment_waste_b","non-library products from well B",MaterialKind.WASTE)
    graph.react("tagment_a","Tagment accessible chromatin in indexed well A",StepKind.TAGMENTATION,
                ((well_a,InputRole.TEMPLATE),(transposomes,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged_a,waste_a))
    graph.react("tagment_b","Tagment accessible chromatin in indexed well B",StepKind.TAGMENTATION,
                ((well_b,InputRole.TEMPLATE),(transposomes,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged_b,waste_b))
    graph.discard(waste_a,"fragments lacking productive A/B ends do not amplify")
    graph.discard(waste_b,"fragments lacking productive A/B ends do not amplify")
    pooled=graph.product("pooled_nuclei","pooled Tn5-pre-indexed nuclei",MaterialKind.MOLECULAR)
    graph.react("pool","Pool the indexed-Tn5 wells",StepKind.POOL,
                ((tagged_a,InputRole.TEMPLATE),(tagged_b,InputRole.TEMPLATE)),(pooled,))
    loading=graph.input("gem_loading","overloaded Chromium ATAC GEM condition",MaterialKind.CONDITION)
    gems=graph.product("productive_gems","GEMs containing pre-indexed nuclei and gel beads",
                       MaterialKind.MOLECULAR)
    excluded=graph.product("excluded_gems","empty or multiplet GEMs",MaterialKind.WASTE)
    graph.react("partition","Partition pooled nuclei into Chromium GEMs",StepKind.PARTITION,
                ((pooled,InputRole.TEMPLATE),(loading,InputRole.CONDITION)),(gems,excluded))
    graph.discard(excluded,"empty and multiplet GEMs are excluded from accepted cells")
    bead=graph.input("bead_oligo","P5-GEM-barcode-s5 gel-bead oligo",MaterialKind.OLIGO,
                     state=tuple(bead_oligo()),inferred=True,source=kit_source)
    short=graph.input("short_sbs","Short SBS B-end primer",MaterialKind.OLIGO)
    polymerase=graph.input("polymerase","GEM and library PCR polymerase",MaterialKind.CATALYST)
    barcoded=graph.product("gem_barcoded","GEM-barcoded A/B tagmentation products",
                           MaterialKind.MOLECULAR,inferred=True,source=kit_source)
    graph.react("gem_pcr","Copy the GEM barcode and amplify A/B fragments in droplets",
                StepKind.PCR,
                ((gems,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
                 (short,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(barcoded,),
                geometry=PCRGeometry.unavailable(final_library(),
                    reason="The stock 10x bead-primer sequence is not published in the defining paper."))
    p7=graph.input("p7_primer","published indexed P7.S701 primer",MaterialKind.OLIGO,
                   state=tuple(p7_primer()))
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final,
                          inferred=True,source=kit_source)
    graph.react("library_pcr","Add the lane index and complete flow-cell ends",StepKind.PCR,
                ((barcoded,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
                 (p7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.unavailable(final,
                    reason="The final inner custom regions are published; the P5/GEM end retains the explicit 10x-kit boundary."))
    graph.sequence(library,SequencingPlan(
        library=final,primers=SEQ_PRIMERS,acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
        read_lengths={"Read 1":51,"Index 1 (i7)":8,"Index 2 (i5)":16,"Read 2":78},
        read_uses={"Read 1":"include","Index 1 (i7)":"demultiplex",
                   "Index 2 (i5)":"extract","Read 2":"extract_then_include"},
        limitation="The stock 10x gel-bead oligo is a marked secondary reconstruction."))
    graph.declare_reviewed_scope((
        "nuclei","tn5_split","indexed_transposomes","tn5","tagment_a","tagment_b",
        "pool","gem_loading","partition","bead_oligo","short_sbs","polymerase",
        "gem_pcr","p7_primer","library_pcr","library"),source)
    return graph.require_publishable()
