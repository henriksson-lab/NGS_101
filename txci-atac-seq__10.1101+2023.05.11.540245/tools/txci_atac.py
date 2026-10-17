"""Molecular construct model for txci-ATAC-seq."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (AnnotationSlot, Highlight, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)


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


def tn5_b_template() -> list[Segment]:
    """Tn5ME-B exactly as Table S4 prints it: the 8 N are the per-well Tn5 barcode."""
    return [seg("partial TruSeq Read 2", SHORT_SBS, "r2"),
            seg("Tn5 barcode (Table S5, one per well)", "N" * 8, "cbc", placeholder=True,
                feature=feature("cell_tn5", "cell_barcode", "combinatorial",
                                group="cell_id", part="Tn5")),
            seg("mosaic end", nx.ME, "me")]


def short_sbs() -> list[Segment]:
    return [seg("Short SBS (partial TruSeq Read 2)", SHORT_SBS, "r2")]


# Published Genome Biology Methods, "txci-ATAC-seq using human lung, mouse lung, and
# mouse liver tissue samples" (the standard protocol), and the authors' protocols.io
# PDF (dx.doi.org/10.17504/protocols.io.dm6gp3o68vzp/v1) steps 1-7 and 50-52.
SEQUENCING_NOTE = (
    "The authors' protocols.io run sheet (step 52) gives Read 1 51, i7 8, i5 16 and Read 2 78 "
    "cycles. The published lung/liver run acquired 10 i7 cycles only to accommodate a co-pooled "
    "bulk ATAC library (the txci library has an 8-bp i7). Read 2 cycles 9-27 read the fixed "
    "mosaic end: the authors spiked 5% bulk ATAC library and recommend PhiX or dark cycles 9-27 "
    "when the library is sequenced alone.")


def protocol_graph() -> ProtocolGraph:
    source="Defining source: https://doi.org/10.1186/s13059-023-03150-1"
    bench=source+"; authors' protocols.io dx.doi.org/10.17504/protocols.io.dm6gp3o68vzp/v1"
    kit_source=source+" (stock 10x gel-bead oligo follows the marked secondary reconstruction)"
    graph=ProtocolGraph("txci-ATAC-seq",authority=source)
    tn5_a_oligo=graph.input("tn5me_a","Tn5ME-A (Table S4)",MaterialKind.OLIGO,state=tuple(tn5_a()))
    tn5_b_oligo=graph.input("tn5me_b","Tn5ME-B, 96 barcoded versions (Tables S4/S5)",
                            MaterialKind.OLIGO,state=tuple(tn5_b_template()))
    tn5_rev=graph.input("tn5me_rev","[phos]Tn5MErev (Table S4)",MaterialKind.OLIGO,
                        state=tuple(tn5_bottom()))
    linkers=graph.product("annealed_linkers","annealed linker A and barcoded linker B duplexes",
                          MaterialKind.ADAPTER)
    graph.react("anneal_linkers","Anneal Tn5ME-A and each Tn5ME-B to Tn5MErev "
                "(95 °C 5 min, -0.1 °C/s to 65 °C, 65 °C 5 min, -0.1 °C/s to 4 °C)",
                StepKind.HYBRIDIZATION,
                ((tn5_a_oligo,InputRole.TEMPLATE),(tn5_b_oligo,InputRole.TEMPLATE),
                 (tn5_rev,InputRole.STRUCTURAL)),(linkers,),source=bench)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    transposomes=graph.product("indexed_transposomes",
                               "iTSM plate: 5 µL 500 nM A/B-loaded Tn5, one Tn5ME-B barcode per well",
                               MaterialKind.ADAPTER)
    graph.react("load_tn5","Load linker A plus one barcoded linker B onto Tn5 in each well "
                "(protocols.io: 23 °C 30 min; the published iTSM plate was pre-loaded by Illumina)",
                StepKind.GENERIC,((linkers,InputRole.TEMPLATE),(tn5,InputRole.CATALYST)),
                (transposomes,),source=bench,changes_molecule=False,
                reason="Non-covalent assembly of the annealed linkers with Tn5 protein.")
    nuclei=graph.input("nuclei","isolated nuclei",MaterialKind.BIOLOGICAL,primary=True)
    well_a=graph.product("tn5_well_a","representative indexed-Tn5 well A nuclei",MaterialKind.BIOLOGICAL)
    well_b=graph.product("tn5_well_b","representative indexed-Tn5 well B nuclei",MaterialKind.BIOLOGICAL)
    graph.react("tn5_split","Split 20,000 nuclei (7 µL PBSB + 13 µL TBS) into each indexed-Tn5 well",
                StepKind.SPLIT,((nuclei,InputRole.TEMPLATE),),(well_a,well_b))
    tagged_a=graph.product("tagged_a","pre-indexed nuclei from well A",MaterialKind.MOLECULAR)
    tagged_b=graph.product("tagged_b","pre-indexed nuclei from well B",MaterialKind.MOLECULAR)
    waste_a=graph.product("tagment_waste_a","A/A and B/B fragments from well A",MaterialKind.WASTE)
    waste_b=graph.product("tagment_waste_b","A/A and B/B fragments from well B",MaterialKind.WASTE)
    tag_label=("Tagment accessible chromatin in indexed well {} (37 °C 1 h, lid 47 °C; "
               "5 min on ice, no stop buffer)")
    graph.react("tagment_a",tag_label.format("A"),StepKind.TAGMENTATION,
                ((well_a,InputRole.TEMPLATE),(transposomes,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged_a,waste_a))
    graph.react("tagment_b",tag_label.format("B"),StepKind.TAGMENTATION,
                ((well_b,InputRole.TEMPLATE),(transposomes,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged_b,waste_b))
    same_end=("A/A and B/B fragments stay in the nuclei and later reactions but lack the "
              "bead-primed A end plus SBS-primed B end; no physical separation is represented")
    graph.discard(waste_a,same_end)
    graph.discard(waste_b,same_end)
    pooled=graph.product("pooled_nuclei","pooled Tn5-pre-indexed nuclei (TMG washes, 40 µm filter)",
                         MaterialKind.MOLECULAR)
    graph.react("pool","Pool the indexed-Tn5 wells",StepKind.POOL,
                ((tagged_a,InputRole.TEMPLATE),(tagged_b,InputRole.TEMPLATE)),(pooled,))
    loading=graph.input("gem_loading",
                        "overloaded Chromium ATAC GEM condition (100,000-200,000 nuclei in 15 µL "
                        "loading buffer with 5 µM short SBS)",MaterialKind.CONDITION)
    gems=graph.product("productive_gems","GEMs containing pre-indexed nuclei and gel beads",
                       MaterialKind.MOLECULAR)
    excluded=graph.product("excluded_gems","empty or multiplet GEMs",MaterialKind.WASTE)
    graph.react("partition","Partition pooled nuclei into Chromium GEMs (CG000209 Rev D step 2)",
                StepKind.PARTITION,
                ((pooled,InputRole.TEMPLATE),(loading,InputRole.CONDITION)),(gems,excluded))
    graph.discard(excluded,"computational exclusion only: empty-GEM and multiplet barcodes are "
                  "filtered in analysis; their DNA is not physically separated")
    bead=graph.input("bead_oligo","P5-GEM-barcode-s5 gel-bead oligo",MaterialKind.OLIGO,
                     state=tuple(bead_oligo()),inferred=True,source=kit_source)
    short=graph.input("short_sbs","Short SBS B-end primer (Table S8)",MaterialKind.OLIGO,
                      state=tuple(short_sbs()))
    polymerase=graph.input("polymerase","GEM and library PCR polymerase",MaterialKind.CATALYST)
    barcoded=graph.product("gem_barcoded","GEM-barcoded A/B tagmentation products",
                           MaterialKind.MOLECULAR,inferred=True,source=kit_source)
    graph.react("gem_pcr","Copy the GEM barcode and amplify A/B fragments in droplets "
                "(72 °C 5 min; 98 °C 30 s; 8 × 98 °C 10 s / 59 °C 30 s / 72 °C 1 min)",
                StepKind.PCR,
                ((gems,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
                 (short,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(barcoded,),
                geometry=PCRGeometry.unavailable(final_library(),
                    reason="The stock 10x bead-primer sequence is not published in the defining paper."),
                source=bench)
    cleaned=graph.product("gem_cleaned","recovered GEM-barcoded products",MaterialKind.MOLECULAR,
                          inferred=True,source=kit_source)
    gem_waste=graph.product("gem_cleanup_waste","emulsion reagents and unrecovered material",
                            MaterialKind.WASTE)
    graph.react("gem_cleanup","Post GEM Incubation Cleanup (CG000209 Rev D step 3)",StepKind.CLEANUP,
                ((barcoded,InputRole.TEMPLATE),),(cleaned,gem_waste))
    graph.discard(gem_waste,"removed by the kit's post-GEM cleanup")
    p7=graph.input("p7_primer","published indexed P7.S701 primer (2.5 µL 25 µM)",MaterialKind.OLIGO,
                   state=tuple(p7_primer()))
    final=final_library()
    amplified=graph.product("amplified_library","sample-indexed library before size selection",
                            MaterialKind.MOLECULAR,inferred=True,source=kit_source)
    graph.react("library_pcr","Add the lane index and complete flow-cell ends "
                "(98 °C 45 s; 5 × 98 °C 20 s / 67 °C 30 s / 72 °C 20 s; 72 °C 1 min)",StepKind.PCR,
                ((cleaned,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
                 (p7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(amplified,),
                geometry=PCRGeometry.unavailable(final,
                    reason="The final inner custom regions are published; the P5/GEM end retains the explicit 10x-kit boundary."),
                source=bench)
    spri=graph.input("size_selection","kit Double Sided Size Selection (CG000209 Rev D step 4.2)",
                     MaterialKind.SELECTION)
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final,
                          inferred=True,source=kit_source)
    size_waste=graph.product("size_excluded","fragments outside the double-sided size window",
                             MaterialKind.WASTE)
    graph.react("size_select","Double-sided size selection",StepKind.SELECTION,
                ((amplified,InputRole.TEMPLATE),(spri,InputRole.SELECTION)),(library,size_waste))
    graph.discard(size_waste,"removed by the double-sided bead size selection")
    graph.sequence(library,SequencingPlan(
        library=final,primers=SEQ_PRIMERS,acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
        read_lengths={"Read 1":51,"Index 1 (i7)":8,"Index 2 (i5)":16,"Read 2":78},
        read_uses={"Read 1":"include","Index 1 (i7)":"demultiplex",
                   "Index 2 (i5)":"extract","Read 2":"extract_then_include"},
        limitation="The stock 10x gel-bead oligo is a marked secondary reconstruction."))
    graph.annotate(Highlight("library","Acquired cycles and base diversity",SEQUENCING_NOTE,
                             slot=AnnotationSlot.SEQUENCING))
    graph.declare_reviewed_scope((
        "tn5me_a","tn5me_b","tn5me_rev","anneal_linkers","tn5","load_tn5",
        "nuclei","tn5_split","tagment_a","tagment_b",
        "pool","gem_loading","partition","bead_oligo","short_sbs","polymerase",
        "gem_pcr","gem_cleanup","p7_primer","library_pcr","size_selection","size_select",
        "library"),source)
    return graph.require_publishable()
