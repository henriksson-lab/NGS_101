"""Molecular construct model for HyDrop-ATAC (De Rop et al., 2022)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind)


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


# Source oligos.  One real member of each barcode/index plate is retained for the
# transcription checks; public molecular diagrams use placeholders for variable bases.
T7_PROMOTER = "TAATACGACTCACTATAGGG"
# The printed oligo has eight consecutive Ts: seven before the promoter's initial T.
ACRYDITE_PRIMER = "T" * 7 + T7_PROMOTER + rt.SMART_HANDLE + "AC"
PLATE1_FIRST = "GCAGTAGCTGTGTAGCAAGTGTACTCTGCG"
PLATE2_FIRST = "AGGGTACTCGTTAGTTGGACGCAGTAGCTG"
PLATE3_FIRST = "CCGAGCCCACGAGACTGACCGTACTAGGGTACTCG"

LINK1 = revcomp(PLATE1_FIRST[:10])  # CAGCTACTGC on the bead
LINK2 = revcomp(PLATE2_FIRST[:10])  # CGAGTACCCT on the bead
BEAD_HYI7_SITE = "GGAAGCAGTGGTATCAACGCAGAGTAC"  # final 27 nt retained in library
HYI7_SPACER = "CTGTCCGC"
I7_OLIGO_INDEX = "CGCTCAGTTC"
I5_OLIGO_INDEX = "TCGTGGAGCG"


def barcode(name: str) -> Segment:
    digit = next(c for c in reversed(name) if c.isdigit())
    symbol = {"1": "B", "2": "C", "3": "D"}[digit]
    return seg(name, symbol * 10, "cbc", placeholder=True,
               feature=feature(f"cell_bc{digit}", "cell_barcode", "combinatorial",
                               group="cell_id", part=f"round {digit}"))


def bead_oligo() -> Construct:
    """Final 117-nt bead oligo after three 96-way extension rounds.

    Linkers and s7 are derived from the reverse complements of the plate templates,
    preventing the extension direction from being represented inconsistently.
    """
    bead = Construct([
        seg("poly(T) before promoter", "T" * 7),
        seg("T7 promoter", T7_PROMOTER),
        seg("SMART handle", rt.SMART_HANDLE, "tso"),
        seg("AC", "AC"),
        barcode("barcode 1'"),
        seg("linker 1", LINK1, "r1"),
        barcode("barcode 2'"),
        seg("linker 2", LINK2, "r2"),
        barcode("barcode 3'"),
        seg("s7 capture", nx.S7, "s7"),
    ], name="HyDrop-ATAC bead oligo")
    if len(bead) != 117:
        raise ValueError(f"HyDrop bead oligo must be 117 nt, got {len(bead)}")
    return bead


def hyi7_primer() -> list[Segment]:
    return [seg("P7", il.P7, "p7"), seg("i7", I7_OLIGO_INDEX, "cbc"),
            seg("spacer", HYI7_SPACER), seg("bead-backbone anneal", BEAD_HYI7_SITE, "tso")]


def hyi5_primer() -> list[Segment]:
    return [seg("P5", il.P5, "p5"), seg("i5", I5_OLIGO_INDEX, "cbc"),
            seg("s5", nx.S5, "s5"), seg("ME prefix", nx.ME[:7], "me")]


def tagmented_scene() -> Scene:
    """One amplifiable s5/s7 Tn5 product before its two 9-nt gaps are filled."""
    gap = nx.TAGMENTATION_GAP
    core = seg("genomic core", "X" * 24, placeholder=True)
    top = [seg("s5", nx.S5, "s5"), seg("left ME", nx.ME, "me"),
           seg("left 9-nt gap", "X" * gap, placeholder=True), core]
    bottom = [seg("s7", nx.S7, "s7"), seg("right ME", nx.ME, "me"),
              seg("right 9-nt gap", "X" * gap, placeholder=True),
              *complement_segments([core])]
    sc = Scene()
    sc.strand("top", top, label="DNA")
    sc.anneal("bottom", bottom, to="top", pair=("genomic core'", "genomic core"),
              label="DNA")
    sc.anneal("left ME bottom", [seg("left ME'", nx.ME_RC, "me")], to="top",
              pair=("left ME'", "left ME"), label="", mod5="p")
    sc.anneal("right ME bottom", [seg("right ME'", nx.ME_RC, "me")], to="bottom",
              pair=("right ME'", "right ME"), label="", mod5="p", above=True)
    sc.mark("top", "left 9-nt gap", "gap on opposite strand")
    sc.mark("bottom", "right 9-nt gap", "gap on opposite strand")
    return sc


def gap_filled_template(insert_nt: int = 28) -> Construct:
    """s5-to-s7' strand available to the released bead primer."""
    return Construct([
        seg("s5", nx.S5, "s5"),
        seg("mosaic end", nx.ME, "me"),
        seg("accessible genomic DNA", "X" * insert_nt, placeholder=True),
        seg("opposite mosaic end", nx.ME_RC, "me"),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
    ], name="gap-filled HyDrop-ATAC fragment")


def capture_scene(insert_nt: int = 28) -> Scene:
    """Released bead oligo annealed to s7' after droplet gap fill."""
    template = gap_filled_template(insert_nt)
    sc = Scene()
    sc.strand("template", list(template), label="template")
    sc.anneal("bead primer", list(bead_oligo()), to="template",
              pair=("s7 capture", "s7 reverse complement"), label="bead primer")
    sc.arrow("bead primer", "linear extension copies the fragment")
    return sc


def bead_primed_product(insert_nt: int = 28) -> Construct:
    """One linear-extension product, 5' bead backbone to copied s5'."""
    return Construct([
        *bead_oligo(),
        seg("copied mosaic end", nx.ME, "me"),
        seg("copied accessible DNA", "X" * insert_nt, placeholder=True),
        seg("copied opposite ME", nx.ME_RC, "me"),
        seg("copied s5", nx.S5_RC, "s5"),
    ], name="cell-barcoded linear product")


SEQ_PRIMERS = tuple(sp.NEXTERA[k] for k in ("R1", "I1", "I2", "R2"))


def final_library(insert_nt: int = 28) -> Construct:
    """Final P5-to-P7' library after HYi5/HYi7 bulk PCR."""
    bead = bead_oligo()
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("i5 sample index", I5_OLIGO_INDEX, "cbc",
            feature=feature("sample_i5", "sample_index", "whitelist", whitelist="published HYi5 primer set")),
        seg("s5", nx.S5, "s5"),
        seg("mosaic end", nx.ME, "me"),
        seg("accessible genomic DNA", "X" * insert_nt, placeholder=True),
        seg("opposite mosaic end", nx.ME_RC, "me"),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
        barcode("barcode 3"),
        seg("linker 2, plate orientation", revcomp(bead.get("linker 2").top), "r2"),
        barcode("barcode 2"),
        seg("linker 1, plate orientation", revcomp(bead.get("linker 1").top), "r1"),
        barcode("barcode 1"),
        seg("bead-backbone tail", revcomp(BEAD_HYI7_SITE), "tso"),
        seg("HYi7 spacer reverse complement", revcomp(HYI7_SPACER)),
        seg("i7 reverse complement", revcomp(I7_OLIGO_INDEX), "cbc",
            feature=feature("sample_i7", "sample_index", "whitelist", whitelist="published HYi7 primer set")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="HyDrop-ATAC sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid final HyDrop-ATAC construct: " + "; ".join(problems))
    return lib


def primer_landings(lib: Construct | None = None) -> list[tuple[sp.SeqPrimer, sp.Landing]]:
    lib = lib or final_library()
    out = []
    for primer in SEQ_PRIMERS:
        hit = sp.locate(lib, primer)
        if hit is None:
            raise ValueError(f"{primer.name} does not land on {lib.name}")
        out.append((primer, hit))
    return out


def index1_layout() -> list[tuple[str, str]]:
    """Published 52-cycle cell-barcode read, derived from final-library segments."""
    lib = final_library()
    names = [
        ("barcode 3", "cell barcode 3"),
        ("linker 2, plate orientation", "linker 2"),
        ("barcode 2", "cell barcode 2"),
        ("linker 1, plate orientation", "linker 1"),
        ("barcode 1", "cell barcode 1"),
    ]
    rows, cycle = [], 1
    for name, label in names:
        end = cycle + len(lib.get(name)) - 1
        rows.append((f"{cycle}&ndash;{end}", label))
        cycle = end + 1
    rows.append((f"{cycle}&ndash;52", "bead-backbone tail"))
    return rows


def protocol_graph() -> ProtocolGraph:
    source="Defining source: https://doi.org/10.7554/eLife.73971"
    graph=ProtocolGraph("HyDrop-ATAC",authority=source)
    precursor=graph.input("bead_precursor","acrydite hydrogel bead primer",
                          MaterialKind.ADAPTER)
    barcode_plates=graph.input("barcode_plates","three 96-way barcode extension plates",
                               MaterialKind.OLIGO)
    bead_pol=graph.input("bead_polymerase","bead-extension polymerase",MaterialKind.CATALYST)
    bead=graph.product("barcoded_bead",bead_oligo().name,MaterialKind.ADAPTER,state=bead_oligo())
    graph.react("bead_synthesis","Build the three-round split-pool-barcoded bead oligo",
                StepKind.EXTENSION,
                ((precursor,InputRole.TEMPLATE),(barcode_plates,InputRole.STRUCTURAL),
                 (bead_pol,InputRole.CATALYST)),(bead,))
    nuclei=graph.input("nuclei","permeabilized nuclei",MaterialKind.BIOLOGICAL,primary=True)
    transposome=graph.input("transposome","s5/s7 Tn5 transposome",MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=graph.product("tagmented_chromatin","s5/s7-tagmented accessible chromatin",
                         MaterialKind.MOLECULAR,state=tagmented_scene())
    tag_waste=graph.product("tagment_waste","non-amplifiable tagmentation products",
                            MaterialKind.WASTE)
    graph.react("tagment","Tagment nuclei in bulk with s5 and s7 adaptors",
                StepKind.TAGMENTATION,
                ((nuclei,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,tag_waste))
    graph.discard(tag_waste,"s5/s5 and s7/s7 products are not the modelled amplifiable fragment")
    droplets=graph.input("droplet_device","HyDrop droplet generator",MaterialKind.CONDITION)
    productive=graph.product("productive_droplets","droplets containing tagmented nucleus and bead",
                              MaterialKind.MOLECULAR,state=tagmented_scene())
    empty=graph.product("unproductive_droplets","droplets without a usable nucleus/bead pair",
                        MaterialKind.WASTE)
    graph.react("partition","Co-encapsulate tagmented nuclei with barcoded beads",
                StepKind.PARTITION,
                ((tagged,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
                 (droplets,InputRole.CONDITION)),(productive,empty))
    graph.discard(empty,"unproductive droplets do not contribute indexed library molecules")
    dtt=graph.input("dtt","DTT hydrogel-solubilization condition",MaterialKind.CONDITION)
    released=graph.product("released_primer","released bead primers with tagmented template",
                           MaterialKind.MOLECULAR,state=bead_oligo())
    matrix=graph.product("bead_matrix","insoluble hydrogel bead remnants",MaterialKind.WASTE)
    graph.react("release_bead","Dissolve the hydrogel bead and release its primers",
                StepKind.SOLUBILIZATION,
                ((productive,InputRole.TEMPLATE),(dtt,InputRole.CONDITION)),(released,matrix))
    graph.discard(matrix,"the hydrogel matrix is not part of the library")
    gap_pol=graph.input("gap_fill_polymerase","gap-fill polymerase",MaterialKind.CATALYST)
    dntps=graph.input("gap_fill_dntps","gap-fill dNTPs",MaterialKind.OLIGO)
    filled=graph.product("gap_filled",gap_filled_template().name,MaterialKind.MOLECULAR,
                         state=gap_filled_template())
    graph.react("gap_fill","Fill the two 9-nt Tn5 gaps",StepKind.EXTENSION,
                ((released,InputRole.TEMPLATE),(dntps,InputRole.STRUCTURAL),
                 (gap_pol,InputRole.CATALYST)),(filled,))
    product=graph.product("linear_product",bead_primed_product().name,MaterialKind.MOLECULAR,
                          state=bead_primed_product())
    graph.react("linear_extension","Prime s7-prime and linearly copy fragments for 13 cycles",
                StepKind.EXTENSION,
                ((filled,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
                 (gap_pol,InputRole.CATALYST)),(product,))
    hyi5=graph.input("hyi5","HYi5 bulk-PCR primer",MaterialKind.OLIGO,
                     state=tuple(hyi5_primer()))
    hyi7=graph.input("hyi7","HYi7 bulk-PCR primer",MaterialKind.OLIGO,
                     state=tuple(hyi7_primer()))
    pcr_pol=graph.input("pcr_polymerase","bulk-PCR polymerase",MaterialKind.CATALYST)
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final)
    graph.react("library_pcr","Add flow-cell adapters and sample indexes with HYi5/HYi7",
                StepKind.PCR,
                ((product,InputRole.TEMPLATE),(hyi5,InputRole.STRUCTURAL),
                 (hyi7,InputRole.STRUCTURAL),(pcr_pol,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.exact(
                    bead_primed_product(),
                    PCRPrimer("HYi7","".join(s.top for s in hyi7_primer()),
                              anneal_nt=len(BEAD_HYI7_SITE)),
                    PCRPrimer("HYi5","".join(s.top for s in hyi5_primer()),
                              anneal_nt=len(nx.S5)+7),
                    statement="Supplementary file 4 prints HYi5/HYi7 and their bead-backbone and s5/ME-prefix landing regions."))
    lengths={"Read 1":50,"Index 1 (i7)":52,"Index 2 (i5)":10,"Read 2":50}
    graph.sequence(library,SequencingPlan(
        library=final,primers=SEQ_PRIMERS,acquired_roles=tuple(lengths),read_lengths=lengths,
        read_uses={"Read 1":"include","Index 1 (i7)":"extract",
                   "Index 2 (i5)":"demultiplex","Read 2":"include"}))
    graph.declare_reviewed_scope((
        "bead_precursor","barcode_plates","bead_polymerase","bead_synthesis","nuclei",
        "transposome","tn5","tagment","droplet_device","partition","dtt","release_bead",
        "gap_fill_polymerase","gap_fill_dntps","gap_fill","linear_extension","hyi5",
        "hyi7","pcr_polymerase","library_pcr","library"),source)
    return graph.require_publishable()
