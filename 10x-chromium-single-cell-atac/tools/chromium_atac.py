"""Molecular construct model for 10x Chromium Single Cell ATAC v1."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments, feature
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal)
from protocol_model import molecule, role
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
from protocol_families.tenx_graph import partition_material

CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="10x-chromium-atac-v1")
SAMPLE_INDEX = feature("sample_index_i7", "sample_index", "fixed")
FEATURES = {"10x cell barcode": CELL_BARCODE, "i7 sample index read": SAMPLE_INDEX}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


# CG000168 Rev A Appendix prints the i7 primer and final library with an 8-N
# "Sample Index N"; each Set A well is a mix of four 8-nt words (10x Set A CSV), so the
# model keeps a placeholder rather than one concrete index.
SAMPLE_INDEX_NT = 8
# CG000168 Rev A Appendix, Step 4.1: SI-PCR Primer B PN-2000128 ("Partial P5") and the
# i7 Sample Index Plate N, Set A PN-3000262 primer P7-N8-partial Read 2N.
SI_PCR_PRIMER_B = il.P5[:22]


def i7_primer_segments() -> list[Segment]:
    return [seg("P7", il.P7, "p7"),
            seg("i7 sample index", "N" * SAMPLE_INDEX_NT, "cbc", placeholder=True,
                feature=SAMPLE_INDEX),
            seg("partial Read 2N", nx.S7, "s7")]


def bead_oligo() -> list[Segment]:
    return [seg("P5", il.P5, "p5"),
            seg("10x cell barcode", "B" * 16, "cbc", placeholder=True),
            seg("s5", nx.S5, "s5")]


def transposed_scene(insert_nt: int = 28) -> Scene:
    gap = nx.TAGMENTATION_GAP
    core = seg("accessible genomic DNA", "X" * insert_nt, placeholder=True)
    top = [seg("s5", nx.S5, "s5"), seg("left ME", nx.ME, "me"),
           seg("left 9-nt gap", "X" * gap, placeholder=True), core]
    bottom = [seg("s7", nx.S7, "s7"), seg("right ME", nx.ME, "me"),
              seg("right 9-nt gap", "X" * gap, placeholder=True),
              *complement_segments([core])]
    sc = Scene()
    sc.strand("top", top, label="DNA")
    sc.anneal("bottom", bottom, to="top", pair=("accessible genomic DNA'", "accessible genomic DNA"), label="DNA")
    sc.anneal("left ME bottom", [seg("left ME'", nx.ME_RC, "me")], to="top",
              pair=("left ME'", "left ME"), label="", mod5="p")
    sc.anneal("right ME bottom", [seg("right ME'", nx.ME_RC, "me")], to="bottom",
              pair=("right ME'", "right ME"), label="", mod5="p", above=True)
    return sc


def gap_filled_template(insert_nt: int = 28) -> Construct:
    return Construct([seg("s7", nx.S7, "s7"), seg("mosaic end", nx.ME, "me"),
                      seg("accessible genomic DNA", "X" * insert_nt, placeholder=True),
                      seg("opposite ME reverse complement", nx.ME_RC, "me"),
                      seg("s5 reverse complement", nx.S5_RC, "s5")],
                     name="gap-filled 10x ATAC fragment")


def capture_scene(insert_nt: int = 28) -> Scene:
    template = gap_filled_template(insert_nt)
    sc = Scene()
    sc.strand("template", list(template), label="template")
    sc.anneal("bead primer", bead_oligo(), to="template", pair=("s5", "s5 reverse complement"),
              label="released gel-bead oligo")
    sc.arrow("bead primer", "linear extension through the accessible fragment")
    return sc


SEQ_PRIMERS = tuple(sp.NEXTERA[k] for k in ("R1", "I1", "I2", "R2"))
READ_LENGTHS = {"Read 1": 50, "Index 1 (i7)": 8,
                "Index 2 (i5)": 16, "Read 2": 50}


def final_library(insert_nt: int = 28) -> Construct:
    lib = Construct([
        *bead_oligo(), seg("mosaic end", nx.ME, "me"),
        seg("accessible genomic DNA", "X" * insert_nt, placeholder=True),
        seg("opposite ME reverse complement", nx.ME_RC, "me"),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 sample index read", "N" * SAMPLE_INDEX_NT, "cbc", placeholder=True),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="10x Chromium Single Cell ATAC library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid 10x ATAC library: " + "; ".join(problems))
    return lib


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")


def protocol_graph():
    authority=('Authoritative source: 10x Genomics <a href="https://www.10xgenomics.com/support/single-cell-atac/documentation/steps/library-prep/chromium-single-cell-atac-reagent-kits-user-guide-v-1-chemistry">'
               'CG000168 Rev A user guide</a>.')
    g=ProtocolGraph("10x Chromium Single Cell ATAC v1",authority=authority)
    nuclei=g.input("nuclei","intact nuclei containing accessible chromatin",
                   MaterialKind.BIOLOGICAL,
                   state=molecule("nuclei",(role("accessible chromatin"),)),primary=True)
    transposome=g.input("transposome","s5/s7 mosaic-end-loaded transposome",
                        MaterialKind.ADAPTER)
    tn5=g.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=g.product("tagged","tagmented accessible chromatin",MaterialKind.MOLECULAR,
                     state=transposed_scene())
    off=g.product("off_products","candidate s5/s5 and s7/s7 tagmentation products",
                  MaterialKind.WASTE,inferred=True,
                  source=authority+" (candidate end classes derived from the two adaptors; not described in CG000168)")
    g.react("tagment","Tagment accessible chromatin in intact nuclei",
            StepKind.TAGMENTATION,((nuclei,InputRole.TEMPLATE),
            (transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,off))
    g.mark_terminal(Terminal(off.id,Disposition.UNRESOLVED,
        "Candidate same-adaptor products lack an s5/s7 end pair; they stay in the nuclei "
        "through partitioning and are not physically removed at tagmentation. Their "
        "recovery through SILANE/SPRIselect cleanups is not assigned."))
    reviewed=["nuclei","transposome","tn5","tagment"]
    occupied=partition_material(g,tagged,reviewed,label="tagmented nuclei",
        unresolved_reason=("CG000168 Rev A Step 2 overview: most GEMs contain no nuclei, and after "
            "incubation the GEMs are broken and pooled; no physical removal of empty or "
            "multiplet GEM products is described. Cell Ranger ATAC assigns barcodes "
            "computationally (Sequencing section); the products' fate is not modelled further."))
    fill_pol=g.input("fill_polymerase","Barcoding Enzyme (gap-fill activity inferred)",
                     MaterialKind.CATALYST)
    filled=g.product("filled","gap-filled accessible-chromatin fragments",
                     MaterialKind.MOLECULAR,state=gap_filled_template(),inferred=True,
                     source=authority+" (GEM incubation starts 72 C 5 min; its gap-fill purpose is inferred, not stated)")
    g.react("gap_fill","INFERRED: fill both nine-nucleotide Tn5 gaps at 72 C before denaturation",StepKind.EXTENSION,
            ((occupied,InputRole.TEMPLATE),(fill_pol,InputRole.CATALYST)),(filled,))
    bead=g.input("bead_primer","P5/cell-barcode/s5 gel-bead primer",MaterialKind.OLIGO)
    extension_pol=g.input("extension_polymerase","Barcoding Enzyme (12-cycle linear amplification)",
                          MaterialKind.CATALYST)
    barcoded=g.product("barcoded","cell-barcoded accessible-DNA fragments",
                       MaterialKind.MOLECULAR,state=capture_scene())
    g.react("barcode","Extend the gel-bead primer through each s5 fragment",
            StepKind.EXTENSION,((filled,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
            (extension_pol,InputRole.CATALYST)),(barcoded,))
    reviewed.extend(("fill_polymerase","gap_fill","bead_primer",
                     "extension_polymerase","barcode"))
    lib_state=final_library()
    primer_b=g.input("index_pcr_forward","SI-PCR Primer B PN-2000128 (partial P5)",
                     MaterialKind.OLIGO,state=Construct([seg("partial P5",SI_PCR_PRIMER_B,"p5")]))
    i7_primer=g.input("index_pcr_reverse","i7 Sample Index Plate N, Set A PN-3000262 primer",
                      MaterialKind.OLIGO,state=Construct(i7_primer_segments()))
    amp=g.input("index_pcr_polymerase","Amp Mix",MaterialKind.CATALYST)
    library=g.product("index_pcr_product","Sample Index PCR product",MaterialKind.LIBRARY,
                      state=lib_state)
    geometry=PCRGeometry.exact(lib_state,PCRPrimer("SI-PCR Primer B",SI_PCR_PRIMER_B),
        PCRPrimer("i7 Sample Index N primer",il.P7+"N"*SAMPLE_INDEX_NT+nx.S7,anneal_nt=len(nx.S7)),
        statement="CG000168 Rev A Appendix Step 4.1 prints both Sample Index PCR primers; "
                  "the i7 primer anneals over partial Read 2N and Primer B over partial P5.")
    g.react("index_pcr","Complete the s7 end and add the i7 sample index",StepKind.PCR,
            ((barcoded,InputRole.TEMPLATE),(primer_b,InputRole.STRUCTURAL),
             (i7_primer,InputRole.STRUCTURAL),(amp,InputRole.CATALYST)),(library,),
            geometry=geometry)
    reviewed.extend(("index_pcr_forward","index_pcr_reverse","index_pcr_polymerase","index_pcr"))
    g.sequence(library,SequencingPlan(library=final_library(),primers=SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS),read_lengths=READ_LENGTHS))
    reviewed.append(library.id)
    g.declare_reviewed_scope(reviewed,authority)
    return g.require_publishable()
