"""Chen et al. single-cell TAPS."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "taps__10.1038+s41587-019-0041-2" / "tools"))
from taps import taps_path
from chemdraw import Construct, MolecularState, Scene, Workflow, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "scTAPS — single-cell 5mC + 5hmC sequencing"
NOTES = "01_sctaps.html"
SOURCE = 'Defining preprint: <a href="https://doi.org/10.1101/2024.02.12.579864">Chen et al. (2024)</a>; authoritative version: <a href="https://doi.org/10.1186/s13059-025-03708-1">Genome Biology (2025)</a>.'
SUMMARY = "Each cell is fragmented with a uniquely barcoded Tn5 transposome. Ninety-six gap-filled libraries are pooled, then TET oxidation and pyridine–borane reduction make both 5mC and 5hmC read as T."
CAVEAT = "Additional file 5 prints the 96 transposome combinations and four custom sequencing primers. This diagram retains a generic Nextera completion shell; its primer sites are not a reconstruction of the custom run."
CONVERTS = ("5mC", "5hmC")
TAPS_PATHS = {base: taps_path(base) for base in ("5mC", "5hmC", "C")}

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"),
    seg("Tn5 i5 cell barcode", "J" * 8, "cbc", placeholder=True,
        feature=feature("cell_barcode_i5", "cell_barcode", "combinatorial", group="cell_barcode", part="i5", whitelist="Additional file 5 Table S3")),
    seg("S5", nx.S5, "s5"), seg("left mosaic end", nx.ME, "me"),
    seg("TAPS-converted genomic insert", "X" * 400, placeholder=True, length_bp=400),
    seg("right mosaic end reverse complement", nx.ME_RC, "me"), seg("S7 reverse complement", nx.S7_RC, "s7"),
    seg("Tn5 i7 cell barcode reverse complement", "I" * 8, "cbc", placeholder=True,
        feature=feature("cell_barcode_i7", "cell_barcode", "combinatorial", group="cell_barcode", part="i7", whitelist="Additional file 5 Table S3")),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="scTAPS Nextera library")
SEQ_PRIMERS = (sp.NEXTERA["R1"], sp.NEXTERA["I1"], sp.NEXTERA["I2"], sp.NEXTERA["R2"])
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS):
    raise ValueError("scTAPS final library: " + "; ".join(problems))
READ_LENGTHS = {"Read 1":120, "Index 1 (i7)":8, "Index 2 (i5)":8, "Read 2":120}
SHOW_READ_CONTENT = True
FINAL_CAPTION = "Custom i5/i7 transposomes encode the cell before 96 wells are pooled. Paired 120-base reads interrogate a genomic insert in which 5mC and 5hmC are read as T."
SEQUENCING_INTRO = "The published run uses four custom NovaSeq primers (Additional file 5, Table S4); the primer binding diagram below shows only a generic Nextera proxy. Two eight-cycle index reads recover the combinatorial Tn5 cell barcode."

_initial = Scene.duplex([seg("single-cell genomic DNA", "D" * 52, placeholder=True)], label="genomic DNA from one FACS-sorted cell")
INITIAL_ROWS=tuple(_initial.rows()); INITIAL_NAME="Single-cell genomic DNA"

def _base_state(title: str, stage: int) -> Scene:
    states = {base: path.after(stage) for base, path in TAPS_PATHS.items()}
    sc=_tagged_scene("cell-barcoded genomic fragment retained through conversion")
    sc.strand("5mC", [seg(f"5mC path: {states['5mC']}", "M", "w1", placeholder=True)], label=f"5mC: {states['5mC']}")
    sc.strand("5hmC", [seg(f"5hmC path: {states['5hmC']}", "H", "w1", placeholder=True)], label=f"5hmC: {states['5hmC']}")
    sc.strand("C", [seg(f"C path: {states['C']}", "C", placeholder=True)], label=f"unmodified C: {states['C']}")
    sc.footer(title, "C"); return sc

def _tagged_scene(label: str) -> Scene:
    sc=Scene.duplex([
        seg("i5 cell barcode", "J"*8, "cbc", placeholder=True,
            feature=feature("cell_barcode_i5", "cell_barcode", "combinatorial",
                            group="cell_barcode", part="i5", whitelist="Additional file 5 Table S3")),
        seg("left Tn5 end", "X"*12, "me", placeholder=True),
        seg("genomic fragment", "D"*40, placeholder=True),
        seg("right Tn5 end", "X"*12, "me", placeholder=True),
        seg("i7 cell barcode", "I"*8, "cbc", placeholder=True,
            feature=feature("cell_barcode_i7", "cell_barcode", "combinatorial",
                            group="cell_barcode", part="i7", whitelist="Additional file 5 Table S3")),
    ],label=label)
    sc.mark("top","i5 cell barcode","well-specific i5/i7 transposome")
    return sc

def workflow():
    tagged=_tagged_scene("barcoded, gap-filled single-cell fragment")
    pooled=_tagged_scene("pooled DNA from 96 uniquely indexed wells plus controls")
    oxidized=_base_state("two mTET1 oxidation rounds", 1)
    reduced=_base_state("pyridine–borane reduction", 2)
    final=Scene.duplex(list(FINAL_LIBRARY), label="PCR readout: modified cytosines become T")
    w=Workflow(MolecularState(INITIAL_NAME, INITIAL_ROWS))
    w.react("Tagment with a well-specific Tn5 and gap-fill", tagged.rows(), name="Cell-barcoded fragments", note="Tn5 acts for 15 min at 50 °C; SDS strips Tn5 and 72 °C extension fills the nine-base gaps.")
    w.react("Pool 96 wells", pooled.rows(), name="Pooled cell-barcoded DNA", note="Filler DNA and modification controls are added before conversion chemistry.")
    w.react("Oxidize 5mC and 5hmC with mTET1", oxidized.rows(), name="5caC-containing DNA", note="Two oxidation rounds convert both modified cytosines to 5caC; unmodified C remains C.")
    w.react("Reduce 5caC with pyridine–borane", reduced.rows(), name="DHU-containing DNA", note="5caC becomes dihydrouracil (DHU).")
    w.react("PCR with Nextera XT index primers", final.rows(), name="Sequencing library", note="Uracil-tolerant PCR copies DHU as T and preserves unmodified C.")
    return w


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    dna=graph.input("single_cell_dna","genomic DNA from one sorted cell",
                    MaterialKind.BIOLOGICAL,state=INITIAL_ROWS,primary=True)
    transposome=graph.input("indexed_transposome","well-specific i5/i7 Tn5 transposome",
                            MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=graph.product("tagged_dna","cell-indexed tagmented DNA",MaterialKind.MOLECULAR)
    transposition_waste=graph.product("transposition_waste","non-library Tn5 products",
                                      MaterialKind.WASTE)
    graph.react("tagment","Tagment with a well-specific transposome",StepKind.TAGMENTATION,
                ((dna,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,transposition_waste))
    graph.discard(transposition_waste,"unusable and non-amplifiable strand-transfer products")
    polymerase=graph.input("gap_fill_polymerase","gap-fill DNA polymerase",MaterialKind.CATALYST)
    gap_filled=graph.product("gap_filled","gap-filled cell-indexed DNA",MaterialKind.MOLECULAR,
                             state=tuple(_tagged_scene("barcoded, gap-filled fragment").rows()))
    graph.react("gap_fill","Fill the nine-base transposition gaps",StepKind.EXTENSION,
                ((tagged,InputRole.TEMPLATE),(polymerase,InputRole.CATALYST)),(gap_filled,))
    controls=graph.input("conversion_controls","filler DNA and modification controls",
                         MaterialKind.MOLECULAR)
    pooled=graph.product("pooled_dna","pooled indexed wells plus conversion controls",
                         MaterialKind.MOLECULAR)
    graph.react("pool_wells","Pool the indexed wells with conversion controls",StepKind.POOL,
                ((gap_filled,InputRole.TEMPLATE),(controls,InputRole.TEMPLATE)),(pooled,),
                changes_molecule=False)
    mtet1=graph.input("mtet1","mTET1 oxidation enzyme",MaterialKind.CATALYST)
    oxidized=graph.product("oxidized_dna","5caC-containing indexed DNA",MaterialKind.MOLECULAR,
                           state=tuple(_base_state("two mTET1 oxidation rounds",1).rows()))
    graph.react("oxidize","Oxidize 5mC and 5hmC to 5caC",StepKind.CHEMICAL_CONVERSION,
                ((pooled,InputRole.TEMPLATE),(mtet1,InputRole.CATALYST)),(oxidized,))
    borane=graph.input("pyridine_borane","pyridine–borane reducing condition",
                       MaterialKind.CONDITION)
    reduced=graph.product("reduced_dna","DHU-containing indexed DNA",MaterialKind.MOLECULAR,
                          state=tuple(_base_state("pyridine–borane reduction",2).rows()))
    graph.react("reduce","Reduce 5caC to dihydrouracil",StepKind.CHEMICAL_CONVERSION,
                ((oxidized,InputRole.TEMPLATE),(borane,InputRole.CONDITION)),(reduced,))
    i5=graph.input("i5_primer","Nextera XT i5 PCR primer",MaterialKind.OLIGO)
    i7=graph.input("i7_primer","Nextera XT i7 PCR primer",MaterialKind.OLIGO)
    pcr_pol=graph.input("pcr_polymerase","uracil-tolerant PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    graph.react("library_pcr","Copy DHU as T and complete the indexed library",StepKind.PCR,
                ((reduced,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
                 (i7,InputRole.STRUCTURAL),(pcr_pol,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The page references the 96-member index set rather than one exact primer pair."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS),
                                          read_lengths=READ_LENGTHS))
    graph.declare_reviewed_scope((
        "single_cell_dna","indexed_transposome","tn5","tagment","gap_fill_polymerase",
        "gap_fill","conversion_controls","pool_wells","mtet1","oxidize",
        "pyridine_borane","reduce","i5_primer","i7_primer","pcr_polymerase",
        "library_pcr","library"),SOURCE)
    return graph.require_publishable()
