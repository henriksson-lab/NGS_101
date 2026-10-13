"""Chen et al. single-cell CAPS+ selective 5hmC library."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from base_conversion import caps_plus_path
from chemdraw import Construct, MolecularState, Scene, Workflow, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)

TITLE = "scCAPS+ — single-cell 5hmC sequencing"
NOTES = "01_sccaps-plus.html"
SOURCE = 'Defining preprint: <a href="https://doi.org/10.1101/2024.02.12.579864">Chen et al. (2024)</a>; authoritative version: <a href="https://doi.org/10.1186/s13059-025-03708-1">Genome Biology (2025)</a>.'
SUMMARY = "Barcoded Tn5 fragments each single-cell genome before pooling. ACT+BF4− and Pinnick oxidation selectively take 5hmC through 5fC to 5caC; pyridine–borane then makes only 5hmC read as T."
CAVEAT = "The 96 published transposome/index oligos are referenced rather than copied here; their two eight-cycle cell-index positions are represented explicitly."
CONVERTS = ("5hmC",)
CAPS_PATHS = {base: caps_plus_path(base) for base in ("5hmC", "5mC", "C")}

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"),
    seg("Tn5 i5 cell barcode", "J"*8, "cbc", placeholder=True,
        feature=feature("cell_barcode_i5", "cell_barcode", "combinatorial", group="cell_barcode", part="i5", whitelist="Additional file 5 Table S3")),
    seg("S5", nx.S5, "s5"), seg("left mosaic end", nx.ME, "me"),
    seg("CAPS+-converted genomic insert", "X"*400, placeholder=True, length_bp=400),
    seg("right mosaic end reverse complement", nx.ME_RC, "me"), seg("S7 reverse complement", nx.S7_RC, "s7"),
    seg("Tn5 i7 cell barcode reverse complement", "I"*8, "cbc", placeholder=True,
        feature=feature("cell_barcode_i7", "cell_barcode", "combinatorial", group="cell_barcode", part="i7", whitelist="Additional file 5 Table S3")),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="scCAPS+ Nextera library")
SEQ_PRIMERS=(sp.NEXTERA["R1"],sp.NEXTERA["I1"],sp.NEXTERA["I2"],sp.NEXTERA["R2"])
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS): raise ValueError("scCAPS+ final library: "+"; ".join(problems))
READ_LENGTHS={"Read 1":120,"Index 1 (i7)":8,"Index 2 (i5)":8,"Read 2":120}
SHOW_READ_CONTENT=True
FINAL_CAPTION="Custom i5/i7 transposomes encode the cell before pooling. Paired reads interrogate genomic DNA in which 5hmC—but not 5mC or C—is read as T."
SEQUENCING_INTRO="Custom NovaSeq primers read the Nextera landing sites; the two eight-cycle index reads recover the combinatorial Tn5 cell barcode."

_initial=Scene.duplex([seg("single-cell genomic DNA","D"*52,placeholder=True)],label="genomic DNA from one FACS-sorted cell")
INITIAL_ROWS=tuple(_initial.rows()); INITIAL_NAME="Single-cell genomic DNA"

def _base_state(title:str,stage:int)->Scene:
    states={base:path.after(stage) for base,path in CAPS_PATHS.items()}
    sc=_tagged_scene("cell-barcoded genomic fragment retained through conversion")
    sc.strand("5hmC",[seg(f"5hmC path: {states['5hmC']}","H","w1",placeholder=True)],label=f"5hmC: {states['5hmC']}")
    sc.strand("5mC",[seg(f"5mC path: {states['5mC']}","M",placeholder=True)],label=f"5mC: {states['5mC']}")
    sc.strand("C",[seg(f"C path: {states['C']}","C",placeholder=True)],label=f"unmodified C: {states['C']}")
    sc.footer(title,"C"); return sc

def _tagged_scene(label:str)->Scene:
    sc=Scene.duplex([
        seg("i5 cell barcode","J"*8,"cbc",placeholder=True,
            feature=feature("cell_barcode_i5","cell_barcode","combinatorial",
                            group="cell_barcode",part="i5",whitelist="Additional file 5 Table S3")),
        seg("left Tn5 end","X"*12,"me",placeholder=True),
        seg("genomic fragment","D"*40,placeholder=True),
        seg("right Tn5 end","X"*12,"me",placeholder=True),
        seg("i7 cell barcode","I"*8,"cbc",placeholder=True,
            feature=feature("cell_barcode_i7","cell_barcode","combinatorial",
                            group="cell_barcode",part="i7",whitelist="Additional file 5 Table S3")),
    ],label=label)
    sc.mark("top","i5 cell barcode","well-specific i5/i7 transposome")
    return sc

def workflow():
    tagged=_tagged_scene("barcoded, gap-filled single-cell fragment")
    pooled=_tagged_scene("pooled DNA from 96 uniquely indexed wells plus controls")
    formyl=_base_state("ACT+BF4− oxidation",1)
    carboxyl=_base_state("Pinnick oxidation with NaClO2",2)
    reduced=_base_state("pyridine–borane reduction",3)
    final=Scene.duplex(list(FINAL_LIBRARY),label="PCR readout: 5hmC becomes T")
    w=Workflow(MolecularState(INITIAL_NAME,INITIAL_ROWS))
    w.react("Tagment with a well-specific Tn5 and gap-fill",tagged.rows(),name="Cell-barcoded fragments",note="Tn5 acts for 15 min at 50 °C; SDS strips Tn5 and 72 °C extension fills the nine-base gaps.")
    w.react("Pool 96 wells",pooled.rows(),name="Pooled cell-barcoded DNA",note="Filler DNA and modification controls are added before conversion chemistry.")
    w.react("Oxidize 5hmC to 5fC with ACT+BF4−",formyl.rows(),name="5fC-containing DNA",note="This first chemical oxidation leaves 5mC and unmodified C unchanged.")
    w.react("Pinnick-oxidize 5fC to 5caC",carboxyl.rows(),name="5caC-containing DNA",note="NaClO2 and 2-methyl-2-butene complete selective oxidation of the original 5hmC.")
    w.react("Reduce 5caC with pyridine–borane",reduced.rows(),name="DHU-containing DNA",note="Only the original 5hmC track reaches DHU.")
    w.react("PCR with Nextera XT index primers",final.rows(),name="Sequencing library",note="Uracil-tolerant PCR copies DHU as T; 5mC and unmodified C remain C.")
    return w


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    genomic=g.input("single_cell_dna","genomic DNA from one FACS-sorted cell",MaterialKind.MOLECULAR,state=_initial,primary=True)
    transposome=g.input("indexed_transposome","well-specific i5/i7 loaded Tn5 adapters",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (96 oligos referenced in Additional file 5 Table S3)")
    tn5=g.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=g.product("tagged","well-indexed single-cell genomic fragments with nine-base gaps",MaterialKind.MOLECULAR,state=_tagged_scene("barcoded single-cell fragments before gap fill"),inferred=True,source=SOURCE+" (96 oligos referenced rather than copied)")
    tagment_waste=g.product("tagment_waste","non-library tagmentation products",MaterialKind.WASTE)
    g.react("tagment","Tagment with a well-specific indexed Tn5",StepKind.TAGMENTATION,((genomic,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,tagment_waste))
    g.discard(tagment_waste,"does not retain a library-eligible heterologous indexed fragment")

    gap_polymerase=g.input("gap_polymerase","nine-base gap-fill polymerase",MaterialKind.CATALYST)
    gap_filled=g.product("gap_filled","barcoded gap-filled single-cell fragments",MaterialKind.MOLECULAR,state=_tagged_scene("barcoded, gap-filled single-cell fragment"),inferred=True,source=SOURCE+" (96 oligos referenced rather than copied)")
    g.react("gap_fill","Fill the nine-base Tn5 gaps",StepKind.EXTENSION,((tagged,InputRole.TEMPLATE),(gap_polymerase,InputRole.CATALYST)),(gap_filled,))

    peer_wells=g.input("peer_wells","corresponding indexed fragments from the other processed wells",MaterialKind.MOLECULAR,state=_tagged_scene("other well-indexed gap-filled fragments"),source=SOURCE)
    pooled=g.product("pooled","pooled cell-barcoded DNA from 96 wells",MaterialKind.MOLECULAR,state=_tagged_scene("pooled DNA from 96 uniquely indexed wells plus controls"),inferred=True,source=SOURCE+" (96 oligos referenced rather than copied)")
    g.react("pool_wells","Pool the 96 independently indexed wells",StepKind.POOL,((gap_filled,InputRole.TEMPLATE),(peer_wells,InputRole.TEMPLATE)),(pooled,),changes_molecule=False)

    act=g.input("act_bf4","ACT+BF4-minus oxidation condition",MaterialKind.CONDITION)
    formyl=g.product("formyl","pooled DNA after selective 5hmC-to-5fC oxidation",MaterialKind.MOLECULAR,state=_base_state("ACT+BF4− oxidation",1))
    g.react("oxidize_formyl","Oxidize 5hmC to 5fC with ACT+BF4-minus",StepKind.CHEMICAL_CONVERSION,((pooled,InputRole.TEMPLATE),(act,InputRole.CONDITION)),(formyl,))
    pinnick=g.input("pinnick","Pinnick oxidation condition with NaClO2",MaterialKind.CONDITION)
    carboxyl=g.product("carboxyl","pooled DNA after 5fC-to-5caC oxidation",MaterialKind.MOLECULAR,state=_base_state("Pinnick oxidation with NaClO2",2))
    g.react("oxidize_carboxyl","Pinnick-oxidize 5fC to 5caC",StepKind.CHEMICAL_CONVERSION,((formyl,InputRole.TEMPLATE),(pinnick,InputRole.CONDITION)),(carboxyl,))
    reducing=g.input("pyridine_borane","pyridine-borane reduction condition",MaterialKind.CONDITION)
    reduced=g.product("reduced","DNA in which original 5hmC has become DHU",MaterialKind.MOLECULAR,state=_base_state("pyridine–borane reduction",3))
    g.react("reduce","Reduce 5caC to DHU with pyridine-borane",StepKind.CHEMICAL_CONVERSION,((carboxyl,InputRole.TEMPLATE),(reducing,InputRole.CONDITION)),(reduced,))

    forward=g.input("pcr_forward","Nextera XT cell-index PCR forward-primer family",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (96 oligos referenced in Additional file 5 Table S3)")
    reverse=g.input("pcr_reverse","Nextera XT cell-index PCR reverse-primer family",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (96 oligos referenced in Additional file 5 Table S3)")
    polymerase=g.input("pcr_polymerase","uracil-tolerant PCR polymerase",MaterialKind.CATALYST)
    library=g.product("library","scCAPS+ sequencing library",MaterialKind.LIBRARY,state=FINAL_LIBRARY,inferred=True,source=SOURCE+" (96 oligos referenced rather than copied)")
    g.react("pcr","PCR with Nextera XT cell-index primers",StepKind.PCR,((reduced,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),(reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),geometry=PCRGeometry.unavailable(FINAL_LIBRARY,reason="The source references the 96 indexed transposome/PCR oligos rather than printing one universal primer pair on this page."))
    g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=SEQ_PRIMERS,acquired_roles=tuple(p.role for p in SEQ_PRIMERS),read_lengths=READ_LENGTHS))
    reviewed=["single_cell_dna","indexed_transposome","tn5","tagment","gap_polymerase","gap_fill","peer_wells","pool_wells","act_bf4","oxidize_formyl","pinnick","oxidize_carboxyl","pyridine_borane","reduce","pcr_forward","pcr_reverse","pcr_polymerase","pcr","library"]
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
