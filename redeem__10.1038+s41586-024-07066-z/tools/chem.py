"""ReDeeM: source-owned declarations and reaction graph.

Source scope (bounded pass 2026-10-17-breadth-first-16):
- Weng et al. 2024 Nature, Supplementary Information (MOESM1) Supplementary Methods
  "ReDeeM Method" (extracted lines 414-478): fixation, GDN tagmentation, CG000338 Rev A
  Steps 2-4/6-7, in-house N70X/SI-PCR-primer-B ATAC PCR, four-pool xGen capture,
  P5/P7 post-capture PCR, NovaSeq read lengths.
- Supplementary Data 1 (MOESM5) "Primers" sheet: Nextera_N701-N706, P5, P7.
- Supplementary ReDeeM protocol 2023-11-20 (MOESM3), compared for the same steps.
- 10x CG000338 Rev A Appendix "Sequences" (cited by Supplementary Methods): ATAC bead
  oligo and final ATAC / gene-expression library layouts.
"""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from protocol_model import molecule,role
from chemdraw import Scene,Segment,complement_segments
from adapter_ligation import AdapterLigationGeometry
seg=Segment
from chemdraw import Construct, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind, Disposition,Terminal,Highlight)

TITLE="ReDeeM — deep mitochondrial mutation profiling"
NOTES="01_redeem.html"
SOURCE='Weng et al., Nature(2024), <a href="https://doi.org/10.1038/s41586-024-07066-z">doi:10.1038/s41586-024-07066-z</a>; SupplementaryMethods/Data1 and explicitly cited originalCG000338RevA.'
SUMMARY="Whole-cell Multiome retains nuclear DNA, mtDNA and RNA as separate populations; source-indexed ATAC aliquots are hybrid-captured in four mtDNA probe reactions and reamplified."
CAVEAT="Source-positive final structures follow the explicitly cited CG000338RevA appendix; this component scope does not establish procedural equivalence to RevB cited by the separate protocol. Native GEM/preamp/capture first-cycle assembly remains unassigned. Barcode plus Tn5start/end defines a computational endogenous UMI, not an added physical UMI sequence. Published probe/revision/run discrepancies remain below."
CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="10x Multiome ATAC gel-bead barcodes (CG000338)")
GEX_BARCODE = feature("gex_cell_barcode", "cell_barcode", "whitelist",
                      whitelist="10x Multiome gene-expression gel-bead barcodes (CG000338)")
UMI = feature("gex_umi", "umi", "random")
N70X_INDEX = feature("nextera_n70x", "sample_index", "whitelist",
                     whitelist="Weng 2024 Supplementary Data 1 Nextera_N701-N706")
GEX_I7 = feature("gex_i7", "sample_index", "unknown")
GEX_I5 = feature("gex_i5", "sample_index", "unknown")
TENX_ATAC_SPACER = "CGCGTCTG"      # CG000338 Rev A Appendix, ATAC gel bead and library

# Supplementary Data 1, Primers sheet.
N70X_INDEXES = {"N701": "TCGCCTTA", "N702": "CTAGTACG", "N703": "TTCTGCCT",
                "N704": "GCTCAGGA", "N705": "AGGAGTCC", "N706": "CATGCCTA"}
REDEEM_P5 = "AATGATACGGCGACCACCGAGA"
REDEEM_P7 = il.P7


def n70x_primer(index: str) -> str:
    return il.P7 + index + nx.S7


def atac_shell(name: str, insert: str) -> Construct:
    """CG000338 Rev A ATAC library with the ReDeeM N70X i7 index."""
    return Construct([
        seg("P5", il.P5, "p5"),
        seg("10x cell barcode", "B"*16, "cbc", placeholder=True, feature=CELL_BARCODE),
        seg("10x spacer", TENX_ATAC_SPACER),
        seg("Read 1N (S5)", nx.S5, "s5"), seg("left mosaic end", nx.ME, "me"),
        seg(insert, "X"*34, placeholder=True),
        seg("right mosaic end reverse complement", nx.ME_RC, "me"),
        seg("Read 2N (S7) reverse complement", nx.S7_RC, "s7"),
        seg("N70X i7 index reverse complement", "I"*8, "cbc", placeholder=True,
            feature=N70X_INDEX),
        seg("P7 reverse complement", il.P7_RC, "p7")], name=name)


ATAC_LIBRARY = atac_shell("ReDeeM ATAC library", "accessible-DNA insert")
MT_LIBRARY = atac_shell("ReDeeM enriched mtDNA library", "mtDNA fragment (Tn5 start–end)")
# CG000338 Rev A prints the GEX P5 arm as AATGATACGGCGACCACCGAGATCT before N10.
GEX_LIBRARY = Construct([
    seg("P5 (Rev A arm)", il.P5[:-4], "p5"),
    seg("i5 sample index", "J"*10, "cbc", placeholder=True, feature=GEX_I5),
    seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
    seg("10x cell barcode", "B"*16, "cbc", placeholder=True, feature=GEX_BARCODE),
    seg("UMI", "U"*12, "umi", placeholder=True, feature=UMI),
    seg("poly(dT)VN", "T"*30 + "VN", placeholder=True),
    seg("cDNA insert", "X"*36, placeholder=True),
    seg("TruSeq Read 2 reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
    seg("i7 sample index reverse complement", "I"*10, "cbc", placeholder=True, feature=GEX_I7),
    seg("P7 reverse complement", il.P7_RC, "p7")], name="ReDeeM RNA (Multiome GEX) library")
ATAC_PRIMERS = tuple(sp.NEXTERA[k] for k in ("R1", "I1", "I2", "R2"))
GEX_PRIMERS = tuple(sp.TRUSEQ[k] for k in ("R1", "I1", "I2", "R2"))
for _lib, _primers in ((ATAC_LIBRARY, ATAC_PRIMERS), (MT_LIBRARY, ATAC_PRIMERS),
                       (GEX_LIBRARY, GEX_PRIMERS)):
    _problems = sp.verify(_lib, _primers, required_roles=tuple(p.role for p in _primers))
    if _problems:
        raise ValueError(f"{_lib.name}: " + "; ".join(_problems))
# Supplementary Methods "Sequencing for ReDeeM libraries": NovaSeq S4 300 cycles,
# Read1 149, i7 8, i5 24, Read2 149.
ATAC_READ_LENGTHS = {"Read 1": 149, "Index 1 (i7)": 8, "Index 2 (i5)": 24, "Read 2": 149}

def oligos():return [
    f"Nextera N70X (Supplementary Data 1): 5′-{il.P7}<b>[8-nt index]</b>{nx.S7}-3′; "
    + ", ".join(f"{k} {v}" for k, v in N70X_INDEXES.items()),
    f"Post-capture P5: 5′-{REDEEM_P5}-3′; P7: 5′-{REDEEM_P7}-3′ (Supplementary Data 1)",
    f"10x ATAC gel-bead oligo (CG000338 Rev A): 5′-{il.P5}<b>[N16]</b>{TENX_ATAC_SPACER}-3′",
]


FINAL_LIBRARIES=(("ATAC library",ATAC_LIBRARY,ATAC_PRIMERS,"SourceRevA structure and localN70X index; complete native precursor geometry bounded.","Published SupplementaryMethods149/8/24/149."),("RNA library",GEX_LIBRARY,GEX_PRIMERS,"SourceRevA GEX structure; published acquisition statements conflict.","RNA run acquisition remains unassigned; both published alternatives shown below."),("Enriched mtDNA library",MT_LIBRARY,ATAC_PRIMERS,"Same indexedATAC architecture enriched by capture; first-copy substrate is bead-bound DNA.","Published SupplementaryMethods149/8/24/149."))
SHOW_READ_CONTENT=False
READ_SPECS={name:tuple(sp.read_spec(read,cycles,"genomic/mtDNA insert" if read in ("Read 1","Read 2") else "N70X sample index" if read=="Index 1 (i7)" else "source8ntspacer and16ntcellbarcode", "include" if read in ("Read 1","Read 2") else "demultiplex" if read=="Index 1 (i7)" else "extract") for read,cycles in ATAC_READ_LENGTHS.items()) for name in ("ATAC library","Enriched mtDNA library")}
READ_SPECS["RNA library"]=tuple(sp.read_spec(read,"unassigned",content+"; conflicting published run descriptions", "undetermined") for read,content in (("Read 1","source16ntbarcode/12ntUMI"),("Index 1 (i7)","physical10ntsampleindex"),("Index 2 (i5)","physical10ntsampleindex"),("Read 2","transcript cDNA")))
READ_SPECS_NOTE="ATAC/mtDNA positive source acquisition retained; RNA conflicting run alternatives do not establish a selected mask. Physical index widths are separate from acquired cycles."

def cell_state(label,*,tagged=False,copied=False,amplified=False):
    sc=Scene()
    for key,description in (("nuclear","nuclear chromatin"),("mitochondrial","mtDNA")):
        parts=[role(("Tn5-tagged " if tagged else "")+description+"; native ends unassigned")]
        if copied:parts=[Segment("cellbarcode","N"*16,placeholder=True,feature=CELL_BARCODE),role(description+" barcoded fragment; native end/gap assembly unassigned")]
        sc.strand(key,parts,label=label+" independent "+description+" population")
    if copied:
        parts=[Segment("sourceRevA partialR1","CTACACGACGCTCTTCCGATCT"),Segment("GEXcellbarcode","N"*16,placeholder=True,feature=GEX_BARCODE),Segment("GEXUMI","N"*12,placeholder=True,feature=UMI),Segment("sourceT30","T"*30),Segment("VN","VN",placeholder=True),role("first-strand cDNA; native copiedTSO/end extent unassigned")]
    else:parts=[role("cellularRNA")]
    sc.strand("RNA",parts,label=label+(" RNA-derived strand population" if copied else " independent RNA population"))
    if amplified:sc.strand("opposite copied strands",[role("complementary preamp strand classes; full native duplex geometry unassigned")],label="independent preamp copies")
    return sc

def rna_fragment_state(*,a_tailed=False,adapted=False):
    """One barcode-retaining target class, with other fragments carried separately."""
    body=[Segment("retained source partialR1","CTACACGACGCTCTTCCGATCT"),Segment("retained GEX cell barcode","N"*16,placeholder=True,feature=GEX_BARCODE),Segment("retained GEX UMI","N"*12,placeholder=True,feature=UMI),role("RNA-derived fragment interval; full native extent unassigned")]
    if adapted:body=[role("source adapter-associated end; native strand/junction unassigned",inferred=True),*body,role("source opposite adapter-associated end; native strand/junction unassigned",inferred=True)]
    opposite=complement_segments(body)
    top=[*body,Segment("top3′dA","A")] if a_tailed else body
    bottom=[*opposite,Segment("bottom3′dA","A")] if a_tailed else opposite
    sc=Scene();sc.strand("RNA target",top,label="barcode-retaining RNA-derived fragment class; native end extent bounded")
    sc.anneal("RNA target opposite",bottom,to="RNA target",pair=(opposite[-1].name,body[0].name),unpaired=tuple(x.name for x in opposite if x.is_role_token()),label="opposite target strand;5′chemistry unassigned")
    if a_tailed:
        sc.overhang("RNA target","top3′dA",polarity="3-prime",owner="substrate",text="source3′dA")
        sc.overhang("RNA target opposite","bottom3′dA",polarity="3-prime",owner="substrate",text="source opposed3′dA")
    sc.strand("carried fragments",[role("other cDNA, nuclearDNA and mtDNA fragment classes; barcode retention/recovery unassigned")],label="independent carried mixed fragments")
    return sc

def capture_state():
    sc=Scene();sc.strand("captured DNA",list(MT_LIBRARY),label="source-denatured indexed DNA strand bound through capture probe")
    sc.strand("probes",[role("biotinylated mtDNA capture probes; exact bound-strand extent unassigned")],label="independent probe oligos")
    sc.strand("background",[role("co-retained off-target indexedDNA; quantitative specificity unassigned")],label="independent background")
    return sc

def protocol_graph():
    source=SOURCE
    graph=ProtocolGraph(TITLE,authority=source)
    cells = graph.input("cells", "whole cells (optionally TotalSeq-A hashed)",
                        MaterialKind.BIOLOGICAL,
                        state=cell_state("whole cell"),
                        primary=True)
    fixperm = graph.input("fix_perm", "0.1% formaldehyde 10 min RT, 0.125 M glycine quench; 0.1% NP40 permeabilization 5 min on ice; low-BSA buffers",
                          MaterialKind.CONDITION)
    permeable = graph.product("fixed_cells", "fixed, permeabilized whole cells retaining mtDNA",
                              MaterialKind.BIOLOGICAL, state=cells.state)
    graph.react("fix_permeabilize", "Mildly fix and permeabilize whole cells",
                StepKind.PERMEABILIZATION,
                ((cells, InputRole.TEMPLATE), (fixperm, InputRole.CONDITION)), (permeable,),
                changes_molecule=False)
    tn5 = graph.input("tn5", "10x ATAC Enzyme B in ATAC buffer B with 0.1% GDN, 37 °C 1 h",
                      MaterialKind.CATALYST)
    tn5_adapters = graph.input("tn5_adapters", "ATACEnzymeB loadedadapter population; native assembly unassigned",
                               MaterialKind.ADAPTER)
    tagmented = graph.product("tagmented_cells", "cells with tagmented nuclear and mitochondrial DNA",
                              MaterialKind.BIOLOGICAL,
                              state=cell_state("tagmented cell",tagged=True))
    graph.react("tagment", "Tagment whole cells with GDN added", StepKind.TAGMENTATION,
                ((permeable, InputRole.TEMPLATE), (tn5_adapters, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagmented,))
    droplets = graph.input("multiome_droplets", "10x Multiome GEM generation, barcoding, cleanup (CG000338 Rev A Steps 2–3)",
                           MaterialKind.CONDITION)
    barcoded = graph.product("barcoded", "matched cell-barcoded RNA and ATAC products",
                             MaterialKind.MOLECULAR,
                             state=cell_state("GEMbarcoded",tagged=True,copied=True))
    failed = graph.product("failed_partitions", "empty or ambiguous droplet outcomes", MaterialKind.MOLECULAR)
    graph.react("multiome", "Run 10x Multiome barcoding on the tagmented whole cells",
                StepKind.MAGIC,
                ((tagmented, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (barcoded, failed),reason="Source RevA GEMgeneration/barcoding/postGEMcleanup scope; RT creates ss cDNA beside independently barcoded DNA. Complete native transfer/recovery geometry remains bounded.")
    graph.mark_terminal(Terminal(failed.id,Disposition.UNRESOLVED,"Droplet occupancy/ambiguity classification is computational, not physical rejection."))
    preamp_f=graph.input("preamp_forward","sourceRevA Pre-AmpPrimersPN2000271 first binding role; species unassigned",MaterialKind.OLIGO)
    preamp_r=graph.input("preamp_reverse","same source commercial mix opposite binding role; species unassigned",MaterialKind.OLIGO)
    preamp_pol=graph.input("preamp_polymerase","sourceRevA pre-amplification PCRmix; exactenzyme/formulation unassigned",MaterialKind.CATALYST)
    preamp=graph.product("preamp","pre-amplified mixedMultiome product;100µlEB",MaterialKind.MOLECULAR,state=cell_state("preamp products",tagged=True,copied=True,amplified=True))
    graph.react("preamplify","Pre-amplify RNA andATAC products before physicalaliquots",StepKind.PCR,((barcoded,InputRole.TEMPLATE),(preamp_f,InputRole.STRUCTURAL),(preamp_r,InputRole.STRUCTURAL),(preamp_pol,InputRole.CATALYST)),(preamp,),geometry=PCRGeometry.unavailable(barcoded.state,reason="Actualss-cDNA beside barcoded nuclear/mitochondrialDNA enters named commercial preampmix; full native first-cycle/primer-to-product geometry unassigned."))
    atac = graph.product("atac_branch", "pre-amp aliquots for ATAC PCR (8 × 5 µl)",
                         MaterialKind.MOLECULAR, state=preamp.state)
    rna = graph.product("rna_branch", "pre-amp aliquot for RNA library (35 µl)",
                        MaterialKind.MOLECULAR, state=preamp.state)
    remainder=graph.product("preamp_remainder","remaining25µl source100µl preamp pool",MaterialKind.MOLECULAR,state=preamp.state)
    graph.react("split_branches", "Split the pre-amplified product into ATAC and RNA aliquots",
                StepKind.SPLIT, ((preamp, InputRole.TEMPLATE),), (atac, rna,remainder),
                changes_molecule=False)
    graph.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,"NamedATAC8×5µl andRNA35µl aliquots leave25µl; subsequent handling unassigned."))
    reviewed = ["cells", "fix_perm", "fix_permeabilize", "tn5", "tn5_adapters", "tagment",
                "multiome_droplets", "multiome", "preamp_forward","preamp_reverse","preamp_polymerase", "preamplify", "split_branches"]

    atac_forward = graph.input("atac_forward", "10x SI PCR Primer B PN2000128; exact species unassigned",
                               MaterialKind.OLIGO)
    atac_reverse = graph.input("atac_reverse", "Nextera N70X primer (Supplementary Data 1)",
                               MaterialKind.OLIGO,
                               state=Construct([seg("P7", il.P7, "p7"),
                                                seg("N70X index", "I"*8, "cbc", placeholder=True,
                                                    feature=N70X_INDEX),
                                                seg("S7", nx.S7, "s7")], name="Nextera N70X"))
    atac_polymerase = graph.input("atac_polymerase", "NEBNext HF PCR master mix, SYBR-monitored, 6–12 cycles",
                                  MaterialKind.CATALYST)
    amplified_atac = graph.product("atac_amplified", "indexed ATAC PCR product",
                                   MaterialKind.MOLECULAR, state=ATAC_LIBRARY)
    graph.react("atac_pcr", "In-house ATAC PCR with N70X and SI PCR Primer B (8 × 50 µl)",
                StepKind.PCR,
                ((atac, InputRole.TEMPLATE), (atac_forward, InputRole.STRUCTURAL),
                 (atac_reverse, InputRole.STRUCTURAL), (atac_polymerase, InputRole.CATALYST)),
                (amplified_atac,), geometry=PCRGeometry.unavailable(
                    atac.state, reason="Actualmixedpreampaliquot retained; SI PCR Primer B is named by part number, so the P5-side primer-to-product geometry is not reconstructed."))
    atac_spri = graph.input("atac_spri", "SPRI double-sided size selection 0.4×–1.55×",
                            MaterialKind.SELECTION)
    indexed_atac = graph.product("atac_indexed", "indexed ATAC library including mtDNA fragments",
                                 MaterialKind.LIBRARY, state=ATAC_LIBRARY)
    atac_rejected = graph.product("atac_rejected_sizes", "out-of-range ATAC fragments",
                                  MaterialKind.WASTE)
    graph.react("atac_size_select", "Pool and size-select the ATAC PCR product", StepKind.SELECTION,
                ((amplified_atac, InputRole.TEMPLATE), (atac_spri, InputRole.SELECTION)),
                (indexed_atac, atac_rejected), changes_molecule=False)
    graph.discard(atac_rejected, "large fragments on the 0.4× beads and the 1.55× supernatant are discarded")
    atac_for_sequencing = graph.product("atac_library", "ATAC library aliquot for sequencing",
                                        MaterialKind.LIBRARY, state=ATAC_LIBRARY)
    atac_for_capture = graph.product("atac_capture_input", "ATAC library for mtDNA capture (4 × 500 ng)",
                                     MaterialKind.LIBRARY, state=ATAC_LIBRARY)
    graph.react("split_atac_library", "Reserve indexed ATAC aliquots for sequencing and mtDNA capture",
                StepKind.SPLIT, ((indexed_atac, InputRole.TEMPLATE),),
                (atac_for_sequencing, atac_for_capture), changes_molecule=False)
    graph.sequence(atac_for_sequencing, SequencingPlan(
        library=ATAC_LIBRARY, primers=ATAC_PRIMERS, read_lengths=ATAC_READ_LENGTHS))
    reviewed.extend(("atac_forward", "atac_reverse", "atac_polymerase", "atac_pcr",
                     "atac_spri", "atac_size_select", "atac_indexed", "atac_library",
                     "split_atac_library"))

    rna_cf=graph.input("rna_cdna_forward","source RevA cDNA-PCR mix first binding role; primer species unassigned",MaterialKind.OLIGO)
    rna_cr=graph.input("rna_cdna_reverse","same historical cDNA-PCR mix opposite binding role; species unassigned",MaterialKind.OLIGO)
    rna_cp=graph.input("rna_cdna_polymerase","source RevA cDNA-amplification reagents; exact enzyme/formulation/cycle count unassigned",MaterialKind.CATALYST)
    rna_amp=graph.product("rna_cdna_amplified","cDNA-amplified RNA-preparation mixture; carried DNA retained",MaterialKind.MOLECULAR,state=cell_state("RNA cDNA-PCR products",tagged=True,copied=True,amplified=True))
    graph.react("rna_cdna_pcr","Amplify cDNA from actual35µl mixed preamp aliquot",StepKind.PCR,((rna,InputRole.TEMPLATE),(rna_cf,InputRole.STRUCTURAL),(rna_cr,InputRole.STRUCTURAL),(rna_cp,InputRole.CATALYST)),(rna_amp,),geometry=PCRGeometry.unavailable(rna.state,reason="Actual mixed preamp aliquot enters cDNA-PCR. Historical primer species and complete native first-cycle/end geometry remain unassigned."))
    frag=graph.input("rna_fragmentation_reagents","source RevA library fragmentation reagents; formulation/conditions unassigned in inspected scope",MaterialKind.CATALYST)
    fragments=graph.product("rna_fragments","native fragmented RNA-library preparation",MaterialKind.MOLECULAR,state=rna_fragment_state())
    graph.react("rna_fragment","Fragment amplified RNA-library preparation",StepKind.FRAGMENTATION,((rna_amp,InputRole.TEMPLATE),(frag,InputRole.CATALYST)),(fragments,))
    repair=graph.input("rna_end_repair_reagents","source historical end-repair reagents; exact formulation unassigned",MaterialKind.CATALYST)
    repaired=graph.product("rna_repaired","end-repaired mixed fragment population",MaterialKind.MOLECULAR,state=rna_fragment_state())
    graph.react("rna_end_repair","End-repair library fragments",StepKind.END_REPAIR,((fragments,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(repaired,))
    tail=graph.input("rna_a_tail_reagents","source dA-tailing reagents; exact enzyme/conditions unassigned",MaterialKind.CATALYST)
    tailed=graph.product("rna_a_tailed","source dA-tailed native fragment classes",MaterialKind.MOLECULAR,state=rna_fragment_state(a_tailed=True))
    graph.react("rna_a_tail","dA-tail repaired fragments",StepKind.TAILING,((repaired,InputRole.TEMPLATE),(tail,InputRole.CATALYST)),(tailed,))
    adapter=graph.input("rna_adapter","source RevA library adapter species; physical strands/ends unassigned",MaterialKind.ADAPTER)
    ligase=graph.input("rna_adapter_ligase","source adapter-ligation reagents; enzyme/formulation unassigned",MaterialKind.CATALYST)
    ligated=graph.product("rna_ligated","actual pre-index adapter-ligated mixed fragments",MaterialKind.MOLECULAR,state=rna_fragment_state(adapted=True))
    graph.react("rna_adapter_ligation","Ligate source-named RNA-library adapters",StepKind.LIGATION,((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(ligated,),geometry=AdapterLigationGeometry.unavailable("Source dA-tailed substrate retained, but full historical adapter strand/compatible-end geometry and sealed junction remain unassigned."))
    index_f=graph.input("rna_index_forward","source RNA sample-index PCR forward binding role; selected species unassigned",MaterialKind.OLIGO)
    index_r=graph.input("rna_index_reverse","source RNA sample-index PCR reverse binding role; selected species unassigned",MaterialKind.OLIGO)
    index_p=graph.input("rna_index_polymerase","source RevA RNA sample-index PCR reagents; exact enzyme/conditions unassigned",MaterialKind.CATALYST)
    rna_library=graph.product("rna_library",GEX_LIBRARY.name,MaterialKind.LIBRARY,state=GEX_LIBRARY)
    graph.react("rna_index_pcr","Index-amplify actual adapter-ligated fragments",StepKind.PCR,((ligated,InputRole.TEMPLATE),(index_f,InputRole.STRUCTURAL),(index_r,InputRole.STRUCTURAL),(index_p,InputRole.CATALYST)),(rna_library,),geometry=PCRGeometry.unavailable(ligated.state,reason="Actual pre-index adapter-ligated fragment pool is the substrate. Source RevA mature structure is a separate reference; unknown primer/first-cycle geometry is not supplied by that endpoint."))
    graph.sequence(rna_library,SequencingPlan(library=GEX_LIBRARY,primers=GEX_PRIMERS,acquisition="unspecified",limitation="Conflicting published RNA masks; no acquired roles/cycles selected."))

    probes = graph.input("capture_probes", "one of four staggered biotinylated 120-nt mtDNA probe sets (mitoV1–V4; Supplementary Data 1)",
                         MaterialKind.SELECTION)
    blockers = graph.input("hyb_condition", "Cot DNA, xGen Universal Blockers NXT, xGen hybridization buffer + enhancer; 95 °C 30 s then 65 °C 16 h; streptavidin Dynabeads and xGen washes",
                           MaterialKind.CONDITION)
    captured = graph.product("captured_mt", "bead-bound captured mtDNA library fragments",
                             MaterialKind.MOLECULAR, state=capture_state())
    rejected = graph.product("rejected_mt", "source unboundDNA/capturewashes; exact mtDNA/off-target recovery unassigned",
                             MaterialKind.WASTE)
    graph.react("capture_mt", "Hybrid-capture mtDNA from the indexed ATAC library (four parallel probe-set reactions)",
                StepKind.AFFINITY_CAPTURE,
                ((atac_for_capture, InputRole.TEMPLATE), (probes, InputRole.SELECTION),
                 (blockers, InputRole.CONDITION)),
                (captured, rejected))
    graph.discard(rejected, "removed with the supernatant and xGen washes")
    mt_forward = graph.input("mt_forward", "P5 primer (Supplementary Data 1)", MaterialKind.OLIGO,
                             state=Construct([seg("P5 primer", REDEEM_P5, "p5")], name="ReDeeM P5"))
    mt_reverse = graph.input("mt_reverse", "P7 primer (Supplementary Data 1)", MaterialKind.OLIGO,
                             state=Construct([seg("P7 primer", REDEEM_P7, "p7")], name="ReDeeM P7"))
    mt_polymerase = graph.input("mt_polymerase", "KAPA HiFi HotStart ReadyMix on beads, 3–5 cycles",
                                MaterialKind.CATALYST)
    mt_amplified = graph.product("mt_amplified", "post-capture amplified mtDNA library",
                                 MaterialKind.MOLECULAR, state=MT_LIBRARY)
    graph.react("mt_pcr", "Post-capture P5/P7 PCR on the washed beads", StepKind.PCR,
                ((captured, InputRole.TEMPLATE), (mt_forward, InputRole.STRUCTURAL),
                 (mt_reverse, InputRole.STRUCTURAL), (mt_polymerase, InputRole.CATALYST)),
                (mt_amplified,), geometry=PCRGeometry.unavailable(
                    captured.state, reason="Actualsingle-strand probe-bound indexedDNA retained. The P5 and P7 primers re-amplify the existing flow-cell ends; this compact model does not draw the bead-bound template duplex."))
    mt_spri = graph.input("mt_spri", "four reactions pooled, 1.6× SPRI", MaterialKind.SELECTION)
    mt_library = graph.product("mt_library", MT_LIBRARY.name, MaterialKind.LIBRARY,
                               state=MT_LIBRARY)
    mt_waste = graph.product("mt_spri_waste", "SPRI supernatant and washes", MaterialKind.WASTE)
    graph.react("mt_cleanup", "Pool the four capture PCRs and purify", StepKind.CLEANUP,
                ((mt_amplified, InputRole.TEMPLATE), (mt_spri, InputRole.SELECTION)),
                (mt_library, mt_waste), changes_molecule=False)
    graph.discard(mt_waste, "bead supernatant and washes are discarded")
    graph.sequence(mt_library, SequencingPlan(
        library=MT_LIBRARY, primers=ATAC_PRIMERS, read_lengths=ATAC_READ_LENGTHS))
    reviewed.extend(("capture_probes", "hyb_condition", "capture_mt", "mt_forward", "mt_reverse",
                     "mt_polymerase", "mt_pcr", "mt_spri", "mt_cleanup", "mt_library"))
    graph.annotate(Highlight("multiome","Native RNA fate across GEM cleanup remains unassigned","The RNA-derived output row tracks copied cDNA. Native RNA templates are not represented as consumed or physically discarded by RT; their persistence, degradation and partition through the source GEM/post-GEM cleanup boundary remain unassigned."))
    graph.declare_reviewed_scope(REVIEWED_TARGETS,source)
    graph.annotate(SourceDiscrepancy(
        "capture_probes", "ReDeeM capture-probe specification changed across published supplements",
        (SourceClaim("Weng et al., Nature (2024), supplementary ReDeeM protocol, extracted lines 64–65",
                     "Four staggered sets of 137 probes are described as 100 nt, with sequences forthcoming."),
         SourceClaim("Same article, Supplementary Methods, extracted lines 456–465",
                     "Each of four sets is described as 138 biotinylated probes of 120 nt."),
         SourceClaim("Same article, Supplementary Data 1, eight MitoProbe sheets",
                     "Individual human and mouse probe bases are supplied, usually 120 nt; sheet counts vary and a few human v1 rows are shorter.")),
        "The model treats Supplementary Data 1 as the sequence inventory for the published version.",
        "The protocol text predates the workbook and disagrees with the final supplementary Methods and sequence inventory.",
        "The workbook lists bases but does not by itself establish each oligo's synthesis modification or a sample-specific pool assignment."))
    graph.annotate(SourceDiscrepancy(
        "multiome", "ReDeeM supplements cite different revisions of the 10x Multiome guide",
        (SourceClaim("Weng et al., Nature (2024), Supplementary Methods, extracted lines 433–440",
                     "GEM generation and RNA library preparation cite CG000338 Rev A."),
         SourceClaim("Same article, supplementary ReDeeM protocol, extracted lines 144–156 and 416",
                     "The corresponding steps cite CG000338 Rev B.")),
        "The model uses the fixed 10x library shells of the Rev A Sequences appendix; the Rev B appendix prints the same shells.",
        "The supplements do not say which revision was used in each experiment.",
        "Procedural differences between the guide revisions remain unresolved."))
    graph.annotate(SourceDiscrepancy(
        "rna_library", "ReDeeM supplements give different read lengths",
        (SourceClaim("Weng et al., Nature (2024), Supplementary Methods, 'Sequencing for ReDeeM libraries'",
                     "All three libraries were sequenced on NovaSeq S4 300-cycle kits as Read1 149, i7 8, i5 24, Read2 149."),
         SourceClaim("Same article, supplementary ReDeeM protocol, steps 26, 35 and 37",
                     "ATAC and mtDNA: R1 150, i7 8, i5 24, R2 150; RNA: R1 28, i7 10, i5 10, R2 90 following the 10x guide.")),
        "ATAC and mtDNA reads use the Supplementary Methods 149/8/24/149. No acquired cycle counts are assigned to the RNA library.",
        "The Supplementary Methods describe the sequencing actually performed; the protocol gives recommendations.",
        "For the RNA library, 8 i7 cycles would read only part of the 10-nt 10x index printed in CG000338, and the protocol recommends 10/10. The cycles actually acquired remain unresolved."))
    return graph.require_publishable()



REVIEWED_TARGETS=("cells","fix_perm","fix_permeabilize","tn5","tn5_adapters","tagment","multiome_droplets","multiome","preamp_forward","preamp_reverse","preamp_polymerase","preamplify","split_branches","atac_forward","atac_reverse","atac_polymerase","atac_pcr","atac_spri","atac_size_select","atac_indexed","atac_library","split_atac_library","rna_cdna_forward","rna_cdna_reverse","rna_cdna_polymerase","rna_cdna_pcr","rna_fragmentation_reagents","rna_fragment","rna_end_repair_reagents","rna_end_repair","rna_a_tail_reagents","rna_a_tail","rna_adapter","rna_adapter_ligase","rna_adapter_ligation","rna_index_forward","rna_index_reverse","rna_index_polymerase","rna_index_pcr","rna_library","capture_probes","hyb_condition","capture_mt","mt_forward","mt_reverse","mt_polymerase","mt_pcr","mt_spri","mt_cleanup","mt_library")
