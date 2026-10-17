"""Single-molecule molecular inversion probe model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, Segment, circle_rows, feature
from circular import CircularizationGeometry, circularize_ssdna
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, MolecularRelation, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role

TITLE="Single-molecule molecular inversion probes (smMIPs)"
NOTES="01_smmip.html"
SOURCE='Hiatt et al. 2013, <a href="https://doi.org/10.1101/gr.147686.112">doi:10.1101/gr.147686.112</a>; later MIP protocol (different 70-mer, 5–8-nt-tag design) supplying the inferred PCR/run primers: <a href="https://doi.org/10.1007/978-1-4939-6442-0_6">Cantsilieris et al. 2017</a>.'
SUMMARY="A 12-nt-tagged inversion probe copies a target gap, circularizes, survives exonuclease and is universally amplified so reads can be collapsed by original molecule."
CAVEAT=("Probe arms, target bases and tag/index values are variable. The original panel's 28-nt invariant probe backbone is specified in unavailable Supplemental Table S9, so the capture circle keeps it symbolic. INFERRED — the drawn PCR handles and custom sequencing primers are transcribed from the later Cantsilieris et al. 2017 MIP protocol (Notes 8 and 18), which describes 70-mer MIPs with 5–8-nt tags rather than this 80-mer N12 panel; the 2013 article assigns its own PCR, sequencing and index primers to the same unavailable Table S9, so their identity with the 2013 primers is not established. "
        "Published reaction-volume discrepancy: Hiatt et al. 2013, Methods ‘Capture and library construction’, gap-fill recipe says to add water to 10 mL and then calls the total reaction 20 µL. Both cannot describe one reaction. This model uses neither volume; the intended water amount remains uncertain.")
MOLECULAR_TAG=feature("molecular_tag","umi","random")
I7_FEATURE=feature("sample_index_i7","sample_index","unknown")

FORWARD_HANDLE="ATACGAGATCCGTAATCGGGAAGCTGAAG"
REVERSE_HANDLE_TOP="ACACTACCGTCGGATCGTGCGTGT"
CAPTURED_GAP_NT=112  # Hiatt et al. 2013, Methods: the panel's fixed gap-fill length
READ1_PRIMER="CATACGAGATCCGTAATCGGGAAGCTGAAG"
READ2_PRIMER="ACACGCACGATCCGACGGTAGTGT"
INDEX1_PRIMER=REVERSE_HANDLE_TOP
UNIVERSAL_FORWARD_PCR=il.P5+FORWARD_HANDLE
INDEXED_REVERSE_PCR=il.P7+"N"*8+READ2_PRIMER

FINAL_LIBRARY=Construct([
 seg("P5",il.P5,"p5"),
 seg("smMIP forward handle",FORWARD_HANDLE,"r1",inferred=True),
 seg("extension targeting arm","X"*18,placeholder=True),
 seg("captured target","N"*CAPTURED_GAP_NT,placeholder=True),
 seg("ligation targeting arm","Y"*22,placeholder=True),
 seg("12-nt single-molecule tag","U"*12,"umi",placeholder=True,feature=MOLECULAR_TAG),
 seg("smMIP reverse handle",REVERSE_HANDLE_TOP,"r2",inferred=True),
 seg("eight-base i7 index reverse complement","I"*8,"cbc",placeholder=True,feature=I7_FEATURE),
 seg("P7 reverse complement",il.P7_RC,"p7"),
],name="single-index smMIP capture product")
SEQ_PRIMERS=(
 sp.custom("Read 1","MIP forward sequencing primer",READ1_PRIMER,
           "Cantsilieris et al. 2017, Note 18 (later protocol; 2013 Table S9 unavailable)",
           "The 5'-terminal C is an unpaired flap; the 31-nt 3' end anneals exactly."),
 sp.custom("Index 1 (i7)","MIP sample-index sequencing primer",INDEX1_PRIMER,
           "Cantsilieris et al. 2017, Note 18 (later protocol; 2013 Table S9 unavailable)"),
 sp.custom("Read 2","MIP reverse sequencing primer",READ2_PRIMER,
           "Cantsilieris et al. 2017, Note 18 (later protocol; 2013 Table S9 unavailable)"),
)
if problems:=sp.verify(FINAL_LIBRARY,SEQ_PRIMERS,
                        required_roles=("Read 1","Index 1 (i7)","Read 2")):
 raise ValueError("smMIP final library: "+"; ".join(problems))

FINAL_CAPTION="INFERRED — handle bases from the later Cantsilieris et al. 2017 MIP primers, not the unavailable 2013 Table S9. Single-index smMIP amplicon (2013 PAGE band ~280 bp). Read 2 encounters the twelve-base molecular tag first; the separate eight-base index read identifies the sample."
SEQUENCING_INTRO=("Hiatt et al. 2013 collected two 101-nt reads and one 8-nt sample-index read on HiSeq 2000 (shown), or two 152-nt reads and one 8-nt read on MiSeq, "
                  "using custom sequencing primers. It has no i5/Index 2 read. The first 12 nt of the reverse read are the molecular tag. "
                  "INFERRED — the primer bases drawn are from the later Cantsilieris et al. 2017 protocol.")


def capture_linear() -> Construct:
 """Symbolic 80-mer probe plus its copied 112-nt gap, before circle closure.

 The original panel uses 40 nt of target arms in total.  The 18+22 split is one
 permitted arm pair; the target-specific bases and 28-nt backbone are not asserted.
 """
 return Construct([
  Segment("ligation targeting arm", "X"*22, placeholder=True),
  Segment("12-nt molecular tag", "U"*12, "umi", placeholder=True, feature=MOLECULAR_TAG),
  Segment("invariant probe backbone", "N"*28, placeholder=True),
  Segment("extension targeting arm", "Y"*18, placeholder=True),
  Segment("copied target gap", "N"*CAPTURED_GAP_NT, placeholder=True),
 ], name="gap-filled smMIP before nick sealing")


def capture_circle():
 return circularize_ssdna(capture_linear(), five_prime_phosphate=True)


def smmip_rows():
 return circle_rows(capture_linear(), "gap fill + ligation")

def sections():
 return [
  ("Hybridize, fill and close one probe",smmip_rows(),
   "Both targeting arms must bind one genomic molecule. Extension across the gap and ligation convert only successful captures into closed circles."),
  ("Digest everything linear",[
   Row(chunks=[("unreacted probe + genomic DNA + incompletely closed products",None,False)]),
   Row(chunks=[("                         ↓ exonuclease",None,False)]),
   Row(chunks=[("closed smMIP circles survive as the selected template pool", "umi", False)]),
  ],"Circular topology, rather than bead affinity, performs the enrichment."),
  ("Amplify from the invariant backbone",[
   Row(chunks=[("circle—[12-nt tag]—captured gap → universal PCR → indexed library", "umi", False)]),
   Row(chunks=[("same tag across PCR descendants → single-molecule consensus", "cbc", False)]),
  ],"The molecular tag is already covalently attached before PCR, allowing family consensus and molecule counting."),
 ]


def protocol_graph():
 graph=ProtocolGraph(TITLE,authority=SOURCE)
 dna=graph.input("genomic_dna","genomic DNA containing target intervals",MaterialKind.BIOLOGICAL,state=molecule("genomic DNA",(role("target gaps flanked by probe-arm sites"),)),primary=True)
 unphosphorylated=graph.input("unphosphorylated_probe","12-nt-tagged molecular inversion probe pool",MaterialKind.OLIGO)
 kinase=graph.input("t4_pnk","T4 polynucleotide kinase (NEB, 50 U, 1× T4 DNA ligase buffer)",MaterialKind.CATALYST)
 probe=graph.product("probe","5′-phosphorylated molecular inversion probe",MaterialKind.MOLECULAR)
 graph.react("phosphorylate_probe","Phosphorylate the pooled probes (37°C 45 min; kinase heat-killed 80°C 20 min)",StepKind.PHOSPHORYLATION,
             ((unphosphorylated,InputRole.TEMPLATE),(kinase,InputRole.CATALYST)),(probe,))
 hybrid=graph.product("hybrid","probe arms hybridized around one target gap",MaterialKind.MOLECULAR,state=dna.state)
 graph.react("hybridize","Hybridize both smMIP targeting arms (500 ng gDNA, 330 fmol probe, 1× Ampligase buffer; 98°C 3 min, 85°C 30 min, 60°C 60 min, 56°C 120 min)",StepKind.HYBRIDIZATION,((dna,InputRole.TEMPLATE),(probe,InputRole.STRUCTURAL)),(hybrid,))
 polymerase=graph.input("gap_polymerase","Phusion DNA polymerase (3.2 U) with dNTPs and betaine",MaterialKind.CATALYST)
 filled=graph.product("filled","gap-filled smMIP capture intermediate",MaterialKind.MOLECULAR,
                      chemical_state=capture_linear())
 graph.react("gap_fill","Copy the 112-nt genomic target gap (combined gap-fill/ligation mix, 56°C 60 min then 72°C 20 min)",StepKind.EXTENSION,((hybrid,InputRole.TEMPLATE),(polymerase,InputRole.CATALYST)),(filled,))
 ligase=graph.input("circle_ligase","Ampligase thermostable DNA ligase (5 U) with NAD+",MaterialKind.CATALYST)
 circle=graph.product("circle","closed smMIP capture circle",MaterialKind.MOLECULAR,
                      state=tuple(smmip_rows()),chemical_state=capture_circle())
 graph.react("circularize","Seal the copied gap to close the probe",StepKind.CIRCULARIZATION,((filled,InputRole.TEMPLATE),(ligase,InputRole.CATALYST)),(circle,),geometry=CircularizationGeometry.known("gap-sealing probe closure",statement="DNA ligase seals the copied target gap between the smMIP probe arms."),relations=(MolecularRelation.CIRCLE_CLOSURE,))
 exonuclease=graph.input("exonuclease","Exonuclease I (20 U) + Exonuclease III (100 U)",MaterialKind.CATALYST)
 selected=graph.product("selected","exonuclease-resistant closed smMIP circles",MaterialKind.MOLECULAR,
                        state=tuple(smmip_rows()),chemical_state=capture_circle())
 graph.react("digest_linear","Digest uncircularized probe and genomic DNA (37°C 45 min; 80°C 20 min inactivation)",StepKind.DIGESTION,((circle,InputRole.TEMPLATE),(exonuclease,InputRole.CATALYST)),(selected,))
 forward=graph.input("universal_forward","invariant-backbone forward PCR primer",MaterialKind.OLIGO,
                     state=UNIVERSAL_FORWARD_PCR,source="Cantsilieris et al. 2017, Note 8 (later protocol; 2013 Table S9 unavailable)",inferred=True)
 reverse=graph.input("indexed_reverse","indexed invariant-backbone reverse PCR primer",MaterialKind.OLIGO,
                     state=INDEXED_REVERSE_PCR,source="Cantsilieris et al. 2017, Note 8 (later protocol; 2013 Table S9 unavailable)",inferred=True)
 pcr_pol=graph.input("pcr_polymerase","Phusion Hot-Start II polymerase (2 U; SYBR green)",MaterialKind.CATALYST)
 amplified=graph.product("amplified","indexed smMIP PCR product with low-molecular-weight artifact",MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
 graph.react("amplify","Universally amplify captured circles (two 50 µL PCRs per capture, 26 cycles: 98°C 15 s, 65°C 15 s, 72°C 45 s)",StepKind.PCR,((selected,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),(reverse,InputRole.STRUCTURAL),(pcr_pol,InputRole.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(FINAL_LIBRARY,reason="The circular template and panel-specific arms are not represented as one exact linear PCR template."))
 beads=graph.input("ampure_1p8x","Ampure XP beads, 1.8×",MaterialKind.SELECTION)
 purified=graph.product("purified","purified per-sample smMIP PCR product",MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
 bead_waste=graph.product("ampure_1p8x_waste","bead supernatant and washes",MaterialKind.WASTE)
 graph.react("purify","Purify each PCR product with 1.8× Ampure XP",StepKind.CLEANUP,((amplified,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(purified,bead_waste))
 graph.discard(bead_waste,"supernatant and ethanol washes removed")
 gel=graph.input("page_gel","6% polyacrylamide TBE gel (SYBR Gold), pooled samples",MaterialKind.SELECTION)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 artifact=graph.product("low_mw_artifact","low-molecular-weight PCR artifact left in the gel",MaterialKind.WASTE)
 graph.react("size_select","Pool and excise the ~280 bp capture-product band (crush/soak, Ampure recovery)",StepKind.SELECTION,((purified,InputRole.TEMPLATE),(gel,InputRole.SELECTION)),(library,artifact))
 graph.discard(artifact,"not excised with the ~280 bp band")
 graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),acquired_roles=("Read 1","Index 1 (i7)","Read 2"),
                read_lengths={"Read 1":101,"Index 1 (i7)":8,"Read 2":101}))
 graph.declare_reviewed_scope(["genomic_dna","unphosphorylated_probe","t4_pnk","phosphorylate_probe","probe","hybridize","gap_polymerase","gap_fill","circle_ligase","circularize","exonuclease","digest_linear","universal_forward","indexed_reverse","pcr_polymerase","amplify","amplified","ampure_1p8x","purify","purified","ampure_1p8x_waste","page_gel","size_select","low_mw_artifact","library"],SOURCE)
 return graph.require_publishable()
