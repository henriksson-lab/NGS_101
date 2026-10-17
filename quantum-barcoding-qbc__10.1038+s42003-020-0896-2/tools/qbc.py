"""Quantum Barcoding (QBC) molecular model from the defining supplement."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from adapter_ligation import reported_adapter_end
from chemdraw import Construct, Scene, Segment, feature
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan, SourceClaim,
                            SourceDiscrepancy, StepKind)
from topology import InterruptedSplintTopology, SplintProductSegment

TITLE="Quantum Barcoding (QBC) — targeted RNA and protein indexing"
SOURCE="Roy et al., Communications Biology (2020), doi:10.1038/s42003-020-0896-2"
def seg(name,seq,tag=None,**kw): return Segment(name,seq,tag,**kw)

LEFT_PARTIAL_R1="CCCTACACGACGCTCTTCCGATCT"
SC3_RIGHT="GCTCACCAAGTCCTAGT"; SC2_LEFT="GCGTTACCC"; SC2_RIGHT="CCCATTCTC"
SCN_LEFT="GTCTTTGCC"; SCN_RIGHT="CCGTTGTGG"; SC1_LEFT="GCGACTTCC"; SC1_RIGHT="CCAACCACC"
ANCHOR="GCTCCCTGTCTGACG"
RIGHT_ARM="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAACTCCATCAGC"   # AHCA 3' arm (Supplementary Data 1/2)
# Methods 'Primer design': RNA-specific 2nd-strand primers carry this Illumina
# annealing sequence on their 5' end; the RNA amplicon therefore ends in its
# reverse complement, without the AHCA-only GAACTCCATCAGC (La4FB-binding) tail.
RNA_SECOND_STRAND_TAIL="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
def _rc(x): return x[::-1].translate(str.maketrans("ACGT","TGCA"))
RNA_RIGHT_ARM=_rc(RNA_SECOND_STRAND_TAIL)
SPLINT6="CGTCAGACAGGGAGCGGTGGTTGGGGAAGTCGCCCACAACGGGGCAAAGACGAGAATGGGGGGTAACGCACTAGGACTTGGTGAGC"

UMI_F=feature("qbc_umi","umi","random")
def bc_feature(n): return feature(f"qbc_cell_bc_{n}","cell_barcode","whitelist",group="cell_barcode",part=f"round {n}",whitelist="QBC 60-well subcode set")
AHCA_F=feature("qbc_antibody_code","feature_barcode","whitelist",whitelist="QBC AHCA set")

def indexed_core(protein=True) -> Construct:
    parts=[seg("10-nt UMI","N"*10,"umi",placeholder=True,feature=UMI_F),
      seg("SC3 cell subcode","B"*7,"cbc",placeholder=True,feature=bc_feature(4)),seg("SC3 right arm",SC3_RIGHT),
      seg("SC2 left arm",SC2_LEFT),seg("SC2 cell subcode","C"*7,"cbc",placeholder=True,feature=bc_feature(3)),seg("SC2 right arm",SC2_RIGHT),
      seg("SCN left arm",SCN_LEFT),seg("SCN cell subcode","D"*7,"cbc",placeholder=True,feature=bc_feature(2)),seg("SCN right arm",SCN_RIGHT),
      seg("SC1 left arm",SC1_LEFT),seg("SC1 cell subcode","E"*7,"cbc",placeholder=True,feature=bc_feature(1)),seg("SC1 right arm",SC1_RIGHT),
      seg("QBC anchor",ANCHOR)]
    if protein: parts.append(seg("9-nt antibody-handle code","H"*9,"cbc",placeholder=True,feature=AHCA_F))
    else: parts.append(seg("target-specific cDNA","X"*40,placeholder=True))
    return Construct(parts, "QBC antibody reporter" if protein else "QBC targeted-RNA reporter")

def final_library(protein=True) -> Construct:
    payload=list(indexed_core(protein))
    return Construct([seg("P5",il.P5,"p5"),seg("i5 sample index","I"*8,"cbc",placeholder=True,
      feature=feature("qbc_i5","sample_index","random")),seg("Read 1 site",il.TRUSEQ_READ1,"r1"),
      *payload,seg("Read 2 site reverse complement",il.TRUSEQ_READ2[::-1].translate(str.maketrans("ACGT","TGCA")),"r2"),
      seg("i7 sample index reverse complement","J"*8,"cbc",placeholder=True,
          feature=feature("qbc_i7","sample_index","random")),seg("P7 reverse complement",il.P7_RC,"p7")],
      "QBC antibody Illumina library" if protein else "QBC targeted-RNA Illumina library")

def pcr_template(protein=True) -> Construct:
    right=(seg("AHCA Read 2-side arm",RIGHT_ARM,"r2") if protein else
           seg("second-strand primer Illumina arm (reverse complement)",RNA_RIGHT_ARM,"r2"))
    return Construct([seg("partial Read 1 arm",LEFT_PARTIAL_R1,"r1"),
                      *list(indexed_core(protein)),
                      right],
                     "QBC antibody pre-PCR template" if protein else "QBC RNA pre-PCR template")

PROTEIN_LIBRARY=final_library(True); RNA_LIBRARY=final_library(False)
PRIMERS=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["I2"],sp.TRUSEQ["R2"])
for _lib in (PROTEIN_LIBRARY,RNA_LIBRARY):
    _problems=sp.verify(_lib,PRIMERS)
    if _problems: raise ValueError(_lib.name+": "+"; ".join(_problems))

def splint_scene() -> Scene:
    # Chunk identities follow the docking arms (first 24 nt pair with the anchor
    # and SC1 right arm, etc.); the Methods ordering string ends in /3SpC3/.
    sc=Scene(); sc.strand("Splint6",[seg("anchor / SC1 docking",SPLINT6[:24]),seg("C3 spacer 1","[C3]",placeholder=True),
      seg("SC1 / SCN docking",SPLINT6[24:42]),seg("C3 spacer 2","[C3]",placeholder=True),
      seg("SCN / SC2 docking",SPLINT6[42:60]),seg("C3 spacer 3","[C3]",placeholder=True),
      seg("SC2 / SC3 docking",SPLINT6[60:]),seg("3′ C3 spacer","[3SpC3]",placeholder=True)],label="Splint6, written 5′→3′")
    sc.note("Splint6","one scaffold co-localizes the four successive barcode junctions")
    return sc


SPLINT6_CHUNKS = (
    ("anchor / SC1 chunk", SPLINT6[:24]),
    ("SC1 / SCN chunk", SPLINT6[24:42]),
    ("SCN / SC2 chunk", SPLINT6[42:60]),
    ("SC2 / SC3 chunk", SPLINT6[60:]),
)
SPLINT6_TOPOLOGY = InterruptedSplintTopology(
    "QBC four-round assembly on Splint6",
    SPLINT6_CHUNKS,
    (
        SplintProductSegment("SC3 PCR handle, UMI and cell subcode",
                             LEFT_PARTIAL_R1 + "N" * 10 + "N" * 7),
        SplintProductSegment("SC3 right arm", SC3_RIGHT, docking=True),
        SplintProductSegment("SC2 left arm", SC2_LEFT, docking=True,
                             junction_before=True),
        SplintProductSegment("SC2 cell subcode", "N" * 7),
        SplintProductSegment("SC2 right arm", SC2_RIGHT, docking=True),
        SplintProductSegment("SCN left arm", SCN_LEFT, docking=True,
                             junction_before=True),
        SplintProductSegment("SCN cell subcode", "N" * 7),
        SplintProductSegment("SCN right arm", SCN_RIGHT, docking=True),
        SplintProductSegment("SC1 left arm", SC1_LEFT, docking=True,
                             junction_before=True),
        SplintProductSegment("SC1 cell subcode", "N" * 7),
        SplintProductSegment("SC1 right arm", SC1_RIGHT, docking=True),
        SplintProductSegment("reporter anchor", ANCHOR, docking=True,
                             junction_before=True),
    ),
)

# Supplementary Data 1 "Generic Versions of Sequences to Order": whole ordered strands,
# with N at the printed variable positions (x7 subcode, x9 AHCA code, N10 UID).
def subcode_oligo(round_no) -> Construct:
    name,body,end=SUBCODE_BODIES[round_no]
    return Construct([seg(f"{name} 5′ part",body),seg(f"{name} 3′ arm",end)],
                     f"{name} generic subcode oligo")
AHCA_OLIGO=Construct([seg("AHCA anchor",ANCHOR),seg("9-nt AHCA code","N"*9),
                      seg("AHCA Read 2-side arm",RIGHT_ARM)],"AHCA generic oligo")

SUBCODE_BODIES={1:("SC1",SC1_LEFT+"N"*7,SC1_RIGHT),2:("SCN",SCN_LEFT+"N"*7,SCN_RIGHT),
                3:("SC2",SC2_LEFT+"N"*7,SC2_RIGHT),4:("SC3",LEFT_PARTIAL_R1+"N"*17,SC3_RIGHT)}
def ligation_geometry(round_no):
    name,body,end=SUBCODE_BODIES[round_no]
    return reported_adapter_end(f"QBC round-{round_no} {name} subcode",kind="splinted",
      polarity="3-prime",terminal_sequence=end,body_sequence=body,
      statement=(f"Supplementary Data 1 'Generic Versions of Sequences to Order' prints {name} with a 3′-OH; "
                 "its 3′ arm is splinted by Splint6 onto the 5′-phosphorylated end of the growing reporter "
                 "(anchor of the phosphorylated RT primer/AHCA oligo, then the previous subcode). "
                 + ("SC1, SCN and SC2 are 5′-phosphorylated." if round_no<4 else
                    "SC3 is printed without a 5′ phosphate and carries the partial Read 1 PCR handle.")),
      substrate_scenes=(splint_scene(),))

def pcr_geometry(template):
    f=PCRPrimer("NEBNext universal primer",il.NEBNEXT_UNIVERSAL_PRIMER,anneal_nt=24)
    r=PCRPrimer("NEBNext indexed primer",il.nebnext_i7_primer("N"*8),anneal_nt=len(il.TRUSEQ_READ2))
    return PCRGeometry.exact(template,f,r,statement="Methods name NEBNext E7600S primers for the targeted-RNA PCR and E7335 for protein-only PCR; Supplementary Data 1 and 'Primer design' print the internal primer-binding arms. Vendor index layouts were not inspected; both branches are drawn with the dual-index form.")

def protocol_graph() -> ProtocolGraph:
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cells=g.input("fixed_cells","fixed cells containing RNA and antibody-bound proteins",MaterialKind.BIOLOGICAL,primary=True)
    rtprimer=g.input("target_rt_primers","5′-phosphorylated anchor–target-specific RT primers (Supplementary Data 5–7)",MaterialKind.OLIGO)
    antibodies=g.input("dna_antibodies","La4FB-conjugated antibodies pre-hybridized to 5′-phosphorylated AHCA oligos",MaterialKind.OLIGO,state=AHCA_OLIGO)
    rt=g.input("reverse_transcriptase","ProtoScript II reverse transcriptase with aminoallyl-dUTP",MaterialKind.CATALYST)
    reporters=g.product("reporters","cells carrying targeted cDNA and antibody DNA reporters",MaterialKind.MOLECULAR,state=(indexed_core(False),indexed_core(True)))
    g.react("capture","In-cell RT of selected RNAs, DTSSP crosslinking of aminoallyl cDNA, then antibody staining and refixation",StepKind.REVERSE_TRANSCRIPTION,
      ((cells,InputRole.TEMPLATE),(rtprimer,InputRole.STRUCTURAL),(antibodies,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(reporters,))
    current=reporters; reviewed=["fixed_cells","target_rt_primers","dna_antibodies","reverse_transcriptase","capture"]
    ligase=g.input("ligase","T4 DNA ligase",MaterialKind.CATALYST); splint=g.input("splint6","Splint6 scaffold",MaterialKind.OLIGO,state=splint_scene())
    reviewed += ["ligase","splint6"]
    for n,name in enumerate(("SC1","SCN","SC2","SC3"),1):
      ol=g.input(f"round_{n}_subcodes",f"{name} 60-member subcode plate ("+("5′-phosphorylated" if name!="SC3" else "no 5′ phosphate; partial Read 1 handle")+")",MaterialKind.ADAPTER,state=subcode_oligo(n))
      out=g.product(f"after_round_{n}",f"reporters after {name} split-pool ligation",MaterialKind.MOLECULAR,state=(indexed_core(False),indexed_core(True)))
      g.react(f"ligate_round_{n}",f"Split into 60 wells, ligate {name}, wash and pool",StepKind.LIGATION,
        ((current,InputRole.TEMPLATE),(ol,InputRole.STRUCTURAL),(splint,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(out,),geometry=ligation_geometry(n))
      current=out; reviewed += [f"round_{n}_subcodes",f"ligate_round_{n}"]
    rna_ss=g.product("rna_released","barcoded single-stranded targeted cDNA",MaterialKind.MOLECULAR,state=indexed_core(False))
    ab_t=g.product("protein_template","barcoded antibody-reporter PCR template",MaterialKind.MOLECULAR,state=pcr_template(True))
    g.react("separate_branches","Lyse cells (RNase H for cDNA samples, proteinase K) and follow targeted-RNA and antibody reporters as separate amplicon classes",StepKind.SPLIT,((current,InputRole.TEMPLATE),),(rna_ss,ab_t),changes_molecule=False)
    ssp=g.input("rna_second_strand_primers","RNA-specific second-strand primers with 5′ Illumina annealing sequence (Supplementary Data 8–10)",MaterialKind.OLIGO)
    hc=g.input("herculase","Herculase II Fusion DNA polymerase",MaterialKind.CATALYST)
    rna_t=g.product("rna_template","barcoded targeted-RNA PCR template",MaterialKind.MOLECULAR,state=pcr_template(False))
    g.react("rna_second_strand","Second-strand synthesis of targeted cDNA (1–20 cycles; optional thermolabile Exo I)",StepKind.EXTENSION,
      ((rna_ss,InputRole.TEMPLATE),(ssp,InputRole.STRUCTURAL),(hc,InputRole.CATALYST)),(rna_t,))
    reviewed += ["separate_branches","rna_second_strand_primers","herculase","rna_second_strand"]
    fp=g.input("universal_pcr_primer","NEBNext universal PCR primer",MaterialKind.OLIGO)
    rp=g.input("indexed_pcr_primer","NEBNext indexed PCR primer",MaterialKind.OLIGO)
    reviewed += ["universal_pcr_primer","indexed_pcr_primer"]
    for key,temp,lib in (("rna",rna_t,RNA_LIBRARY),("protein",ab_t,PROTEIN_LIBRARY)):
      out=g.product(f"{key}_library",lib.name,MaterialKind.LIBRARY,state=lib)
      g.react(f"{key}_pcr",f"Indexed library PCR for the {key} branch",StepKind.PCR,
        ((temp,InputRole.TEMPLATE),(fp,InputRole.STRUCTURAL),(rp,InputRole.STRUCTURAL)),(out,),geometry=pcr_geometry(pcr_template(key=="protein")))
      roles=tuple(p.role for p in PRIMERS)
      g.sequence(out,SequencingPlan(library=lib,primers=PRIMERS,acquired_roles=roles,
        read_uses={"Read 1":"include","Index 1 (i7)":"extract","Index 2 (i5)":"extract","Read 2":"include"}))
      reviewed.extend((f"{key}_pcr",f"{key}_library"))
    g.annotate(SourceDiscrepancy(
        "dna_antibodies","Published AHCA oligo anchor differs by one 5′ base",
        (SourceClaim("Methods, 'Cell staining with oligonucleotide-conjugated antibodies' (PMC7205613 JATS line 234)",
                     "5′-(phos)-CTCCCTGTCTGACG(xxxxxxxxx)AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAACTCCATCAGC-3′"),
         SourceClaim("Supplementary Data 1 (AHCA row) and Supplementary Data 2 (all 60 AHCA sequences); Methods 'Primer design' anchor",
                     "GCTCCCTGTCTGACG xxxxxxxxx AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC GAACTCCATCAGC")),
        "The model uses the 15-nt anchor GCTCCCTGTCTGACG for the AHCA oligo.",
        "Supplementary Data 1/2 and the RT-primer anchor in 'Primer design' print the 5′ G; the Splint6 5′ docking arm CGTCAGACAGGGAGC is its exact reverse complement, so the 14-nt Methods form would leave a one-base gap at the SC1 junction (derived).",
        "Which form was ordered is not stated beyond these printed strings."))
    g.annotate(SourceDiscrepancy(
        "splint6","Published Splint6 5′ modification differs",
        (SourceClaim("Methods, 'Quantum barcoding by split-pool synthesis' (PMC7205613 JATS line 238)",
                     "5′-CGTCAGACAGGGAGCGGTGGTTGG/iSpC3/GGAAGTCGCCCACAACGG/iSpC3/GGCAAAGACGAGAATGGG-/iSpC3/GGGTAACGCACTAGGACTTGGTGAGC/3SpC3/-3′ (no 5′ modification)"),
         SourceClaim("Supplementary Data 1, 'Generic Versions of Sequences to Order'",
                     "p-CGTCAGACAGGGAGCGGTGGTTGG/iSpC3/GGAAGTCGCCCACAACGG/iSpC3/GGCAAAGAC GAGAATGGG/iSpC3/GGGTAACGC ACTAGGACTTGGTGAGC/iSpC3/")),
        "Splint6 is drawn with its 86 bases, three internal C3 spacers and a 3′ C3 spacer; no 5′ phosphate is drawn.",
        "Bases and spacer positions agree; only the 5′ phosphate differs between the two printed forms.",
        "Whether the ordered splint carried a 5′ phosphate, and whether it matters to any reaction, is not stated."))
    g.declare_reviewed_scope(tuple(reviewed),SOURCE)
    return g.require_publishable()
