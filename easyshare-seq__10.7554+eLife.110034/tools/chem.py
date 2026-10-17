"""easySHARE-seq dual-modality, three-part combinatorial indexing.

Authority: the published eLife article (Soltys et al., eLife 2026 version of record, eLife.110034,
PMC13384510) Methods, and its Supplementary file 1 (md5-identical to the bioRxiv
media-6 workbook; sheets "PCR oligos", "BC oligos R1", "BC oligos R2").
"""
import illumina as il
import nextera as nx
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Row, Scene, Workflow, feature
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SourceClaim,
                            SourceDiscrepancy, StepKind)
from protocol_model import molecule, role

TITLE="easySHARE-seq — joint RNA and chromatin accessibility"
NOTES="01_easyshare-seq.html"
SOURCE=('Defining source: <a href="https://doi.org/10.7554/eLife.110034">Soltys et al., <i>eLife</i> (2026, version of record)</a> '
        'and its Supplementary file 1; first posted as <a href="https://doi.org/10.1101/2024.02.26.581705">bioRxiv 2024.02.26.581705</a>.')
SUMMARY=("Fixed nuclei are tagmented with a mixed Tn5-A/Tn5-B2S transposome and reverse-transcribed "
         "with a biotinylated RT primer; both handles receive two 7-nt ligation barcodes. "
         "Streptavidin separates cDNA from ATAC fragments, and an N5XX PCR index completes each library.")
CAVEAT=("Representative barcodes Round1_001, Round2_001 and N501 are shown; each well or sub-library "
        "substitutes another published family member. Fragment-end classes and the 72 °C gap fill "
        "are derived from the published oligos, not drawn by the source.")
R1_001="GGTTAGCGTCTCGT"
R1_001_SCAFFOLD="CCGAGCCCACGAGACGCTAACCTCT"
R2_001="CAAGCAGAAGACGGCATACGAGATTCTGGCAAGA"
R2_001_SCAFFOLD="TGCCAGAATCTCGTATGCCGTCTTCTGCTTG"
TN5_A="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
TN5_REVERSE="CTGTCTCTTATACACATCT"
TN5_B2S="GGGCTCGGAGATGTGTATAAGAGACAG"
RT_PRIMER="GGGCTCGGAGATGTGTATAAGAGACAGNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTVN"  # T/iBiodT/T24
TSO="UUUTCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGG"
KRT_B2S="CAGGGCATTTGCACCCCATGCC"   # eLife Key resources table value (see discrepancy)
I7_TRUSEQ_LONG=il.P7
N501=nx.n5xx_primer("AGCGCTAG")
LINKER="TCT"  # Supplementary file 1: bases 8-10 of Index 1
READ_LENGTHS={"Read 1":150,"Index 1 (i7)":17,"Index 2 (i5)":8,"Read 2":150}

def _bc1(): return seg("round-1 barcode (Index 1 bases 1–7)","N"*7,"cbc",placeholder=True,feature=feature("cell_bc1","cell_barcode","whitelist",group="cell_barcode",part="round 1",whitelist="Supplementary file 1 BC oligos R1 (192)"))
def _bc2(): return seg("round-2 barcode (Index 1 bases 11–17)","N"*7,"cbc",placeholder=True,feature=feature("cell_bc2","cell_barcode","whitelist",group="cell_barcode",part="round 2",whitelist="Supplementary file 1 BC oligos R2 (192)"))
def _n5(): return seg("N5XX index (round 3 / sub-library)","AGCGCTAG","cbc",feature=feature("cell_bc3","cell_barcode","whitelist",group="cell_barcode",part="round 3",whitelist="Supplementary file 1 Nextera N5XX set"))
def _umi(): return seg("10-nt UMI","N"*10,"umi",placeholder=True,feature=feature("easyshare_umi","umi","random"))
def _barcode_end(): return [seg("ME + S7 (Tn5-B2S/RT handle + BC1 tail)",nx.INDEX1_PRIMER,"r2"),_bc1(),seg("TCT linker",LINKER),_bc2(),seg("P7 reverse complement",il.P7_RC,"p7")]
def _polya(): return [seg("VN anchor (NB on this strand)","NB",placeholder=True),seg("RT-primer T26 tail (A26 on this strand)","A"*26)]

def round1_geometry():
    """Round-1 duplex meets the 5′-phosphorylated GGGCTCGG overhang of a B2S/RT end."""
    sc=Scene()
    sc.strand("round-1 oligo",[seg("Round1_001 barcode bases",R1_001[:7],"cbc"),seg("Round1_001 S7 tail",R1_001[7:],"r2")],label="/5Phos/Round1_001")
    sc.anneal("round-1 scaffold",[seg("scaffold B2S/RT arm",R1_001_SCAFFOLD[:8]),seg("scaffold barcode arm",R1_001_SCAFFOLD[8:-3]),seg("TCT 3′ overhang",R1_001_SCAFFOLD[-3:])],
              to="round-1 oligo",pair=("scaffold barcode arm","Round1_001 barcode bases"),label="Round1_001 scaffold")
    sc.anneal("B2S end",[seg("B2S/RT 5′ arm",TN5_B2S[:8],"r2"),seg("ME",nx.ME,"me"),seg("insert or UMI-T tail","X"*10,placeholder=True)],
              to="round-1 scaffold",pair=("B2S/RT 5′ arm","scaffold B2S/RT arm"),above=True,label="/5Phos/ Tn5-B2S or RT-primer end",unpaired=("ME","insert or UMI-T tail"))
    sc.ligation_end("round-1 scaffold","scaffold B2S/RT arm",kind="cohesive",polarity="5-prime",owner="adapter",text="8-nt scaffold arm")
    sc.ligation_end("B2S end","B2S/RT 5′ arm",kind="cohesive",polarity="5-prime",owner="substrate",text="5′-phosphorylated GGGCTCGG")
    return AdapterLigationGeometry.known("cohesive",sc,statement="Methods 'Preparing the barcoding oligonucleotides': the round-1 duplex has two single-stranded linkers, one complementary to the 5′ end of the cDNA or transposed molecule (RT primer or Tn5 adapter). Appendix 1 says this initial overhang is 8 bp and that the Read 2 primer site is split between barcode oligo and Tn5/RT primer. Supplementary file 1 prints Round1_001 and its scaffold; Scene checks that the scaffold pairs with both the barcode oligo and the GGGCTCGG arm.")

def round2_geometry():
    sc=Scene()
    sc.strand("round-2 oligo",[seg("Round2_001 P7",R2_001[:24],"p7"),seg("Round2_001 barcode bases",R2_001[24:31],"cbc"),seg("AGA 3′ overhang",R2_001[31:])],label="Round2_001")
    sc.anneal("round-1 scaffold",[seg("scaffold B2S/RT arm",R1_001_SCAFFOLD[:8]),seg("scaffold barcode arm",R1_001_SCAFFOLD[8:-3]),seg("TCT 3′ overhang",R1_001_SCAFFOLD[-3:])],
              to="round-2 oligo",pair=("TCT 3′ overhang","AGA 3′ overhang"),label="Round1_001 scaffold (in ligated round-1 product)")
    sc.anneal("round-1 oligo",[seg("Round1_001 barcode bases",R1_001[:7],"cbc"),seg("Round1_001 S7 tail",R1_001[7:],"r2")],
              to="round-1 scaffold",pair=("Round1_001 barcode bases","scaffold barcode arm"),above=True,label="/5Phos/Round1_001")
    sc.anneal("round-2 scaffold",[seg("Round2_001 scaffold",R2_001_SCAFFOLD)],to="round-2 oligo",pair=("Round2_001 scaffold","Round2_001 P7"),label="Round2_001 scaffold")
    sc.overhang("round-2 oligo","AGA 3′ overhang",polarity="3-prime",owner="adapter",text="Round2_001 3′-AGA")
    sc.ligation_end("round-1 scaffold","TCT 3′ overhang",kind="cohesive",polarity="3-prime",owner="substrate",text="round-1 scaffold 3′-TCT")
    return AdapterLigationGeometry.known("cohesive",sc,statement="Methods: the second 3-bp overhang of the round-1 duplex ligates the round-2 barcodes; round-1 oligos are 5′-phosphorylated. Supplementary file 1 prints Round2_001 and its scaffold; Scene checks the AGA/TCT pairing.")

def _library(name,payload_parts):
    lib=Construct([seg("P5",il.P5,"p5"),_n5(),seg("Nextera Read 1 arm (S5 + ME)",nx.READ1_PRIMER,"r1"),*payload_parts,*_barcode_end()],name=name)
    primers=sp.standard_primers(lib)
    problems=sp.verify(lib,primers)
    if problems: raise ValueError(name+": "+"; ".join(problems))
    return lib,primers
RNA,RP=_library("easySHARE-seq RNA library",[seg("cDNA insert","X"*42,placeholder=True),*_polya(),_umi()])
ATAC,AP=_library("easySHARE-seq ATAC library",[seg("accessible chromatin insert","X"*42,placeholder=True)])
ATAC_TEMPLATE=Construct([seg("Tn5-A end (S5 + ME)",nx.READ1_PRIMER,"r1"),seg("accessible chromatin insert","X"*42,placeholder=True),*_barcode_end()],name="released ATAC fragment after 72 °C fill-in")
CDNA_TEMPLATE=Construct([seg("TSO-derived S5 + ME + GGG",TSO[3:].replace("U","T"),"r1"),seg("full-length cDNA","X"*42,placeholder=True),*_polya(),_umi(),*_barcode_end()],name="template-switched barcoded cDNA")
RNA_TEMPLATE=Construct([seg("Tn5-A end (S5 + ME)",nx.READ1_PRIMER,"r1"),seg("cDNA insert","X"*42,placeholder=True),*_polya(),_umi(),*_barcode_end()],name="Tn5-A-tagmented cDNA after 72 °C fill-in")
FINAL_LIBRARIES=(("RNA library",RNA,RP,"Tn5-A-fragmented cDNA: Read 2 reads the UMI and poly(T); Index 1 reads BC1–TCT–BC2.","Published run: Read 1 150, Index 1 17, Index 2 8, Read 2 150 cycles (NovaSeq 6000 S4)."),
                 ("ATAC library",ATAC,AP,"Tn5-A/Tn5-B2S fragment carrying the same three-part cellular identity.","Published run: Read 1 150, Index 1 17, Index 2 8, Read 2 150 cycles (NovaSeq 6000 S4)."))
def _o(name,seq,mod=""): return f"<b>{name}</b> — 5′-{mod}{seq}-3′"
def oligos(): return [
    _o("Tn5-A",TN5_A),
    _o("Tn5-reverse",TN5_REVERSE,"/5Phos/"),
    _o("Tn5-B2S (Supplementary file 1)",TN5_B2S,"/5Phos/"),
    "<b>Reverse-transcription primer</b> — 5′-/5Phos/GGGCTCGGAGATGTGTATAAGAGACAGNNNNNNNNNNT/iBiodT/TTTTTTTTTTTTTTTTTTTTTTTTVN-3′",
    "<b>Template-switch oligo</b> — 5′-rUrUrUTCGTCGGCAGCGTCAGATGTGTATAAGAGACArGrGrG-3′",
    _o("Round1_001 barcode oligo",R1_001,"/5Phos/")+" with scaffold 5′-"+R1_001_SCAFFOLD+"-3′",
    _o("Round2_001 barcode oligo",R2_001)+" with scaffold 5′-"+R2_001_SCAFFOLD+"-3′",
    _o("i7-Tru-Seq-long primer",I7_TRUSEQ_LONG),
    _o("Nextera N501 (representative N5XX)",N501),
]
def sections(): return [
 ("Tag accessible DNA and copy RNA in fixed nuclei",[Row(chunks=[("chromatin ** Tn5-A / Tn5-B2S     RNA → biotin RT primer","me",False)])],"Both handles end in a 5′-phosphorylated GGGCTCGG arm."),
 ("Ligate round-1 and round-2 barcodes",round1_geometry().rows()+[Row()]+round2_geometry().rows(),"BC1–TCT–BC2 becomes Index 1."),
]

def workflow():
    w=Workflow(MolecularState("RNA and ATAC handles",tuple(sections()[0][1])))
    w.ligate_adapter("Ligate round-1 barcodes",round1_geometry(),[Row(chunks=[("BC1 ** B2S/RT handle","cbc",False)])],note="Round-1 duplexes ligate to the handle ends.")
    w.ligate_adapter("Ligate round-2 barcodes",round2_geometry(),[Row(chunks=[("BC2 ** TCT ** BC1 ** handle","cbc",False)])],note="Round-2 duplexes ligate through the TCT/AGA overhang.")
    return w

def _exact(name,seq): return Construct([seg(name,seq)],name=name)

def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    roles=lambda label,*r: molecule(label,tuple(role(x) for x in r))
    nuclei = g.input("nuclei", "formaldehyde-fixed nuclei or cells (0.35% liver, 0.1% cell lines)",
        MaterialKind.BIOLOGICAL, state=roles("fixed nuclei","poly(A) RNA","chromatin"), primary=True)
    tn5_a = g.input("tn5_a_oligo","Tn5-A oligo",MaterialKind.OLIGO,state=_exact("Tn5-A",TN5_A))
    tn5_b = g.input("tn5_b2s_oligo","/5Phos/ Tn5-B2S oligo (Supplementary file 1)",MaterialKind.OLIGO,state=_exact("Tn5-B2S",TN5_B2S))
    tn5_r = g.input("tn5_reverse","/5Phos/ Tn5-reverse oligo",MaterialKind.OLIGO,state=_exact("Tn5-reverse",TN5_REVERSE))
    tn5 = g.input("tn5", "in-house Tn5 loaded with pooled A and B2S duplexes", MaterialKind.CATALYST)
    tagged = g.product("tagged", "nuclei with Tn5-A/B2S-tagged chromatin fragments",
        MaterialKind.MOLECULAR, state=roles("tagged nuclei","Tn5-A or B2S end","accessible DNA","Tn5-A or B2S end","poly(A) RNA"))
    g.react("tagment", "Tagment 10,000 nuclei per reaction, 37 °C 30 min", StepKind.TAGMENTATION,
        ((nuclei, InputRole.TEMPLATE), (tn5_a, InputRole.STRUCTURAL), (tn5_b, InputRole.STRUCTURAL),
         (tn5_r, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)), (tagged,))
    # Printed /5Phos/...NNNNNNNNNNT/iBiodT/T24VN: the biotin-dT is one of the 26 T.
    rt_primer = g.input("rt_primer", "/5Phos/ biotin-dT RT primer with 10-nt UMI", MaterialKind.OLIGO,
        state=Construct([seg("B2S/RT handle",RT_PRIMER[:27]),_umi(),
                         seg("T/iBiodT/T24 (26 T)",RT_PRIMER[37:63]),seg("VN anchor",RT_PRIMER[63:],placeholder=True)],
                        name="RT primer"))
    rt = g.input("rt", "Maxima H Minus reverse transcriptase (with PEG)", MaterialKind.CATALYST)
    copied = g.product("copied", "nuclei with tagged chromatin and biotinylated first-strand cDNA",
        MaterialKind.MOLECULAR, state=roles("tagged and copied nuclei","B2S/RT handle","cDNA:RNA","Tn5-tagged DNA"))
    g.react("reverse_transcribe", "In-nucleus RT, 52 °C 12 min, two low-temperature ramps, 52 °C 5 min",
        StepKind.REVERSE_TRANSCRIPTION, ((tagged, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
         (rt, InputRole.CATALYST)), (copied,))
    r1 = g.input("round1_barcodes", "annealed round-1 barcode/scaffold duplexes (192, two plates)",
        MaterialKind.ADAPTER, state=round1_geometry().scenes[0])
    ligase = g.input("ligase", "T4 DNA ligase", MaterialKind.CATALYST)
    round1 = g.product("round1", "handles carrying round-1 barcodes", MaterialKind.MOLECULAR,
        state=roles("round-1 product","BC1","B2S/RT handle","insert"))
    g.react("ligate_round1", "Round-1 split ligation, 25 °C 30 min", StepKind.LIGATION,
        ((copied, InputRole.TEMPLATE), (r1, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)),
        (round1,), geometry=round1_geometry())
    pool = g.input("pool_condition", "pool, wash, resuspend in fresh ligation mix and redistribute", MaterialKind.CONDITION)
    pooled = g.product("pooled", "round-1 nuclei redistributed into round-2 plates", MaterialKind.MOLECULAR, state=round1.state)
    lost = g.product("redistribution_losses", "nuclei lost in pooling/washing", MaterialKind.WASTE)
    g.react("redistribute", "Pool and redistribute nuclei", StepKind.PARTITION,
        ((round1, InputRole.TEMPLATE), (pool, InputRole.CONDITION)), (pooled, lost), changes_molecule=False)
    g.discard(lost, "removed with supernatants or left in plates")
    r2 = g.input("round2_barcodes", "annealed round-2 barcode/scaffold duplexes (192, two plates)",
        MaterialKind.ADAPTER, state=round2_geometry().scenes[0])
    indexed = g.product("indexed", "RNA and ATAC molecules carrying BC2–TCT–BC1", MaterialKind.MOLECULAR,
        state=roles("round-2 product","P7","BC2","TCT","BC1","B2S/RT handle","insert"))
    g.react("ligate_round2", "Round-2 split ligation, 25 °C 45 min", StepKind.LIGATION,
        ((pooled, InputRole.TEMPLATE), (r2, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)),
        (indexed,), geometry=round2_geometry())
    reversal = g.input("reversal_condition", "3500-cell sub-library; SDS/NaCl, proteinase K 62 °C 2 h, 80 °C 10 min, Tween-20 quench",
        MaterialKind.CONDITION)
    released = g.product("released", "released barcoded cDNA and ATAC molecules", MaterialKind.MOLECULAR, state=indexed.state)
    g.react("reverse_crosslinks", "Split into sub-libraries and reverse crosslinks",
        StepKind.DEPROTEINIZATION, ((indexed, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)), (released,))
    beads = g.input("streptavidin", "M280 streptavidin beads", MaterialKind.SELECTION)
    bound = g.product("cdna_beads", "biotinylated cDNA on beads", MaterialKind.MOLECULAR,
                      state=roles("bead-bound cDNA","barcode end","cDNA:RNA"))
    supernatant = g.product("atac_supernatant", "transposed DNA in the supernatant", MaterialKind.MOLECULAR,
                            state=roles("ATAC supernatant","Tn5 end","accessible DNA","barcode end"))
    g.react("streptavidin_split", "Streptavidin pull-down separates cDNA from ATAC fragments",
        StepKind.AFFINITY_CAPTURE, ((released, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (bound, supernatant))
    minelute = g.input("minelute", "QIAGEN MinElute cleanup", MaterialKind.SELECTION)
    atac = g.product("atac_branch", "purified ATAC fragments", MaterialKind.MOLECULAR, state=ATAC_TEMPLATE)
    atac_waste = g.product("atac_cleanup_waste", "MinElute flow-through", MaterialKind.WASTE)
    g.react("atac_cleanup", "MinElute cleanup of the supernatant", StepKind.CLEANUP,
        ((supernatant, InputRole.TEMPLATE), (minelute, InputRole.SELECTION)), (atac, atac_waste))
    g.discard(atac_waste, "column flow-through")
    tso = g.input("tso", "template-switch oligo (rU×3 … rG×3)", MaterialKind.OLIGO, state=TSO)
    rt2 = g.input("ts_rt", "Maxima H Minus reverse transcriptase (template switch)", MaterialKind.CATALYST)
    switched = g.product("switched", "template-switched cDNA on beads", MaterialKind.MOLECULAR, state=CDNA_TEMPLATE)
    g.react("template_switch", "On-bead template switching, 25 °C 30 min then 42 °C 90 min",
        StepKind.TEMPLATE_SWITCHING, ((bound, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL), (rt2, InputRole.CATALYST)), (switched,))
    n5_c = g.input("cdna_n5", "Nextera N5XX indexing primer (cDNA PCR)", MaterialKind.OLIGO, state=_exact("Nextera N501",N501))
    i7_c = g.input("cdna_i7", "i7-Tru-Seq-long primer (cDNA PCR)", MaterialKind.OLIGO, state=_exact("i7-Tru-Seq-long",I7_TRUSEQ_LONG))
    q5 = g.input("cdna_q5", "Q5 polymerase", MaterialKind.CATALYST)
    cdna = g.product("cdna_amplified", "amplified full-length barcoded cDNA (mean >700 bp)", MaterialKind.MOLECULAR,
        state=roles("amplified cDNA","P5 + N5XX","S5 + ME + GGG","full-length cDNA","poly(A)","UMI","barcode end"))
    g.react("cdna_pcr", "cDNA PCR, 13–15 cycles; 0.7× bead cleanup", StepKind.PCR,
        ((switched, InputRole.TEMPLATE), (n5_c, InputRole.STRUCTURAL), (i7_c, InputRole.STRUCTURAL), (q5, InputRole.CATALYST)), (cdna,),
        geometry=PCRGeometry.exact(CDNA_TEMPLATE, PCRPrimer("Nextera N501", N501, anneal_nt=14),
            PCRPrimer("i7-Tru-Seq-long", I7_TRUSEQ_LONG),
            statement="N5XX binds the TSO-derived S5 arm; i7-Tru-Seq-long binds the round-2 P7 arm."))
    a_duplex = g.input("tn5_a_only_adapter", "Tn5-A/Tn5-reverse duplex (10 µM / 10.5 µM)", MaterialKind.ADAPTER)
    tn5_a_only = g.input("tn5_a_only", "in-house Tn5 loaded only with the Tn5-A duplex", MaterialKind.CATALYST)
    rna = g.product("rna_branch", "Tn5-A-fragmented cDNA retaining the barcode end", MaterialKind.MOLECULAR, state=RNA_TEMPLATE)
    rna_other = g.product("rna_tn5_other", "candidate fragments lacking the barcode end", MaterialKind.WASTE)
    g.react("tagment_cdna", "Tagment 25 ng cDNA with Tn5-A-only, 55 °C 10 min; SDS stop; 1.3× beads",
        StepKind.TAGMENTATION, ((cdna, InputRole.TEMPLATE), (a_duplex, InputRole.STRUCTURAL), (tn5_a_only, InputRole.CATALYST)), (rna, rna_other))
    g.discard(rna_other, "candidate only: no P7 barcode end for the i7 primer; not separated physically")
    reviewed = [
        "nuclei", "tn5_a_oligo", "tn5_b2s_oligo", "tn5_reverse", "tn5", "tagment",
        "rt_primer", "rt", "reverse_transcribe", "round1_barcodes", "ligase", "ligate_round1",
        "pool_condition", "redistribute", "round2_barcodes", "ligate_round2",
        "reversal_condition", "reverse_crosslinks", "streptavidin", "streptavidin_split",
        "minelute", "atac_cleanup", "tso", "ts_rt", "template_switch",
        "cdna_n5", "cdna_i7", "cdna_q5", "cdna_pcr", "tn5_a_only_adapter", "tn5_a_only", "tagment_cdna",
    ]
    finish_pcr_branch(
        g, "rna", rna, RNA, RP, reviewed, label="easySHARE-seq RNA", read_lengths=READ_LENGTHS,
        geometry=PCRGeometry.exact(
            RNA_TEMPLATE, PCRPrimer("Nextera N501", N501, anneal_nt=14),
            PCRPrimer("i7-Tru-Seq-long", I7_TRUSEQ_LONG),
            statement="Methods: 72 °C 6 min then 5 cycles with the same N5XX as the cDNA PCR and i7-Tru-Seq-long; the 72 °C step is taken as Tn5 gap fill (derived)."))
    finish_pcr_branch(
        g, "atac", atac, ATAC, AP, reviewed, label="easySHARE-seq ATAC", read_lengths=READ_LENGTHS,
        geometry=PCRGeometry.exact(
            ATAC_TEMPLATE, PCRPrimer("Nextera N501", N501, anneal_nt=14),
            PCRPrimer("i7-Tru-Seq-long", I7_TRUSEQ_LONG),
            statement="Methods: 72 °C 6 min, then 15–17 cycles with an N5XX primer and i7-Tru-Seq-long; only Tn5-A/barcode-end fragments have both sites (derived)."))
    g.annotate(SourceDiscrepancy(
        "tn5_b2s_oligo", "Tn5-B2S sequence differs between the eLife Key resources table and Supplementary file 1",
        (SourceClaim("eLife Methods, Key resources table (eLife v1 XML line 257; PMC13384510 XML lines 286–289)", KRT_B2S),
         SourceClaim("eLife Supplementary file 1, sheet 'PCR oligos' row Tn5-B2S", "/5Phos/"+TN5_B2S)),
        "The model uses the Supplementary file 1 sequence.",
        "The table value has no Tn5 mosaic end and does not match the RT primer's shared handle or the round-1 scaffold arm; the workbook value does both.",
        "The origin of the table value is not stated; no correction notice was found in this bounded pass."))
    g.declare_reviewed_scope(reviewed, SOURCE)
    return g.require_publishable()
