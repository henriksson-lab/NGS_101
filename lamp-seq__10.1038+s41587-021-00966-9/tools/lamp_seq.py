"""Source-backed models for named isothermal-amplification sequencing protocols."""
from __future__ import annotations

import illumina as il
from chemdraw import Construct, Scene, Segment, feature, oligo, revcomp, complement_segments
from pcr import PCRGeometry
from protocol_graph import (InputRole as I, MaterialKind as M, ProtocolGraph,
                            SequencingPlan, StepKind as K, Terminal, Disposition)
from protocol_model import role
import seqprimers as sp


LAMP_F3 = "A+AC+AC+AA+GC+TTTCGGCAG"
LAMP_B3 = "G+AA+AT+TT+GG+ATCTTTGTCATCC"
LAMP_FIP_LEFT = "TGCGGCCAATGTTTGTAATCAG"
LAMP_FIP_BARCODE = "AGCAAACTAC"  # first validated Tier-1 barcode
LAMP_FIP_RIGHT = "CCAAGGAAATTTTGGGGAC"
LAMP_FIP = LAMP_FIP_LEFT + LAMP_FIP_BARCODE + LAMP_FIP_RIGHT
LAMP_BIP = "CGCATTGGCATGGAAGTCACTTTGATGGCACCTGTGTAG"
LAMP_LF = "TTCCTTGTCTGATTAGTTC"
LAMP_LB = "ACCTTCGGGAACGTGGTT"
LAMP_PCR_FWD_SITE = "AACGCTGAAGCGCTGGGGGCAAA"
LAMP_PCR_REV_SITE = "GTTTGTAATCAGTTCCTTGTCTG"
LAMP_PCR_FWD = il.TRUSEQ_READ1 + LAMP_PCR_FWD_SITE
# The source prints this first-stage tail with TGACT rather than canonical GTGACT.
LAMP_PCR_REV_ARM = "TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
LAMP_PCR_REV = LAMP_PCR_REV_ARM + LAMP_PCR_REV_SITE
LAMP_POOL_I5 = "TCGTGGAGCG"
LAMP_POOL_I7 = "CGCTCAGTTC"
# Supplementary Data 1 "Pool barcoding primers" (row Pool_fwd_1/Pool_rev_1) and the
# Methods two-step PCR print 23- and 24-nt partial sequencing-primer tails, not full
# TruSeq arms: the first PCR supplies the remaining Read 1 / Read 2 bases.
LAMP_POOL_FWD_TAIL = il.TRUSEQ_READ1[:23]   # ACACTCTTTCCCTACACGACGCT
LAMP_POOL_REV_TAIL = il.TRUSEQ_READ2[:24]   # GTGACTGGAGTTCAGACGTGTGCT
LAMP_POOL_FWD = il.P5 + LAMP_POOL_I5 + LAMP_POOL_FWD_TAIL
LAMP_POOL_REV = il.P7 + LAMP_POOL_I7 + LAMP_POOL_REV_TAIL
# The 5' G of the pool reverse tail completes canonical GTGACT on the TGACT first-PCR arm.
if LAMP_POOL_REV_TAIL[1:] != LAMP_PCR_REV_ARM[:len(LAMP_POOL_REV_TAIL) - 1]:
    raise ValueError("pool reverse tail does not overlap the first-PCR reverse arm")

LAMP_BARCODE_FEATURE = feature(
    "lamp_seq_fip_barcode", "sample_index", "whitelist",
    whitelist="LAMP-Seq Supplementary Data 1 Tier-1 or Tier-2 FIP barcodes")

ACTB_BARCODE_FEATURE = feature(
    "lamp_seq_actb_barcode", "sample_index", "unknown", group="ACTB_control",
    note="Source barcoded ACTB-FIP control; native barcode width/member and primer architecture unassigned")

# A LAMP reaction contains several branched and inverted isomers. This construct is only
# the source-established PCR-selected interval: exact landing sites and the embedded FIP
# barcode are retained, while their intervening order within a selected concatemer is not
# guessed from the six-primer reaction.
LAMP_SELECTED_INSERT = [
    Segment("PCR-N forward landing region", LAMP_PCR_FWD_SITE),
    role("LAMP-derived N-gene sequence before the recovered FIP barcode"),
    Segment("FIP sample barcode (representative)", LAMP_FIP_BARCODE, "cbc",
            feature=LAMP_BARCODE_FEATURE),
    role("LAMP-derived N-gene sequence after the recovered FIP barcode"),
    Segment("PCR-N reverse landing region reverse complement",
            revcomp(LAMP_PCR_REV_SITE)),
]
LAMP_PREAMP = Construct([
    Segment("TruSeq Read 1 arm", il.TRUSEQ_READ1, "r1"),
    *LAMP_SELECTED_INSERT,
    Segment("first-PCR reverse arm reverse complement", revcomp(LAMP_PCR_REV_ARM), "r2"),
], name="first-stage LAMP-Seq amplicon")
# Both physical10nt pool indexes come from the explicitly selected published pair1.
LAMP_LIBRARY = Construct([
    Segment("sourceP5",il.P5),
    Segment("physicali5 pool1",LAMP_POOL_I5,feature=feature("lamp_pool_i5","sample_index","fixed",note="sourcePool_fwd_1; physical10nt")),
    Segment("sourceRead1",il.TRUSEQ_READ1), *LAMP_SELECTED_INSERT,
    Segment("sourceRead2 rebuilt by overlapping PCR tails",revcomp("G"+LAMP_PCR_REV_ARM)),
    Segment("physicali7 pool1 reverse complement",revcomp(LAMP_POOL_I7),feature=feature("lamp_pool_i7","sample_index","fixed",note="sourcePool_rev_1; physical10nt")),
    Segment("sourceP7rc",il.P7_RC)],"printed-primer-derived LAMP-Seq read reference")
LAMP_SEQ_PRIMERS=tuple(sp.TRUSEQ[r] for r in ("R1","I1","I2","R2"))



def _lamp_oligos():
    feature_seg = Segment("FIP sample barcode", LAMP_FIP_BARCODE, "cbc",
                          feature=LAMP_BARCODE_FEATURE)
    yield oligo("N-F3-LNA", [Segment("LNA-modified N-F3", LAMP_F3, placeholder=True)])
    yield oligo("N-B3-LNA", [Segment("LNA-modified N-B3", LAMP_B3, placeholder=True)])
    yield oligo("N-FIP-Barcode (representative Tier-1 barcode)", [
        Segment("F1c", LAMP_FIP_LEFT), feature_seg, Segment("F2", LAMP_FIP_RIGHT)])
    yield oligo("N-BIP", [Segment("B1c + B2", LAMP_BIP)])
    yield oligo("N-LF", [Segment("loop forward", LAMP_LF)])
    yield oligo("N-LB", [Segment("loop backward", LAMP_LB)])
    yield oligo("PCR-N-fwd-10", [Segment("Read 1 tail + N target", LAMP_PCR_FWD)])
    yield oligo("PCR-N-rev-11", [Segment("Read 2 tail + N target", LAMP_PCR_REV)])
    yield oligo("Pool_fwd_1", [Segment("P5", il.P5, "p5"),
                               Segment("i5 pool index", LAMP_POOL_I5, "cbc",feature=feature("lamp_pool_i5","sample_index","fixed")),
                               Segment("partial Read 1 tail", LAMP_POOL_FWD_TAIL, "r1")])
    yield oligo("Pool_rev_1", [Segment("P7", il.P7, "p7"),
                               Segment("i7 pool index", LAMP_POOL_I7, "cbc",feature=feature("lamp_pool_i7","sample_index","fixed")),
                               Segment("partial Read 2 tail", LAMP_POOL_REV_TAIL, "r2")])


def native(stage):
    scene=Scene()
    scene.strand("original_RNA",[role("possible original viral RNA and human RNA; heat/processing recovery unresolved")],label="independent original RNA lineage")
    scene.strand("carrier_DNA",[role("pX330 carrier plasmid DNA and possible host DNA; not viral RNA or synthesized cDNA")])
    if stage not in ("sample","lysate"):
        scene.strand("LAMP_DNA",[role("RT-derived DNA plus heterogeneous branched/stem-loop/inverted LAMP DNA with embedded FIP10 sample identifier",molecular_feature=LAMP_BARCODE_FEATURE)],label="native heterogeneous LAMP population; selected PCR interval is not substituted for it")
        scene.strand("control_DNA",[role("ACTB control RT-LAMP populations from independent hostRNA; native control barcode width/member unassigned; no viral diagnostic result inferred",molecular_feature=ACTB_BARCODE_FEATURE)])
        scene.strand("native_primers",[role("possible original DNA/LNA primers and off-target products; recovery unresolved")],label="Original outer-primer LNA chemistry stays on original molecules; copied DNA does not inherit LNA sugars")
    if stage in ("preamp","indexed"):
        con=LAMP_PREAMP if stage=="preamp" else LAMP_LIBRARY
        scene.strand("PCR_one",list(con),label="PCR daughter component reference; concatemer selection/firstcycle assembly unassigned")
        scene.strand("PCR_two",complement_segments(list(con)),label="independent complementary PCR daughter DNA")
        scene.strand("PCR_native_residuals",[role("possible original LAMP templates, original PCR primers and off-products; no cleanup purity assumed")])
        if stage=="indexed":scene.strand("other_pooled_lineages",[role("other pooled libraries have independently assigned poolindexes; pair1 reference does not assign their identities")])
    return scene


def strand(name,seq):
    return Construct([Segment(name,seq)],name)


def _lamp_graph(source):
    g=ProtocolGraph("LAMP-Seq",authority=source)
    def inp(i,l,k=M.CONDITION,state=None):return g.input(i,l,k,state=state)
    def take(i,l,x):
        y=g.product(i+"_used",l,M.MOLECULAR,state=x.state)
        rem=g.product(i+"_remainder","untransferred preparation",M.MOLECULAR,state=x.state)
        c=inp(i+"_allocation",l)
        g.react(i,l,K.SPLIT,((x,I.TEMPLATE),(c,I.CONDITION)),(y,rem))
        g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,"Physical unsampled remainder; fate unassigned, not assumed discarded."))
        return y
    def clean(i,l,x,kind=K.CLEANUP):
        c=inp(i+"_selection",l,M.SELECTION)
        y=g.product(i+"_kept",l+" retained mixture",x.kind,state=x.state)
        w=g.product(i+"_waste",l+" nonretained fraction",M.WASTE)
        g.react(i,l,kind,((x,I.TEMPLATE),(c,I.SELECTION)),(y,w));g.discard(w,"Actual source physical fraction; class-specific recovery/purity unresolved.")
        return y
    sample=g.input("sample","swab QuickExtract with pX330 carrierDNA; commercial branch",M.BIOLOGICAL,state=native("sample"),primary=True)
    heat=inp("heat","source commercial pretreatment: ≥10 min at room temperature, then 95 °C for 5 min")
    lysate=g.product("heated_lysate","heat-inactivated mixed lysate",M.MOLECULAR,state=native("lysate"))
    g.react("inactivate","Heat-inactivate mixed swab preparation",K.THERMAL_TREATMENT,((sample,I.TEMPLATE),(heat,I.CONDITION)),(lysate,))
    x=take("pretreat_aliquot","source100µl lysate aliquot",lysate)
    beads=inp("carbon_ion_beads","source activatedcarbon/ionexchange treatment",M.SELECTION)
    treated=g.product("bead_treated","bead-treated mixed lysate",M.MOLECULAR,state=x.state)
    g.react("bead_treatment","Incubate with carbon/ion-exchange adsorbents",K.BINDING,((x,I.TEMPLATE),(beads,I.SELECTION)),(treated,))
    x=clean("bead_transfer","source90µl transfer into fresh tube; residual bead mixture left behind",treated)
    acid=inp("hcl","source HCl pH adjustment")
    adjusted=g.product("adjusted_lysate","pH-adjusted mixed lysate",M.MOLECULAR,state=x.state)
    g.react("adjust_ph","Adjust lysate pH",K.CHEMICAL_CONVERSION,((x,I.TEMPLATE),(acid,I.CONDITION)),(adjusted,))
    x=take("lamp_aliquot","source swab-lysate allocation into RT-LAMP well",adjusted)
    fip=inp("lamp_fip","representativeTier1 C-FIP DNA barcodeAGCAAACTAC",M.OLIGO,Construct([Segment("F1c",LAMP_FIP_LEFT),Segment("physical10nt samplebarcode",LAMP_FIP_BARCODE,feature=LAMP_BARCODE_FEATURE),Segment("F2",LAMP_FIP_RIGHT)]))
    outer=Scene()
    for n,q in (("C-F3-LNA",LAMP_F3),("C-B3-LNA",LAMP_B3)):
        outer.strand(n,[Segment(n,q,placeholder=True,note="source+ marks LNA at following base; chemical sugar substitutions, not extra nucleotides")])
    lna=inp("lamp_outer","source LNA-modified outer primers",M.OLIGO,outer)
    other=Scene()
    for n,q in (("C-BIP",LAMP_BIP),("C-LF",LAMP_LF),("C-LB",LAMP_LB)):other.strand(n,[Segment(n,q)])
    others=inp("lamp_bip_loops","source separate BIP/LF/LB DNA primers",M.OLIGO,other)
    control_state=Scene()
    control_state.strand("ACTB_FIP",[role("physical barcoded ACTB-FIP primer; ordered strand and barcode width/member unassigned",molecular_feature=ACTB_BARCODE_FEATURE)],label="independent control primer lineage")
    control_state.strand("other_ACTB_primers",[role("other independent ACTB control primers; ordered strands unassigned")])
    control=inp("actb_control_primers","commercial master mix includes independent ACTB primers and barcoded ACTB-FIP; native architecture unexpanded",M.OLIGO,control_state)
    carrierstate=Scene();carrierstate.strand("mastermix_carrierDNA",[role("additional sourcepX330 carrierDNA; not RNA or RT-derived DNA")])
    carrier=inp("mastermix_carrier","additionalpX330 DNA in source commercialmastermix",M.MOLECULAR,carrierstate)
    bst=inp("bst","source commercial LAMP mix+Bst3; 65 °C for 80 min",M.CATALYST)
    lamp=g.product("lamp_product","mixed native RT-LAMP DNA/original polymers",M.MOLECULAR,state=native("lamp"))
    g.react("rt_lamp","Source commercial RT-LAMP of mixed lysate",K.EXTENSION,((x,I.TEMPLATE),(fip,I.STRUCTURAL),(lna,I.STRUCTURAL),(others,I.STRUCTURAL),(control,I.STRUCTURAL),(carrier,I.STRUCTURAL),(bst,I.CATALYST)),(lamp,))
    peers=inp("other_lamp_reactions","other independently barcoded RT-LAMP wells in source pool; exact membership unassigned",M.MOLECULAR,native("lamp"))
    pooled=g.product("pooled_lamp","physical multiwell RT-LAMP pool",M.MOLECULAR,state=lamp.state)
    g.react("pool_lamp","Physically pool independently barcoded RT-LAMP reactions",K.POOL,((lamp,I.TEMPLATE),(peers,I.TEMPLATE)),(pooled,))
    dilution=inp("dilution","source commercial pre-PCR dilution 1:100,000; no physical class purification")
    diluted=g.product("diluted_lamp","diluted native RT-LAMP pool",M.MOLECULAR,state=pooled.state)
    g.react("dilute","Dilute pooled mixture",K.GENERIC,((pooled,I.TEMPLATE),(dilution,I.CONDITION)),(diluted,),changes_molecule=False)
    x=take("pcr1_aliquot","source5µl diluted mixture into targetPCR",diluted)
    f=inp("pcr_f","sourcePCR-N-fwd-10 / PCR-C-fwd",M.OLIGO,strand("PCRforward",LAMP_PCR_FWD))
    r=inp("pcr_r","sourcePCR-N-rev-11 / PCR-C-rev",M.OLIGO,strand("PCRreverse",LAMP_PCR_REV))
    pol=inp("pcr_polymerase","NEBNext 2× mix; each source PCR uses 18 cycles of 98 °C 20 s, 65 °C 20 s, 72 °C 30 s",M.CATALYST)
    pre=g.product("preamp","targetPCR daughters plus original templates/primers",M.MOLECULAR,state=native("preamp"))
    g.react("target_pcr","TargetPCR from actual diluted heterogeneous LAMP population",K.PCR,((x,I.TEMPLATE),(f,I.STRUCTURAL),(r,I.STRUCTURAL),(pol,I.CATALYST)),(pre,),geometry=PCRGeometry.unavailable(x.state,reason="Source PCR target primers known; actual native branched/inverted LAMP mixture enters PCR. Selected linear interval is only a read-reference, not replacement template or proved firstcycle product."))
    x=take("pcr2_aliquot","source5µl priorPCR into poolindexPCR",pre)
    pf=inp("pool_f","sourcePool_fwd_1 physical10nt index",M.OLIGO,Construct([Segment("P5",il.P5),Segment("i5pool1",LAMP_POOL_I5,feature=feature("lamp_pool_i5","sample_index","fixed")),Segment("partialRead1",LAMP_POOL_FWD_TAIL)]))
    pr=inp("pool_r","sourcePool_rev_1 physical10nt index",M.OLIGO,Construct([Segment("P7",il.P7),Segment("i7pool1",LAMP_POOL_I7,feature=feature("lamp_pool_i7","sample_index","fixed")),Segment("partialRead2",LAMP_POOL_REV_TAIL)]))
    indexed=g.product("indexed_pcr","indexed PCR daughters plus native residuals",M.LIBRARY,state=native("indexed"))
    g.react("index_pcr","Pool-indexPCR from actual firstPCR aliquot",K.PCR,((x,I.TEMPLATE),(pf,I.STRUCTURAL),(pr,I.STRUCTURAL),(pol,I.CATALYST)),(indexed,),geometry=PCRGeometry.unavailable(x.state,reason="Actual unpurified firstPCR aliquot enters secondPCR; source primer overlaps establish read-reference arms, not class-pure input or complete native product distribution."))
    otherstate=Scene();otherstate.strand("independent_libraries",[role("other independently indexed PCR mixtures; index identities/native polymers unassigned; not all pair1")])
    otherlib=inp("other_indexed_libraries","other separately pool-indexed PCR preparations; membership and ratios unassigned",M.MOLECULAR,otherstate)
    total=g.product("pooled_indexed","physically pooled indexed PCR mixtures",M.LIBRARY,state=indexed.state)
    g.react("pool_indexed","Pool indexed PCR products before gel",K.POOL,((indexed,I.TEMPLATE),(otherlib,I.TEMPLATE)),(total,))
    x=clean("gel_band","source2% Egel library-band excision",total,K.SELECTION)
    for n in (1,2):x=clean("silica"+str(n),"sourceQiagen silica column round"+str(n),x)
    g.sequence(x,SequencingPlan(library=LAMP_LIBRARY,primers=tuple(LAMP_SEQ_PRIMERS),acquisition="unspecified",limitation="Commercial twoPCR route names MiSeq/iSeq but does not uniquely assign its run recipe. The separate named sequencing sections print100Read1+8i7+8i5, reading8 of each physical10nt index; noRead2 in that recipe. Exact acquired recipe for this selected preparation and pool membership remain unassigned. Correct viral sequence excluding primer-covered bases plus barcode cooccurrence is computational diagnostic analysis; this molecular graph establishes no detection threshold or clinical performance."))
    g.declare_reviewed_scope(REVIEWED,source)
    return g.require_publishable()


source=('<a href="https://doi.org/10.1038/s41587-021-00966-9">Ludwig et al.,Nature Biotechnology2021</a>; published PMC8678193 commercialtwoPCR Methods and SupplementaryData1.')
TITLE="LAMP-Seq — barcoded RT-LAMP with Illumina readout"
NOTES="01_lamp-seq.html"
SOURCE=source
SUMMARY="Selected published commercial twoPCR branch with representative Tier1 FIP barcode and pool-primer pair1. Native RNA, carrierDNA, modified original primers, LAMP concatemer populations and PCR daughters remain distinct."
CAVEAT="Independent ACTB control populations are retained but not expanded. Branched LAMP intermediates and PCR firstcycle assembly remain unresolved. No single linear selected interval substitutes for native incoming PCR material. Open-source singlePCR and N100/PhiX run-library admixtures are outside this graph."
oligos=_lamp_oligos
FINAL_LIBRARY=LAMP_LIBRARY
FINAL_CAPTION="Printed-primer-derived component/read reference with representative sample barcode and actual10nt pool1 indexes. Target-selected interval order inside native LAMP concatemers remains unresolved."
SEQ_PRIMERS=LAMP_SEQ_PRIMERS
SEQUENCING_INTRO="The source named MiSeq/iSeq recipes acquire100Read1+8i7+8i5, so each index acquisition covers only8 of10physical bases. The selected commercial branch's actual run linkage is not assigned."

def protocol_graph():return _lamp_graph(source)

REVIEWED=("sample","heat","inactivate","pretreat_aliquot_allocation","pretreat_aliquot","carbon_ion_beads","bead_treatment","bead_transfer_selection","bead_transfer","hcl","adjust_ph","lamp_aliquot_allocation","lamp_aliquot","lamp_fip","lamp_outer","lamp_bip_loops","actb_control_primers","mastermix_carrier","bst","rt_lamp","other_lamp_reactions","pool_lamp","dilution","dilute","pcr1_aliquot_allocation","pcr1_aliquot","pcr_f","pcr_r","pcr_polymerase","target_pcr","pcr2_aliquot_allocation","pcr2_aliquot","pool_f","pool_r","index_pcr","other_indexed_libraries","pool_indexed","gel_band_selection","gel_band","silica1_selection","silica1","silica2_selection","silica2","silica2_kept")
