"""10x Chromium 3' v3 Feature Barcode chemistry for direct sgRNA capture."""
from __future__ import annotations
import illumina as il
import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct,Scene,Segment,feature,revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind, Terminal)
from protocol_model import molecule, role
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
from protocol_families.tenx_graph import partition_material

CS1 = "TTGCTAGGACCGGCCTTAAAGC"
CS2 = "CCTTAGCCGCTAATAGGTGAGC"
CS1_RC = revcomp(CS1)
CS2_RC = revcomp(CS2)
CR1 = "GTTTAAGAGCTAAGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
PARTIAL_R1N = nx.READ1_PRIMER[-22:]
TSO_DNA = rt.SMART_HANDLE + "ACAT"
FEATURE_CDNA_FORWARD = nx.READ1_PRIMER[-27:]
FEATURE_CDNA_REVERSE = rt.SMART_HANDLE[:22]
FEATURE_SI_FORWARD = il.P5 + nx.READ1_PRIMER
FEATURE_SI_REVERSE = il.TRUSEQ_READ2 + FEATURE_CDNA_REVERSE
CELL_BARCODE=feature("cell_barcode","cell_barcode","whitelist",whitelist="10x-chromium-3prime-v3")
UMI_FEATURE=feature("umi","umi","random")
GUIDE_FEATURE=feature("guide","guide_barcode","unknown")
I7_FEATURE=feature("sample_index_i7","sample_index","unknown")
FEATURES={"cell barcode":CELL_BARCODE,"UMI":UMI_FEATURE,"protospacer":GUIDE_FEATURE,
          "protospacer complement":GUIDE_FEATURE,"i7 reverse complement":I7_FEATURE}
def seg(name,top,tag=None,**kw):
    kw.setdefault("feature",FEATURES.get(name)); return Segment(name=name,top=top,tag=tag,**kw)
def bead_primer(): return [seg("Partial Read 1N",PARTIAL_R1N,"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*12,"umi",placeholder=True),seg("Capture Sequence 1",CS1,"tso")]
def tso(): return [seg("SMART handle",rt.SMART_HANDLE,"tso"),seg("ACAT","ACAT","tso"),seg("rGrGrG","GGG","tso")]
def sgrna():
    return [seg("protospacer", "N" * 20, "cbc", placeholder=True),
            seg("CR1 scaffold", CR1, "insert"),
            seg("Capture Sequence 1 rc", CS1_RC, "tso"),
            seg("Pol III terminator", "U" * 7, "insert")]
def copied_sgrna(): return [seg("CR1 complement",revcomp(CR1),"insert"),seg("protospacer complement","n"*20,"cbc",placeholder=True)]
def feature_cdna_primers(): return ([seg("Read 1N completion",FEATURE_CDNA_FORWARD,"r1")],[seg("partial TSO",FEATURE_CDNA_REVERSE,"tso")])
def feature_si_primers(): return ([seg("P5 + Read 1N",FEATURE_SI_FORWARD,"p5+r1")],[seg("Read 2 + partial TSO",FEATURE_SI_REVERSE,"r2+tso")])
def capture_scene():
    sc=Scene(); sc.strand("sgRNA",sgrna(),label="capture-sequence sgRNA")
    sc.anneal("bead primer",bead_primer(),to="sgRNA",pair=("Capture Sequence 1","Capture Sequence 1 rc"),label="gel-bead primer",mod5="bead",unpaired=("UMI",))
    sc.mark("sgRNA","protospacer","Feature Barcode / guide identity")
    sc.mark("sgRNA","Capture Sequence 1 rc","engineered capture sequence")
    sc.arrow("bead primer","reverse transcription"); return sc
def first_strand(): return Construct([*bead_primer(),*copied_sgrna(),seg("CCC",rt.UNTEMPLATED_TAIL,"tso")],name="captured sgRNA first strand")
def switch_scene():
    first=first_strand(); sc=Scene(); sc.strand("first",list(first),label="first strand")
    sc.anneal("TSO",tso(),to="first",pair=("rGrGrG","CCC"),label="TSO",above=True); sc.arrow("first","RT switches template"); return sc
def tso_complement(): return rt.UNTEMPLATED_TAIL+revcomp(TSO_DNA)
def feature_cdna(): return Construct([*bead_primer(),*copied_sgrna(),seg("TSO complement",tso_complement(),"tso")],name="CRISPR feature cDNA")
SEQ_PRIMERS=(sp.NEXTERA["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["R2"])
def final_library():
    lib=Construct([seg("P5",il.P5,"p5"),seg("Nextera Read 1",nx.READ1_PRIMER,"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*12,"umi",placeholder=True),seg("Capture Sequence 1",CS1,"tso"),*copied_sgrna(),seg("TSO complement",tso_complement(),"tso"),seg("Read 2 leading A","A","r2"),seg("Index 1 / Read 2 arm",il.INDEX1_PRIMER,"r2"),seg("i7 reverse complement","I"*8,"cbc",placeholder=True),seg("P7 reverse complement",il.P7_RC,"p7")],name="10x 3' v3 CRISPR Screening library")
    problems=sp.verify(lib,SEQ_PRIMERS,required_roles=("Read 1","Index 1 (i7)","Read 2"))
    if problems: raise ValueError("invalid CRISPR Screening library: "+"; ".join(problems))
    return lib
def read_layout():
    spans=sp.feature_spans(final_library(),SEQ_PRIMERS,READ_LENGTHS)
    rows=lambda role:[(f"{x.read}, {x.cycles}",x.feature.label) for x in spans if x.read==role]
    return [*rows("Read 1"),("Read 2, 1–30","TSO-derived sequence"),
            *rows("Read 2"),("Read 2, 51–91","CR1 scaffold (first 41 bases)")]


def amplified_cdna():
    """CG000184 Rev A appendix Step 2.2: Partial Read 1N primer (27 nt) extends the bead handle."""
    return Construct([seg("Partial Read 1N",FEATURE_CDNA_FORWARD,"r1"),*list(feature_cdna())[1:]],name="amplified CRISPR feature cDNA")
def feature_pcr_product():
    """CG000184 Rev A appendix Step 4.2 Feature PCR product (P5/Read 1N ... TSO/Read 2)."""
    return Construct([seg("P5",il.P5,"p5"),seg("Nextera Read 1",nx.READ1_PRIMER,"r1"),*list(feature_cdna())[1:],
                      seg("Read 2 reverse complement",revcomp(il.TRUSEQ_READ2),"r2")],name="CRISPR Feature PCR product")
I7_PRIMER=il.P7+"N"*8+il.TRUSEQ_READ2[:21]   # Chromium i7 Sample Index PN-220103 (4-oligo well mix)
READ_LENGTHS={"Read 1":28,"Index 1 (i7)":8,"Read 2":91}   # CG000184 Rev A step 5 run parameters
READ2_MINIMUM=70                                          # same table footnote


def protocol_graph():
    authority=('Commercial protocol: <a href="https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/CG000184_ChromiumSingleCellSingleCell3v3_FeatureBarcodingtechnology_CRISPR_RevA.pdf">'
               '10x Genomics CG000184 Rev A</a> (vendor PDF mirrored by scg_lib_structs); guide design authority CG000197 Rev A.')
    g=ProtocolGraph("10x 3-prime Feature Barcode direct guide capture",authority=authority)
    cells=g.input("cells","cells expressing Capture Sequence 1 guide RNA",
                  MaterialKind.BIOLOGICAL,
                  state=molecule("programmed cells",(role("CS1 guide RNA"),role("poly(A) mRNA"))),
                  primary=True)
    reviewed=["cells"]
    occupied=partition_material(g,cells,reviewed,label="guide-expressing cells")
    bead=g.input("bead_primer","Capture Sequence 1 cell-barcode/UMI bead primer",
                 MaterialKind.OLIGO,state=Construct(bead_primer(),name="Capture Sequence 1 gel-bead primer"))
    rt_enzyme=g.input("rt","RT Enzyme C (GEM-RT 53 °C 45 min, 85 °C 5 min)",MaterialKind.CATALYST)
    first=g.product("first","cell-barcoded first-strand guide cDNA",MaterialKind.MOLECULAR,
                    state=first_strand())
    g.react("reverse_transcribe","Reverse-transcribe the engineered guide from CS1",
            StepKind.REVERSE_TRANSCRIPTION,
            ((occupied,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
             (rt_enzyme,InputRole.CATALYST)),(first,))
    tso_oligo=g.input("tso","Template Switch Oligo",MaterialKind.OLIGO,state=Construct(tso(),name="Template Switch Oligo"))
    switched=g.product("switched","template-switched guide feature cDNA",
                       MaterialKind.MOLECULAR,state=feature_cdna())
    g.react("template_switch","Template-switch at the copied guide end",
            StepKind.TEMPLATE_SWITCHING,
            ((first,InputRole.TEMPLATE),(tso_oligo,InputRole.STRUCTURAL),
             (rt_enzyme,InputRole.CATALYST)),(switched,))
    recovery=g.input("gem_recovery","Recovery Agent GEM break and Dynabeads MyOne SILANE cleanup",
                     MaterialKind.SELECTION)
    recovered=g.product("recovered","recovered GEM-RT cDNA",MaterialKind.MOLECULAR,state=feature_cdna())
    recovery_waste=g.product("recovery_waste","oil phase and bead-cleanup supernatant",MaterialKind.WASTE)
    g.react("recover","Break GEMs and purify GEM-RT products (step 2.1)",StepKind.SELECTION,
            ((switched,InputRole.TEMPLATE),(recovery,InputRole.SELECTION)),(recovered,recovery_waste),
            changes_molecule=False)
    g.discard(recovery_waste,"removed during post GEM-RT cleanup")
    reviewed.extend(("bead_primer","rt","reverse_transcribe","tso","template_switch","gem_recovery","recover"))

    def oligo_input(identifier,label,parts):
        return g.input(identifier,label,MaterialKind.OLIGO,state=Construct(parts,name=label))
    def pcr(prefix,label,template,fwd,rev,product,geometry):
        pol=g.input(prefix+"_polymerase","Amp Mix (PN-2000047)",MaterialKind.CATALYST)
        g.react(prefix,label,StepKind.PCR,((template,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),
                (rev,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(product,),geometry=geometry)
        reviewed.extend((fwd.id,rev.id,pol.id,prefix))
        return product
    def cleanup(prefix,label,template,state,kept_label):
        beads=g.input(prefix+"_beads",label,MaterialKind.SELECTION)
        kept=g.product(prefix+"_kept",kept_label,MaterialKind.MOLECULAR,state=state)
        waste=g.product(prefix+"_waste",label+" unbound/excluded fraction and washes",MaterialKind.WASTE)
        g.react(prefix,label,StepKind.SELECTION,((template,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),
                (kept,waste),changes_molecule=False)
        g.discard(waste,"removed by "+label)
        reviewed.extend((beads.id,prefix))
        return kept

    f1=oligo_input("feature_cdna_pcr_forward","Feature cDNA Primers 1 Partial Read 1N (PN-2000096)",feature_cdna_primers()[0])
    r1=oligo_input("feature_cdna_pcr_reverse","Feature cDNA Primers 1 Partial TSO (PN-2000096)",feature_cdna_primers()[1])
    enriched=pcr("feature_cdna_pcr","Amplify cDNA with Feature cDNA Primers 1 (11–13 cycles)",recovered,f1,r1,
        g.product("feature_cdna_pcr_product","amplified sgRNA cDNA and mRNA cDNA",MaterialKind.MOLECULAR,state=amplified_cdna()),
        PCRGeometry.exact(feature_cdna(),PCRPrimer("Feature cDNA forward",FEATURE_CDNA_FORWARD),
            PCRPrimer("Feature cDNA reverse",FEATURE_CDNA_REVERSE),
            statement="CG000184 Rev A appendix Step 2.2 prints the sgRNA primer pair; the mRNA pair in the same mix is not drawn."))
    spri=g.input("size_selection","0.6X SPRIselect",MaterialKind.SELECTION)
    supernatant=g.product("selected","transferred 150 µl supernatant with sgRNA cDNA",MaterialKind.MOLECULAR,state=amplified_cdna())
    gex=g.product("gex_pellet","bead-bound gene-expression cDNA (pellet)",MaterialKind.MOLECULAR)
    g.react("size_select","Separate short sgRNA cDNA (supernatant) from gene-expression cDNA (pellet)",
            StepKind.SELECTION,((enriched,InputRole.TEMPLATE),(spri,InputRole.SELECTION)),(supernatant,gex),
            changes_molecule=False)
    g.mark_terminal(Terminal(gex.id,Disposition.STORED,
        "pellet cleanup (2.3A) feeds the separate 3-prime Gene Expression library, step 3; not drawn"))
    reviewed.extend(("size_selection","size_select"))
    sup=cleanup("supernatant_cleanup","1.2X SPRIselect (2.3B)",supernatant,amplified_cdna(),"purified sgRNA cDNA (50 µl)")
    guide=cleanup("guide_cdna_cleanup","1.0X SPRIselect Guide RNA cDNA Cleanup (4.1)",sup,amplified_cdna(),"guide RNA cDNA (5 µl of 30 µl used)")
    f2=oligo_input("feature_pcr_forward","Feature SI Primers 1 P5/Read 1N primer (PN-2000098)",feature_si_primers()[0])
    r2=oligo_input("feature_pcr_reverse","Feature SI Primers 1 Read 2/partial TSO primer (PN-2000098)",feature_si_primers()[1])
    featured=pcr("feature_pcr","Feature PCR: add P5 and the Read 2 arm (15 cycles)",guide,f2,r2,
        g.product("feature_pcr_product","Feature PCR product",MaterialKind.MOLECULAR,state=feature_pcr_product()),
        PCRGeometry.exact(amplified_cdna(),PCRPrimer("Feature SI forward",FEATURE_SI_FORWARD),
            PCRPrimer("Feature SI reverse",FEATURE_SI_REVERSE),
            statement="CG000184 Rev A appendix Step 4.2 prints both primers (labelled PN-2000099 there; protocol step 4.2 uses Feature SI Primers 1 PN-2000098)."))
    post=cleanup("feature_pcr_cleanup","1.0X SPRIselect Post Feature PCR Cleanup (4.3)",featured,feature_pcr_product(),"Feature PCR product (5 µl of 30 µl used)")
    lib=final_library()
    f3=oligo_input("index_pcr_forward","Feature SI Primers 2 P5/Read 1N primer (PN-2000099)",feature_si_primers()[0])
    r3=g.input("index_pcr_reverse","Chromium i7 Sample Index (PN-220103)",MaterialKind.OLIGO,
               state=Construct([seg("P7",il.P7,"p7"),seg("i7","I"*8,"cbc",placeholder=True,feature=I7_FEATURE),
                                seg("Partial Read 2",il.TRUSEQ_READ2[:21],"r2")],name="Chromium i7 Sample Index"))
    indexed=pcr("index_pcr","Sample Index PCR: add P7 and the i7 sample index (5 cycles)",post,f3,r3,
        g.product("index_pcr_product","indexed CRISPR Screening library",MaterialKind.MOLECULAR,state=lib),
        PCRGeometry.exact(feature_pcr_product(),PCRPrimer("Feature SI Primers 2",FEATURE_SI_FORWARD),
            PCRPrimer("Chromium i7 Sample Index",I7_PRIMER),
            statement="CG000184 Rev A appendix Step 4.4 prints both primers; each i7 well is a mix of four 8-nt index oligos."))
    beads=g.input("double_sided_selection","0.7X/1.0X double-sided SPRIselect (4.5)",MaterialKind.SELECTION)
    library=g.product("library","10x 3' v3 CRISPR Screening library",MaterialKind.LIBRARY,state=lib)
    large=g.product("large_fraction","0.7X bead-bound large fragments",MaterialKind.WASTE)
    small=g.product("small_fraction","1.0X supernatant (primers, short products)",MaterialKind.WASTE)
    g.react("double_sided","Post Sample Index PCR double-sided size selection",StepKind.SELECTION,
            ((indexed,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(library,large,small),changes_molecule=False)
    g.discard(large,"0.7X pellet discarded; supernatant carried forward")
    g.discard(small,"1.0X supernatant removed")
    reviewed.extend(("double_sided_selection","double_sided"))
    g.annotate(SourceDiscrepancy(target="feature_pcr",title="Feature PCR primer catalogue number",
        claims=(SourceClaim("CG000184 Rev A step 4.2 Feature PCR Mix table (p.51)","Feature SI Primers 1, PN-2000098"),
                SourceClaim("CG000184 Rev A appendix 'Protocol Step 4.2 – Feature PCR' (p.70)","labelled Feature SI Primers 2, PN-2000099, but prints the P5/Read 1N and Read 2/partial-TSO pair")),
        model_choice="The model uses the printed primer pair and names it Feature SI Primers 1 (PN-2000098).",
        evidence="The step 4.2 reaction table and the kit contents list Feature SI Primers 1 for Feature PCR; step 4.4 uses Feature SI Primers 2 with the Chromium i7 index, and only the Read 2/partial-TSO primer converts the TSO end.",
        uncertainty="The appendix label may be a typographic error; the physical composition of PN-2000098 versus PN-2000099 beyond the printed sequences is not disclosed."))
    g.annotate(SourceDiscrepancy(target="index_pcr",title="TSO segment in the final library drawing",
        claims=(SourceClaim("CG000184 Rev A appendix Step 4.2 Feature PCR product (p.70)","…sgRNA_Feature_Barcode-CCCATGTACTCTGCGTTGATACCACTGCTT-AGATCGGAAGAGC…"),
                SourceClaim("CG000184 Rev A appendix Step 4.4 Sample Index PCR product (p.70)","…sgRNA_Feature_Barcode-AGATCGGAAGAGC… (no TSO bases, though the label row lists TS0)")),
        model_choice="The final library keeps the 30 TSO-derived bases between the protospacer complement and TruSeq Read 2.",
        evidence="Neither Sample Index PCR primer anneals within or removes the TSO-derived segment; the i7 primer lands on Read 2 and the Feature SI Primers 2 primer on Read 1N.",
        uncertainty="The step 4.4 sequence is read as a drawing abbreviation; no source states that the TSO segment is removed."))
    g.sequence(library,SequencingPlan(library=lib,primers=tuple(SEQ_PRIMERS),
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS),read_lengths=dict(READ_LENGTHS)))
    reviewed.append(library.id)
    g.declare_reviewed_scope(reviewed,authority)
    return g.require_publishable()
