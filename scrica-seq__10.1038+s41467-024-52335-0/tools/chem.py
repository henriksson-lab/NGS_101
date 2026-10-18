"""scRICA-seq circular cDNA isoform reconstruction plus ATAC."""
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Row, Scene, Segment, duplex_rows, feature, oligo, revcomp
from circular import CircularizationGeometry, ClosedCircle, rolling_circle
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (AnnotationSlot, Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind)
TITLE="scRICA-seq — RNA isoforms and chromatin accessibility"
NOTES="01_scrica-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-024-52335-0">Jiang et al., <i>Nature Communications</i> (2024)</a>.'
SUMMARY="Full-length cDNA receives end UMIs, is circularized and amplified so short Tn5 fragments can be grouped back into isoforms; a matched ATAC branch reports chromatin accessibility from the same cell."
CAVEAT="Supplementary Data 2 prints the low-throughput N6 cell barcode, end-label, linker, RCA, TSS/TES enrichment and terminal index-primer oligos, and Supplementary Fig. 3 draws their products; this page models that branch. INFERRED: the TTE Mix Tn5 adapters are not printed (only their S5/ME end is drawn), and the low-throughput ATAC library follows a cited method (ref. 11) whose native assembly remains unassigned here; a generic Nextera shell is not substituted. Index-read cycles are not stated (PE150). The separate high-throughput 10x Multiome variant uses a 16-base cell barcode and different downstream primers and is not substituted here."
ATAC=Construct([Segment("low-throughput ATAC endpoint; native assembly unassigned","[accessibility library]",placeholder=True)],name="scRICA low-throughput ATAC boundary")
AP=()

DT_HANDLE="AAGCAGTGGTATCAACGCAGAGT"
TSO="GCAATGAAGTCGCAGGGTTG"
TAG1="AGTTAACGCTATACGGCACGCG"
LINKER=revcomp(DT_HANDLE)+TSO
TSO_ISPCR=TSO
OLIGO_PRIMER=DT_HANDLE
PHI29_PRIMER=TSO
TN5_I5=nx.S5+nx.ME[:6]
TN5_I7=nx.S7+nx.ME[:5]
ANTI_TAG1=revcomp(TAG1)
ANTI_TSO=revcomp(TSO)
I7_ANTI_TAG1=nx.ADAPTOR_S7+ANTI_TAG1
I7_ANTI_TSO=nx.ADAPTOR_S7+ANTI_TSO
INDEX_N7=il.P7+"N"*8+nx.ADAPTOR_S7
INDEX_N5=il.P5+"N"*8+nx.ADAPTOR_S5
POLYT="T"*30+"VN"
PE150={"Read 1":150,"Read 2":150}
CELL=feature("cell_barcode","cell_barcode","random")
UMI_F=feature("umi","umi","random")

def _cell(): return Segment("cell barcode","B"*6,"cbc",placeholder=True,feature=CELL)
def _umi(): return Segment("UMI","U"*12,"umi",placeholder=True,feature=UMI_F)
def _tso_site(): return Segment("TSO complement",revcomp(TSO))

# Supplementary Fig. 3 step 1: the amplified cDNA strand begins with the oligo(dT)
# primer (handle, N6, N12, tag1, T30VN) and ends with the rG-templated CCC and the
# TSO complement.  The opposite strand begins with the TSO.
LABELLED_CDNA=Construct([
    Segment("oligo(dT) handle",DT_HANDLE),_cell(),_umi(),Segment("tag1",TAG1),
    Segment("oligo(dT)VN",POLYT,placeholder=True),
    Segment("full-length cDNA","X"*30,placeholder=True),
    Segment("template-switch CCC","CCC"),_tso_site(),
],name="SMART-amplified full-length cDNA (oligo(dT)-primer strand)")
# NEBuilder joins the TSO-complement end to the oligo(dT)-handle end: the printed
# linker must span exactly that junction (Supplementary Fig. 3 step 2).
if LINKER!=revcomp(LABELLED_CDNA.get("TSO complement").top+LABELLED_CDNA.get("oligo(dT) handle").top):
    raise AssertionError("printed linker does not bridge the TSO-complement/oligo(dT)-handle junction")
CIRCLE=ClosedCircle(LABELLED_CDNA,"TSO complement","oligo(dT) handle")
RCA=rolling_circle(CIRCLE,PHI29_PRIMER,copies=3)

# Supplementary Fig. 3 steps 4-5: anti-tag1 copies from tag1 through the UMI, cell
# barcode, handle and TSO into the 5' (TSS) end; anti-TSO copies from the TSO through
# the handle, barcode, UMI, tag1 and T30 into the 3' (TES) end.  Libraries are drawn
# with the Tn5 S5/ME end on the left; the right arm is primer-introduced.
def _tss_insert():
    return [Segment("TSS-proximal cDNA (antisense)","X"*30,placeholder=True),
            Segment("template-switch CCC","CCC"),_tso_site(),
            Segment("oligo(dT) handle",DT_HANDLE),_cell(),_umi(),Segment("tag1",TAG1)]
def _tes_insert():
    return [Segment("TES-proximal cDNA (sense)","X"*30,placeholder=True),
            Segment("poly(A) / T30VN complement","NB"+"A"*30,placeholder=True),
            Segment("tag1 complement",revcomp(TAG1)),_umi(),_cell(),
            Segment("oligo(dT) handle complement",revcomp(DT_HANDLE)),Segment("TSO",TSO)]
def _s5_arm(): return [Segment("S5",nx.S5,"s5"),Segment("left mosaic end",nx.ME,"me")]
def _s7_arm(): return [Segment("right mosaic end reverse complement",nx.ME_RC,"me"),Segment("S7 reverse complement",nx.S7_RC,"s7")]
TSS_PCR1=Construct([*_s5_arm(),*_tss_insert()],name="TSS PCR1 product")
TES_PCR1=Construct([*_s5_arm(),*_tes_insert()],name="TES PCR1 product")
TSS_PCR2=Construct([*_s5_arm(),*_tss_insert(),*_s7_arm()],name="TSS PCR2 product")
TES_PCR2=Construct([*_s5_arm(),*_tes_insert(),*_s7_arm()],name="TES PCR2 product")
I5=feature("scrica_i5","sample_index","unknown",note="printed N5 pool N8; member unassigned")
I7=feature("scrica_i7","sample_index","unknown",note="printed N7 pool N8; member unassigned")

def _indexed(insert,name):
    lib=Construct([Segment("P5",il.P5,"p5"),Segment("N5 index","N"*8,"cbc",placeholder=True,feature=I5),
                   *_s5_arm(),*insert,*_s7_arm(),Segment("N7 index complement","N"*8,"cbc",placeholder=True,feature=I7),Segment("P7 complement",il.P7_RC,"p7")],name=name)
    primers=tuple(sp.NEXTERA[k] for k in ("R1","I1","I2","R2"))
    errors=sp.verify(lib,primers)
    if errors: raise ValueError("; ".join(errors))
    return lib,primers
TSS_LONG,TSS_PRIMERS=_indexed(_tss_insert(),"scRICA TSS long source-primer endpoint")
TSS_SHORT,_=_indexed(_tss_insert(),"scRICA TSS short source-primer endpoint")
TES_LONG,TES_PRIMERS=_indexed(_tes_insert(),"scRICA TES long source-primer endpoint")
TES_SHORT,_=_indexed(_tes_insert(),"scRICA TES short source-primer endpoint")

FINAL_LIBRARIES=(
    ("TSS long-fragment library",TSS_LONG,TSS_PRIMERS,"Tn5 cDNA end on the Read 1 side; on the Read 2 side the primer-introduced arm is followed by tag1, UMI, cell barcode, oligo(dT) handle, TSO and the transcript 5' end.","Nextera primers on the TSS long-fragment branch."),
    ("TSS short-fragment library",TSS_SHORT,TSS_PRIMERS,"Short TSS fragments are a separately retained 0.55x-0.9x size fraction of the same indexed TSS pool.","Nextera primers on the TSS short-fragment branch."),
    ("TES long-fragment library",TES_LONG,TES_PRIMERS,"Tn5 cDNA end on the Read 1 side; Read 2 starts in the TSO and continues through the oligo(dT) handle, cell barcode, UMI, tag1 and T30VN into the transcript 3' end.","Nextera primers on the TES long-fragment branch."),
    ("TES short-fragment library",TES_SHORT,TES_PRIMERS,"Short TES fragments are a separately retained 0.55x-0.9x size fraction of the same indexed TES pool.","Nextera primers on the TES short-fragment branch."),
    ("Accessibility library",ATAC,AP,"The low-throughput ATAC fraction is paired by source cell. Its native adapter/primer assembly is an explicit boundary for the cited predecessor route.","Accessibility sequencing-primer geometry remains unassigned."))

def oligos():
    yield oligo("Oligo(dT) primer",[Segment("SMART handle",DT_HANDLE),Segment("cell barcode","B"*6,"cbc",placeholder=True,feature=CELL),Segment("UMI","U"*12,"umi",placeholder=True,feature=UMI_F),Segment("tag1",TAG1),Segment("oligo(dT)VN",POLYT,placeholder=True)])
    yield oligo("TSO primer",[Segment("TSO",TSO)],mods="/rG//rG//rG/ at 3-prime end")
    for name,seq in (("TSO ISPCR primer",TSO_ISPCR),("linker primer",LINKER),("oligo primer",OLIGO_PRIMER),("Phi29 primer",PHI29_PRIMER),("Tn5-i5 primer",TN5_I5),("Tn5-i7 primer",TN5_I7),("anti-tag1 primer",ANTI_TAG1),("anti-TSO primer",ANTI_TSO),("i7-anti-tag1 primer",I7_ANTI_TAG1),("i7-anti-TSO primer",I7_ANTI_TSO),("Index 1 primer (N7)",INDEX_N7),("Index 2 primer (N5)",INDEX_N5)):
        yield oligo(name,[Segment(name,seq,placeholder="N" in seq)])
def _single(parts,label):
 sc=Scene();sc.strand("strand",parts,label=label);return sc


def _fragments():
 sc=Scene()
 sc.strand("tss",TSS_PCR1.segments,label="INFERRED — TSS-eligible native fragment arm; opposite end/gap unassigned")
 sc.strand("tes",TES_PCR1.segments,label="INFERRED — TES-eligible native fragment arm; opposite end/gap unassigned")
 sc.strand("other",[Segment("other RCA fragments","[other fragment classes]",placeholder=True)],label="independent off-target fragments; amplification/recovery unassigned")
 return sc

def _amplified(target):
 sc=Scene()
 sc.strand("target",target.segments,label="source-primer-derived target reference; native product extent bounded")
 sc.strand("other",[Segment("co-carried RCA and off-target products","[other fragment classes]",placeholder=True)],label="other molecular populations; PCR and bead-class recovery unassigned")
 return sc

SHOW_READ_CONTENT=False
READ_SPECS={name: (
 sp.read_spec("Read 1",150,"TSS/TES cDNA-derived insert; conditional source-primer endpoint","include"),
 sp.read_spec("Read 2",150,"Tag/UMI/N6 barcode and cDNA-derived sequence","extract_then_include"),
 sp.read_spec("Index 1 (i7)","unassigned","Printed N8 physical sample index; acquisition cycles unassigned","undetermined"),
 sp.read_spec("Index 2 (i5)","unassigned","Printed N8 physical sample index; acquisition cycles unassigned","undetermined")) if name!="Accessibility library" else tuple(
 sp.read_spec(role,"unassigned","Low-throughput accessibility run and primer geometry unassigned in this scoped reference","undetermined") for role in ("Read 1","Read 2","Index 1 (i7)","Index 2 (i5)")) for name,*_ in FINAL_LIBRARIES}
READ_SPECS_NOTE="Published low-throughput isoform library sequencing is NovaSeq6000 PE150. Physical N8 index windows do not assign index acquisition/cycles; the high-throughput Multiome/MGI route is separate."


def _cleanup(g, identifier, template, label, state):
    reagent=g.input(identifier+"_beads",label,MaterialKind.SELECTION)
    recovered=g.product(identifier+"_recovered","bead-recovered material",MaterialKind.MOLECULAR,state=state)
    waste=g.product(identifier+"_waste","bead supernatant/washes",MaterialKind.WASTE)
    g.react(identifier,label,StepKind.CLEANUP,((template,InputRole.TEMPLATE),(reagent,InputRole.SELECTION)),(recovered,waste),changes_molecule=False)
    g.discard(waste,"physical supernatant/washes; molecular-class recovery not inferred")
    return recovered


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cell=g.input("cell","one cell containing RNA and chromatin",MaterialKind.BIOLOGICAL,primary=True)
    rna=g.product("rna","full-length RNA branch",MaterialKind.MOLECULAR,state=_single([Segment("RNA","[RNA]",placeholder=True)],"physically separated RNA fraction; source ref11 boundary"))
    chromatin=g.product("chromatin","matched accessible-chromatin branch",MaterialKind.MOLECULAR,state=_single([Segment("chromatin","[chromatin]",placeholder=True)],"matched nuclear fraction; source ref11 boundary"))
    g.react("separate","Capture RNA and accessible DNA from the same cell",StepKind.SPLIT,((cell,InputRole.TEMPLATE),),(rna,chromatin),changes_molecule=False)

    rt_primer=g.input("rt_primer","printed Oligo(dT) primer: SMART handle, N6 cell barcode, N12 UMI, tag1, T30VN",MaterialKind.OLIGO,state=Construct(LABELLED_CDNA.segments[:5],name="source low-throughput Oligo(dT) primer"))
    tso_state=Scene()
    tso_state.strand("TSO",[Segment("TSO DNA portion",TSO),Segment("three terminal riboguanosines","GGG",note="Source /rG//rG//rG/: these three residues are RNA, not DNA G")],label="printed TSO; DNA20 + RNA rGrGrG")
    tso=g.input("tso","printed TSO primer ending in three riboguanosines",MaterialKind.OLIGO,state=tso_state)
    rt=g.input("reverse_transcriptase","Maxima H Minus reverse transcriptase",MaterialKind.CATALYST)
    cdna=g.product("cdna","template-switched first-strand cDNA labelled at both ends; RNA fate unassigned",MaterialKind.MOLECULAR,state=_single(LABELLED_CDNA.segments,"switched ss cDNA; native first-cycle extent unassigned"))
    g.react("reverse_transcribe","Lyse (72 °C 10 min), reverse-transcribe with Oligo(dT) primer and template-switch onto the TSO (42 °C 90 min, 10 × 50/42 °C cycles, 70 °C 15 min)",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
    exo1=g.input("exo1","Exonuclease I",MaterialKind.CATALYST)
    rt_beads=g.input("rt_beads","1x AMPure XP beads",MaterialKind.SELECTION)
    rt_purified=g.product("rt_purified","1x bead-recovered RT material",MaterialKind.MOLECULAR,state=cdna.state)
    rt_waste=g.product("rt_waste","first RT cleanup supernatant/washes",MaterialKind.WASTE)
    g.react("purify_rt","Purify RT product with 1x AMPure XP",StepKind.CLEANUP,((cdna,InputRole.TEMPLATE),(rt_beads,InputRole.SELECTION)),(rt_purified,rt_waste),changes_molecule=False)
    g.discard(rt_waste,"physical AMPure supernatant/washes; RNA class-specific recovery unassigned")
    exo_treated=g.product("exo_treated","ExoI-treated RT material",MaterialKind.MOLECULAR,state=cdna.state)
    exo_debris=g.product("exo_debris","ExoI hydrolysates",MaterialKind.MOLECULAR)
    g.react("exo1_rt","Digest with ExoI (37 °C 30 min, 80 °C 10 min)",StepKind.DIGESTION,((rt_purified,InputRole.TEMPLATE),(exo1,InputRole.CATALYST)),(exo_treated,exo_debris))
    g.mark_terminal(Terminal(exo_debris.id,Disposition.UNRESOLVED,"Chemical hydrolysis is not physical removal; class-specific recovery at subsequent beads unassigned."))
    cdna_treated=g.product("cdna_treated","second 1x bead-recovered ss cDNA",MaterialKind.MOLECULAR,state=cdna.state)
    rt_second_waste=g.product("rt_second_waste","second RT cleanup supernatant/washes",MaterialKind.WASTE)
    g.react("purify_exo_rt","Purify ExoI-treated material with 1x AMPure XP",StepKind.CLEANUP,((exo_treated,InputRole.TEMPLATE),(rt_beads,InputRole.SELECTION)),(cdna_treated,rt_second_waste),changes_molecule=False)
    g.discard(rt_second_waste,"physical bead supernatant/washes")
    oligo_primer=g.input("oligo_primer","printed oligo primer",MaterialKind.OLIGO,state=OLIGO_PRIMER)
    ispcr_primer=g.input("ispcr_primer","printed TSO ISPCR primer",MaterialKind.OLIGO,state=TSO_ISPCR)
    cdna_polymerase=g.input("cdna_polymerase","KAPA HiFi HotStart ReadyMix",MaterialKind.CATALYST)
    labelled=g.product("labelled","amplified full-length cDNA labelled at both ends",MaterialKind.MOLECULAR,state=LABELLED_CDNA)
    g.react("amplify_cdna","Amplify cDNA with oligo primer and TSO ISPCR primer (18-25 cycles)",StepKind.PCR,((cdna_treated,InputRole.TEMPLATE),(oligo_primer,InputRole.STRUCTURAL),(ispcr_primer,InputRole.STRUCTURAL),(cdna_polymerase,InputRole.CATALYST)),(labelled,),geometry=PCRGeometry.unavailable(cdna_treated.state,reason="Actual switched first-strand cDNA is the PCR precursor; conditional duplex reference is its product, not its template. First-cycle completeness is unproved."))
    linker=g.input("linker","printed linker primer (reverse complement of the oligo(dT) handle followed by the TSO)",MaterialKind.OLIGO,state=LINKER)
    assembly_mix=g.input("assembly_mix","NEBuilder HiFi DNA Assembly Master Mix",MaterialKind.CATALYST)
    assembled=g.product("assembled","assembly mixture containing closed cDNA circles",MaterialKind.MOLECULAR,state=CIRCLE)
    g.react("circularize","Circularize 200 ng cDNA with 0.04 µM linker in NEBuilder HiFi (55 °C 60 min)",StepKind.CIRCULARIZATION,((labelled,InputRole.TEMPLATE),(linker,InputRole.STRUCTURAL),(assembly_mix,InputRole.CATALYST)),(assembled,),geometry=CircularizationGeometry.known("overlap assembly",statement="The printed linker is the reverse complement of the TSO-complement/oligo(dT)-handle junction; the Methods specify NEBuilder HiFi assembly."))
    linear_exo=g.input("linear_exo","Exonuclease I plus Exonuclease III:Lambda exonuclease (1:10) in NEB Buffer 1",MaterialKind.CATALYST)
    circle_beads=g.input("circle_beads","1.8x AMPure XP beads",MaterialKind.SELECTION)
    digested_circle=g.product("digested_circle","circle-containing exonuclease-treated material",MaterialKind.MOLECULAR,state=CIRCLE)
    linear_debris=g.product("linear_debris","hydrolysed linear cDNA/linker products; residual linear species unresolved",MaterialKind.MOLECULAR)
    g.react("digest_linear","Digest linear material with ExoI and ExoIII:Lambda (37 °C 30 min, 80 °C 10 min)",StepKind.DIGESTION,((assembled,InputRole.TEMPLATE),(linear_exo,InputRole.CATALYST)),(digested_circle,linear_debris))
    g.mark_terminal(Terminal(linear_debris.id,Disposition.UNRESOLVED,"Digestion does not physically discard polymers; hydrolysates and residual linear-class recovery unassigned."))
    circle=g.product("circle","bead-recovered circular cDNA",MaterialKind.MOLECULAR,state=CIRCLE)
    circle_waste=g.product("circle_waste","1.8x AMPure supernatant/washes",MaterialKind.WASTE)
    g.react("purify_circle","Purify digested circle material with 1.8x AMPure XP",StepKind.CLEANUP,((digested_circle,InputRole.TEMPLATE),(circle_beads,InputRole.SELECTION)),(circle,circle_waste),changes_molecule=False)
    g.discard(circle_waste,"physical bead supernatant/washes; residual linear-class recovery not presumed")
    rca_polymerase=g.input("rca_polymerase","phi29 DNA polymerase",MaterialKind.CATALYST)
    rca_primer=g.input("rca_primer","printed Phi29 primer",MaterialKind.OLIGO,state=PHI29_PRIMER)
    concatemer=g.product("concatemer","rolling-circle copies of full-length cDNA",MaterialKind.MOLECULAR,state=RCA)
    g.react("rolling_circle","Denature with Phi29 primer (95 °C 5 min, ice), extend with phi29 polymerase (37 °C 16 h, 65 °C 10 min)",StepKind.EXTENSION,((circle,InputRole.TEMPLATE),(rca_primer,InputRole.STRUCTURAL),(rca_polymerase,InputRole.CATALYST)),(concatemer,))

    concatemer=_cleanup(g,"purify_rca",concatemer,"0.6x AMPure XP after phi29 RCA",RCA)

    rna_transposome=g.input("rna_transposome","TTE Mix loaded Tn5 adapters (S5/ME end as drawn in Supplementary Fig. 3; adapter oligos not printed)",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (TTE Mix adapter oligos not printed; Supplementary Data 3 lists Vazyme Hyperactive pG-Tn5 S602-01)")
    rna_tn5=g.input("rna_tn5","Tn5 in 5x Tagment Buffer L (55 °C 10 min, 50 ng)",MaterialKind.CATALYST)
    rna_template=g.product("rna_template","tagmented rolling-circle cDNA fragments",MaterialKind.MOLECULAR,state=_fragments())
    rna_waste=g.product("rna_tagment_waste","non-library rolling-circle tagmentation products",MaterialKind.MOLECULAR)
    g.react("tagment_rna","Tagment the rolling-circle cDNA copies",StepKind.TAGMENTATION,((concatemer,InputRole.TEMPLATE),(rna_transposome,InputRole.STRUCTURAL),(rna_tn5,InputRole.CATALYST)),(rna_template,rna_waste));g.mark_terminal(Terminal(rna_waste.id,Disposition.UNRESOLVED,"Non-library end combinations do not establish physical removal or zero amplification; subsequent bead-class recovery unresolved."))
    rna_template=_cleanup(g,"purify_fragments",rna_template,"1x AMPure XP after TTE fragmentation",rna_template.state)
    tss_fragments=g.product("tss_fragments","10 µL mixed RCA fragments for TSS-directed PCR",MaterialKind.MOLECULAR,state=rna_template.state)
    tes_fragments=g.product("tes_fragments","10 µL mixed RCA fragments for TES-directed PCR",MaterialKind.MOLECULAR,state=rna_template.state)
    g.react("split_ends","Split tagmented cDNA for separate TSS and TES library construction",StepKind.SPLIT,((rna_template,InputRole.TEMPLATE),),(tss_fragments,tes_fragments),changes_molecule=False)
    tn5_i5=g.input("tn5_i5","published Tn5-i5 enrichment primer",MaterialKind.OLIGO,state=TN5_I5)
    anti_tag1=g.input("anti_tag1","published anti-tag1 TSS primer",MaterialKind.OLIGO,state=ANTI_TAG1)
    anti_tso=g.input("anti_tso","published anti-TSO TES primer",MaterialKind.OLIGO,state=ANTI_TSO)
    i7_anti_tag1=g.input("i7_anti_tag1","published i7-anti-tag1 TSS primer",MaterialKind.OLIGO,state=I7_ANTI_TAG1)
    i7_anti_tso=g.input("i7_anti_tso","published i7-anti-TSO TES primer",MaterialKind.OLIGO,state=I7_ANTI_TSO)
    enrich_polymerase=g.input("enrich_polymerase","KAPA HiFi HotStart ReadyMix",MaterialKind.CATALYST)
    tss_first=g.product("tss_first","TSS-enriched PCR1 fragments",MaterialKind.MOLECULAR,state=_amplified(TSS_PCR1))
    tes_first=g.product("tes_first","TES-enriched PCR1 fragments",MaterialKind.MOLECULAR,state=_amplified(TES_PCR1))
    tss_second=g.product("tss_second","TSS-enriched PCR2 fragments",MaterialKind.MOLECULAR,state=_amplified(TSS_PCR2))
    tes_second=g.product("tes_second","TES-enriched PCR2 fragments",MaterialKind.MOLECULAR,state=_amplified(TES_PCR2))
    g.react("tss_pcr1","Enrich TSS fragments with Tn5-i5 and anti-tag1 (6 cycles)",StepKind.PCR,((tss_fragments,InputRole.TEMPLATE),(tn5_i5,InputRole.STRUCTURAL),(anti_tag1,InputRole.STRUCTURAL),(enrich_polymerase,InputRole.CATALYST)),(tss_first,),geometry=PCRGeometry.unavailable(tss_fragments.state,reason="Actual preceding native material is the template. The source-primer-derived endpoint is separate; native gap/end and first-cycle transformations remain bounded."))
    g.react("tes_pcr1","Enrich TES fragments with Tn5-i5 and anti-TSO",StepKind.PCR,((tes_fragments,InputRole.TEMPLATE),(tn5_i5,InputRole.STRUCTURAL),(anti_tso,InputRole.STRUCTURAL),(enrich_polymerase,InputRole.CATALYST)),(tes_first,),geometry=PCRGeometry.unavailable(tes_fragments.state,reason="Actual preceding native material is the template. The source-primer-derived endpoint is separate; native gap/end and first-cycle transformations remain bounded."))
    tss_first=_cleanup(g,"purify_tss_pcr1",tss_first,"1x AMPure XP after TSS PCR1",tss_first.state)
    tes_first=_cleanup(g,"purify_tes_pcr1",tes_first,"1x AMPure XP after TES PCR1",tes_first.state)
    g.react("tss_pcr2","Add the TSS i7-anti-tag1 arm (5 cycles)",StepKind.PCR,((tss_first,InputRole.TEMPLATE),(tn5_i5,InputRole.STRUCTURAL),(i7_anti_tag1,InputRole.STRUCTURAL),(enrich_polymerase,InputRole.CATALYST)),(tss_second,),geometry=PCRGeometry.unavailable(tss_first.state,reason="Actual preceding native material is the template. The source-primer-derived endpoint is separate; native gap/end and first-cycle transformations remain bounded."))
    g.react("tes_pcr2","Add the TES i7-anti-TSO arm",StepKind.PCR,((tes_first,InputRole.TEMPLATE),(tn5_i5,InputRole.STRUCTURAL),(i7_anti_tso,InputRole.STRUCTURAL),(enrich_polymerase,InputRole.CATALYST)),(tes_second,),geometry=PCRGeometry.unavailable(tes_first.state,reason="Actual preceding native material is the template. The source-primer-derived endpoint is separate; native gap/end and first-cycle transformations remain bounded."))
    tss_second=_cleanup(g,"purify_tss_pcr2",tss_second,"1x AMPure XP after TSS PCR2",tss_second.state)
    tes_second=_cleanup(g,"purify_tes_pcr2",tes_second,"1x AMPure XP after TES PCR2",tes_second.state)
    n5_index=g.input("n5_index","published N5 index-primer pool",MaterialKind.OLIGO,state=Construct([Segment("P5",il.P5),Segment("N5 index","N"*8,"cbc",placeholder=True,feature=I5),*_s5_arm()],name="source N5 pool"))
    n7_index=g.input("n7_index","published N7 index-primer pool",MaterialKind.OLIGO,state=Construct([Segment("P7",il.P7),Segment("N7 index","N"*8,"cbc",placeholder=True,feature=I7),Segment("S7/ME",nx.ADAPTOR_S7)],name="source N7 pool"))
    index_polymerase=g.input("index_polymerase","KAPA HiFi index-PCR polymerase",MaterialKind.CATALYST)
    tss_indexed=g.product("tss_indexed","indexed TSS fragments before size separation",MaterialKind.MOLECULAR,state=TSS_LONG)
    tes_indexed=g.product("tes_indexed","indexed TES fragments before size separation",MaterialKind.MOLECULAR,state=TES_LONG)
    g.react("tss_index_pcr","Index the TSS library with N5/N7 primers (4-10 cycles)",StepKind.PCR,((tss_second,InputRole.TEMPLATE),(n5_index,InputRole.STRUCTURAL),(n7_index,InputRole.STRUCTURAL),(index_polymerase,InputRole.CATALYST)),(tss_indexed,),geometry=PCRGeometry.unavailable(tss_second.state,reason="Actual preceding native material is the template. The source-primer-derived endpoint is separate; native gap/end and first-cycle transformations remain bounded."))
    g.react("tes_index_pcr","Index the TES library with N5/N7 primers (4-10 cycles)",StepKind.PCR,((tes_second,InputRole.TEMPLATE),(n5_index,InputRole.STRUCTURAL),(n7_index,InputRole.STRUCTURAL),(index_polymerase,InputRole.CATALYST)),(tes_indexed,),geometry=PCRGeometry.unavailable(tes_second.state,reason="Actual preceding native material is the template. The source-primer-derived endpoint is separate; native gap/end and first-cycle transformations remain bounded."))
    size_beads=g.input("size_beads","AMPure XP: 0.55x retains long fragments; 0.55x-0.9x retains short fragments",MaterialKind.SELECTION)
    tss_long=g.product("tss_long","long-fragment TSS library",MaterialKind.LIBRARY,state=TSS_LONG)
    tss_short=g.product("tss_short","short-fragment TSS library",MaterialKind.LIBRARY,state=TSS_SHORT)
    tes_long=g.product("tes_long","long-fragment TES library",MaterialKind.LIBRARY,state=TES_LONG)
    tes_short=g.product("tes_short","short-fragment TES library",MaterialKind.LIBRARY,state=TES_SHORT)
    g.react("select_tss_sizes","Retain long and short TSS size fractions",StepKind.SELECTION,((tss_indexed,InputRole.TEMPLATE),(size_beads,InputRole.SELECTION)),(tss_long,tss_short),changes_molecule=False)
    g.react("select_tes_sizes","Retain long and short TES size fractions",StepKind.SELECTION,((tes_indexed,InputRole.TEMPLATE),(size_beads,InputRole.SELECTION)),(tes_long,tes_short),changes_molecule=False)
    for material,construct,primers in ((tss_long,TSS_LONG,TSS_PRIMERS),(tss_short,TSS_SHORT,TSS_PRIMERS),(tes_long,TES_LONG,TES_PRIMERS),(tes_short,TES_SHORT,TES_PRIMERS)):
        g.sequence(material,SequencingPlan(library=construct,primers=tuple(primers),acquired_roles=(),acquisition="unspecified",read_lengths=PE150))  # Methods: NovaSeq 6000 PE150; index cycles not stated

    atac_method=g.input("atac_method","Low-throughput accessibility preparation delegated to ref11; native reagents/assembly not assigned in this source cap",MaterialKind.CONDITION)
    atac_library=g.product("atac_library","paired low-throughput accessibility endpoint; assembly unassigned",MaterialKind.LIBRARY,state=ATAC)
    g.react("atac_prepare","Cited low-throughput ATAC preparation — native reaction details deferred",StepKind.MAGIC,((chromatin,InputRole.TEMPLATE),(atac_method,InputRole.CONDITION)),(atac_library,),reason="The inspected defining low-throughput route delegates matched ATAC preparation to predecessor ref11; this capped pass does not assign its native reactions, adapters or run identity.")
    g.sequence(atac_library,SequencingPlan(library=ATAC,primers=(),acquired_roles=(),acquisition="unspecified",read_lengths={},limitation="Low-throughput ATAC native assembly and acquisition are deferred; no generic shell or high-throughput run is substituted."))
    g.annotate(SourceDiscrepancy(
        target="rt_primer", title="Low-throughput barcode length in Supplementary Data 2",
        claims=(
            SourceClaim("Jiang et al. 2024, Supplementary Data 2, Oligo(dT) primer row",
                        "The oligo is printed as (N)6(N)12 before tag1."),
            SourceClaim("Jiang et al. 2024, Supplementary Data 2, table footnote",
                        "The footnote calls (N)6 a 16bp cell barcode."),
            SourceClaim("Jiang et al. 2024, Supplementary Fig. 3, low-throughput chemistry",
                        "The low-throughput drawing labels N6 as a 6bp cell barcode; Supplementary Fig. 2 separately labels high-throughput N16 as 16bp.")),
        model_choice="The low-throughput oligo and isoform library use a six-base cell barcode; the paired low-throughput ATAC library has no molecular cell-barcode segment.",
        evidence="The printed oligo length and low-throughput molecular drawing agree on N6. Cells are paired by physical separation in this route.",
        uncertainty="The supplementary table footnote appears to carry over the high-throughput N16 length; the correction status is unassigned in this bounded pass.",
        slot=AnnotationSlot.INTRO))
    g.declare_reviewed_scope(REVIEWED_TARGETS,SOURCE)
    return g.require_publishable()

REVIEWED_TARGETS=("cell","separate","rt_primer","tso","reverse_transcriptase","reverse_transcribe","exo1","rt_beads","purify_rt","exo1_rt","purify_exo_rt","oligo_primer","ispcr_primer","cdna_polymerase","amplify_cdna","linker","assembly_mix","circularize","linear_exo","circle_beads","digest_linear","purify_circle","rca_primer","rca_polymerase","rolling_circle","rna_transposome","rna_tn5","tagment_rna","split_ends","tn5_i5","anti_tag1","anti_tso","i7_anti_tag1","i7_anti_tso","enrich_polymerase","tss_pcr1","tes_pcr1","tss_pcr2","tes_pcr2","n5_index","n7_index","index_polymerase","tss_index_pcr","tes_index_pcr","size_beads","select_tss_sizes","select_tes_sizes","tss_long","tss_short","tes_long","tes_short","atac_method","atac_prepare","atac_library","purify_rca_beads","purify_rca","purify_fragments_beads","purify_fragments","purify_tss_pcr1_beads","purify_tss_pcr1","purify_tes_pcr1_beads","purify_tes_pcr1","purify_tss_pcr2_beads","purify_tss_pcr2","purify_tes_pcr2_beads","purify_tes_pcr2")
