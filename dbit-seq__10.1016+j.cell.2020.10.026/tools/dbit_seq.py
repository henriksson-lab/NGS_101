"""Original Liu et al. DBiT-seq architecture."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import nextera as nx
import seqprimers as sp
from dataclasses import replace
from protocol_model import molecule, role
from adapter_ligation import AdapterLigationGeometry
from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow,
                      complement_segments, feature, oligo, revcomp)
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

seg=Segment

TITLE = "DBiT-seq — deterministic barcoding in tissue"
NOTES = "01_dbit-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1016/j.cell.2020.10.026">'
          'Liu et al., <i>Cell</i> (2020)</a>; exact oligos: '
          '<a href="https://doi.org/10.1016/j.xpro.2021.100532">'
          'Su et al., <i>STAR Protocols</i> (2021)</a>.')
SUMMARY = ("Two perpendicular 50-channel microfluidic flows install an 8-nt A barcode "
           "during reverse transcription and ligate an 8-nt B barcode plus UMI at each "
           "channel intersection, creating up to 2,500 addressable tissue pixels.")
CAVEAT = ("The exact A and B whitelist entries vary by channel. This page uses A-1 and "
          "B-1 printed in the authors' 2021 STAR Protocols Table S2 to show the splinted "
          "ligation geometry. The defining 2020 Cell article points to its own Table S2; "
          "that original workbook remains an unresolved source target, so exact 2020/2021 "
          "oligo-table identity remains under source review. "
          "<b>Model correction (2026-10-17):</b> earlier versions of this page placed the "
          "spatial address at the Read 1 end and drew the Primer 2 handle in the wrong "
          "polarity. The published indexed P7 oligos end in the Barcode B handle, and the "
          "Cell analysis extracts UMI/Barcode A/Barcode B from Read 2, so the address is now "
          "read by Read 2 and the first-strand 3′ end carries the template-switch-oligo "
          "complement. The P5 oligo in Table S2 carries a fixed i5 (TAGATCGC), drawn here "
          "as its source-defined fixed index. The detailed native graph selects the 2021 implementation reference (including its 0.8× cDNA cleanup), not an asserted reconstruction of the 2020 Table S2. The 2020 Cell Methods instead use 0.6× cDNA fractionation.  The antibody-tag (ADT) branch (0.6× supernatant, 2× "
          "AMPure, 10-cycle PCR with P7 and P5-cite primers, pooled 1:9) is not drawn.")

BARCODE_A = feature("spatial_barcode_a", "spatial_barcode", "whitelist",
                    group="spatial_barcode", part="flow A",
                    whitelist="DBiT-seq barcode A table")
BARCODE_B = feature("spatial_barcode_b", "spatial_barcode", "whitelist",
                    group="spatial_barcode", part="flow B",
                    whitelist="DBiT-seq barcode B table")
UMI_FEATURE = feature("umi", "umi", "random")

# STAR Protocols Table S2, written 5′→3′.
A_LINKER = "AGGCCAGAGCATTCG"
A1_BARCODE = "AACGTGAT"
A_POLYT_NT=15 # literal2021 TableS2 A1 tract after its8-base barcode
B_LINKER = "ATCCACGTGCTTGAG"
B1_BARCODE = "AACGTGAT"
B_PCR_HANDLE = "CAAGCGTTGGCTTCTCGCATCT"
PRIMER_2 = "AAGCAGTGGTATCAACGCAGAGT"
TEMPLATE_SWITCH_PRIMER = "AAGCAGTGGTATCAACGCAGAGTGAATrGrG+G"
PRIMER_3_CITE_SEQ = "CCTTGGCACCCGAGAATT*C*C"
LIGATION_SPLINT = "CGAATGCTCTGGCCTCTCAAGCACGTGGAT"

# DNA-equivalent of the template-switch oligo (rGrG+G written as GGG): the first strand
# copies it, so its 3' end carries the reverse complement, which Primer 2 then primes.
TSO_DNA = PRIMER_2 + "GAATGGG"
if not TEMPLATE_SWITCH_PRIMER.startswith(PRIMER_2 + "GAAT"):
    raise ValueError("Primer 2 must be the 5' part of the template-switch oligo")

if LIGATION_SPLINT != revcomp(B_LINKER + A_LINKER):
    raise ValueError("published DBiT splint must bridge the B and A linker arms")


def a_oligo() -> list[Segment]:
    return [seg("A linker", A_LINKER, "me"),
            seg("barcode A", A1_BARCODE, "cbc", feature=BARCODE_A),
            seg("poly(dT)15", "T" * A_POLYT_NT), seg("V", "V", placeholder=True),
            seg("N", "N", placeholder=True)]


def b_oligo() -> list[Segment]:
    return [seg("B PCR handle", B_PCR_HANDLE, "r1"),
            seg("UMI", "N" * 10, "umi", placeholder=True, feature=UMI_FEATURE),
            seg("barcode B", B1_BARCODE, "cbc", feature=BARCODE_B),
            seg("B linker", B_LINKER, "me")]


def rt_scene() -> Scene:
    mrna = [seg("transcript", "R" * 34, placeholder=True),
            seg("poly(A)", "A" * 18, placeholder=True)]
    sc = Scene(); sc.strand("mRNA", mrna, label="polyadenylated RNA")
    sc.anneal("A oligo", a_oligo(), to="mRNA", pair=("poly(dT)15", "poly(A)"),
              label="Barcode A-1", mod5="p", unpaired=("A linker", "barcode A", "V", "N"))
    sc.mark("A oligo", "A linker", "5′-phosphorylated A linker remains unpaired")
    sc.arrow("A oligo", "reverse transcriptase")
    return sc


def ligation_ready_scene() -> Scene:
    """The two oligos annealed side-by-side before the central nick is sealed."""
    a = [*a_oligo(), seg("cDNA", "X" * 34, placeholder=True)]
    b = b_oligo()
    splint = [seg("A splint arm", revcomp(A_LINKER), "me"),
              seg("B splint arm", revcomp(B_LINKER), "me")]
    if "".join(s.top for s in splint) != LIGATION_SPLINT:
        raise ValueError("derived DBiT splint differs from the published oligo")
    sc = Scene(); sc.strand("ligation splint", splint, label="universal ligation linker")
    sc.anneal("B oligo", b, to="ligation splint", pair=("B linker", "B splint arm"),
              label="Barcode B-1", mod5="biotin",
              unpaired=tuple(s.name for s in b if s.name != "B linker"))
    sc.anneal("A oligo/cDNA", a, to="ligation splint", pair=("A linker", "A splint arm"),
              label="Barcode A-1/cDNA", mod5="p",
              unpaired=tuple(s.name for s in a if s.name != "A linker"))
    sc.ligation_end("B oligo", "B linker", kind="splinted", polarity="3-prime",
                    owner="adapter", text="B 3′-OH at nick")
    sc.ligation_end("A oligo/cDNA", "A linker", kind="splinted", polarity="5-prime",
                    owner="adapter", text="A 5′-phosphate at nick")
    sc.footer("The barcode/UMI/primer bodies project beyond opposite splint ends",
              "ligation splint", "A splint arm")
    return sc


def ligation_geometry() -> AdapterLigationGeometry:
    return AdapterLigationGeometry.known(
        "splinted", ligation_ready_scene(),
        statement="Table S2 prints both 15-nt linker arms and their 30-nt universal splint.")


def ligation_scene() -> Scene:
    """Show the exact sealed B–A nick and both long, unpaired oligo bodies."""
    a_body = [s for s in a_oligo() if s.name != "A linker"]
    joined = [*b_oligo(), seg("A linker", A_LINKER, "me"), *a_body,
              seg("cDNA", "X" * 34, placeholder=True)]
    # The published splint is constructed from the arms it must pair with; the explicit
    # equality above prevents either source transcription from drifting independently.
    arms = [seg("B linker arm", B_LINKER, "me"), seg("A linker arm", A_LINKER, "me")]
    splint = complement_segments(arms)
    if "".join(s.top for s in splint) != LIGATION_SPLINT:
        raise ValueError("derived DBiT splint differs from the published oligo")

    unpaired = tuple(s.name for s in joined if s.name not in ("B linker", "A linker"))
    sc = Scene(); sc.strand("joined strand", joined, label="B-1 + A-1/cDNA", mod3="OH")
    sc.anneal("ligation splint", splint, to="joined strand",
              pair=("B linker arm'", "B linker"), label="universal ligation linker",
              unpaired=unpaired)
    sc.component("joined strand", "B PCR handle", "B linker", "Barcode B oligo")
    sc.component("joined strand", "A linker", "cDNA", "Barcode A oligo + cDNA")
    sc.junction("joined strand", "B linker", "A linker", "T4 ligation")
    sc.footer("B oligo body and A/cDNA body are overhangs beyond the 30-nt splint",
              "joined strand", "B PCR handle")
    return sc


def amplified_scene() -> Scene:
    first = [*b_oligo(), seg("A linker", A_LINKER, "me"),
             *[s for s in a_oligo() if s.name != "A linker"],
             seg("cDNA", "X" * 34, placeholder=True),
             seg("template-switch oligo complement", revcomp(TSO_DNA), "r2")]
    return Scene.duplex(first, label="spatially barcoded cDNA amplicon",
                        unpaired=("V'", "N'", "cDNA'"))


# STAR Protocols Table S2: the P5 oligo ends in the Nextera S5/mosaic end (Tn5 side) and
# every indexed P7 oligo ends ...AGACAG + CAAGCGTTGGCTTCTCGCATCT (the Barcode B handle).
# The B-handle end is therefore the P7/Read 2 end; the Tn5-cut transcript end is Read 1,
# matching Cell STAR Methods (Read 2 -> UMI, Barcode A, Barcode B; Read 1 -> genome).
P7_701 = ("CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
          + B_PCR_HANDLE)
FINAL_INSERT = [seg("cDNA insert (transcript sense)", "X" * 34, placeholder=True),
                seg("N complement", "N", placeholder=True),
                seg("V complement", "B", placeholder=True),
                seg("poly(A)15", "A" * A_POLYT_NT),
                seg("barcode A (rc)", revcomp(A1_BARCODE), "cbc", feature=BARCODE_A),
                seg("A linker (rc)", revcomp(A_LINKER), "me"),
                seg("B linker (rc)", revcomp(B_LINKER), "me"),
                seg("barcode B (rc)", revcomp(B1_BARCODE), "cbc", feature=BARCODE_B),
                seg("UMI (rc)", "N" * 10, "umi", placeholder=True, feature=UMI_FEATURE),
                seg("B PCR handle (rc)", revcomp(B_PCR_HANDLE), "r1")]
P5_OLIGO = il.P5 + "TAGATCGC" + nx.ADAPTOR_S5
I5_FEATURE=feature("sample_index_i5","sample_index","fixed",note="STAR2021 TableS2 P5 fixed TAGATCGC")
I7_FEATURE=feature("sample_index_i7","sample_index","whitelist",whitelist="STAR2021 TableS2 P7 701–706")
# Source-owned physical primers install the ends; no endpoint recipe chooses indexes.
FINAL_LIBRARY=Construct([
    Segment("P5",il.P5,"p5"),Segment("fixed i5","TAGATCGC",feature=I5_FEATURE),
    Segment("S5",nx.S5,"s5"),Segment("left mosaic end",nx.ME,"me"),*FINAL_INSERT,
    Segment("right mosaic end reverse complement",nx.ME_RC,"me"),
    Segment("S7 reverse complement",nx.S7_RC,"s7"),
    Segment("P7-701 index reverse complement",revcomp(P7_701[len(il.P7):len(il.P7)+8]),feature=I7_FEATURE),
    Segment("P7 reverse complement",il.P7_RC,"p7")],name="2021 TableS2 DBiT mRNA library reference (A1/B1/P7-701)")
SEQ_PRIMERS=(sp.NEXTERA["R1"],sp.NEXTERA["R2"])
errors=sp.verify(FINAL_LIBRARY,SEQ_PRIMERS,required_roles=("Read 1","Read 2"))
if errors:raise ValueError("DBiT source endpoint: "+"; ".join(errors))

if "".join(s.top for s in FINAL_LIBRARY if s.name in (
        "B PCR handle (rc)", "right mosaic end reverse complement")) != revcomp(
        P7_701[P7_701.index("AGATGTGTATAAGAGACAG"):]):
    raise ValueError("B-handle end differs from the published indexed P7 oligo")
FINAL_CAPTION = ("mRNA library: Nextera XT Amplicon Tagment Mix cuts the amplified cDNA; "
                 "PCR with the P5 oligo (Tn5 S5 end) and an indexed P7 oligo whose 3′ end is "
                 "the Barcode B handle keeps the transcript fragment joined to the A–B "
                 "address (STAR Protocols steps 60–66, Table S2).")
SEQUENCING_INTRO = ("Geometry-only 2021 TableS2 reference: Read1 can enter the Tn5-cut transcript end, and Read2 can cross the B handle into UMI/B/A identifiers. No acquired reads or cycle spans are assigned to this later implementation. The defining 2020 Cell article separately reports HiSeq4000 paired100/100; that historical run is not transferred to this reference.")
READ_LENGTHS={}
READ_USES={"Read 1":"include","Read 2":"extract"}
SEQUENCING_ENDING="2021 physical index values and compatible primer sites do not establish acquisition, run length or physical sequencing reagent identity. The source2020 run is a separate historical reference."


def workflow() -> Workflow:
    start = Scene(); start.strand("mRNA", [seg("transcript", "R" * 34, placeholder=True),
                                           seg("poly(A)", "A" * 18, placeholder=True)],
                                    label="polyadenylated RNA in fixed tissue")
    w = Workflow(MolecularState("Fixed tissue RNA", tuple(start.rows())))
    w.react("Flow phosphorylated Barcode A oligo-dT and reverse-transcribe",
            rt_scene().rows(), name="A-coordinate cDNA")
    w.ligate_adapter(
        "T4 DNA ligase seals the B–A nick", ligation_geometry(),
        ligation_scene().rows(), ready_name="B and A oligos annealed on the splint",
        name="sealed B–A ligation",
        note="The B 3′-OH is joined to the A 5′-phosphate at the central nick.")
    w.react("Recover biotinylated material, template-switch and PCR-amplify",
            amplified_scene().rows(), name="spatially barcoded cDNA")
    w.react("Tagment with Nextera XT and index by PCR",
            Scene.duplex(list(FINAL_LIBRARY), label=FINAL_LIBRARY.name,
                         unpaired=("cDNA insert (transcript sense)'", "N complement'",
                                   "V complement'", "UMI (rc)'")).rows(),
            name="sequencing library")
    return w


def oligos():
    return [
        oligo("Barcode A-1 (5′-phosphorylated)", a_oligo()),
        oligo("Barcode B-1 (5′ biotin)", b_oligo()),
        oligo("universal ligation linker", [seg("splint", LIGATION_SPLINT)]),
        oligo("Primer 1", [seg("B PCR handle", B_PCR_HANDLE)]),
        oligo("Primer 2", [seg("template-switch/PCR primer", PRIMER_2)]),
        f"<p>Template Switch Primer: 5′- {TEMPLATE_SWITCH_PRIMER} -3′</p>",
        f"<p>Primer 3 (optional antibody-tag branch): 5′- {PRIMER_3_CITE_SEQ} -3′</p>",
    ]


def first_cdna_parts(*, switched=False):
    parts=[*a_oligo(),role("transcript-derived cDNA; variable length")]
    if switched:parts.append(Segment("template-switch oligo complement",revcomp(TSO_DNA)))
    return parts


def joined_cdna_parts(*, switched=False):
    return [*b_oligo(),*first_cdna_parts(switched=switched)]


def protocol_graph():
    detailed="Su et al. STAR Protocols2021 PMC8132129 detailed implementation; exact TableS2 examples; not established2020 TableS2 identity"
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    tissue=graph.input("tissue_rna","polyadenylated RNA in fixed tissue (mRNA route)",MaterialKind.BIOLOGICAL,
        state=molecule("tissue RNA",(role("transcript"),role("poly(A)"))),primary=True)
    a=graph.input("a_oligo","2021 TableS2 5′phosphorylated BarcodeA1",MaterialKind.OLIGO,
        state=molecule("BarcodeA1",tuple(a_oligo())),source=detailed)
    rt=graph.input("reverse_transcriptase","Maxima H Minus reverse transcriptase",MaterialKind.CATALYST,source=detailed)
    first=graph.product("a_cdna","single-stranded A-coordinate cDNA",MaterialKind.MOLECULAR,
        state=molecule("first-strand A-cDNA",tuple(first_cdna_parts())))
    graph.react("reverse_transcribe","2021 reference — flow A oligo and reverse-transcribe",StepKind.REVERSE_TRANSCRIPTION,
        ((tissue,InputRole.TEMPLATE),(a,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,),source=detailed)
    b=graph.input("b_oligo","2021 TableS2 5′biotin BarcodeB1",MaterialKind.ADAPTER,state=molecule("BarcodeB1",tuple(b_oligo())),source=detailed)
    splint=graph.input("ligation_splint","universal B–A linker",MaterialKind.OLIGO,
        state=molecule("splint",(Segment("linker",LIGATION_SPLINT),)),source=detailed)
    ligase=graph.input("ligase","T4 DNA ligase",MaterialKind.CATALYST)
    sealed=graph.product("sealed","covalently sealed B–A spatially barcoded cDNA",MaterialKind.MOLECULAR,state=ligation_scene())
    graph.react("ligate_ba","Seal B3′OH to A5′phosphate at the splinted nick",StepKind.LIGATION,
        ((first,InputRole.TEMPLATE),(b,InputRole.STRUCTURAL),(splint,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(sealed,),geometry=ligation_geometry(),source=detailed)
    beads=graph.input("streptavidin","Dynabeads MyOne StreptavidinC1; twoB&W washes and Tris/Tween rinse",MaterialKind.SELECTION)
    lysis=graph.input("lysis","2mg/mL ProteinaseK,2%SDS,55°C2h; PMSF quench",MaterialKind.CONDITION)
    recovered=graph.product("recovered","bead-retained biotinylated ss cDNA",MaterialKind.MOLECULAR,
        state=molecule("recovered B–A cDNA",tuple(joined_cdna_parts())))
    waste=graph.product("unrecovered","removed lysate/wash fractions",MaterialKind.WASTE)
    graph.react("recover","Lyse tissue and capture biotinylated spatial cDNA",StepKind.AFFINITY_CAPTURE,
        ((sealed,InputRole.TEMPLATE),(beads,InputRole.SELECTION),(lysis,InputRole.CONDITION)),(recovered,waste),source=detailed)
    graph.discard(waste,"Source capture/wash transfers remove supernatants; component recovery is not quantitatively assigned")
    tso=graph.input("tso","TableS2 TSO rGrG+G; Maxima buffer/Ficoll/dNTP",MaterialKind.OLIGO,
        state=molecule("TSO",(Segment("DNA part",PRIMER_2+"GAAT"),Segment("modified guanosines","GGG",note="rGrG+G"))),source=detailed)
    switched=graph.product("switched","single-stranded switched B–A cDNA",MaterialKind.MOLECULAR,
        state=molecule("switched ss cDNA",tuple(joined_cdna_parts(switched=True))))
    graph.react("template_switch","Template-switch recovered cDNA (RT30min,42°C90min)",StepKind.TEMPLATE_SWITCHING,
        ((recovered,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(switched,),source=detailed)
    p1=graph.input("primer1","TableS2 Primer1",MaterialKind.OLIGO,state=Construct([Segment("Primer1",B_PCR_HANDLE)]),source=detailed)
    p2=graph.input("primer2","TableS2 Primer2",MaterialKind.OLIGO,state=Construct([Segment("Primer2",PRIMER_2)]),source=detailed)
    pol=graph.input("preamp_polymerase","KAPA HiFi HotStart; five cycles on beads, then monitored amplification",MaterialKind.CATALYST)
    amplicon=graph.product("amplified","duplex B–A cDNA amplicon",MaterialKind.MOLECULAR,state=amplified_scene())
    graph.react("preamp","Amplify with Primer1/Primer2",StepKind.PCR,
        ((switched,InputRole.TEMPLATE),(p1,InputRole.STRUCTURAL),(p2,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(amplicon,),
        geometry=PCRGeometry.unavailable(switched.state,reason="Actual incoming cDNA is single stranded; native first-cycle primer pairing is not implemented. Known final amplicon is not substituted as its template."),source=detailed)
    clean_beads=graph.input("cdna_beads","2021 steps45–59:0.8×KAPA Pure,15µl elution",MaterialKind.SELECTION)
    clean=graph.product("clean_cdna","purified cDNA for tagmentation",MaterialKind.MOLECULAR,state=amplicon.state)
    clean_waste=graph.product("cdna_waste","0.8×supernatant and ethanol washes",MaterialKind.WASTE)
    graph.react("cdna_cleanup","2021 implementation: purify cDNA with0.8×beads",StepKind.CLEANUP,
        ((amplicon,InputRole.TEMPLATE),(clean_beads,InputRole.SELECTION)),(clean,clean_waste),changes_molecule=False,source=detailed)
    graph.discard(clean_waste,"2021 steps49–53 physically discard supernatant/ethanol; not substituted for2020 Cell0.6×modality split")
    transposome=graph.input("transposome","Nextera XT Amplicon Tagment Mix;0.75–1ng,55°C5min; NT buffer5minRT",MaterialKind.ADAPTER,source=detailed)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=graph.product("tagged","tagmented cDNA population; native ends unassigned",MaterialKind.MOLECULAR,
        state=molecule("actual pre-index tagmented population",(role("tagmented cDNA with possible B-handle retention; native end geometry unassigned"),)))
    off=graph.product("tagment_other","candidate fragments lacking the B handle",MaterialKind.MOLECULAR)
    graph.react("tagment","Tagment the amplified cDNA",StepKind.TAGMENTATION,
        ((clean,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,off),source=detailed)
    from protocol_graph import Disposition,Terminal
    graph.mark_terminal(Terminal(off.id,Disposition.UNRESOLVED,"Lack of exponential amplification is not physical removal; recovery through later cleanup unassigned."))
    ip5=graph.input("i5_primer","2021 P5 fixed i5TAGATCGC",MaterialKind.OLIGO,state=Construct([Segment("P5 oligo",P5_OLIGO)]),source=detailed)
    ip7=graph.input("i7_primer","2021 P7-701 representative",MaterialKind.OLIGO,state=Construct([Segment("P7-701",P7_701)]),source=detailed)
    ipol=graph.input("index_polymerase","PCR master mix,12cycles per2021step65",MaterialKind.CATALYST)
    raw=graph.product("index_pcr_product","2021 primer-derived indexed mRNA library",MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
    graph.react("index_pcr","2021 reference — amplify with P5/P7-701",StepKind.PCR,
        ((tagged,InputRole.TEMPLATE),(ip5,InputRole.STRUCTURAL),(ip7,InputRole.STRUCTURAL),(ipol,InputRole.CATALYST)),(raw,),
        geometry=PCRGeometry.unavailable(tagged.state,reason="Actual tagmented incoming population and first-cycle ends are unassigned. The source-primer-derived final reference is not a pre-PCR duplex."),source=detailed)
    final_beads=graph.input("final_beads","2021 step66:0.7×KAPA Pure cleanup",MaterialKind.SELECTION)
    final=graph.product("library","purified2021 reference mRNA sequencing library",MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    final_waste=graph.product("final_waste","0.7×supernatant and ethanol wash material",MaterialKind.WASTE)
    graph.react("index_cleanup","Purify indexed library with0.7×KAPA Pure",StepKind.CLEANUP,
        ((raw,InputRole.TEMPLATE),(final_beads,InputRole.SELECTION)),(final,final_waste),changes_molecule=False,source=detailed)
    graph.discard(final_waste,"2021step66 follows physical supernatant/wash removals")
    graph.sequence(final,SequencingPlan(library=FINAL_LIBRARY,primers=SEQ_PRIMERS,acquisition="unspecified",read_uses=READ_USES,
        limitation="Geometry-only2021 TableS2 endpoint reference; acquired reads/cycles unassigned. Published2020PE100/100 remains a separate historical narrative, not a run declaration for this later implementation. Exact2020oligo-table identity unresolved."))
    graph.declare_reviewed_scope(("tissue_rna","a_oligo","reverse_transcriptase","reverse_transcribe","b_oligo","ligation_splint","ligase","ligate_ba","streptavidin","lysis","recover","tso","template_switch","primer1","primer2","preamp_polymerase","preamp","cdna_beads","cdna_cleanup","transposome","tn5","tagment","i5_primer","i7_primer","index_polymerase","index_pcr","final_beads","index_cleanup","library"),detailed)
    return graph.require_publishable()
