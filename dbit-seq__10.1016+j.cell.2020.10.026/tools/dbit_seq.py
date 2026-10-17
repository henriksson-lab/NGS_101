"""Original Liu et al. DBiT-seq architecture."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from batch_ngs import nextera_library, seg
from adapter_ligation import AdapterLigationGeometry
from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow,
                      complement_segments, feature, oligo, revcomp)
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)


TITLE = "DBiT-seq — deterministic barcoding in tissue"
NOTES = "01_dbit-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1016/j.cell.2020.10.026">'
          'Liu et al., <i>Cell</i> (2020)</a>; exact oligos: '
          '<a href="https://doi.org/10.1016/j.xpro.2021.100532">'
          'Liu et al., <i>STAR Protocols</i> (2021)</a>.')
SUMMARY = ("Two perpendicular 50-channel microfluidic flows install an 8-nt A barcode "
           "during reverse transcription and ligate an 8-nt B barcode plus UMI at each "
           "channel intersection, creating up to 2,500 addressable tissue pixels.")
CAVEAT = ("The exact A and B whitelist entries vary by channel. This page uses A-1 and "
          "B-1 printed in the authors' 2021 STAR Protocols Table S2 to show the splinted "
          "ligation geometry. The defining 2020 Cell article points to its own Table S2; "
          "that original workbook has not yet been recovered here, so exact 2020/2021 "
          "oligo-table identity remains under source review. "
          "<b>Model correction (2026-10-17):</b> earlier versions of this page placed the "
          "spatial address at the Read 1 end and drew the Primer 2 handle in the wrong "
          "polarity. The published indexed P7 oligos end in the Barcode B handle, and the "
          "Cell analysis extracts UMI/Barcode A/Barcode B from Read 2, so the address is now "
          "read by Read 2 and the first-strand 3′ end carries the template-switch-oligo "
          "complement. The P5 oligo in Table S2 carries a fixed i5 (TAGATCGC), drawn here "
          "as a generic i5 window. The antibody-tag (ADT) branch (0.6× supernatant, 2× "
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
            seg("poly(dT)16", "T" * 16), seg("V", "V", placeholder=True),
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
    sc.anneal("A oligo", a_oligo(), to="mRNA", pair=("poly(dT)16", "poly(A)"),
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
                seg("poly(A)16", "A" * 16),
                seg("barcode A (rc)", revcomp(A1_BARCODE), "cbc", feature=BARCODE_A),
                seg("A linker (rc)", revcomp(A_LINKER), "me"),
                seg("B linker (rc)", revcomp(B_LINKER), "me"),
                seg("barcode B (rc)", revcomp(B1_BARCODE), "cbc", feature=BARCODE_B),
                seg("UMI (rc)", "N" * 10, "umi", placeholder=True, feature=UMI_FEATURE),
                seg("B PCR handle (rc)", revcomp(B_PCR_HANDLE), "r1")]
FINAL_LIBRARY, SEQ_PRIMERS = nextera_library(FINAL_INSERT, "DBiT-seq mRNA library")
if "".join(s.top for s in FINAL_LIBRARY if s.name in (
        "B PCR handle (rc)", "right mosaic end reverse complement")) != revcomp(
        P7_701[P7_701.index("AGATGTGTATAAGAGACAG"):]):
    raise ValueError("B-handle end differs from the published indexed P7 oligo")
FINAL_CAPTION = ("mRNA library: Nextera XT Amplicon Tagment Mix cuts the amplified cDNA; "
                 "PCR with the P5 oligo (Tn5 S5 end) and an indexed P7 oligo whose 3′ end is "
                 "the Barcode B handle keeps the transcript fragment joined to the A–B "
                 "address (STAR Protocols steps 60–66, Table S2).")
SEQUENCING_INTRO = ("HiSeq 4000 paired-end 100×100 (Cell STAR Methods). Read 1 enters the "
                    "Tn5-cut transcript end; Read 2 crosses the B handle and reads UMI, "
                    "Barcode B and Barcode A, as the Cell analysis extracts them from Read 2.")
READ_LENGTHS = {"Read 1": 100, "Index 1 (i7)": 8,
                "Index 2 (i5)": 8, "Read 2": 100}
READ_USES = {"Read 1": "include", "Index 1 (i7)": "demultiplex",
             "Index 2 (i5)": "demultiplex", "Read 2": "extract"}


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


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    tissue=graph.input("tissue_rna","polyadenylated RNA in fixed tissue",
                       MaterialKind.BIOLOGICAL,primary=True)
    a_oligo_input=graph.input("a_oligo","phosphorylated Barcode A oligo-dT",
                              MaterialKind.OLIGO,state=tuple(rt_scene().rows()))
    rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
    a_cdna=graph.product("a_cdna","A-coordinate first-strand cDNA",
                         MaterialKind.MOLECULAR,state=tuple(rt_scene().rows()))
    graph.react("reverse_transcribe","Flow Barcode A and reverse-transcribe tissue RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((tissue,InputRole.TEMPLATE),(a_oligo_input,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(a_cdna,))
    b_oligo_input=graph.input("b_oligo","biotinylated Barcode B oligo",
                              MaterialKind.ADAPTER,state=tuple(b_oligo()))
    splint=graph.input("ligation_splint","universal B-A ligation splint",
                       MaterialKind.ADAPTER,state=(Segment("universal ligation linker",LIGATION_SPLINT),))
    ligase=graph.input("ligase","T4 DNA ligase",MaterialKind.CATALYST)
    sealed=graph.product("sealed","sealed B-A spatial barcode junction",
                         MaterialKind.MOLECULAR,state=tuple(ligation_scene().rows()))
    graph.react("ligate_ba","Seal the Barcode B-to-A nick",StepKind.LIGATION,
                ((a_cdna,InputRole.TEMPLATE),(b_oligo_input,InputRole.STRUCTURAL),
                 (splint,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),
                (sealed,),geometry=ligation_geometry())
    streptavidin=graph.input("streptavidin","Dynabeads MyOne Streptavidin C1",
                            MaterialKind.SELECTION)
    lysis=graph.input("lysis","Proteinase K lysis (2 mg/mL, 2% SDS, 55 °C 2 h), PMSF quench",
                      MaterialKind.CONDITION)
    recovered=graph.product("recovered","recovered biotinylated spatial cDNA",
                            MaterialKind.MOLECULAR)
    unrecovered=graph.product("unrecovered","unbiotinylated tissue material",
                              MaterialKind.WASTE)
    graph.react("recover","Recover biotinylated spatial cDNA",StepKind.AFFINITY_CAPTURE,
                ((sealed,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION),
                 (lysis,InputRole.CONDITION)),
                (recovered,unrecovered))
    graph.discard(unrecovered,"material without the biotinylated Barcode B oligo is washed away")
    tso=graph.input("tso","template-switch oligo: "+TEMPLATE_SWITCH_PRIMER,MaterialKind.OLIGO,
                    state=(Segment("template-switch oligo",PRIMER_2+"GAATGGG"),),
                    source=SOURCE+" (STAR Protocols Supplementary Table S2)")
    switched=graph.product("switched","template-switched spatially barcoded cDNA",
                           MaterialKind.MOLECULAR,state=tuple(amplified_scene().rows()))
    graph.react("template_switch","Install the distal cDNA handle by template switching",
                StepKind.TEMPLATE_SWITCHING,
                ((recovered,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(switched,))
    p1=graph.input("primer1","Primer 1 at the Barcode B handle",MaterialKind.OLIGO,
                   state=(Segment("Primer 1",B_PCR_HANDLE),))
    p2=graph.input("primer2","Primer 2 at the template-switch handle",MaterialKind.OLIGO,
                   state=(Segment("Primer 2",PRIMER_2),))
    preamp_pol=graph.input("preamp_polymerase","Kapa HiFi HotStart (5 cycles on beads, then qPCR-monitored to plateau)",
                           MaterialKind.CATALYST)
    amplified=graph.product("amplified","amplified spatially barcoded cDNA",
                            MaterialKind.MOLECULAR,state=tuple(amplified_scene().rows()))
    graph.react("preamp","Amplify the spatially barcoded cDNA",StepKind.PCR,
                ((switched,InputRole.TEMPLATE),(p1,InputRole.STRUCTURAL),
                 (p2,InputRole.STRUCTURAL),(preamp_pol,InputRole.CATALYST)),
                (amplified,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The exact Primer 1/2 sequences are known, but the compact graph does not reconstruct both strands of the preamplification substrate."))
    transposome=graph.input("transposome","Nextera XT Amplicon Tagment Mix (0.75–1 ng cDNA, 55 °C 5 min)",
                            MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=graph.product("tagged","Nextera-tagged spatial cDNA",MaterialKind.MOLECULAR)
    tagment_waste=graph.product("tagment_waste","candidate fragments lacking the B-handle end",
                                MaterialKind.WASTE)
    graph.react("tagment","Tagment amplified spatial cDNA",StepKind.TAGMENTATION,
                ((amplified,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,tagment_waste))
    graph.discard(tagment_waste,"not exponentially amplified by the P5 + B-handle P7 pair; "
                  "not separately removed, only diluted and size-trimmed by the 0.7× KAPA Pure cleanup")
    i5=graph.input("i5_primer","P5 oligo (Table S2; S5/mosaic-end 3′ end, fixed i5)",MaterialKind.OLIGO)
    i7=graph.input("i7_primer","indexed P7 oligo 701–706 (Table S2; 3′ end = Barcode B handle)",
                   MaterialKind.OLIGO,state=(Segment("P7 oligo (701)",P7_701),))
    index_pol=graph.input("index_polymerase","index PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY)
    graph.react("index_pcr","Complete and index the Nextera library",StepKind.PCR,
                ((tagged,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
                 (i7,InputRole.STRUCTURAL),(index_pol,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="Index-member identities vary; the P7 oligo primes the Barcode B handle while the P5 oligo primes the Tn5 S5 end (STAR Protocols step 64, 12 cycles)."))
    graph.sequence(library,SequencingPlan(
        library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),read_lengths=READ_LENGTHS,
        read_uses=READ_USES))
    graph.declare_reviewed_scope((
        "tissue_rna","a_oligo","reverse_transcriptase","reverse_transcribe",
        "b_oligo","ligation_splint","ligase","ligate_ba","streptavidin","lysis",
        "recover","tso","template_switch","primer1","primer2",
        "preamp_polymerase","preamp","transposome","tn5","tagment",
        "i5_primer","i7_primer","index_polymerase","index_pcr","library"),SOURCE)
    return graph.require_publishable()
