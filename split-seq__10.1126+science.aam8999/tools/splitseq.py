"""Construct model shared by the individual SPLiT-seq and microSPLiT pages."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, MolecularRelation,
                            ProtocolGraph, SequencingPlan, StepKind)
import illumina as il
import seqprimers as sp

R1_HANDLE = "ACTGTGG"
R1_BARCODE = "ACTCGTAA"
SPLIT_R1_HANDLE = "AGGCCAGAGCATTCG"
SPLIT_R1_BARCODE = "AACGTGAT"  # Supplementary Table S12, Round1_01
R2_HEAD = "CATCGGCGTACGACT"
R2_BARCODE = "AACGTGAT"
R2_TAIL = "ATCCACGTGCTTGAG"
R3_HEAD = "CAGACGTGTGCTCTTCCGATCT"
R3_UMI_NT = 10
R3_BARCODE = "AACGTGAT"
R3_TAIL = "GTGGCCGATGTTTCG"
LINK2 = "CCACAGTCTCAAGCACG"
LINK3 = "TACGCCGATGCGAAACATCG"
TSO_BODY = "AAGCAGTGGTATCAACGCAGAGTGAAT"
TSO = TSO_BODY + "GGG"  # chemical 3' tail is rGrG+G
MICRO_PREAMP_FWD = R3_HEAD                    # BC_0062
MICRO_PREAMP_REV = "AAGCAGTGGTATCAACGCAGAGT"  # BC_0108
MICRO_LIBRARY_FWD = il.TRUSEQ_P5_FULL         # BC_0027
MICRO_I7 = "GATCTG"                          # BC_0076 / TSBC07
MICRO_LIBRARY_REV = il.P7 + MICRO_I7 + il.TRUSEQ_READ2
SPLIT_LIBRARY_FWD = "AATGATACGGCGACCACCGAGATCTACACTAGATCGCTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"  # BC_0118
SPLIT_LIBRARY_REV = "CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"  # BC_0076
SPLIT_LINK2 = "CGAATGCTCTGGCCTCTCAAGCACGTGGAT"  # BC_0215
SPLIT_LINK3 = "AGTCGTACGCCGATGCGAAACATCGGCCAC"  # BC_0060
MICRO_SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])

BC1_FEATURE=feature("cell_barcode_round1","cell_barcode","combinatorial",group="cell_barcode",part="RT round")
BC2_FEATURE=feature("cell_barcode_round2","cell_barcode","combinatorial",group="cell_barcode",part="ligation round 2")
BC3_FEATURE=feature("cell_barcode_round3","cell_barcode","combinatorial",group="cell_barcode",part="ligation round 3")
UMI_FEATURE=feature("umi","umi","random")
FEATURES={"barcode 1":BC1_FEATURE,"barcode 2":BC2_FEATURE,"barcode 3":BC3_FEATURE,"UMI":UMI_FEATURE}

def _s(name, top, tag=None, **kw):
    kw.setdefault("feature",FEATURES.get(name)); return Segment(name, top, tag, **kw)

def round1(protocol: str):
    if protocol == "microSPLiT":
        return [_s("round-1 handle", R1_HANDLE), _s("barcode 1", R1_BARCODE, "cbc"),
                _s("oligo-dT", "T" * 15), _s("VN", "VN", placeholder=True)]
    if protocol == "SPLiT-seq":
        return [_s("round-1 handle", SPLIT_R1_HANDLE),
                _s("barcode 1", SPLIT_R1_BARCODE, "cbc"),
                _s("oligo-dT", "T" * 15),
                _s("VN", "VN", placeholder=True)]
    raise ValueError(protocol)

def rt_scene(protocol: str):
    mrna = [_s("body", "XXXXXXXXXXXX", placeholder=True), _s("junction", "BB", placeholder=True),
            _s("poly(A)", "A" * 15)]
    sc = Scene(); sc.strand("mRNA", mrna, label="mRNA")
    sc.anneal("RT primer", round1(protocol), to="mRNA", pair=("oligo-dT", "poly(A)"),
              label="RT primer")
    sc.arrow("RT primer", "reverse transcriptase")
    return sc

def round2():
    return [_s("round-2 head", R2_HEAD), _s("barcode 2", R2_BARCODE, "cbc"),
            _s("round-2 tail", R2_TAIL)]

def round3():
    return [_s("round-3 handle", R3_HEAD, "r2"),
            _s("UMI", "N" * R3_UMI_NT, "umi", placeholder=True),
            _s("barcode 3", R3_BARCODE, "cbc"), _s("round-3 tail", R3_TAIL)]

def barcoded_first_strand(protocol: str, through: int) -> Construct:
    r1 = round1(protocol)
    segs = [*r1, _s("cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True)]
    if through >= 2:
        segs = [*round2(), *segs]
    if through >= 3:
        segs = [*round3(), *segs]
    return Construct(segs, name=f"{protocol} after barcode round {through}")

def splint_scene(round_no: int) -> Scene:
    if round_no == 2:
        left, right, source = round2(), [_s("round-1 handle", R1_HANDLE)], LINK2
        arms = [_s("round-2 tail arm", R2_TAIL[-10:]), right[0]]
    elif round_no == 3:
        left, right, source = round3(), [_s("round-2 head", R2_HEAD)], LINK3
        arms = [_s("round-3 tail arm", R3_TAIL[-10:]),
                _s("round-2 head arm", R2_HEAD[:10])]
    else: raise ValueError(round_no)
    junction = [*left, *right]
    # The source splint must equal the derived complement across the nick.
    derived = complement_segments(arms)
    if "".join(s.top for s in derived) != source:
        raise ValueError(f"round-{round_no} linker no longer bridges its two arms")
    sc = Scene(); sc.strand("joined strand", junction, label="barcoded cDNA")
    sc.anneal("linker", derived, to="joined strand",
              pair=(derived[0].name, right[0].name), label="splint")
    sc.junction("joined strand", left[-1].name, right[0].name)
    return sc

def splint_ligation_scene(round_no: int) -> Scene:
    """The two still-separate oligos aligned at the source-defined nick."""
    if round_no == 2:
        adapter, substrate, source = round2(), [_s("round-1 handle", R1_HANDLE)], LINK2
        arms = [_s("round-2 tail arm", R2_TAIL[-10:]), substrate[0]]
    elif round_no == 3:
        adapter, substrate, source = round3(), [
            _s("round-2 terminal arm", R2_HEAD[:10])], LINK3
        arms = [_s("round-3 tail arm", R3_TAIL[-10:]),
                _s("round-2 head arm", substrate[0].top)]
    else:
        raise ValueError(round_no)
    linker = complement_segments(arms)
    if "".join(s.top for s in linker) != source:
        raise ValueError(f"round-{round_no} linker no longer bridges its two arms")
    sc = Scene(); sc.strand("linker", linker, label=f"round-{round_no} splint")
    sc.anneal("substrate", substrate, to="linker",
              pair=(substrate[0].name, linker[0].name),
              label="previous-round cDNA", mod5="p")
    sc.anneal("adapter", adapter, to="linker",
              pair=(adapter[-1].name, linker[-1].name),
              label=f"round-{round_no} barcode oligo")
    sc.ligation_end("substrate", substrate[0].name, kind="splinted",
                    polarity="5-prime", owner="substrate",
                    text="5′-phosphate at the ligation nick")
    sc.ligation_end("adapter", adapter[-1].name, kind="splinted",
                    polarity="3-prime", owner="adapter",
                    text="3′-OH at the ligation nick")
    return sc

def final_boundary(protocol: str) -> Construct:
    inferred = protocol == "SPLiT-seq"
    inner = [_s("cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
             _s("barcode 1", "B" * 8, "cbc", placeholder=True, inferred=inferred),
             _s("barcode 2", "B" * 8, "cbc", placeholder=True),
             _s("barcode 3", "B" * 8, "cbc", placeholder=True),
             _s("UMI", "N" * R3_UMI_NT, "umi", placeholder=True)]
    if protocol == "microSPLiT":
        final = Construct([
            _s("P5", il.P5, "p5"),
            _s("TruSeq Read 1 remainder", il.TRUSEQ_READ1[4:], "r1"),
            *inner,
            _s("TruSeq Read 2 reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
            _s("i7 reverse complement", revcomp(MICRO_I7), "cbc",
               feature=feature("microsplit_sublibrary_i7", "sample_index", "whitelist",
                               whitelist="TSBC07-TSBC14")),
            _s("P7 reverse complement", il.P7_RC, "p7"),
        ], name="microSPLiT final library (representative TSBC07)")
        problems = sp.verify(final, MICRO_SEQ_PRIMERS,
                             required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
        if problems:
            raise ValueError("invalid microSPLiT library: " + "; ".join(problems))
        return final
    return Construct([_s("unpublished left library arm", "[LEFT ARM]", inferred=True,
                         placeholder=True), *inner,
                      _s("unpublished indexed right arm", "[RIGHT ARM + i7]", inferred=True,
                         placeholder=True)], name=f"{protocol} final library boundary")

def _validate():
    rt_scene("SPLiT-seq").rows(); rt_scene("microSPLiT").rows()
    splint_scene(2).rows(); splint_scene(3).rows()
    splint_ligation_scene(2).rows(); splint_ligation_scene(3).rows()
_validate()


def graph_for(protocol: str):
    if protocol not in ("SPLiT-seq","microSPLiT"):
        raise ValueError(protocol)
    doi=("10.1126/science.aam8999" if protocol=="SPLiT-seq"
         else "10.1126/science.aba5257")
    source=f"Defining source: https://doi.org/{doi}"
    graph=ProtocolGraph(protocol+" split-pool RNA library",authority=source)
    cells=graph.input("fixed_cells","fixed permeabilized cells",MaterialKind.BIOLOGICAL,
                      primary=True)
    rt_oligo=graph.input("round1_primer","round-1 barcoded oligo-dT primer",
                         MaterialKind.OLIGO,state=tuple(round1(protocol)),
                         source=source+(" Supplementary Table S12" if protocol=="SPLiT-seq" else ""))
    rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
    first=graph.product("round1_cdna","round-1 barcoded first-strand cDNA",
                        MaterialKind.MOLECULAR,state=barcoded_first_strand(protocol,1),
                        source=source+(" Supplementary Table S12" if protocol=="SPLiT-seq" else ""))
    graph.react("reverse_transcribe","Install barcode 1 during in-cell reverse transcription",
                StepKind.REVERSE_TRANSCRIPTION,
                ((cells,InputRole.TEMPLATE),(rt_oligo,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(first,))
    ligase=graph.input("dna_ligase","DNA ligase",MaterialKind.CATALYST)
    current=first
    for n,adapter_parts in ((2,round2()),(3,round3())):
        adapter=graph.input(f"round{n}_adapter",f"round-{n} split-pool barcode adapter",
                            MaterialKind.ADAPTER,state=tuple(adapter_parts),
                            source=source+(" Supplementary Table S12" if protocol=="SPLiT-seq" else ""))
        splint=graph.input(f"round{n}_splint",f"round-{n} ligation splint",
                          MaterialKind.ADAPTER,
                          state=Construct([_s("splint",SPLIT_LINK2 if n==2 else SPLIT_LINK3)])
                          if protocol=="SPLiT-seq" else None,
                          source=source+(" Supplementary Table S12" if protocol=="SPLiT-seq" else ""))
        product=graph.product(f"round{n}_cdna",f"cDNA after barcode round {n}",
                              MaterialKind.MOLECULAR,
                              state=barcoded_first_strand(protocol,n),
                              source=source+(" Supplementary Table S12" if protocol=="SPLiT-seq" else ""))
        geometry=(AdapterLigationGeometry.known(
                    "splinted",splint_ligation_scene(n),splint_scene(n),
                    statement=f"The microSPLiT oligo table establishes the round-{n} splinted nick.")
                  if protocol=="microSPLiT" else
                  AdapterLigationGeometry.unavailable(
                    f"The published SPLiT-seq Table S12 prints the oligos, but this model has not yet reconstructed the round-{n} annealing scene."))
        graph.react(f"ligate_round{n}",f"Ligate split-pool barcode round {n}",
                    StepKind.LIGATION,
                    ((current,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
                     (splint,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),
                    (product,),geometry=geometry)
        current=product
    beads=graph.input("streptavidin","streptavidin capture beads",MaterialKind.SELECTION)
    captured=graph.product("captured_cdna","captured round-3-biotinylated cDNA",
                           MaterialKind.MOLECULAR,
                           state=barcoded_first_strand(protocol,3))
    uncaptured=graph.product("uncaptured","uncaptured cellular material",MaterialKind.WASTE)
    graph.react("capture","Capture biotinylated barcoded cDNA",StepKind.AFFINITY_CAPTURE,
                ((current,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),
                (captured,uncaptured),changes_molecule=False,
                relations=(MolecularRelation.STRAND_ORDER,))
    graph.discard(uncaptured,"material not linked to the biotinylated round-3 oligo is washed away")
    tso=graph.input("tso","template-switch oligo",MaterialKind.OLIGO,
                    state=Construct([_s("BC_0127",TSO_BODY+"GGG")]),
                    source=source+(" Supplementary Table S12; 3' rGrG+G" if protocol=="SPLiT-seq" else " Supplementary Table S3; 3' rGrG+G"))
    switched=graph.product("switched_cdna","template-switched three-barcode cDNA",
                           MaterialKind.MOLECULAR)
    graph.react("template_switch","Install the distal amplification handle",
                StepKind.TEMPLATE_SWITCHING,
                ((captured,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(switched,))
    preamp_fwd=graph.input(
        "preamp_forward","round-3-handle cDNA amplification primer",MaterialKind.OLIGO,
        state=Construct([_s("BC_0062", MICRO_PREAMP_FWD)]),
        source=source+(" Supplementary Table S3" if protocol=="microSPLiT" else " Supplementary Table S12"))
    preamp_rev=graph.input(
        "preamp_reverse","template-switch-handle cDNA amplification primer",MaterialKind.OLIGO,
        state=Construct([_s("BC_0108", MICRO_PREAMP_REV)]),
        source=source+(" Supplementary Table S3" if protocol=="microSPLiT" else " Supplementary Table S12"))
    preamp_pol=graph.input("preamp_polymerase","cDNA PCR polymerase",MaterialKind.CATALYST)
    amplified=graph.product("amplified_cdna","amplified three-barcode cDNA",
                            MaterialKind.MOLECULAR,state=barcoded_first_strand(protocol,3))
    graph.react("cdna_pcr","Amplify captured template-switched cDNA",StepKind.PCR,
                ((switched,InputRole.TEMPLATE),(preamp_fwd,InputRole.STRUCTURAL),
                 (preamp_rev,InputRole.STRUCTURAL),(preamp_pol,InputRole.CATALYST)),
                (amplified,),geometry=PCRGeometry.unavailable(
                    barcoded_first_strand(protocol,3),
                    reason="The available source record does not establish both primer-binding scenes."))
    if protocol=="SPLiT-seq":
        transposome=graph.input("transposome","Nextera transposome",MaterialKind.ADAPTER)
        tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
        pcr_template=graph.product("tagged_cdna","tagmented three-barcode cDNA",
                                   MaterialKind.MOLECULAR)
        tag_waste=graph.product("tagment_waste","non-amplifiable tagmentation products",
                                MaterialKind.WASTE)
        graph.react("tagment","Tagment amplified barcoded cDNA",StepKind.TAGMENTATION,
                    ((amplified,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                     (tn5,InputRole.CATALYST)),(pcr_template,tag_waste))
        graph.discard(tag_waste,"products without the required two library primer sites")
        downstream_scope=("transposome","tn5","tagment")
    else:
        fragmentation=graph.input("fragmentation_mix","WGS enzymatic fragmentation mix",
                                  MaterialKind.CATALYST)
        fragmented=graph.product("fragmented_cdna","enzymatically fragmented barcoded cDNA",
                                 MaterialKind.MOLECULAR)
        graph.react("fragment","Fragment amplified cDNA",StepKind.FRAGMENTATION,
                    ((amplified,InputRole.TEMPLATE),(fragmentation,InputRole.CATALYST)),
                    (fragmented,))
        prepared=graph.product("prepared_fragments","end-repaired, dA-tailed cDNA fragments",
                               MaterialKind.MOLECULAR,inferred=True,
                               source=source+" (dA tail inferred from the disclosed 3′-T adapter)")
        graph.react("end_prepare","Prepare fragment ends for the 3-prime-T adapter",
                    StepKind.END_REPAIR,
                    ((fragmented,InputRole.TEMPLATE),(fragmentation,InputRole.CATALYST)),
                    (prepared,),source=source+" (end repair/dA tail inferred from adapter compatibility)")
        adapter=graph.input("library_adapter","pre-annealed BC_0243/BC_0244 adapter",
                            MaterialKind.ADAPTER)
        library_ligase=graph.input("library_ligase","WGS DNA ligase",MaterialKind.CATALYST)
        pcr_template=graph.product("adapter_ligated","adapter-ligated barcoded cDNA fragments",
                                   MaterialKind.MOLECULAR)
        graph.react("adapter_ligation","Ligate the 3-prime-T library adapter",
                    StepKind.LIGATION,
                    ((prepared,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
                     (library_ligase,InputRole.CATALYST)),(pcr_template,),
                    geometry=reported_adapter_end(
                        "BC_0243/BC_0244 adapter",kind="cohesive",polarity="3-prime",
                        terminal_sequence="T",
                        statement="The microSPLiT oligo table and fragmentation method disclose a duplex adapter with a terminal 3′ T."))
        downstream_scope=("fragmentation_mix","fragment","end_prepare","library_adapter",
                          "library_ligase","adapter_ligation")
    exact_micro = protocol == "microSPLiT"
    fwd_seq = MICRO_LIBRARY_FWD if exact_micro else SPLIT_LIBRARY_FWD
    rev_seq = MICRO_LIBRARY_REV if exact_micro else SPLIT_LIBRARY_REV
    fwd=graph.input("pcr_forward","library PCR forward primer",MaterialKind.OLIGO,
                    state=Construct([_s("BC_0027" if exact_micro else "BC_0118",fwd_seq)]),
                    source=source+(" Supplementary Table S3" if exact_micro
                                   else " Supplementary Table S12"))
    rev=graph.input("pcr_reverse","indexed library PCR reverse primer",MaterialKind.OLIGO,
                    state=Construct([_s("BC_0076 / TSBC07",rev_seq)]),
                    source=source+(" Supplementary Table S3" if exact_micro
                                   else " Supplementary Table S12"))
    pol=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    final=final_boundary(protocol)
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final,
                          inferred=not exact_micro,
                          source=source+(" Supplementary Table S3" if exact_micro
                                        else " (full geometry not yet reconstructed)"))
    graph.react("library_pcr","Build the indexed sequencing library",StepKind.PCR,
                ((pcr_template,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),
                 (rev,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(library,),
                geometry=(PCRGeometry.exact(
                    final, PCRPrimer("BC_0027", MICRO_LIBRARY_FWD),
                    PCRPrimer("BC_0076 / TSBC07", MICRO_LIBRARY_REV),
                    statement="Supplementary Table S3 prints the complete terminal PCR pair.")
                    if exact_micro else PCRGeometry.unavailable(
                        final,reason="The published primers are known; their sites on the full SPLiT-seq library have not yet been reconstructed.")))
    lengths=({"Read 1":74,"Read 2":86,"Index 1 (i7)":6}
             if protocol=="microSPLiT" else {})
    graph.sequence(library,SequencingPlan(
        library=final if exact_micro else None,
        primers=MICRO_SEQ_PRIMERS if exact_micro else (),
        acquired_roles=("Read 1","Read 2","Index 1 (i7)"),read_lengths=lengths,
        read_uses={"Read 1":"include","Read 2":"extract","Index 1 (i7)":"demultiplex"},
        limitation=("" if exact_micro else
                    "Published PCR primers are known, but the full library's sequencing-primer sites have not yet been reconstructed.")))
    graph.declare_reviewed_scope(("fixed_cells","round1_primer","reverse_transcriptase",
      "reverse_transcribe","dna_ligase","round2_adapter","round2_splint","ligate_round2",
      "round3_adapter","round3_splint","ligate_round3","streptavidin","capture","tso",
      "template_switch","preamp_forward","preamp_reverse","preamp_polymerase","cdna_pcr",
      *downstream_scope,"pcr_forward","pcr_reverse","pcr_polymerase","library_pcr",
      "library"),source)
    return graph.require_publishable()


def protocol_graph() -> ProtocolGraph:
    """The defining SPLiT-seq chemistry for source and molecular audits."""
    return graph_for("SPLiT-seq")


def start_split_pool(graph: ProtocolGraph, *, protocol: str, source: str):
    cells = graph.input("fixed_cells", "fixed permeabilized cells", MaterialKind.BIOLOGICAL,
                        primary=True)
    rt_oligo = graph.input("round1_primer", "round-1 barcoded oligo-dT primer",
                           MaterialKind.OLIGO, state=tuple(round1(protocol)), source=source)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    first = graph.product("round1_cdna", "round-1 barcoded first-strand cDNA",
                          MaterialKind.MOLECULAR, state=barcoded_first_strand(protocol, 1),
                          source=source)
    graph.react("reverse_transcribe", "Install barcode 1 during in-cell reverse transcription",
                StepKind.REVERSE_TRANSCRIPTION,
                ((cells, InputRole.TEMPLATE), (rt_oligo, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (first,))
    ligase = graph.input("dna_ligase", "DNA ligase", MaterialKind.CATALYST)
    return first, rt, ligase


def ligate_split_pool_round(graph: ProtocolGraph, current, ligase, *, n: int,
                            protocol: str, source: str):
    adapter = graph.input(f"round{n}_adapter", f"round-{n} split-pool barcode adapter",
                          MaterialKind.ADAPTER,
                          state=tuple(round2() if n == 2 else round3()), source=source)
    splint = graph.input(f"round{n}_splint", f"round-{n} ligation splint",
                         MaterialKind.ADAPTER, source=source)
    product = graph.product(f"round{n}_cdna", f"cDNA after barcode round {n}",
                            MaterialKind.MOLECULAR,
                            state=barcoded_first_strand(protocol, n), source=source)
    graph.react(f"ligate_round{n}", f"Ligate split-pool barcode round {n}",
                StepKind.LIGATION,
                ((current, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (splint, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)),
                (product,), geometry=AdapterLigationGeometry.known(
                    "splinted", splint_ligation_scene(n), splint_scene(n),
                    statement=f"The microSPLiT oligo table establishes the round-{n} splinted nick."))
    return product


def capture_split_pool_cdna(graph: ProtocolGraph, current, *, protocol: str):
    beads = graph.input("streptavidin", "streptavidin capture beads", MaterialKind.SELECTION)
    captured = graph.product("captured_cdna", "captured round-3-biotinylated cDNA",
                             MaterialKind.MOLECULAR,
                             state=barcoded_first_strand(protocol, 3))
    waste = graph.product("uncaptured", "uncaptured cellular material", MaterialKind.WASTE)
    graph.react("capture", "Capture biotinylated barcoded cDNA", StepKind.AFFINITY_CAPTURE,
                ((current, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (captured, waste), changes_molecule=False,
                relations=(MolecularRelation.STRAND_ORDER,))
    graph.discard(waste, "material not linked to the biotinylated round-3 oligo is washed away")
    return captured


def switch_and_amplify_split_pool(graph: ProtocolGraph, captured, rt, *,
                                  protocol: str, source: str):
    tso = graph.input("tso", "template-switch oligo", MaterialKind.OLIGO,
                      state=Construct([_s("BC_0127", TSO_BODY + "GGG")]),
                      source=source + " Supplementary Table S3; 3' rGrG+G")
    switched = graph.product("switched_cdna", "template-switched three-barcode cDNA",
                             MaterialKind.MOLECULAR)
    graph.react("template_switch", "Install the distal amplification handle",
                StepKind.TEMPLATE_SWITCHING,
                ((captured, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (switched,))
    forward = graph.input("preamp_forward", "round-3-handle cDNA amplification primer",
                          MaterialKind.OLIGO,
                          state=Construct([_s("BC_0062", MICRO_PREAMP_FWD)]),
                          source=source + " Supplementary Table S3")
    reverse = graph.input("preamp_reverse", "template-switch-handle cDNA amplification primer",
                          MaterialKind.OLIGO,
                          state=Construct([_s("BC_0108", MICRO_PREAMP_REV)]),
                          source=source + " Supplementary Table S3")
    polymerase = graph.input("preamp_polymerase", "cDNA PCR polymerase", MaterialKind.CATALYST)
    amplified = graph.product("amplified_cdna", "amplified three-barcode cDNA",
                              MaterialKind.MOLECULAR,
                              state=barcoded_first_strand(protocol, 3))
    graph.react("cdna_pcr", "Amplify captured template-switched cDNA", StepKind.PCR,
                ((switched, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(
                    barcoded_first_strand(protocol, 3),
                    reason="The available source record does not establish both primer-binding scenes."))
    return amplified
