"""Human NET-seq library model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import seqprimers as sp
from adapter_ligation import reported_adapter_end
from chemdraw import Construct, Row, Scene, Segment, feature, revcomp
from circular import CircularizationGeometry, circularize_ssdna
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "Human NET-seq — nascent RNA 3′ ends"
NOTES = "01_human-net-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1016/j.cell.2015.03.010">Mayer et al., <i>Cell</i> (2015)</a>; '
          '<a href="https://www.cecerelab.com/_files/ugd/594009_59c7fec32bcc42cf8c824804675fe319.pdf">original Extended Experimental Procedures</a>; '
          '<a href="https://doi.org/10.1038/nprot.2016.086">detailed protocol</a>.')
SUMMARY = ("Ligate a pre-adenylated, six-base-barcoded DNA linker to nascent RNA 3′ ends, fragment and select the linked RNA, reverse-transcribe, circularize the cDNA, deplete abundant mature-RNA products, and amplify for single-read sequencing.")
CAVEAT = ("The 2015 Extended Experimental Procedures print the N6 linker, modified RT primer and custom read primer. The detailed protocol names a commercial Illumina Index forward primer without its sequence; that outer PCR arm is a role token. The 2015 source describes depletion of 19 abundant RNAs, whereas the 2016 protocol lists twenty oligos. The original Table S1 remains unavailable, so the model does not assign the later full probe list to the 2015 experiment.")

LINKER = "N" * 6 + "CTGTAGGCACCATCAAT"
RT_PRIMER_5 = "ATCTCGTATGCCGTCTTCTGCTTG"
RT_PRIMER_MID = "CACTCA"
RT_PRIMER_3 = "TCCGACGATCATTGATGGTGCCTACAG"
RT_PRIMER_WRITTEN = ("/5Phos/" + RT_PRIMER_5 + "/iSp18/" + RT_PRIMER_MID
                     + "/iSp18/" + RT_PRIMER_3)
REVERSE_PCR = "CAAGCAGAAGACGGCATACGA"
SEQUENCING = "TCCGACGATCATTGATGGTGCCTACAG"
MOLECULAR_BARCODE = feature("molecular_barcode", "umi", "random")

def circular_cdna():
    linear = Construct([
        Segment("RT-primer 5-prime tract", RT_PRIMER_5),
        Segment("RT-primer internal spacer 1", "[iSp18]", placeholder=True),
        Segment("RT-primer middle tract", RT_PRIMER_MID),
        Segment("RT-primer internal spacer 2", "[iSp18]", placeholder=True),
        Segment("custom sequencing-primer site", SEQUENCING, "r1"),
        Segment("six-base molecular barcode", "N" * 6, "umi", placeholder=True,
                feature=MOLECULAR_BARCODE),
        Segment("nascent-RNA 3′-end cDNA", "X" * 36, placeholder=True),
    ], name="human NET-seq first-strand cDNA")
    return circularize_ssdna(linear, five_prime_phosphate=True, three_prime_oh=True)

CIRCLE = circular_cdna()

def final_library() -> tuple[Construct, tuple]:
    lib = Construct([
        Segment("Illumina P5/index PCR arm", "[Illumina P5 + index]", "p5+cbc", placeholder=True,
                feature=feature("sample_index_region", "sample_index", "unknown",
                                note="index sub-boundary is not printed")),
        Segment("custom sequencing-primer site", SEQUENCING, "r1"),
        Segment("six-base molecular barcode", "N" * 6, "umi", placeholder=True,
                feature=MOLECULAR_BARCODE),
        Segment("nascent-RNA 3′-end cDNA", "X" * 36, placeholder=True),
        Segment("P7-side reverse-primer complement", revcomp(REVERSE_PCR), "p7"),
    ], name="human NET-seq library")
    primer = sp.custom("Read 1", "oLSC006", SEQUENCING, "Mayer and Churchman 2016, Table 1")
    problems = sp.verify(lib, (primer,), required_roles=("Read 1",))
    if problems: raise ValueError("invalid human NET-seq endpoint: " + "; ".join(problems))
    return lib, (primer,)

FINAL_LIBRARY, SEQ_PRIMERS = final_library()
FINAL_CAPTION = ("The custom primer reads the six-base molecular barcode first, followed immediately by cDNA copied from the nascent RNA 3′ end.")
SEQUENCING_INTRO = ("oLSC006 is the published custom sequencing primer. The run is single-read; its first genomic base marks the engaged polymerase.")

def sections():
    lig = Construct([
        Segment("nascent RNA", "X" * 28, placeholder=True),
        Segment("ligation junction", "N", placeholder=True),
        Segment("six-base molecular barcode", "N" * 6, "umi", placeholder=True,
                feature=MOLECULAR_BARCODE),
        Segment("linker common region", LINKER[6:]),
    ], name="linker-ligated nascent RNA")
    sc = Scene.duplex(list(lig), label="schematic RNA/linker product", unpaired=tuple(s.name + "'" for s in lig))
    sc.junction("top", "nascent RNA", "ligation junction", "RNA–DNA linker ligation")
    sc.labels("top")
    return [
        ("Capture the nascent 3′ end", sc.rows(), "A 5′ pre-adenylated DNA linker with N6 and a blocked 3′ ddC is ligated to the RNA 3′-OH. ** marks the ligation."),
        ("Fragment, reverse-transcribe and circularize", [Row(chunks=[("partially hydrolyze linked RNA → gel select 35–100 nt RNA → RT → gel select 85–160 nt cDNA → CircLigase closes 5′-phosphate to 3′-OH", None, False)]), *CIRCLE.rows("CircLigase closure")], "Fragmentation follows linker ligation, preserving the native nascent-RNA 3′ end next to the N6 barcode."),
        ("Deplete and amplify", [Row(chunks=[("biotinylated depletion oligos remove abundant mature-RNA cDNAs → 6–12 cycle indexed PCR", None, False)])], "The minimal productive PCR cycle count is used; the informative product is about 150 nt."),
    ]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    rna=graph.input("nascent_rna","isolated nascent RNA",MaterialKind.MOLECULAR,
                    primary=True)
    linker=graph.input("linker","pre-adenylated N6 DNA linker with 3-prime ddC block",
                       MaterialKind.ADAPTER)
    ligase=graph.input("linker_ligase","RNA 3-prime linker ligase",MaterialKind.CATALYST)
    linked=graph.product("linked","N6-linker-ligated nascent RNA",MaterialKind.MOLECULAR)
    geometry=reported_adapter_end(
        "NET-seq linker",kind="single-strand",polarity="5-prime",
        terminal_sequence=LINKER,
        statement="The detailed protocol prints the pre-adenylated N6 linker and its blocked 3-prime end; its 5-prime end joins the nascent-RNA 3-prime hydroxyl.")
    graph.react("ligate_linker","Ligate the molecular-barcode linker to RNA 3-prime ends",
                StepKind.LIGATION,
                ((rna,InputRole.TEMPLATE),(linker,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(linked,),geometry=geometry)
    alkali=graph.input("alkaline_hydrolysis","partial alkaline hydrolysis",MaterialKind.CONDITION)
    fragmented=graph.product("fragmented_rna","fragmented linker-bearing RNA",MaterialKind.MOLECULAR)
    graph.react("fragment_rna","Partially hydrolyze linker-ligated RNA",StepKind.FRAGMENTATION,
                ((linked,InputRole.TEMPLATE),(alkali,InputRole.CONDITION)),(fragmented,))
    rna_size_rule=graph.input("rna_size_rule","35-to-100-nt RNA gel-selection rule",
                              MaterialKind.SELECTION)
    selected_rna=graph.product("selected_rna","35-to-100-nt linker-bearing RNA",MaterialKind.MOLECULAR)
    rejected_rna=graph.product("rejected_rna_sizes","off-size RNA fragments",MaterialKind.WASTE)
    graph.react("select_rna","Gel-select linker-bearing RNA fragments",StepKind.SELECTION,
                ((fragmented,InputRole.TEMPLATE),(rna_size_rule,InputRole.SELECTION)),
                (selected_rna,rejected_rna))
    graph.discard(rejected_rna,"RNA outside the 35-to-100-nt interval is excised away")
    rt_primer=graph.input("rt_primer","published NET-seq reverse-transcription primer",
                          MaterialKind.OLIGO)
    rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
    cdna=graph.product("cdna","linker-primed first-strand cDNA",MaterialKind.MOLECULAR,
                       state=CIRCLE.linear)
    graph.react("reverse_transcribe","Copy linker-ligated nascent RNA",StepKind.REVERSE_TRANSCRIPTION,
                ((selected_rna,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(cdna,))
    size_rule=graph.input("size_rule","85-to-160-nt gel-selection rule",
                          MaterialKind.SELECTION)
    selected=graph.product("selected","size-selected first-strand cDNA",MaterialKind.MOLECULAR,
                           state=CIRCLE.linear)
    rejected=graph.product("rejected_sizes","off-size cDNA",MaterialKind.WASTE)
    graph.react("size_select","Gel-select informative cDNA sizes",StepKind.SELECTION,
                ((cdna,InputRole.TEMPLATE),(size_rule,InputRole.SELECTION)),
                (selected,rejected))
    graph.discard(rejected,"cDNA outside the 85-to-160-nt interval is excised away")
    circligase=graph.input("circligase","CircLigase",MaterialKind.CATALYST)
    circle=graph.product("circle","circularized NET-seq cDNA",MaterialKind.MOLECULAR,
                         state=CIRCLE)
    graph.react("circularize","Circularize size-selected single-stranded cDNA",
                StepKind.CIRCULARIZATION,
                ((selected,InputRole.TEMPLATE),(circligase,InputRole.CATALYST)),
                (circle,),geometry=CircularizationGeometry.known(
                    "single-strand intramolecular ligation",
                    statement="CircLigase joins the phosphorylated cDNA 5-prime end to its 3-prime hydroxyl."))
    depletion=graph.input("depletion_oligos","biotinylated mature-RNA depletion oligos",
                          MaterialKind.SELECTION)
    depleted=graph.product("depleted","mature-RNA-depleted cDNA circles",
                           MaterialKind.MOLECULAR,state=CIRCLE)
    abundant=graph.product("abundant_cdna","captured abundant mature-RNA cDNAs",
                           MaterialKind.WASTE)
    graph.react("deplete","Remove abundant mature-RNA cDNA circles",StepKind.AFFINITY_CAPTURE,
                ((circle,InputRole.TEMPLATE),(depletion,InputRole.SELECTION)),
                (depleted,abundant))
    graph.discard(abundant,"biotinylated depletion probes remove abundant mature-RNA products")
    forward=graph.input("index_forward","commercial Illumina Index forward primer",
                        MaterialKind.OLIGO,inferred=True,
                        source=SOURCE+" (primer sequence not printed)")
    reverse=graph.input("reverse_primer","published NET-seq reverse PCR primer",
                        MaterialKind.OLIGO)
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY,inferred=True,
                          source=SOURCE+" (commercial forward PCR arm not printed)")
    graph.react("library_pcr","Amplify and index the depleted circles",StepKind.PCR,
                ((depleted,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The reverse primer is published, but the commercial Illumina Index forward-primer sequence is not."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "nascent_rna","linker","linker_ligase","ligate_linker",
        "alkaline_hydrolysis","fragment_rna","rna_size_rule","select_rna","rt_primer",
        "reverse_transcriptase","reverse_transcribe","size_rule","size_select",
        "circligase","circularize","depletion_oligos","deplete","index_forward",
        "reverse_primer","pcr_polymerase","library_pcr","library"),SOURCE)
    return graph.require_publishable()
