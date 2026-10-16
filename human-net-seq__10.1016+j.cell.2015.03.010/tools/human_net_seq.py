"""Human NET-seq library model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import seqprimers as sp
from adapter_ligation import reported_adapter_end
from chemdraw import Construct, Scene, Segment, feature, revcomp
from circular import CircularizationGeometry, circularize_ssdna
from pcr import PCRGeometry
from protocol_model import role
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind, Highlight)

TITLE = "Human NET-seq — nascent RNA 3′ ends"
NOTES = "01_human-net-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1016/j.cell.2015.03.010">Mayer et al., <i>Cell</i> (2015)</a>; '
          '<a href="https://www.cecerelab.com/_files/ugd/594009_59c7fec32bcc42cf8c824804675fe319.pdf">original Extended Experimental Procedures</a>; '
          '<a href="https://doi.org/10.1002/0471142727.mb0414s98">inherited 2012 procedure</a>; '
          '<a href="https://doi.org/10.1038/nprot.2016.047">2016 successor protocol</a>.')
SUMMARY = ("Ligate a pre-adenylated, six-base-barcoded DNA linker to nascent RNA 3′ ends, fragment and select the linked RNA, reverse-transcribe, circularize the cDNA, deplete abundant mature-RNA products, and amplify for single-read sequencing.")
CAVEAT = ("The selected 2015 implementation prints its N6 linker, modified RT primer and custom read primer. "
          "The ordered outer PCR arm and selected sample index remain unassigned model data; no exact commercial cap is supplied. "
          "The separate 2015 Table S1 depletion-probe original is still route-blocked. The twenty probes in the 2016 successor "
          "are not substituted for the nineteen targets in the 2015 experiment. Variable biological regions have no fixed example length. Possible primer-carryover products are represented by roles; their exact 2015 structure and residual abundance are unassigned.")
READ_LENGTHS = {"Read 1": 75}
ACQUIRED_ROLES = ("Read 1",)
READ_SPECS = (sp.read_spec("Read 1", 75, "N6 molecular barcode, then nascent-RNA-end cDNA", "include"),
              sp.read_spec("Index 1 (i7)", None, "acquisition and selected index unassigned", "undetermined"),
              sp.read_spec("Index 2 (i5)", None, "acquisition and selected index unassigned", "undetermined"))

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
        role("variable nascent-RNA-end cDNA"),
    ], name="human NET-seq first-strand cDNA")
    return circularize_ssdna(linear, five_prime_phosphate=True, three_prime_oh=True)

CIRCLE = circular_cdna()

def final_library() -> tuple[Construct, tuple]:
    lib = Construct([
        role("unassigned outer PCR cap/index region", molecular_feature=feature("sample_index_region", "sample_index", "unknown", note="selected ordered arm and index boundary unassigned")),
        Segment("custom sequencing-primer site", SEQUENCING, "r1"),
        Segment("six-base molecular barcode", "N" * 6, "umi", placeholder=True,
                feature=MOLECULAR_BARCODE),
        role("variable nascent-RNA-end cDNA"),
        Segment("retained RT-primer circle-junction tract", RT_PRIMER_5[:RT_PRIMER_5.index(revcomp(REVERSE_PCR))]),
        Segment("reverse-primer complement", revcomp(REVERSE_PCR), "p7"),
    ], name="human NET-seq library")
    primer = sp.custom("Read 1", "oLSC006", SEQUENCING, "Mayer 2015 Extended Experimental Procedures, oLSC006")
    problems = sp.verify(lib, (primer,), required_roles=("Read 1",))
    if problems: raise ValueError("invalid human NET-seq endpoint: " + "; ".join(problems))
    return lib, (primer,)

FINAL_LIBRARY, SEQ_PRIMERS = final_library()
FINAL_CAPTION = ("The custom primer reads the six-base molecular barcode first, followed immediately by cDNA copied from the nascent RNA 3′ end.")
SEQUENCING_INTRO = ("oLSC006 is the published custom sequencing primer. The run is single-read; its first genomic base marks the engaged polymerase.")

def linker_state():
    scene = Scene()
    scene.strand("linker", [Segment("six-base molecular barcode", "N" * 6, "umi", placeholder=True,
                                 feature=MOLECULAR_BARCODE), Segment("linker common region", LINKER[6:])],
                 label="source DNA linker", mod5="App", mod3="ddC")
    return scene


def rt_state():
    scene = Scene()
    scene.strand("RT primer", list(CIRCLE.linear)[:-2], label="oLSC007", mod5="Phos")
    return scene


def oligos():
    return ["2015 linker (5′→3′): /5App/" + LINKER + "/3ddC/",
            "oLSC007 (5′→3′): " + RT_PRIMER_WRITTEN,
            "oLSC006 (5′→3′): " + SEQUENCING]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    rna=graph.input("nascent_rna","isolated nascent RNA",MaterialKind.MOLECULAR,
                    primary=True, state=Construct([role("variable nascent RNA")], name="isolated nascent RNA"))
    linker=graph.input("linker","pre-adenylated N6 DNA linker with 3-prime ddC block",
                       MaterialKind.ADAPTER, state=linker_state())
    ligase=graph.input("linker_ligase","truncated T4 RNA ligase 2 (inherited 2012 method)",MaterialKind.CATALYST)
    linked_state=Construct([role("variable nascent RNA"), Segment("six-base molecular barcode", "N"*6, "umi", placeholder=True, feature=MOLECULAR_BARCODE), Segment("linker common region", LINKER[6:])], name="linked nascent RNA")
    linked=graph.product("linked","N6-linker-ligated nascent RNA",MaterialKind.MOLECULAR, state=linked_state)
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
                          MaterialKind.OLIGO, state=rt_state())
    rt=graph.input("reverse_transcriptase","SuperScript III (inherited 2012 method)",MaterialKind.CATALYST)
    cdna=graph.product("cdna","linker-primed first-strand cDNA",MaterialKind.MOLECULAR,
                       state=CIRCLE.linear)
    graph.react("reverse_transcribe","Copy linker-ligated nascent RNA",StepKind.REVERSE_TRANSCRIPTION,
                ((selected_rna,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(cdna,))
    rna_removal=graph.input("rna_removal", "NaOH hydrolysis after RT, then neutralization", MaterialKind.CONDITION)
    clean_cdna=graph.product("clean_cdna", "cDNA after removal of its RNA template", MaterialKind.MOLECULAR, state=CIRCLE.linear)
    graph.react("remove_rna", "Hydrolyze RNA after reverse transcription", StepKind.CHEMICAL_CONVERSION,
                ((cdna,InputRole.TEMPLATE),(rna_removal,InputRole.CONDITION)), (clean_cdna,))
    size_rule=graph.input("size_rule","inherited cDNA gel-selection rule; numerical 2015 bounds unassigned",
                          MaterialKind.SELECTION)
    selected=graph.product("selected","size-selected first-strand cDNA",MaterialKind.MOLECULAR,
                           state=CIRCLE.linear)
    rejected=graph.product("rejected_sizes","off-size cDNA",MaterialKind.WASTE)
    graph.react("size_select","Gel-select informative cDNA sizes",StepKind.SELECTION,
                ((clean_cdna,InputRole.TEMPLATE),(size_rule,InputRole.SELECTION)),
                (selected,rejected))
    graph.discard(rejected,"cDNA outside the chosen informative cDNA gel band is excised away")
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
    forward=graph.input("index_forward","index-forward PCR strand — actual 2015 ordered member unassigned",
                        MaterialKind.OLIGO, source=SOURCE+" (ordered strand unassigned in this model)")
    reverse=graph.input("reverse_primer","published NET-seq reverse PCR primer",
                        MaterialKind.OLIGO, state=Construct([Segment("reverse PCR primer", REVERSE_PCR)], name="oNTI231"))
    polymerase=graph.input("pcr_polymerase","Phusion DNA polymerase (inherited 2012 method)",MaterialKind.CATALYST)
    amplified=graph.product("amplified", "amplified NET-seq target", MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    background=graph.product("primer_background", "possible primer-carryover PCR products; exact 2015 state unassigned", MaterialKind.MOLECULAR,
                             state=Construct([role("possible primer-derived PCR background")], name="primer background"))
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY, source=SOURCE+" (outer PCR arm unassigned)")
    graph.react("library_pcr","Amplify and index the depleted circles",StepKind.PCR,
                ((depleted,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (amplified,background),geometry=PCRGeometry.unavailable(
                    CIRCLE.linear, reason="The ordered forward PCR strand and its binding site have not been assigned from the inherited original procedure; the outer arm remains a role. No publication absence is accepted."))
    final_gel=graph.input("final_gel", "post-PCR PAGE selection in the inherited 2012 procedure", MaterialKind.SELECTION)
    rejected_pcr=graph.product("rejected_pcr", "gel-excluded PCR products, including primer carryover when present", MaterialKind.WASTE)
    graph.react("select_library", "Excise the informative PCR band; avoid primer-carryover bands", StepKind.SELECTION,
                ((amplified,InputRole.TEMPLATE),(background,InputRole.TEMPLATE),(final_gel,InputRole.SELECTION)), (library,rejected_pcr))
    graph.discard(rejected_pcr, "The inherited 2012 procedure physically excludes primer-carryover and undesired PCR gel bands; no complete removal is assumed.")
    graph.annotate(Highlight("library_pcr", "Predecessor primers and background are version-specific",
        "The 2012 original prints oNTI230/oNTI231 and indexed oMHL001/002 alternatives. This does not identify the selected 2015 forward member. "
        "The inherited procedure warns about primer carryover in PCR; the 2016 successor explicitly describes empty-circle products near 120 nt. "
        "That successor size is not assigned to the 2015 branch. The retained ATC before the reverse-primer complement comes from the 2015 RT primer across circle closure."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS), acquired_roles=ACQUIRED_ROLES, read_lengths=READ_LENGTHS))
    graph.annotate(Highlight("library", "Source-selected run and computational filtering",
        "This view selects NextSeq 500 single-read 75 cycles; the same 2015 methods also report HiSeq 2000 single-read 50 cycles. "
        "After removing N6, the first aligned base reports the nascent RNA 3′ end. Mispriming, PCR duplicates and exact splice-junction intermediate positions are filtered computationally; those filters do not discard physical library molecules."))
    graph.declare_reviewed_scope((
        "nascent_rna","linker","linker_ligase","ligate_linker",
        "alkaline_hydrolysis","fragment_rna","rna_size_rule","select_rna","rt_primer",
        "reverse_transcriptase","reverse_transcribe","rna_removal","remove_rna","size_rule","size_select",
        "circligase","circularize","depletion_oligos","deplete","index_forward",
        "reverse_primer","pcr_polymerase","library_pcr","final_gel","select_library","library"),SOURCE)
    return graph.require_publishable()
