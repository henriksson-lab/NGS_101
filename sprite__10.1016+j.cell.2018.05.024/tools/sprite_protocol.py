from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct, Row, Segment, feature, revcomp
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from interaction_capture import FixedMolecularComplex
from protocol_graph import (AnnotationSlot, Disposition, Highlight, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, SourceClaim,
                            SourceDiscrepancy, StepKind, Terminal)
from sprite import CompositeIdentifier, ligation_round_rows

def tag(name: str, part: int):
    return Segment(name, "["+name+"]", "cbc", placeholder=True,
        feature=feature(name.lower().replace(" ", "_"), "inline_barcode", "combinatorial",
                        group="sprite_cluster", part=f"round {part}",
                        whitelist="original 2018 Table S5 tag families"))


def bulk_barcode() -> tuple[list[Segment], CompositeIdentifier]:
    parts=[tag(name,i) for i,name in enumerate(
        ("DPM tag", "Odd tag 1", "Even tag", "Odd tag 2", "Terminal tag"),1)]
    return parts, CompositeIdentifier("sprite_cluster", "spatial complex", tuple(x.name for x in parts))


# Original 2018 Methods, Final Library Amplification Primers. The source calls
# these Read1/Read2 primers, but they are library PCR oligos, not loaded primers.
PCR_FORWARD=il.TRUSEQ_P5_FULL
ARTICLE_PRIMER_INDEX="GCCTAGCC"
PCR_REVERSE=il.P7+ARTICLE_PRIMER_INDEX+il.TRUSEQ_READ2
DPM_QC_FORWARD="TACACGACGCTCTTCCGATCT"
DPM_QC_REVERSE="TGACTTGTCATGTCTTCCGATCT"
AUTHOR_GUIDE="Guttman lab DNA SPRITE guide, last edited August 9, 2018, sections 3.5 and 5.4 (pp. 9, 20–22)"
ORDERING_OLIGOS=(
 ("DPM Top Plate6-A1", "/5Phos/AAACACCCAAGATCGGAAGAGCGTCGTGTA/3SpC3/"),
 ("DPM Bottom DPM6bot1", "/5Phos/TGACTTGTCATGTCTTCCGATCTTGGGTGTTTT"),
 ("Odd Top Odd2Top1", "/5Phos/CAAGTCAAGCTAGATTCCACGAAG"),
 ("Odd Bottom Odd2Bo1", "/5Phos/GACAACTCTTCGTGGAATCTAGCT"),
 ("Even Top Even2Top1", "/5Phos/AGTTGTCACGTCAGCCGCAGTATC"),
 ("Even Bottom Even2Bo1", "/5Phos/TGACTTGGATACTGCGGCTGACGT"),
 ("Terminal Top NYTop1_Stg", "/5Phos/AGTTGTCACCATAATAAGATCGGAAGA"),
 ("Terminal Bottom NYBot1_Stg", "CAGACGTGTGCTCTTCCGATCTTATTATGGT"),
)


def oligos():
    return ["Representative first-well ordering oligos from original Table S5; these do not assign one well to every library.",
            *[name+": "+seq for name,seq in ORDERING_OLIGOS],
            "Library PCR forward: "+PCR_FORWARD,
            "Library PCR reverse (article example; guide name 2P_barcoded_85): "+PCR_REVERSE,
            "DPMQCprimerF (historical author-guide QC branch): "+DPM_QC_FORWARD,
            "DPMQCprimerR (historical author-guide QC branch): "+DPM_QC_REVERSE]


TITLE="SPRITE — split-pool tagging of spatial complexes"
NOTES="01_sprite.html"
SOURCE='Defining source: <a href="https://doi.org/10.1016/j.cell.2018.05.024">Quinodoz et al., <i>Cell</i> (2018)</a>.'
SUMMARY="Crosslinked multi-molecule complexes are immobilized and repeatedly split across 96 tag wells; molecules that remain in one complex acquire the same five-part barcode without proximity-ligating the genomic fragments to one another."
CAVEAT="This page models the DNA branch. Barcode-tag ligations are covalent; genomic fragments in a cluster remain separate. The displayed PCR arms follow the original 2018 printed primers, with one i7 and no i5. Internal tag roles and unassembled junctions are bracketed rather than assigned physical barcode widths."
TAGS,CLUSTER_ID=bulk_barcode()
# Methods final-DNA diagram establishes DPM at both genomic ends. Retain
# this distinction; source tag roles are not stand-ins for all intervening bases.
FINAL_LIBRARY=Construct([
    Segment("PCR forward arm",PCR_FORWARD,"p5"), TAGS[0],
    Segment("DPM/genome junction", "[junction]",placeholder=True),
    Segment("genomic DNA", "[genomic DNA]",placeholder=True),
    Segment("opposite DPM and junction", "[opposite DPM/junction]",placeholder=True),
    TAGS[1],Segment("Odd/Even junction","[junction]",placeholder=True),
    TAGS[2],Segment("Even/Odd junction","[junction]",placeholder=True),
    TAGS[3],Segment("Odd/Terminal junction","[junction]",placeholder=True),
    TAGS[4],Segment("Terminal/PCR junction","[junction]",placeholder=True),
    Segment("Read2 arm reverse complement",revcomp(il.TRUSEQ_READ2),"r2"),
    Segment("i7 reverse complement",revcomp(ARTICLE_PRIMER_INDEX),"cbc",
            feature=feature("sample_index_i7","sample_index","fixed")),
    Segment("P7 reverse complement",il.P7_RC,"p7")],name="SPRITE DNA endpoint: known PCR arms, logical tag roles")
SEQ_PRIMERS=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["R2"])
CLUSTER_ID.bind(FINAL_LIBRARY)
FINAL_CAPTION="Known PCR arms flank a schematic DNA/tag chain. DPM is present at both genomic ends; five logical tag roles define the cluster identifier. Bracketed junctions and tag roles do not assign physical lengths or exact barcode-cycle coordinates."
SEQUENCING_INTRO="The original article specifies paired-end reads of at least 115 bases in Read 1 and 100 in Read 2 on HiSeq 2500 or NextSeq 500. Read 1 supplies DPM and genomic information; Read 2 supplies the remaining tag chain. The model does not substitute these minima for an exact acquired cycle count."
READ_MINIMUMS={"Read 1":115,"Read 2":100}
READ_LENGTHS={}
READ_SPECS=(sp.read_spec("Read 1","at least 115","DPM at beginning, then genomic information; exact modeled spans unresolved","extract_then_include"),
            sp.read_spec("Read 2","at least 100","remaining split-pool tags; physical junction spans not reconstructed","extract"),
            sp.read_spec("Index 1 (i7)","unassigned run; "+str(len(ARTICLE_PRIMER_INDEX))+"-base index","sample index from printed PCR primer example","demultiplex"))
SHOW_READ_CONTENT=False
SEQUENCING_UNAVAILABLE="The full physical tag and junction chain is not reconstructed here. Known PCR arms establish standard sequencing-primer sites, but bracketed interior roles cannot establish exact barcode-cycle coordinates. This is a model limitation, not a claim of absent source sequences."

COMPLEX=FixedMolecularComplex(("DNA fragment A","DNA fragment B","RNA/protein-associated nucleic acid"),"DSG/formaldehyde",label="spatial complex")
INITIAL_ROWS=tuple(COMPLEX.complex().rows())
INITIAL_NAME="Native multi-molecule chromatin complex"

def sections(): return [
 ("Crosslink, isolate and fragment chromatin",COMPLEX.complex().rows(),"DSG and formaldehyde preserve separate molecules in one multiway complex; sonication and DNase make sequenceable fragments."),
 ("Immobilize intact complexes and repair DNA ends",COMPLEX.complex().rows(),"The drawn members remain separate molecules in one bead-bound complex. Low complex loading and detergent washes reduce accidental co-barcoding; repair creates a common ligation substrate."),
 ("Perform five split-pool tag ligations",ligation_round_rows(TAGS),"The published order is DPM, Odd, Even, Odd, Terminal. Original Table S5 supplies the actual ordering oligos; barcode-column lengths and invariant flanks must be distinguished from complete tag lengths."),
 ("Reverse crosslinks and PCR-complete the library",[Row(chunks=[("PCR arm → DPM → genomic DNA → opposite DPM → Odd → Even → Odd → Terminal → PCR arm",None,False)])],"Terminal and DPM arms provide the two library-PCR entry sites; molecules are grouped after sequencing by the full tag combination."),
]


def _dpm_qc_branch(graph, tagged):
    """The author guide tests a separate aliquot before the four remaining tags."""
    remaining=graph.product("dpm_remaining", "DPM-tagged complexes continuing to split-pool tagging",
                            MaterialKind.MOLECULAR, state=tagged.state)
    aliquot=graph.product("dpm_qc_aliquot", "author-guide DPM QC aliquot", MaterialKind.MOLECULAR,
                          state=tagged.state)
    graph.react("split_dpm_qc", "Author guide: separate DPM-ligation QC aliquot", StepKind.SPLIT,
                ((tagged,InputRole.TEMPLATE),), (remaining,aliquot),
                source=AUTHOR_GUIDE, changes_molecule=False)
    release=graph.input("dpm_qc_release", "crosslink reversal and Proteinase K for QC aliquot",
                        MaterialKind.CONDITION, source=AUTHOR_GUIDE)
    template=Construct([Segment("DPM-tagged DNA", "[DPM-tagged DNA]",placeholder=True)],
                       name="unassembled DPM QC template")
    released=graph.product("dpm_qc_released", "released DPM-tagged DNA for QC",
                           MaterialKind.MOLECULAR, state=template)
    graph.react("release_dpm_qc", "Release DPM-tagged DNA from the QC aliquot", StepKind.DEPROTEINIZATION,
                ((aliquot,InputRole.TEMPLATE),(release,InputRole.CONDITION)), (released,),source=AUTHOR_GUIDE)
    primers=[]
    for identifier,sequence in (("dpm_qc_forward",DPM_QC_FORWARD),("dpm_qc_reverse",DPM_QC_REVERSE)):
        primers.append(graph.input(identifier,identifier.replace("_"," "),MaterialKind.OLIGO,
                                   state=Segment(identifier,sequence),source=AUTHOR_GUIDE))
    polymerase=graph.input("dpm_qc_polymerase", "Q5 Hot Start Master Mix for QC", MaterialKind.CATALYST,
                           source=AUTHOR_GUIDE)
    product=graph.product("dpm_qc_amplicons", "DPM QC amplicons for concentration/size assessment",
                          MaterialKind.QC)
    graph.react("dpm_qc_pcr", "Author guide: amplify DPM-ligation QC aliquot", StepKind.PCR,
                ((released,InputRole.TEMPLATE),(primers[0],InputRole.STRUCTURAL),
                 (primers[1],InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(product,),
                source=AUTHOR_GUIDE, geometry=PCRGeometry.unavailable(template,
                    reason="The guide prints both QC primers; this model has not reconstructed their precursor duplex and amplicon boundaries."))
    graph.mark_terminal(Terminal(product.id,Disposition.QC,
        "The author guide assesses DPM-ligation yield and size; this aliquot is not the sequenced SPRITE library."))
    targets=("split_dpm_qc","dpm_qc_aliquot","dpm_qc_release","release_dpm_qc",
             "dpm_qc_forward","dpm_qc_reverse","dpm_qc_polymerase","dpm_qc_pcr","dpm_qc_amplicons")
    graph.declare_reviewed_scope(targets,AUTHOR_GUIDE)
    return remaining


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    chromatin=graph.input("chromatin_complex","native multi-molecule chromatin complex",
                          MaterialKind.BIOLOGICAL,state=INITIAL_ROWS,primary=True)
    crosslinker=graph.input("crosslinkers","DSG and formaldehyde",MaterialKind.CONDITION)
    fixed=graph.product("fixed_complex","crosslinked spatial complex",MaterialKind.MOLECULAR,state=COMPLEX.complex())
    graph.react("crosslink","Preserve multiway chromatin complexes",StepKind.CROSS_LINKING,
                ((chromatin,InputRole.TEMPLATE),(crosslinker,InputRole.CONDITION)),(fixed,))
    fragmentation=graph.input("fragmentation_condition","sonication and DNase fragmentation",
                              MaterialKind.CONDITION)
    fragmented=graph.product("fragmented_complex","fragmented but intact spatial complex",
                              MaterialKind.MOLECULAR,state=COMPLEX.complex())
    graph.react("fragment","Fragment nucleic acids within each fixed complex",
                StepKind.FRAGMENTATION,
                ((fixed,InputRole.TEMPLATE),(fragmentation,InputRole.CONDITION)),(fragmented,))
    nhs=graph.input("nhs_beads","NHS-functionalized beads",MaterialKind.SELECTION)
    immobilized=graph.product("immobilized_complex","bead-immobilized spatial complex",
                              MaterialKind.MOLECULAR,state=COMPLEX.complex())
    unbound=graph.product("unbound_complex","unbound spatial complexes",MaterialKind.WASTE)
    graph.react("immobilize","Covalently immobilize intact complexes on NHS beads",
                StepKind.AFFINITY_CAPTURE,
                ((fragmented,InputRole.TEMPLATE),(nhs,InputRole.SELECTION)),
                (immobilized,unbound))
    graph.discard(unbound,"complexes not immobilized on the sparsely loaded beads are washed away")
    end_enzyme=graph.input("end_repair_enzymes","DNA end-repair and dA-tailing enzymes",
                           MaterialKind.CATALYST)
    repaired=graph.product("repaired_complex","immobilized complex with ligatable DNA ends",
                           MaterialKind.MOLECULAR)
    graph.react("end_repair","Repair, phosphorylate and dA-tail DNA ends",StepKind.END_REPAIR,
                ((immobilized,InputRole.TEMPLATE),(end_enzyme,InputRole.CATALYST)),(repaired,))
    ligase=graph.input("dna_ligase","DNA ligase",MaterialKind.CATALYST)
    current=repaired
    round_ids=("dpm","odd_1","even","odd_2","terminal")
    round_labels=("DPM","Odd 1","Even","Odd 2","Terminal")
    for number,(identifier,label,tag) in enumerate(zip(round_ids,round_labels,TAGS),1):
        adapter=graph.input(identifier+"_adapter",label+" split-pool tag adapter",
                            MaterialKind.ADAPTER,state=tag)
        product=graph.product(identifier+"_tagged",label+"-tagged bead complex",
                              MaterialKind.MOLECULAR,
                              state=tuple(ligation_round_rows(TAGS[:number])))
        graph.react("ligate_"+identifier,"Ligate the "+label+" split-pool tag",
                    StepKind.LIGATION,
                    ((current,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
                     (ligase,InputRole.CATALYST)),(product,),
                    geometry=AdapterLigationGeometry.unavailable(
                        "The source establishes the ordered covalent tag ligation but the compact module does not encode the junction-end geometry."))
        current=product
        if identifier=="dpm":
            current=_dpm_qc_branch(graph,current)
    reversal=graph.input("crosslink_reversal","crosslink-reversal and protein-digestion condition",
                         MaterialKind.CONDITION)
    released=graph.product("released_tagged_dna","released five-part-tagged DNA",
                           MaterialKind.MOLECULAR,
                           state=tuple(ligation_round_rows(TAGS)))
    graph.react("reverse_crosslinks","Reverse crosslinks and release tagged DNA",
                StepKind.DEPROTEINIZATION,
                ((current,InputRole.TEMPLATE),(reversal,InputRole.CONDITION)),(released,))
    fwd=graph.input("dpm_pcr_primer","DPM-arm library PCR primer",MaterialKind.OLIGO, state=Segment("PCR forward",PCR_FORWARD))
    rev=graph.input("terminal_pcr_primer","Terminal-arm library PCR primer",MaterialKind.OLIGO, state=Segment("PCR reverse",PCR_REVERSE))
    pol=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    graph.react("library_pcr","PCR-complete the SPRITE library",StepKind.PCR,
                ((released,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),
                 (rev,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The compact model does not reconstruct both PCR-primer binding states."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS),
                                          read_lengths=READ_LENGTHS))
    graph.annotate(SourceDiscrepancy(
        target="dpm_adapter", title="DPM sequencing-read labels disagree",
        claims=(
            SourceClaim("2018 Methods, DNA Phosphate Modified tag; JATS lines 1123–1126",
                        "Says the barcode is read only in Read 2."),
            SourceClaim("2018 Methods, Barcode identification and Alignment to genome; JATS lines 1222–1250",
                        "Parses DPM at the beginning of Read 1 and trims 11 bases from Read 1 before alignment.")),
        model_choice="The displayed read roles follow the explicit Read 1 parsing and alignment instructions.",
        evidence="Those processing sections agree with the preceding library-read description and accepted manuscript.",
        uncertainty="The earlier Read 2 label is an apparent source error; no correction is established here.",
        slot=AnnotationSlot.INTRO))
    graph.annotate(Highlight(
        target="terminal_pcr_primer", title="The displayed index is the article's primer example",
        explanation="The article's PCR reverse primer contains GCCTAGCC. Table S5 primer 2P_657 contains GGCTAGGC, with GCCTAGCC in its Barcode column. The article does not name its example as 2P_657. These reverse-complement strings are retained as distinct source examples; they do not establish an error or the index used by every library.",
        slot=AnnotationSlot.INTRO))
    graph.annotate(SourceDiscrepancy(
        target="library_pcr",title="Historical author-guide primer lengths and PCR cycles",
        claims=(
            SourceClaim("April 5 and August 9, 2018 DNA guides, section 3.4 (pp. 10 and 9)",
                        "Describe both full primers as around 30 nucleotides and five further cycles after the initial four."),
            SourceClaim("The same guides, printed 2P_universal/2P_barcoded_85 strands in section 3.4; PCR program in section 6.2",
                        f"The literal full strands are {len(PCR_FORWARD)} and {len(PCR_REVERSE)} bases; the program specifies seven further cycles after the initial four.")),
        model_choice="Keep the literal full primer sequences; leave the final PCR cycle count unassigned.",
        evidence="The named full strands agree with the original article's PCR examples. Their lengths are computed from those strands.",
        uncertainty="The shorter length wording may refer to a primer subregion, but the guide does not identify one. No source correction is established.",
        slot=AnnotationSlot.INTRO))
    graph.annotate(Highlight(
        target="end_repair",title="Historical guide overhang wording remains qualified",
        explanation="Both 2018 DNA guides, section 3.1, call the target end a 5-prime A overhang, while their diagram places A at the genomic 3-prime ends. This page records the stated dA-tailing operation without resolving the physical DPM junction or silently adopting the conflicting polarity label.",
        slot=AnnotationSlot.INTRO))
    scope=["chromatin_complex","crosslinkers","crosslink","fragmentation_condition",
           "fragment","nhs_beads","immobilize","end_repair_enzymes","end_repair",
           "dna_ligase"]
    for identifier in round_ids:
        scope.extend((identifier+"_adapter","ligate_"+identifier))
    scope.extend(("crosslink_reversal","reverse_crosslinks","dpm_pcr_primer",
                  "terminal_pcr_primer","pcr_polymerase","library_pcr","library"))
    graph.declare_reviewed_scope(scope,SOURCE)
    return graph.require_publishable()
