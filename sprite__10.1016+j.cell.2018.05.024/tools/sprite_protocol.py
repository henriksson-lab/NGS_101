from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row, Segment
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from interaction_capture import FixedMolecularComplex
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from sprite import CompositeIdentifier, final_library, ligation_round_rows, tag as _tag

def tag(name: str, length: int, part: int, *, group: str, role: str):
    return _tag(name, length, part, group=group, role=role,
                whitelist="published 96-well SPRITE tag set")


def bulk_barcode() -> tuple[list[Segment], CompositeIdentifier]:
    parts = [tag("DPM tag", 9, 1, group="sprite_cluster", role="inline_barcode"),
             tag("Odd tag 1", 17, 2, group="sprite_cluster", role="inline_barcode"),
             tag("Even tag", 17, 3, group="sprite_cluster", role="inline_barcode"),
             tag("Odd tag 2", 17, 4, group="sprite_cluster", role="inline_barcode"),
             tag("Terminal tag", 9, 5, group="sprite_cluster", role="inline_barcode")]
    return parts, CompositeIdentifier("sprite_cluster", "spatial complex", tuple(x.name for x in parts))


TITLE="SPRITE — split-pool tagging of spatial complexes"
NOTES="01_sprite.html"
SOURCE='Defining source: <a href="https://doi.org/10.1016/j.cell.2018.05.024">Quinodoz et al., <i>Cell</i> (2018)</a>.'
SUMMARY="Crosslinked multi-molecule complexes are immobilized and repeatedly split across 96 tag wells; molecules that remain in one complex acquire the same five-part barcode without proximity-ligating the genomic fragments to one another."
CAVEAT="Barcode-tag ligations are covalent, but genomic fragments in a cluster are not ligated to each other. Shared tag combinations encode co-membership of a crosslinked spatial complex."
TAGS,CLUSTER_ID=bulk_barcode()
FINAL_LIBRARY,SEQ_PRIMERS=final_library(TAGS,"bulk SPRITE library")
CLUSTER_ID.bind(FINAL_LIBRARY)
FINAL_CAPTION="Read 1 reports the DPM tag and genomic DNA; Read 2 reports the remaining split-pool tag chain."
SEQUENCING_INTRO="Paired-end sequencing reads genomic sequence from the DPM end and reconstructs the cluster identifier from DPM plus the four tag parts read from the opposite end."
READ_LENGTHS={"Read 1":100,"Read 2":75}
SHOW_READ_CONTENT=True

COMPLEX=FixedMolecularComplex(("DNA fragment A","DNA fragment B","RNA/protein-associated nucleic acid"),"DSG/formaldehyde",label="spatial complex")
INITIAL_ROWS=tuple(COMPLEX.complex().rows())
INITIAL_NAME="Native multi-molecule chromatin complex"

def sections(): return [
 ("Crosslink, isolate and fragment chromatin",COMPLEX.complex().rows(),"DSG and formaldehyde preserve separate molecules in one multiway complex; sonication and DNase make sequenceable fragments."),
 ("Immobilize intact complexes and repair DNA ends",COMPLEX.complex().rows(),"The drawn members remain separate molecules in one bead-bound complex. Low complex loading and detergent washes reduce accidental co-barcoding; repair creates a common ligation substrate."),
 ("Perform five split-pool tag ligations",ligation_round_rows(TAGS),"The published order is DPM, Odd, Even, Odd, Terminal. DPM/Terminal have 9-nt identifiers; Odd/Even tags have 17-nt identifiers."),
 ("Reverse crosslinks and PCR-complete the library",[Row(chunks=[("P5 -- DPM -- genomic DNA -- Terminal -- Odd -- Even -- Odd -- P7",None,False)])],"Terminal and DPM arms provide the two library-PCR entry sites; molecules are grouped after sequencing by the full tag combination."),
]


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
    reversal=graph.input("crosslink_reversal","crosslink-reversal and protein-digestion condition",
                         MaterialKind.CONDITION)
    released=graph.product("released_tagged_dna","released five-part-tagged DNA",
                           MaterialKind.MOLECULAR,
                           state=tuple(ligation_round_rows(TAGS)))
    graph.react("reverse_crosslinks","Reverse crosslinks and release tagged DNA",
                StepKind.DEPROTEINIZATION,
                ((current,InputRole.TEMPLATE),(reversal,InputRole.CONDITION)),(released,))
    fwd=graph.input("dpm_pcr_primer","DPM-arm library PCR primer",MaterialKind.OLIGO)
    rev=graph.input("terminal_pcr_primer","Terminal-arm library PCR primer",MaterialKind.OLIGO)
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
    scope=["chromatin_complex","crosslinkers","crosslink","fragmentation_condition",
           "fragment","nhs_beads","immobilize","end_repair_enzymes","end_repair",
           "dna_ligase"]
    for identifier in round_ids:
        scope.extend((identifier+"_adapter","ligate_"+identifier))
    scope.extend(("crosslink_reversal","reverse_crosslinks","dpm_pcr_primer",
                  "terminal_pcr_primer","pcr_polymerase","library_pcr","library"))
    graph.declare_reviewed_scope(scope,SOURCE)
    return graph.require_publishable()
