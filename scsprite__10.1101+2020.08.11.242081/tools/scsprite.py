from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct, Row, Segment, feature
import seqprimers as sp
from sprite import CompositeIdentifier, ligation_round_rows
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)
from interaction_capture import FixedMolecularComplex

def tag(name: str, part: int, *, group: str, role: str):
    # The biological identifier role/order is known; physical tag lengths and
    # junctions are not assembled from the scSPRITE ordering oligos here.
    return Segment(name, "[" + name + "]", "cbc", placeholder=True,
        feature=feature(name.lower().replace(" ", "_"), role, "combinatorial",
                        group=group, part=f"round {part}",
                        whitelist="source-linked 96-well SPRITE/scSPRITE tag families"))


def single_cell_barcode() -> tuple[list[Segment], tuple[CompositeIdentifier, ...]]:
    names = ("DPM tag", "Odd cell tag", "Even cell tag", "Odd spatial tag",
             "Even spatial tag", "Y-even terminal tag")
    parts=[]
    for i,name in enumerate(names,1):
        role="cell_barcode" if i <= 3 else "inline_barcode"
        group="cell_id" if i <= 3 else "spatial_extension"
        parts.append(tag(name,i,group=group,role=role))
    ids=(CompositeIdentifier("cell_id","cell",tuple(x.name for x in parts[:3])),
         CompositeIdentifier("spatial_cluster","spatial complex",tuple(x.name for x in parts)))
    return parts,ids


TITLE="scSPRITE — cell and spatial-complex barcoding"
NOTES="01_scsprite.html"
SOURCE='Defining preprint: <a href="https://doi.org/10.1101/2020.08.11.242081">Arrastia et al. (2020)</a>; detailed published method: <a href="https://doi.org/10.1038/s41587-021-00998-1"><i>Nature Biotechnology</i> (2022)</a>.'
SUMMARY="Three split-pool ligations barcode intact nuclei; after selecting individual nuclei and fragmenting chromatin, three more rounds extend the same molecules with a spatial-cluster barcode."
CAVEAT="The first three tag parts identify the cell and participate again in the full six-part spatial-cluster identifier. The original SPRITE Supplementary Table S5 prints the DPM, Odd, Even and Terminal oligo families; the scSPRITE authors' read parser lists all 96 Y barcode strings, which match the older Terminal Bottom barcodes. The complete scSPRITE Y-even ligation oligo and PCR primer pair are not assembled here. The displayed endpoint is a logical tag-order model with unknown physical lengths, adapters and sequencing entry; no TruSeq shell or exact barcode-cycle coordinates are inferred. Source-scope discrepancy: the author pipeline's embedded tag_layout.png shows a five-part DPM–Odd–Even–Odd–Y chain, whereas the published Methods specify six scSPRITE parts. This model follows the scSPRITE Methods; the embedded image's experimental scope remains unassigned."
TAGS,IDENTIFIERS=single_cell_barcode()
# This is a logical identifier model, not a reconstructed physical library.
# Do not acquire a TruSeq shell or exact cycle spans from a family helper.
FINAL_LIBRARY=Construct([Segment("genomic payload", "[genomic DNA]", placeholder=True),
                         *TAGS], name="scSPRITE payload and logical identifier order")
SEQ_PRIMERS=()
PARSER_PREFIX="TCATGTCTTCCGATCT"
for x in IDENTIFIERS: x.bind(FINAL_LIBRARY)
FINAL_CAPTION="Logical identifier order DPM–Odd–Even–Odd–Even–Y-even: first three parts identify the nucleus and all six identify a spatial cluster. The panel does not specify physical adapter placement, tag lengths or sequencing-primer sites."
SEQUENCING_INTRO="The published article specifies at least 120 bases in Read 1 and 95 in Read 2. Read 1 supplies DPM and genomic information; Read 2 supplies the remaining tag information. Biological identifier order is distinct from parser matching order."
SEQUENCING_ENDING=SEQUENCING_INTRO
SEQUENCING_UNAVAILABLE="The exact scSPRITE ordering oligos, completed amplicon and loaded sequencing-primer sites have not been assembled in this model. The author parser READ2 pattern is Y–Even–Odd–Even–Odd–DPM; this parsing declaration alone does not establish physical per-cycle spans. No source-absence conclusion is asserted."
READ_SPECS=(sp.read_spec("Read 1","at least 120","DPM at the beginning and genomic information; exact modeled spans unresolved","include"),
            sp.read_spec("Read 2","at least 95","remaining barcode information; exact physical cycle spans unresolved","extract"))
READ_MINIMUMS={"Read 1":120,"Read 2":95}
READ_LENGTHS={}
SHOW_READ_CONTENT=False
COMPLEX=FixedMolecularComplex(("cell-barcoded DNA fragment A","cell-barcoded DNA fragment B","associated RNA/protein"),"formaldehyde",label="spatial complex")

def sections(): return [
 ("Crosslink, permeabilize, repair and dA-tail nuclei",[Row(chunks=[("intact nucleus: repaired genomic fragments remain compartmentalized",None,False)])],"Keeping nuclei intact makes the first barcode rounds cell-specific."),
 ("Add three in-nucleus split-pool tags",ligation_round_rows(TAGS[:3]),"DPM, Odd and Even across three 96-well rounds form the cell barcode."),
 ("Select nuclei, lyse and immobilize spatial complexes",COMPLEX.complex().rows(),"The members drawn separately share the first three-part cell barcode and remain together as one bead-bound complex. About 1,500 filtered nuclei were taken forward."),
 ("Extend each complex with three spatial tags",ligation_round_rows(TAGS[3:]),"Odd, Even and Y-even tags distinguish complexes while retaining the first three-part cell identity."),
 ("Reverse crosslinks and PCR-complete the library",[Row(chunks=[("genomic payload + identifier: DPM → Odd → Even → Odd → Even → Y-even",None,False)])],"Only reads with the complete six-tag order are usable for cell and cluster assignment."),
]

def protocol_graph():
 source="Defining preprint: https://doi.org/10.1101/2020.08.11.242081; detailed method: https://doi.org/10.1038/s41587-021-00998-1"
 graph=ProtocolGraph("scSPRITE",authority=source)
 nuclei=graph.input("nuclei","intact nuclei with compartmentalized chromatin",
                    MaterialKind.BIOLOGICAL,primary=True)
 crosslinker=graph.input("crosslinker","chromatin cross-linking condition",MaterialKind.CONDITION)
 crosslinked=graph.product("crosslinked_nuclei","cross-linked intact nuclei",MaterialKind.BIOLOGICAL)
 graph.react("crosslink","Cross-link and permeabilize nuclei",StepKind.CROSS_LINKING,
             ((nuclei,InputRole.TEMPLATE),(crosslinker,InputRole.CONDITION)),(crosslinked,))
 repair=graph.input("repair_enzymes","DNA end-repair enzymes",MaterialKind.CATALYST)
 repaired=graph.product("repaired_nuclei","nuclei containing repaired genomic fragments",
                        MaterialKind.BIOLOGICAL)
 graph.react("end_repair","Repair genomic fragment ends in intact nuclei",StepKind.END_REPAIR,
             ((crosslinked,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(repaired,))
 tailing=graph.input("tailing_polymerase","dA-tailing polymerase",MaterialKind.CATALYST)
 tailed=graph.product("tailed_nuclei","nuclei containing dA-tailed genomic fragments",
                      MaterialKind.BIOLOGICAL)
 graph.react("da_tail","Add 3′ dA overhangs in nuclei",StepKind.TAILING,
             ((repaired,InputRole.TEMPLATE),(tailing,InputRole.CATALYST)),(tailed,))
 cell_tags=graph.input("cell_tags","first three split-pool cell-tag adapters",MaterialKind.ADAPTER)
 ligase=graph.input("ligase","split-pool DNA ligase",MaterialKind.CATALYST)
 cell_barcoded=graph.product("cell_barcoded_nuclei","nuclei carrying three-part cell barcodes",
                             MaterialKind.BIOLOGICAL)
 graph.react("cell_tag_ligation","Add three in-nucleus split-pool tags",StepKind.LIGATION,
             ((tailed,InputRole.TEMPLATE),(cell_tags,InputRole.STRUCTURAL),
              (ligase,InputRole.CATALYST)),(cell_barcoded,),
             geometry=AdapterLigationGeometry.unavailable(
                 "The source establishes three in-nucleus tag rounds; their combined ligation-ready ends are not assembled in this model."))
 selection=graph.input("nucleus_selection","single-nucleus selection rule",MaterialKind.SELECTION)
 selected=graph.product("selected_nuclei","selected cell-barcoded nuclei",MaterialKind.BIOLOGICAL)
 excluded=graph.product("excluded_nuclei","nuclei excluded from the spatial-complex branch",
                        MaterialKind.WASTE)
 graph.react("select_nuclei","Select nuclei for spatial-complex indexing",StepKind.SAMPLING,
             ((cell_barcoded,InputRole.TEMPLATE),(selection,InputRole.SELECTION)),
             (selected,excluded),changes_molecule=False)
 graph.discard(excluded,"nuclei outside the selected set are not carried into the library")
 lysis=graph.input("lysis_condition","nuclear lysis condition",MaterialKind.CONDITION)
 lysate=graph.product("chromatin_lysate","released cross-linked chromatin complexes",
                      MaterialKind.MOLECULAR)
 graph.react("lysis","Lyse selected nuclei",StepKind.EXTRACTION,
             ((selected,InputRole.TEMPLATE),(lysis,InputRole.CONDITION)),(lysate,))
 sonication=graph.input("sonication","sonication condition",MaterialKind.CONDITION)
 complexes=graph.product("spatial_complexes","sonicated cross-linked chromatin complexes",
                         MaterialKind.MOLECULAR,state=COMPLEX.complex())
 graph.react("sonicate","Sonicate the cross-linked chromatin",StepKind.FRAGMENTATION,
             ((lysate,InputRole.TEMPLATE),(sonication,InputRole.CONDITION)),(complexes,))
 beads=graph.input("beads","spatial-complex immobilization beads",MaterialKind.SELECTION)
 immobilized=graph.product("immobilized_complexes","bead-bound spatial complexes",
                           MaterialKind.MOLECULAR,state=COMPLEX.complex())
 uncaptured=graph.product("uncaptured_complexes","uncaptured chromatin complexes",
                          MaterialKind.WASTE)
 graph.react("immobilize","Immobilize spatial complexes on beads",StepKind.AFFINITY_CAPTURE,
             ((complexes,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),
             (immobilized,uncaptured))
 graph.discard(uncaptured,"uncaptured complexes are removed during bead washing")
 spatial_tags=graph.input("spatial_tags","three spatial-cluster tag adapters",MaterialKind.ADAPTER)
 six_tagged=graph.product("six_tagged_complexes","six-part cell and spatial tagged complexes",
                          MaterialKind.MOLECULAR)
 graph.react("spatial_tag_ligation","Extend each complex with three spatial tags",
             StepKind.LIGATION,
             ((immobilized,InputRole.TEMPLATE),(spatial_tags,InputRole.STRUCTURAL),
              (ligase,InputRole.CATALYST)),(six_tagged,),
             geometry=AdapterLigationGeometry.unavailable(
                 "The source establishes three spatial tag rounds; their combined ligation-ready ends are not assembled in this model."))
 reversal=graph.input("crosslink_reversal","cross-link reversal condition",MaterialKind.CONDITION)
 purified=graph.product("purified_tagged_dna","de-cross-linked six-tag DNA",MaterialKind.MOLECULAR)
 graph.react("reverse_crosslinks","Reverse cross-links and recover tagged DNA",
             StepKind.DEPROTEINIZATION,
             ((six_tagged,InputRole.TEMPLATE),(reversal,InputRole.CONDITION)),(purified,))
 forward=graph.input("pcr_forward","unassigned indexed Illumina forward primer",MaterialKind.OLIGO)
 reverse=graph.input("pcr_reverse","unassigned indexed Illumina reverse primer",MaterialKind.OLIGO)
 polymerase=graph.input("polymerase","library PCR polymerase",MaterialKind.CATALYST)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 graph.react("library_pcr","PCR-complete the six-tag library",StepKind.PCR,
             ((purified,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),(reverse,InputRole.STRUCTURAL),
              (polymerase,InputRole.CATALYST)),(library,),
             geometry=PCRGeometry.unavailable(
                 FINAL_LIBRARY,reason="The original article names indexed Illumina PCR primers; a selected physical primer pair and its precursor binding sites are not assembled in this model."))
 graph.sequence(library,SequencingPlan(
     library=FINAL_LIBRARY,primers=(),
     acquired_roles=("Read 1","Read 2"),read_lengths=READ_LENGTHS,limitation=SEQUENCING_UNAVAILABLE,
     read_uses={"Read 1":"include","Read 2":"extract_then_include"}))
 graph.declare_reviewed_scope((
     "nuclei","crosslinker","crosslink","repair_enzymes","end_repair",
     "tailing_polymerase","da_tail","cell_tags","ligase","cell_tag_ligation",
     "nucleus_selection","select_nuclei","lysis_condition","lysis","sonication",
     "sonicate","beads","immobilize","spatial_tags","spatial_tag_ligation",
     "crosslink_reversal","reverse_crosslinks","pcr_forward","pcr_reverse","polymerase",
     "library_pcr","library"),source)
 return graph.require_publishable()
