from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row
from sprite import final_library, ligation_round_rows, single_cell_barcode
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)
from interaction_capture import FixedMolecularComplex

TITLE="scSPRITE — cell and spatial-complex barcoding"
NOTES="01_scsprite.html"
SOURCE='Defining preprint: <a href="https://doi.org/10.1101/2020.08.11.242081">Arrastia et al. (2020)</a>; detailed published method: <a href="https://doi.org/10.1038/s41587-021-00998-1"><i>Nature Biotechnology</i> (2022)</a>.'
SUMMARY="Three split-pool ligations barcode intact nuclei; after selecting individual nuclei and fragmenting chromatin, three more rounds extend the same molecules with a spatial-cluster barcode."
CAVEAT="The first three tag parts identify the cell and participate again in the full six-part spatial-cluster identifier. The data model records both overlapping composite meanings explicitly."
TAGS,IDENTIFIERS=single_cell_barcode()
FINAL_LIBRARY,SEQ_PRIMERS=final_library(TAGS,"scSPRITE library")
for x in IDENTIFIERS: x.bind(FINAL_LIBRARY)
FINAL_CAPTION="Read 1 reports genomic DNA plus DPM; Read 2 traverses the other five tags. The first three tags identify the nucleus and all six identify a spatial cluster."
SEQUENCING_INTRO="Read 1 is at least 120 bp and begins with genomic information plus DPM; Read 2 is at least 95 bp and reads Odd–Even–Odd–Even–Y-even."
READ_LENGTHS={"Read 1":120,"Read 2":95}
SHOW_READ_CONTENT=True
COMPLEX=FixedMolecularComplex(("cell-barcoded DNA fragment A","cell-barcoded DNA fragment B","associated RNA/protein"),"formaldehyde",label="spatial complex")

def sections(): return [
 ("Crosslink, permeabilize, repair and dA-tail nuclei",[Row(chunks=[("intact nucleus: repaired genomic fragments remain compartmentalized",None,False)])],"Keeping nuclei intact makes the first barcode rounds cell-specific."),
 ("Add three in-nucleus split-pool tags",ligation_round_rows(TAGS[:3]),"DPM, Odd and Even across three 96-well rounds form the cell barcode."),
 ("Select nuclei, lyse and immobilize spatial complexes",COMPLEX.complex().rows(),"The members drawn separately share the first three-part cell barcode and remain together as one bead-bound complex. About 1,500 filtered nuclei were taken forward."),
 ("Extend each complex with three spatial tags",ligation_round_rows(TAGS[3:]),"Odd, Even and Y-even tags distinguish complexes while retaining the first three-part cell identity."),
 ("Reverse crosslinks and PCR-complete the library",[Row(chunks=[("P5 -- DPM -- genomic DNA -- Y-even -- Even -- Odd -- Even -- Odd -- P7",None,False)])],"Only reads with the complete six-tag order are usable for cell and cluster assignment."),
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
                 "The shared SPRITE tag model preserves tag order but does not expose one combined ligating-end scene for the three rounds."))
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
                 "The shared SPRITE tag model preserves the three-round order but does not expose one combined ligating-end scene."))
 reversal=graph.input("crosslink_reversal","cross-link reversal condition",MaterialKind.CONDITION)
 purified=graph.product("purified_tagged_dna","de-cross-linked six-tag DNA",MaterialKind.MOLECULAR)
 graph.react("reverse_crosslinks","Reverse cross-links and recover tagged DNA",
             StepKind.DEPROTEINIZATION,
             ((six_tagged,InputRole.TEMPLATE),(reversal,InputRole.CONDITION)),(purified,))
 primers=graph.input("pcr_primers","scSPRITE library PCR-primer pair",MaterialKind.OLIGO)
 polymerase=graph.input("polymerase","library PCR polymerase",MaterialKind.CATALYST)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 graph.react("library_pcr","PCR-complete the six-tag library",StepKind.PCR,
             ((purified,InputRole.TEMPLATE),(primers,InputRole.STRUCTURAL),
              (polymerase,InputRole.CATALYST)),(library,),
             geometry=PCRGeometry.unavailable(
                 FINAL_LIBRARY,reason="The shared SPRITE model supplies the completed library rather than separate PCR landing objects.",symmetric_species=True))
 graph.sequence(library,SequencingPlan(
     library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
     acquired_roles=tuple(p.role for p in SEQ_PRIMERS),read_lengths=READ_LENGTHS,
     read_uses={"Read 1":"include","Read 2":"extract_then_include"}))
 graph.declare_reviewed_scope((
     "nuclei","crosslinker","crosslink","repair_enzymes","end_repair",
     "tailing_polymerase","da_tail","cell_tags","ligase","cell_tag_ligation",
     "nucleus_selection","select_nuclei","lysis_condition","lysis","sonication",
     "sonicate","beads","immobilize","spatial_tags","spatial_tag_ligation",
     "crosslink_reversal","reverse_crosslinks","pcr_primers","polymerase",
     "library_pcr","library"),source)
 return graph.require_publishable()
