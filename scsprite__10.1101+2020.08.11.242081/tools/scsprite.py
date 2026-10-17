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
SOURCE='Defining source (published version): Arrastia et al., <a href="https://doi.org/10.1038/s41587-021-00998-1"><i>Nature Biotechnology</i> (2022)</a>, Methods; the catalogue directory is named for the superseded <a href="https://doi.org/10.1101/2020.08.11.242081">bioRxiv preprint (2020)</a>.'
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
 ("Crosslink, isolate nuclei, HpyCH4V-digest and dA-tail",[Row(chunks=[("intact nucleus: blunt HpyCH4V fragments, dA-tailed, remain compartmentalized",None,False)])],"Keeping nuclei intact makes the first barcode rounds cell-specific."),
 ("Add three in-nucleus split-pool tags",ligation_round_rows(TAGS[:3]),"DPM, Odd and Even across three 96-well rounds form the cell barcode."),
 ("Select nuclei, lyse and immobilize spatial complexes",COMPLEX.complex().rows(),"The members drawn separately share the first three-part cell barcode and remain together as one bead-bound complex. About 1,500 filtered nuclei were taken forward."),
 ("Extend each complex with three spatial tags",ligation_round_rows(TAGS[3:]),"Odd, Even and Y-even tags distinguish complexes while retaining the first three-part cell identity."),
 ("Reverse crosslinks and PCR-complete the library",[Row(chunks=[("genomic payload + identifier: DPM → Odd → Even → Odd → Even → Y-even",None,False)])],"Only reads with the complete six-tag order are usable for cell and cluster assignment."),
]

PUBLISHED="Arrastia et al., Nature Biotechnology 40:64–73 (2022), doi:10.1038/s41587-021-00998-1 (PMC11588347 author manuscript), Methods"

def protocol_graph():
 source=PUBLISHED
 graph=ProtocolGraph("scSPRITE",authority=source)
 cells=graph.input("cells","trypsinized mESC single-cell suspension",
                   MaterialKind.BIOLOGICAL,primary=True)
 dsg=graph.input("dsg","2 mM DSG in PBS, 45 min RT, then glycine quench",MaterialKind.CONDITION)
 formaldehyde=graph.input("formaldehyde","1% formaldehyde in PBS, 10 min RT, then glycine quench",
                          MaterialKind.CONDITION)
 crosslinked=graph.product("crosslinked_cells","DSG/formaldehyde-crosslinked cells",MaterialKind.BIOLOGICAL)
 graph.react("crosslink","Dual-crosslink cells with DSG then formaldehyde",StepKind.CROSS_LINKING,
             ((cells,InputRole.TEMPLATE),(dsg,InputRole.CONDITION),(formaldehyde,InputRole.CONDITION)),
             (crosslinked,))
 nuclear=graph.input("nuclei_isolation","lysis buffers #1 and #2, 1.2× CutSmart, 0.3% SDS 37 °C 60 min, "
                     "Triton X-100 quench 37 °C 60 min",MaterialKind.CONDITION)
 nuclei=graph.product("nuclei","isolated permeabilized crosslinked nuclei",MaterialKind.BIOLOGICAL)
 graph.react("isolate_nuclei","Lyse cells and isolate permeabilized nuclei",StepKind.PERMEABILIZATION,
             ((crosslinked,InputRole.TEMPLATE),(nuclear,InputRole.CONDITION)),(nuclei,))
 hpy=graph.input("hpych4v","HpyCH4V (NEB R0620L), in nuclei 4 h 37 °C",MaterialKind.CATALYST)
 digested=graph.product("digested_nuclei","nuclei containing blunt HpyCH4V fragments (TGCA sites), mean 823 bp",
                        MaterialKind.BIOLOGICAL)
 graph.react("digest","Digest chromatin in nuclei with blunt-cutting HpyCH4V",StepKind.DIGESTION,
             ((nuclei,InputRole.TEMPLATE),(hpy,InputRole.CATALYST)),(digested,))
 tailing=graph.input("tailing_polymerase","NEBNext dA-Tailing Module (E6053L) Klenow fragment, 37 °C 90 min; EDTA/EGTA stop",
                     MaterialKind.CATALYST)
 tailed=graph.product("tailed_nuclei","nuclei containing dA-tailed genomic fragments",
                      MaterialKind.BIOLOGICAL)
 graph.react("da_tail","dA-tail the blunt fragments in nuclei",StepKind.TAILING,
             ((digested,InputRole.TEMPLATE),(tailing,InputRole.CATALYST)),(tailed,))
 cell_tags=graph.input("cell_tags","DPM, Odd and Even tag plates (2.4 µl of 45 µM per well; tags of the original SPRITE paper)",
                       MaterialKind.ADAPTER)
 ligase=graph.input("ligase","Instant Sticky-end Ligase Master Mix + Quick Ligase Buffer + 1,2-propanediol, 20 °C",
                    MaterialKind.CATALYST)
 cell_barcoded=graph.product("cell_barcoded_nuclei","nuclei carrying three-part cell barcodes",
                             MaterialKind.BIOLOGICAL)
 graph.react("cell_tag_ligation","Three in-nucleus split-pool ligations (3 h each; EDTA/EGTA stop, pool, wash)",
             StepKind.LIGATION,
             ((tailed,InputRole.TEMPLATE),(cell_tags,InputRole.STRUCTURAL),
              (ligase,InputRole.CATALYST)),(cell_barcoded,),
             geometry=AdapterLigationGeometry.unavailable(
                 "The Methods cite the original SPRITE paper for the DPM/Odd/Even tags; their ligation-ready ends are not assembled in this model."))
 mesh=graph.input("mesh_filter","10 µm pluriStrainer mesh",MaterialKind.SELECTION)
 filtered=graph.product("filtered_nuclei","filtered single nuclei",MaterialKind.BIOLOGICAL)
 clumps=graph.product("nucleus_clumps","nucleus aggregates retained by the mesh",MaterialKind.WASTE)
 graph.react("filter_nuclei","Filter nuclei through a 10 µm mesh to remove aggregates",StepKind.SELECTION,
             ((cell_barcoded,InputRole.TEMPLATE),(mesh,InputRole.SELECTION)),(filtered,clumps),
             changes_molecule=False)
 graph.discard(clumps,"retained on the mesh filter")
 selection=graph.input("nucleus_selection","count and withdraw 1,500 nuclei after lysis buffer #3 washes",
                       MaterialKind.SELECTION)
 selected=graph.product("selected_nuclei","1,500 selected cell-barcoded nuclei",MaterialKind.BIOLOGICAL)
 excluded=graph.product("excluded_nuclei","nuclei not withdrawn for the spatial branch",
                        MaterialKind.WASTE)
 graph.react("select_nuclei","Withdraw 1,500 nuclei for spatial-complex indexing",StepKind.SAMPLING,
             ((filtered,InputRole.TEMPLATE),(selection,InputRole.SELECTION)),
             (selected,excluded),changes_molecule=False)
 graph.discard(excluded,"nuclei outside the withdrawn 1,500 are not carried into the library")
 sonication=graph.input("sonication","Covaris M220, 15 µl lysis buffer #3, 2 min (30 W, duty 3.3, 6 °C)",
                        MaterialKind.CONDITION)
 complexes=graph.product("spatial_complexes","released crosslinked chromatin complexes",
                         MaterialKind.MOLECULAR,state=COMPLEX.complex())
 graph.react("sonicate","Sonicate all selected nuclei to release chromatin complexes",StepKind.FRAGMENTATION,
             ((selected,InputRole.TEMPLATE),(sonication,InputRole.CONDITION)),(complexes,))
 beads=graph.input("beads","NHS-activated magnetic beads (HCl-activated), overnight 4 °C in PBS/0.1% SDS; Tris quench",
                   MaterialKind.SELECTION)
 immobilized=graph.product("immobilized_complexes","covalently bead-coupled spatial complexes",
                           MaterialKind.MOLECULAR,state=COMPLEX.complex())
 uncaptured=graph.product("uncaptured_complexes","coupling flow-through and washes (RLT2+, M2)",
                          MaterialKind.WASTE)
 graph.react("immobilize","Couple complexes to NHS beads, quench and wash",StepKind.AFFINITY_CAPTURE,
             ((complexes,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),
             (immobilized,uncaptured))
 graph.discard(uncaptured,"removed on the magnet with the flow-through and washes")
 spatial_tags=graph.input("spatial_tags","Odd, Even and Y-even tag plates (2.4 µl of 4.5 µM per well)",
                          MaterialKind.ADAPTER)
 six_tagged=graph.product("six_tagged_complexes","six-part cell and spatial tagged complexes",
                          MaterialKind.MOLECULAR)
 graph.react("spatial_tag_ligation","Three on-bead split-pool ligations (60 min each; RLT2+ stop, pool, M2 washes)",
             StepKind.LIGATION,
             ((immobilized,InputRole.TEMPLATE),(spatial_tags,InputRole.STRUCTURAL),
              (ligase,InputRole.CATALYST)),(six_tagged,),
             geometry=AdapterLigationGeometry.unavailable(
                 "The Methods cite the original SPRITE protocol for spatial tagging; the Odd/Even/Y-even ligation-ready ends are not assembled in this model."))
 reversal=graph.input("crosslink_reversal","proteinase K in MyK buffer, 60 °C overnight (ten tubes)",
                      MaterialKind.CONDITION)
 released=graph.product("released_tagged_dna","de-crosslinked six-tag DNA released from beads",MaterialKind.MOLECULAR)
 graph.react("reverse_crosslinks","Reverse crosslinks and release tagged DNA",
             StepKind.DEPROTEINIZATION,
             ((six_tagged,InputRole.TEMPLATE),(reversal,InputRole.CONDITION)),(released,))
 column=graph.input("zymo_column","Zymo Clean & Concentrator-5 (5× binding buffer), elute 12 µl",MaterialKind.SELECTION)
 purified=graph.product("purified_tagged_dna","column-purified six-tag DNA",MaterialKind.MOLECULAR)
 column_waste=graph.product("column_waste","column flow-through and washes",MaterialKind.WASTE)
 graph.react("column_cleanup","Purify DNA on Zymo columns",StepKind.CLEANUP,
             ((released,InputRole.TEMPLATE),(column,InputRole.SELECTION)),(purified,column_waste))
 graph.discard(column_waste,"removed by the spin column")
 forward=graph.input("pcr_forward","indexed Illumina primer (sequence not assembled)",MaterialKind.OLIGO)
 reverse=graph.input("pcr_reverse","indexed Illumina primer partner (sequence not assembled)",MaterialKind.OLIGO)
 polymerase=graph.input("polymerase","Q5 Hot Start High-Fidelity 2× Master Mix + EvaGreen, real-time PCR 14 cycles",
                        MaterialKind.CATALYST)
 amplified=graph.product("amplified","amplified six-tag library",MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 graph.react("library_pcr","PCR-complete the six-tag library",StepKind.PCR,
             ((purified,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),(reverse,InputRole.STRUCTURAL),
              (polymerase,InputRole.CATALYST)),(amplified,),
             geometry=PCRGeometry.unavailable(
                 FINAL_LIBRARY,reason="The Methods name 25 µM indexed Illumina primers; a selected physical primer pair and its precursor binding sites are not assembled in this model."))
 gel=graph.input("gel_selection","1% E-gel, excise 300–1,000 bp; Zymo gel purification, elute 20 µl",
                 MaterialKind.SELECTION)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 gel_waste=graph.product("gel_waste","primer dimers, small non-specific and long amplicons outside 300–1,000 bp",
                         MaterialKind.WASTE)
 graph.react("gel_select","Gel size-select 300–1,000 bp and purify",StepKind.SELECTION,
             ((amplified,InputRole.TEMPLATE),(gel,InputRole.SELECTION)),(library,gel_waste))
 graph.discard(gel_waste,"left in the agarose outside the excised 300–1,000 bp window")
 graph.sequence(library,SequencingPlan(
     library=FINAL_LIBRARY,primers=(),
     acquired_roles=("Read 1","Read 2"),read_lengths=READ_LENGTHS,limitation=SEQUENCING_UNAVAILABLE,
     read_uses={"Read 1":"include","Read 2":"extract_then_include"}))
 graph.declare_reviewed_scope(tuple(graph.materials)+tuple(graph.steps),source)
 return graph.require_publishable()
