"""Stergachis et al. Fiber-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Row, Scene, Segment
from dumbbell import smrtbell_circuit, smrtbell_hairpin_geometry, smrtbell_sequencing_plan
from endprep import dA_tailed_scene, repair_and_dA_tail
from longread_footprinting import conventional_smrtbell, marked_fiber, pacbio_entry, smrtbell_rows
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind

TITLE = "Fiber-seq — long-read chromatin fibers"
NOTES = "01_fiber-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1126/science.aaz1646">Stergachis et al., <i>Science</i> (2020)</a>.'
SUMMARY = "Hia5 methylates exposed adenines in intact nuclei; PCR-free PacBio sequencing preserves that m6dA stencil together with the underlying long DNA molecule."
CAVEAT = "The defining paper establishes Hia5 marking and PacBio single-molecule readout. The page uses the current standard PCR-free SMRTbell geometry; hairpin bases are proprietary and are not invented."
SMRTBELL = conventional_smrtbell()
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
SEQUENCING_ENDING = pacbio_entry(SMRTBELL)

def sections():
    return [
        ("Stencil accessible DNA with Hia5", marked_fiber("Hia5"),
         "SAM-dependent m6dA is deposited where chromatin does not protect adenines."),
        ("Extract high-molecular-weight marked DNA", [Row(chunks=[("one native, m6dA-bearing DNA fiber", "w1", False)])],
         "No amplification is introduced: the modification pattern must remain on the sequenced template."),
        ("Shear, repair and add one 3-prime dA to each end",
         Scene.duplex([Segment("3-prime dA", "A"), Segment("m6dA-marked long DNA", "X" * 48, placeholder=True), Segment("3-prime dA", "A")], label="dA-tailed marked DNA").rows(),
         "The long marked insert is shown before either SMRTbell hairpin is attached."),
        ("Ligate PacBio hairpins", smrtbell_rows(SMRTBELL),
         "Both hairpins close the duplex into a nuclease-resistant SMRTbell; the drawing marks both sealed insert-adapter junctions."),
    ]


def protocol_graph():
    authority=SOURCE+" Current PCR-free SMRTbell endpoint: PacBio SMRTbell prep kit 3.0."
    graph=ProtocolGraph(TITLE,authority=authority)
    nuclei=graph.input("nuclei","intact nuclei containing native chromatin",
                       MaterialKind.BIOLOGICAL,primary=True)
    hia5=graph.input("hia5","Hia5 adenine methyltransferase",MaterialKind.CATALYST)
    sam=graph.input("sam","S-adenosyl methionine methyl donor",MaterialKind.MOLECULAR)
    marked=graph.product("marked_chromatin","Hia5-stencilled native chromatin",
                         MaterialKind.MOLECULAR,state=tuple(marked_fiber("Hia5")))
    graph.react("mark_accessible","Methylate solvent-accessible adenines with Hia5",
                StepKind.CHEMICAL_CONVERSION,
                ((nuclei,InputRole.TEMPLATE),(hia5,InputRole.CATALYST),
                 (sam,InputRole.STRUCTURAL)),(marked,))
    extraction=graph.input("extraction","high-molecular-weight DNA extraction condition",
                           MaterialKind.CONDITION)
    dna=graph.product("marked_dna","native high-molecular-weight m6dA-marked DNA",
                      MaterialKind.MOLECULAR,state=SMRTBELL.insert)
    graph.react("extract","Extract native marked DNA without amplification",
                StepKind.EXTRACTION,
                ((marked,InputRole.TEMPLATE),(extraction,InputRole.CONDITION)),(dna,))
    shear=graph.input("shear_condition","long-insert shearing condition",
                      MaterialKind.CONDITION)
    sheared=graph.product("sheared","long m6dA-marked DNA inserts",MaterialKind.MOLECULAR,
                          state=SMRTBELL.insert)
    graph.react("shear","Shear marked DNA to long inserts",StepKind.FRAGMENTATION,
                ((dna,InputRole.TEMPLATE),(shear,InputRole.CONDITION)),(sheared,))
    end_mix=graph.input("end_mix","damage-repair, end-repair and dA-tailing mix",
                        MaterialKind.CATALYST)
    end_product=repair_and_dA_tail(SMRTBELL.insert,inferred=True)
    tailed=graph.product("tailed","repaired dA-tailed marked DNA inserts",
                         MaterialKind.MOLECULAR,state=dA_tailed_scene(end_product))
    graph.react("end_prep","Repair and dA-tail both insert ends",StepKind.END_REPAIR,
                ((sheared,InputRole.TEMPLATE),(end_mix,InputRole.CATALYST)),(tailed,))
    hairpins=graph.input("hairpins","T-overhang PacBio SMRTbell hairpin adapters",
                         MaterialKind.ADAPTER,inferred=True,
                         source=authority+" (hairpin-body bases proprietary)")
    ligase=graph.input("ligase","SMRTbell ligase",MaterialKind.CATALYST)
    closed=graph.product("closed","closed m6dA-marked SMRTbell molecules",
                         MaterialKind.MOLECULAR,state=tuple(smrtbell_rows(SMRTBELL)),
                         inferred=True,source=authority+" (hairpin-body bases proprietary)")
    graph.react("ligate_hairpins","Ligate hairpins to both marked-DNA ends",
                StepKind.LIGATION,
                ((tailed,InputRole.TEMPLATE),(hairpins,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(closed,),
                geometry=smrtbell_hairpin_geometry(
                    SMRTBELL,statement="Current PCR-free Fiber-seq preparation ligates T-overhang SMRTbell hairpins to dA-tailed marked DNA."))
    nuclease=graph.input("nuclease","free-end nuclease",MaterialKind.CATALYST)
    library=graph.product("library",SMRTBELL.name,MaterialKind.LIBRARY,
                          state=tuple(smrtbell_rows(SMRTBELL)),inferred=True,
                          source=closed.source)
    open_waste=graph.product("open_molecules","molecules retaining exposed DNA ends",
                             MaterialKind.WASTE)
    graph.react("remove_open","Digest molecules retaining exposed ends",StepKind.DIGESTION,
                ((closed,InputRole.TEMPLATE),(nuclease,InputRole.CATALYST)),
                (library,open_waste),changes_molecule=False)
    graph.discard(open_waste,"does not form a completely closed SMRTbell")
    graph.sequence(library,smrtbell_sequencing_plan(
        SMRTBELL,smrtbell_circuit(SMRTBELL),
        limitation="The current PacBio hairpin and sequencing-primer bases are proprietary."))
    graph.declare_reviewed_scope((
        "nuclei","hia5","sam","mark_accessible","extraction","extract",
        "shear_condition","shear","end_mix","end_prep","hairpins","ligase",
        "ligate_hairpins","nuclease","remove_open","library"),authority)
    return graph.require_publishable()
