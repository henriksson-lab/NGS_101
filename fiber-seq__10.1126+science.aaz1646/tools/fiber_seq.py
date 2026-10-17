"""Stergachis et al. Fiber-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Row, Scene, Segment
from dumbbell import smrtbell_circuit, smrtbell_hairpin_geometry, smrtbell_sequencing_plan
from endprep import dA_tailed_scene, repair_and_dA_tail
from longread_footprinting import conventional_smrtbell, marked_fiber, pacbio_entry, smrtbell_rows
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SourceClaim,
                            SourceDiscrepancy, StepKind)

TITLE = "Fiber-seq — long-read chromatin fibers"
NOTES = "01_fiber-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1126/science.aaz1646">Stergachis et al., <i>Science</i> (2020)</a>.'
SUMMARY = "Hia5 methylates exposed adenines in intact nuclei; PCR-free PacBio sequencing preserves that m6dA stencil together with the underlying long DNA molecule."
CAVEAT = ("The 2020 Supplementary Materials (Materials and Methods, 'Fiber-seq' and 'AdMTase-seq' paragraphs) give the original reaction: nuclei from 2.5 million S2 or 300,000 K562 cells (0.05% IGEPAL in Buffer A on ice, 10 min) are marked in 60 &micro;l Buffer A with in-house purified Hia5 and 0.8 mM SAM for 10 min at 37 &deg;C, stopped with SDS to 1%, treated with RNase A and proteinase K, phenol-chloroform extracted and ethanol precipitated twice, AMPure-purified, sheared on a Megaruptor to 10&ndash;15 kb and made into a PCR-free PacBio library for Sequel II (System 2.0; half of one K562 sample on System 1.0; the 20 U replicates on Sequel). The supplement names no SMRTbell prep kit, so the repair, dA-tailing, hairpin-ligation and nuclease steps follow the authors' current bench protocol (SMRTbell prep kit 3.0) as a later implementation. The current protocol marks at 25 &deg;C with commercial Hia5; this is a version difference, not a source error. Hairpin bases are proprietary and are not invented.")
SMRTBELL = conventional_smrtbell("m6dA-marked genomic insert")
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
SEQUENCING_ENDING = pacbio_entry(SMRTBELL)

def sections():
    return [
        ("Stencil accessible DNA with Hia5", marked_fiber("Hia5", target="accessible adenines"),
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
    authority=SOURCE+(' 2020 reaction, stop, extraction and shearing: Science Supplementary Materials, Materials and Methods '
                      '(\'AdMTase-seq\' and \'Fiber-seq\'); sample records: '
                      '<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146942">GEO GSE146942</a>.'
                      ' SMRTbell steps (kit not named in 2020): '
                      '<a href="https://fiberseq.github.io/protocol/bench-protocol.html">authors\' current Fiber-seq bench protocol</a> (SMRTbell prep kit 3.0).')
    graph=ProtocolGraph(TITLE,authority=authority)
    nuclei=graph.input("nuclei","intact nuclei (2020 SM: 2.5 million S2 or 300,000 K562 cells; 0.05% IGEPAL in Buffer A on ice, 10 min) containing native chromatin",
                       MaterialKind.BIOLOGICAL,primary=True)
    hia5=graph.input("hia5","Hia5 adenine methyltransferase (2020: in-house purified; SM 200 U, GEO 500 U, or 20 U low dose; 10 min at 37 °C)",MaterialKind.CATALYST)
    sam=graph.input("sam","S-adenosyl methionine (0.8 mM final; SM)",MaterialKind.MOLECULAR)
    marked=graph.product("marked_chromatin","Hia5-stencilled native chromatin",
                         MaterialKind.MOLECULAR,state=tuple(marked_fiber("Hia5", target="accessible adenines")))
    graph.react("mark_accessible","Methylate solvent-accessible adenines with Hia5",
                StepKind.CHEMICAL_CONVERSION,
                ((nuclei,InputRole.TEMPLATE),(hia5,InputRole.CATALYST),
                 (sam,InputRole.STRUCTURAL)),(marked,))
    graph.annotate(SourceDiscrepancy(
        "mark_accessible","Published Hia5 amount for the main Fiber-seq samples differs",
        (SourceClaim("Science Supplementary Materials, Materials and Methods, 'Fiber-seq'",
                     "MTase reaction 'with the exception of using 200 units (2.5 μl) of Hia5', 10 min; "
                     "500 U only for the deproteinized S2 control"),
         SourceClaim("GEO GSE146942 sample records (e.g. GSM4411218, GSM4411222-GSM4411224)",
                     "titles 'Fiber-seq (500 U Hia5)'; characteristics '500 U Hia5 for 10min at 37C'")),
        "The model does not choose a dose: the Hia5 input lists 200 U (SM) and 500 U (GEO); "
        "the 20 U replicates agree in both sources.",
        "Both are author-published records of the same 2020 samples.",
        "Which amount was used for the main S2 and K562 Fiber-seq libraries is unresolved."))
    # 2020 SM 'Fiber-seq': stopped with 3 µl 20% SDS (1% final); current bench protocol agrees.
    sds=graph.input("sds","SDS to 1% final (2020 SM; same in the current protocol)",MaterialKind.CONDITION)
    quenched=graph.product("quenched","denatured m6dA-marked chromatin lysate",
                           MaterialKind.MOLECULAR,state=marked.state)
    graph.react("quench","Stop methylation by denaturing Hia5",StepKind.QUENCHING,
                ((marked,InputRole.TEMPLATE),(sds,InputRole.CONDITION)),(quenched,),
                changes_molecule=False)
    extraction=graph.input("extraction","RNase A + proteinase K (50 °C, 1 h), two phenol-chloroform extractions and ethanol precipitations, AMPure (2020 SM; current protocol uses MagAttract or Wizard HMW kits)",
                           MaterialKind.CONDITION)
    dna=graph.product("marked_dna","native high-molecular-weight m6dA-marked DNA",
                      MaterialKind.MOLECULAR,state=SMRTBELL.insert)
    graph.react("extract","Extract native marked DNA without amplification",
                StepKind.EXTRACTION,
                ((quenched,InputRole.TEMPLATE),(extraction,InputRole.CONDITION)),(dna,))
    shear=graph.input("shear_condition","Megaruptor shearing to 10-15 kb (2020 SM; Megaruptor 3 or g-TUBE in the current protocol)",
                      MaterialKind.CONDITION)
    sheared=graph.product("sheared","long m6dA-marked DNA inserts",MaterialKind.MOLECULAR,
                          state=SMRTBELL.insert)
    graph.react("shear","Shear marked DNA to 10-15 kb inserts",StepKind.FRAGMENTATION,
                ((dna,InputRole.TEMPLATE),(shear,InputRole.CONDITION)),(sheared,))
    end_mix=graph.input("end_mix","damage-repair, end-repair and dA-tailing mix (kit not named in 2020 SM; SMRTbell prep kit 3.0 in the current protocol)",
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
        "nuclei","hia5","sam","mark_accessible","sds","quench","extraction","extract",
        "shear_condition","shear","end_mix","end_prep","hairpins","ligase",
        "ligate_hairpins","nuclease","remove_open","library"),authority)
    return graph.require_publishable()
