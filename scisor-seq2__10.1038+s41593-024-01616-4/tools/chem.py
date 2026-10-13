"""ScISOr-seq2 barcoded 10x cDNA enrichment for PacBio and ONT."""
from batch_ngs import seg
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, MolecularState, Row, Scene, Workflow, feature, revcomp
from dumbbell import smrtbell_circuit, smrtbell_hairpin_geometry, smrtbell_sequencing_plan
from endprep import dA_tailed_scene, repair_and_dA_tail
from longread_footprinting import conventional_smrtbell, pacbio_entry, smrtbell_rows
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
TITLE="ScISOr-seq2 — enriched single-cell isoforms by long reads"
NOTES="01_scisor-seq2.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41593-024-01616-4">Joglekar et al., <i>Nature Neuroscience</i> (2024)</a>.'
SUMMARY="10x 3-prime cDNA is strand-selectively re-amplified to retain cell barcodes, exon-capture probes enrich spliced molecules, and the same enriched cDNA can enter PacBio HiFi or ONT ligation sequencing."
INSERT=Construct([seg("Partial Read 1 handle","CTACACGACGCTCTTCCGATCT","r1"),seg("cell barcode","B"*16,"cbc",placeholder=True,feature=feature("cell_barcode","cell_barcode","whitelist",whitelist="10x Chromium 3-prime v3 whitelist")),seg("UMI","U"*12,"umi",placeholder=True,feature=feature("umi","umi","random")),seg("full-length spliced cDNA","X"*42,placeholder=True),seg("Partial TSO handle","AAGCAGTGGTATCAACGCAGAGTACAT","tso")],name="LAP-CAP enriched barcoded cDNA")
SMRTBELL=conventional_smrtbell(INSERT)
FINAL_LIBRARY=None; SEQ_PRIMERS=(); SEQUENCING_ENDING=pacbio_entry(SMRTBELL)+" The same enriched cDNA may instead receive an ONT motor adapter, in which case no synthesis sequencing primer is used."
def oligos():
    return ["<b>Partial Read 1</b> — 5′-CTACACGACGCTCTTCCGATCT-3′","<b>Partial TSO</b> — 5′-AAGCAGTGGTATCAACGCAGAGTACAT-3′"]
def workflow():
    w=Workflow(MolecularState("10x 3′ cDNA",(Row(chunks=[("Read1 handle—cell barcode—UMI—transcript cDNA—TSO","cbc",False)]),)))
    w.react("Perform linear/asymmetric then exponential amplification",[Row(chunks=[("Partial Read1-only cycles → strand bias; + Partial TSO → full amplicon","r1",False)])],name="LAP amplicons",note="LAP removes non-barcoded products while retaining both terminal handles.")
    w.react("Capture spliced cDNA with exon probes",[Row(chunks=[("blocked cDNA + biotin exon probes ● streptavidin → CAP enrichment","w1",False)])],name="CAP-enriched cDNA",note="Hybrid capture enriches multi-exonic transcripts.")
    tailed=repair_and_dA_tail(Construct([seg("ScISOr-seq2 enriched cDNA","X"*42,placeholder=True)]))
    w.react("Repair and add 3′ dA",dA_tailed_scene(tailed).rows(),name="dA-tailed enriched cDNA",note="PacBio end prep creates the two hairpin-ligation substrates.")
    adapter=Scene(); adapter.strand("hairpin end",[seg("proprietary PacBio hairpin core","[proprietary hairpin]",placeholder=True),seg("3′ T","T")],label="PacBio hairpin ligating end"); adapter.overhang("hairpin end","3′ T",polarity="3-prime",owner="adapter")
    geometry=AdapterLigationGeometry.known("cohesive",dA_tailed_scene(tailed),adapter,statement="The PacBio SMRTbell workflow uses proprietary hairpin sequences with a disclosed 3′-T ligating overhang complementary to each insert 3′-dA.")
    w.ligate_adapter("Ligate PacBio hairpins",geometry,smrtbell_rows(SMRTBELL),name="closed SMRTbell",note="The ONT alternative instead uses a proprietary motor-bearing ligation adapter whose end bases are not disclosed here.")
    return w


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    cdna=graph.input("tenx_cdna","10x 3-prime full-length barcoded cDNA",
                     MaterialKind.MOLECULAR,state=INSERT,primary=True)
    partial_r1=graph.input("partial_r1","published Partial Read 1 primer",
                           MaterialKind.OLIGO)
    polymerase=graph.input("polymerase","LAP amplification polymerase",
                           MaterialKind.CATALYST)
    biased=graph.product("biased","barcode-bearing strand-biased LAP products",
                         MaterialKind.MOLECULAR,state=INSERT)
    graph.react("linear_amplification","Preferentially copy the barcode-bearing strand",
                StepKind.EXTENSION,
                ((cdna,InputRole.TEMPLATE),(partial_r1,InputRole.STRUCTURAL),
                 (polymerase,InputRole.CATALYST)),(biased,))
    partial_tso=graph.input("partial_tso","published Partial TSO primer",MaterialKind.OLIGO)
    amplified=graph.product("amplified","full LAP amplicon pool",MaterialKind.MOLECULAR,
                            state=INSERT)
    graph.react("exponential_pcr","Complete exponential amplification with Partial TSO",
                StepKind.PCR,
                ((biased,InputRole.TEMPLATE),(partial_r1,InputRole.STRUCTURAL),
                 (partial_tso,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (amplified,),geometry=PCRGeometry.unavailable(
                    INSERT,reason="Both exact partial primers are published, but the role-level cDNA construct does not assert a universal transcript-end coordinate."))
    probes=graph.input("exon_probes","biotinylated exon-capture probe pool",
                       MaterialKind.SELECTION)
    enriched=graph.product("enriched","CAP-enriched spliced barcoded cDNA",
                           MaterialKind.MOLECULAR,state=INSERT)
    rejected=graph.product("uncaptured","uncaptured cDNA",MaterialKind.WASTE)
    graph.react("capture","Capture spliced cDNA with biotinylated exon probes",
                StepKind.AFFINITY_CAPTURE,
                ((amplified,InputRole.TEMPLATE),(probes,InputRole.SELECTION)),
                (enriched,rejected))
    graph.discard(rejected,"does not hybridize to the selected exon-probe pool")
    pacbio_input=graph.product("pacbio_input","PacBio-directed enriched cDNA aliquot",
                               MaterialKind.MOLECULAR,state=INSERT)
    ont_input=graph.product("ont_input","ONT-directed enriched cDNA aliquot",
                            MaterialKind.MOLECULAR,state=INSERT)
    graph.react("split_platforms","Split enriched cDNA between PacBio and ONT workflows",
                StepKind.SPLIT,((enriched,InputRole.TEMPLATE),),(pacbio_input,ont_input),
                changes_molecule=False)

    end_mix=graph.input("pacbio_end_mix","PacBio repair and dA-tailing mix",
                        MaterialKind.CATALYST)
    tailed_state=repair_and_dA_tail(INSERT,inferred=True)
    tailed=graph.product("pacbio_tailed","repaired dA-tailed enriched cDNA",
                         MaterialKind.MOLECULAR,state=dA_tailed_scene(tailed_state))
    graph.react("pacbio_end_prep","Repair and dA-tail the PacBio aliquot",
                StepKind.END_REPAIR,
                ((pacbio_input,InputRole.TEMPLATE),(end_mix,InputRole.CATALYST)),(tailed,))
    hairpins=graph.input("hairpins","T-overhang PacBio SMRTbell hairpins",
                         MaterialKind.ADAPTER,inferred=True,
                         source=SOURCE+" (hairpin-body bases proprietary)")
    smrt_ligase=graph.input("smrt_ligase","SMRTbell ligase",MaterialKind.CATALYST)
    pacbio_library=graph.product("pacbio_library",SMRTBELL.name,MaterialKind.LIBRARY,
                                 state=tuple(smrtbell_rows(SMRTBELL)),inferred=True,
                                 source=SOURCE+" (hairpin-body bases proprietary)")
    graph.react("ligate_hairpins","Ligate PacBio hairpins",StepKind.LIGATION,
                ((tailed,InputRole.TEMPLATE),(hairpins,InputRole.STRUCTURAL),
                 (smrt_ligase,InputRole.CATALYST)),(pacbio_library,),
                geometry=smrtbell_hairpin_geometry(
                    SMRTBELL,statement="PacBio SMRTbell preparation ligates T-overhang hairpins to dA-tailed cDNA."))
    graph.sequence(pacbio_library,smrtbell_sequencing_plan(
        SMRTBELL,smrtbell_circuit(SMRTBELL),
        limitation="The proprietary PacBio sequencing-primer site lies in the hairpin adapter."))

    ont_adapter=graph.input("ont_adapter","motor-bearing ONT ligation adapter",
                            MaterialKind.ADAPTER,inferred=True,
                            source=SOURCE+" (motor-adapter terminal bases proprietary)")
    ont_ligase=graph.input("ont_ligase","ONT adapter ligase",MaterialKind.CATALYST)
    ont_library=graph.product("ont_library","motor-loaded ScISOr-seq2 ONT library",
                              MaterialKind.LIBRARY,inferred=True,
                              source=SOURCE+" (motor-adapter terminal bases proprietary)")
    graph.react("ligate_ont","Attach the motor-bearing ONT adapter",StepKind.LIGATION,
                ((ont_input,InputRole.TEMPLATE),(ont_adapter,InputRole.STRUCTURAL),
                 (ont_ligase,InputRole.CATALYST)),(ont_library,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The source names the ONT ligation workflow but does not disclose the motor-adapter terminal bases."))
    graph.sequence(ont_library,SequencingPlan(
        acquired_roles=("Sequencing output",),
        limitation="Nanopore entry is motor-mediated and primer-free; adapter bases are proprietary."))
    graph.declare_reviewed_scope((
        "tenx_cdna","partial_r1","polymerase","linear_amplification","partial_tso",
        "exponential_pcr","exon_probes","capture","split_platforms","pacbio_end_mix",
        "pacbio_end_prep","hairpins","smrt_ligase","ligate_hairpins","pacbio_library",
        "ont_adapter","ont_ligase","ligate_ont","ont_library"),SOURCE)
    return graph.require_publishable()
