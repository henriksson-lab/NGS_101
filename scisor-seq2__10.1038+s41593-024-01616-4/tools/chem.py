"""ScISOr-seq2 barcoded 10x cDNA enrichment for PacBio and ONT."""
from batch_ngs import seg
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, MolecularState, Row, Scene, Workflow, feature, revcomp
from dumbbell import smrtbell_circuit, smrtbell_hairpin_geometry, smrtbell_sequencing_plan
from endprep import dA_tailed_scene, repair_and_dA_tail
from longread_footprinting import conventional_smrtbell, pacbio_entry, smrtbell_rows
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
TITLE="ScISOr-seq2 — enriched single-cell isoforms by long reads"
NOTES="01_scisor-seq2.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41593-024-01616-4">Joglekar et al., <i>Nature Neuroscience</i> (2024)</a>.'
SUMMARY="10x 3-prime cDNA is strand-selectively re-amplified to retain cell barcodes, exon-capture probes enrich spliced molecules, and the same enriched cDNA can enter PacBio HiFi or ONT ligation sequencing."
PARTIAL_READ1 = "CTACACGACGCTCTTCCGATCT"
PARTIAL_TSO = "AAGCAGTGGTATCAACGCAGAGTACAT"
METHODS = "Joglekar et al. 2024 Methods (PMC11156538)"
INSERT=Construct([seg("Partial Read 1 handle",PARTIAL_READ1,"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True,feature=feature("cell_barcode","cell_barcode","whitelist",whitelist="10x Chromium 3-prime v3 whitelist")),seg("UMI","U"*12,"umi",placeholder=True,feature=feature("umi","umi","random")),seg("full-length spliced cDNA","X"*42,placeholder=True),seg("Partial TSO handle reverse complement",revcomp(PARTIAL_TSO),"tso")],name="LAP-CAP enriched barcoded cDNA")
SMRTBELL=conventional_smrtbell(INSERT)
FINAL_LIBRARY=None; SEQ_PRIMERS=(); SEQUENCING_ENDING=pacbio_entry(SMRTBELL)+" The same enriched cDNA may instead receive an ONT motor adapter, in which case no synthesis sequencing primer is used."
def oligos():
    return [f"<b>Partial Read 1</b> — 5′-{PARTIAL_READ1}-3′",f"<b>Partial TSO</b> — 5′-{PARTIAL_TSO}-3′"]
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


def _spri(graph,prefix,material,label):
    kept=graph.product(prefix+"_purified",material.label+" after SPRI",
                       MaterialKind.MOLECULAR,state=material.state)
    waste=graph.product(prefix+"_spri_waste","SPRI supernatant and washes",MaterialKind.WASTE)
    graph.react(prefix+"_spri",label,StepKind.CLEANUP,((material,InputRole.TEMPLATE),),
                (kept,waste),changes_molecule=False)
    graph.discard(waste,"left in the SPRI supernatant")
    return kept


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    cdna=graph.input("tenx_cdna","30 ng 10x Chromium Single Cell 3-prime v3 barcoded cDNA",
                     MaterialKind.MOLECULAR,state=INSERT,primary=True)
    partial_r1=graph.input("partial_r1","published Partial Read 1 primer",
                           MaterialKind.OLIGO,state=Construct([seg("Partial Read1",PARTIAL_READ1)]),
                           source=METHODS+", LAP")
    polymerase=graph.input("polymerase","KAPA HiFi HotStart ReadyMix (2x)",
                           MaterialKind.CATALYST)
    biased=graph.product("biased","barcode-bearing strand-biased LAP products",
                         MaterialKind.MOLECULAR,state=INSERT)
    graph.react("linear_amplification","12 linear/asymmetric cycles with Partial Read1 only (95 °C 3 min; 98 °C 20 s, 64 °C 30 s, 72 °C 60 s)",
                StepKind.EXTENSION,
                ((cdna,InputRole.TEMPLATE),(partial_r1,InputRole.STRUCTURAL),
                 (polymerase,InputRole.CATALYST)),(biased,))
    biased=_spri(graph,"linear",biased,"0.8x SPRIselect, two 80% ethanol washes")
    partial_tso=graph.input("partial_tso","published Partial TSO primer",MaterialKind.OLIGO,
                            state=Construct([seg("Partial TSO",PARTIAL_TSO)]),
                            source=METHODS+", LAP")
    amplified=graph.product("amplified","full LAP amplicon pool",MaterialKind.MOLECULAR,
                            state=INSERT)
    graph.react("exponential_pcr","Six exponential cycles with Partial Read1 + Partial TSO (same conditions)",
                StepKind.PCR,
                ((biased,InputRole.TEMPLATE),(partial_r1,InputRole.STRUCTURAL),
                 (partial_tso,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (amplified,),geometry=PCRGeometry.exact(
                    INSERT,
                    PCRPrimer("Partial Read 1", PARTIAL_READ1),
                    PCRPrimer("Partial TSO", PARTIAL_TSO),
                    statement=("The defining Methods print both complete primers; their "
                               "terminal 10x handles flank every LAP amplicon independently "
                               "of the transcript's internal coordinates.")))
    amplified=_spri(graph,"exponential",amplified,"0.6x SPRIselect, ethanol washes, elute 30 µl EB")
    blockers=graph.input("block_oligos","block oligo mix: Partial Read1 + Partial TSO (100 ng/µl)",
                         MaterialKind.OLIGO,state=(PARTIAL_READ1,PARTIAL_TSO),
                         source=METHODS+", CAP block oligo mix")
    probes=graph.input("exon_probes","biotinylated Agilent SureSelect All Exon probes (Human V8 or Mouse) with SureSelectXT HSQ reagents",
                       MaterialKind.SELECTION)
    beads=graph.input("streptavidin_beads","M-270 Streptavidin Dynabeads; WB1 RT and WB2 65 °C washes",
                      MaterialKind.SELECTION)
    enriched=graph.product("enriched","CAP-enriched spliced barcoded cDNA",
                           MaterialKind.MOLECULAR,state=INSERT)
    rejected=graph.product("uncaptured","uncaptured cDNA",MaterialKind.WASTE)
    graph.react("capture","Block handles (95 °C 5 min, 65 °C), hybridize exon probes 24 h at 65 °C, capture on streptavidin beads and wash",
                StepKind.AFFINITY_CAPTURE,
                ((amplified,InputRole.TEMPLATE),(blockers,InputRole.STRUCTURAL),
                 (probes,InputRole.SELECTION),(beads,InputRole.SELECTION)),
                (enriched,rejected))
    graph.discard(rejected,"removed with the bead washes; not hybridized to the exon-probe pool")
    post=graph.product("post_capture","amplified captured spliced cDNA",
                       MaterialKind.MOLECULAR,state=INSERT)
    graph.react("post_capture_pcr","12 on-bead cycles with Partial Read1 + Partial TSO (98 °C 20 s, 64 °C 60 s, 72 °C 3 min); keep supernatant",
                StepKind.PCR,
                ((enriched,InputRole.TEMPLATE),(partial_r1,InputRole.STRUCTURAL),
                 (partial_tso,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (post,),geometry=PCRGeometry.exact(
                    INSERT,
                    PCRPrimer("Partial Read 1", PARTIAL_READ1),
                    PCRPrimer("Partial TSO", PARTIAL_TSO),
                    statement="Methods name the same two printed primers for post-capture amplification."))
    enriched=_spri(graph,"post_capture",post,"0.8x SPRIselect")
    pacbio_input=graph.product("pacbio_input","PacBio-directed enriched cDNA aliquot",
                               MaterialKind.MOLECULAR,state=INSERT)
    ont_input=graph.product("ont_input","ONT-directed enriched cDNA aliquot",
                            MaterialKind.MOLECULAR,state=INSERT)
    graph.react("split_platforms","Split enriched cDNA between PacBio and ONT workflows",
                StepKind.SPLIT,((enriched,InputRole.TEMPLATE),),(pacbio_input,ont_input),
                changes_molecule=False)

    damage_mix=graph.input("pacbio_damage_mix","SMRTbell Express 2.0 DNA damage repair mix",
                           MaterialKind.CATALYST)
    repaired=graph.product("pacbio_repaired","damage-repaired enriched cDNA (~500 ng)",
                           MaterialKind.MOLECULAR,state=INSERT)
    graph.react("pacbio_damage_repair","DNA damage repair, 37 °C 30 min",StepKind.END_REPAIR,
                ((pacbio_input,InputRole.TEMPLATE),(damage_mix,InputRole.CATALYST)),(repaired,))
    end_mix=graph.input("pacbio_end_mix","SMRTbell Express 2.0 end-repair/A-tailing mix",
                        MaterialKind.CATALYST)
    tailed_state=repair_and_dA_tail(INSERT,inferred=True)
    tailed=graph.product("pacbio_tailed","repaired dA-tailed enriched cDNA",
                         MaterialKind.MOLECULAR,state=dA_tailed_scene(tailed_state))
    graph.react("pacbio_end_prep","End-repair/A-tailing, 20 °C 30 min then 65 °C 30 min",
                StepKind.END_REPAIR,
                ((repaired,InputRole.TEMPLATE),(end_mix,InputRole.CATALYST)),(tailed,))
    hairpins=graph.input("hairpins","T-overhang PacBio SMRTbell hairpins",
                         MaterialKind.ADAPTER,inferred=True,
                         source=SOURCE+" (hairpin-body bases proprietary)")
    smrt_ligase=graph.input("smrt_ligase","SMRTbell ligase",MaterialKind.CATALYST)
    ligated=graph.product("pacbio_ligated","ligated SMRTbell reaction",MaterialKind.MOLECULAR,
                          state=tuple(smrtbell_rows(SMRTBELL)),inferred=True,
                          source=SOURCE+" (hairpin-body bases proprietary)")
    pacbio_library=graph.product("pacbio_library",SMRTBELL.name,MaterialKind.LIBRARY,
                                 state=tuple(smrtbell_rows(SMRTBELL)),inferred=True,
                                 source=SOURCE+" (hairpin-body bases proprietary)")
    graph.react("ligate_hairpins","Ligate PacBio hairpins, 20 °C 60 min",StepKind.LIGATION,
                ((tailed,InputRole.TEMPLATE),(hairpins,InputRole.STRUCTURAL),
                 (smrt_ligase,InputRole.CATALYST)),(ligated,),
                geometry=smrtbell_hairpin_geometry(
                    SMRTBELL,statement="PacBio SMRTbell preparation ligates T-overhang hairpins to dA-tailed cDNA."))
    pacbio_waste=graph.product("pacbio_spri_waste","unligated adapters and short species",
                               MaterialKind.WASTE)
    graph.react("pacbio_cleanup","0.6x SPRIselect purification",StepKind.CLEANUP,
                ((ligated,InputRole.TEMPLATE),),(pacbio_library,pacbio_waste),
                changes_molecule=False)
    graph.discard(pacbio_waste,"left in the SPRI supernatant")
    graph.sequence(pacbio_library,smrtbell_sequencing_plan(
        SMRTBELL,smrtbell_circuit(SMRTBELL),
        limitation="The proprietary PacBio sequencing-primer site lies in the hairpin adapter."))

    ont_adapter=graph.input("ont_adapter","motor-bearing ONT ligation adapter",
                            MaterialKind.ADAPTER,inferred=True,
                            source=SOURCE+" (motor-adapter terminal bases proprietary)")
    ont_ligase=graph.input("ont_ligase","SQK-LSK110 end-prep and adapter-ligation reagents (~75 fmol cDNA)",MaterialKind.CATALYST)
    ont_library=graph.product("ont_library","motor-loaded ScISOr-seq2 ONT library",
                              MaterialKind.LIBRARY,inferred=True,
                              source=SOURCE+" (motor-adapter terminal bases proprietary)")
    graph.react("ligate_ont","Prepare and ligate the motor-bearing ONT adapter per protocol ACDE_9110_v110_revC_10Nov2020 (vendor steps not inspected)",StepKind.LIGATION,
                ((ont_input,InputRole.TEMPLATE),(ont_adapter,InputRole.STRUCTURAL),
                 (ont_ligase,InputRole.CATALYST)),(ont_library,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The source names the ONT ligation workflow but does not disclose the motor-adapter terminal bases."))
    graph.sequence(ont_library,SequencingPlan(
        acquired_roles=("Sequencing output",),
        limitation="PromethION FLO-PRO002, 72 h; Guppy basecalls with quality >7. Nanopore entry is motor-mediated and primer-free; adapter bases are proprietary."))
    graph.declare_reviewed_scope(tuple(graph.materials)+tuple(graph.steps),SOURCE)
    return graph.require_publishable()
