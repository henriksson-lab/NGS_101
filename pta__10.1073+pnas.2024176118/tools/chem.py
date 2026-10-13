"""Primary template-directed amplification and its defining library route."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import pta
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import seg, truseq_library
from chemdraw import Construct, MolecularState, Scene, Workflow
from endprep import dA_tailed_scene, repair_and_dA_tail
from pcr import PCRGeometry
from protocol_graph import (InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind)


TITLE = "PTA — primary template-directed whole-genome amplification"
NOTES = "01_pta.html"
SOURCE = ('<a href="https://doi.org/10.1073/pnas.2024176118">Gonzalez-Pena et al., '
          '<i>PNAS</i> (2021)</a>, Methods and Supplementary Figure S1.')
SUMMARY = ("Random-primed phi29 synthesis is interrupted by exonuclease-resistant "
           "alpha-thio-ddNTPs. Shorter products undergo less daughter-on-daughter "
           "amplification, biasing yield toward repeated copying of the primary genome.")
CAVEAT = ("The defining paper does not publish the random-primer or IDT UDI adapter "
          "sequences. Their roles are preserved; canonical TruSeq adapter bases are "
          "shown as inferred only to locate standard sequencing primers.")

FINAL_LIBRARY, SEQ_PRIMERS = truseq_library(
    [seg("PTA-amplified genomic insert", "X" * 48, placeholder=True)],
    "PTA KAPA HyperPlus library", dual_index=True, inferred_adapters=True)
READ_LENGTHS = {"Index 1 (i7)": 8, "Index 2 (i5)": 8}
SHOW_READ_CONTENT = True
FINAL_CAPTION = ("INFERRED adapter bases: the published KAPA HyperPlus plus IDT UDI route "
                 "fixes a standard dual-index Illumina topology but does not print its "
                 "oligos. The genomic insert was already shortened by PTA, so HyperPlus "
                 "fragmentation is omitted.")
SEQUENCING_INTRO = ("The study used NovaSeq 6000 but did not report genomic read lengths "
                    "or exact UDI oligos. Canonical primer sites follow the stated KAPA "
                    "HyperPlus/IDT dual-index library geometry.")


def oligos():
    return ()


def _initial_rows():
    sc = Scene.duplex([seg("single-cell genomic DNA", "X" * 48, placeholder=True)],
                      label="single-cell genome")
    return sc.rows()


def workflow():
    wf = Workflow(MolecularState("Denatured single-cell genome", tuple(_initial_rows())))
    wf.react("Anneal the undisclosed exonuclease-resistant random primer",
             pta.priming_rows(), name="Random-primed primary template",
             note="The paper gives 500 µM primer and 15 minutes at 65 °C, but does not disclose its bases, length or modification pattern.")
    wf.react("Extend with phi29 and strand displacement at 30 °C for 8 hours",
             pta.product_rows(), name="Terminated primary-template copies",
             note="Equal alpha-thio-ddATP, -ddCTP, -ddGTP and -ddTTP stop synthesis irreversibly; shorter displaced products are poor secondary templates, producing quasi-linear amplification.")
    body = seg("extendable PTA product", "X" * 48, placeholder=True)
    prepared = repair_and_dA_tail(Construct([body], name="PTA product"))
    wf.react("End-repair and dA-tail the purified PTA products without fragmentation",
             dA_tailed_scene(prepared, label="ligation-ready PTA product").rows(),
             name="End-repaired PTA product",
             note="The defining KAPA HyperPlus workflow explicitly omits fragmentation because PTA already creates shorter products.")
    sc = Scene.duplex(list(FINAL_LIBRARY), label="INFERRED adapter-completed library")
    sc.labels("top")
    wf.ligate_adapter("INFERRED — ligate IDT UDI adapters and amplify 15 PCR cycles",
             AdapterLigationGeometry.unavailable("The defining PTA paper identifies IDT UDI adapters used with KAPA HyperPlus but does not print their oligos or state their ligation-ready end geometry."), sc.rows(), name="Dual-index Illumina library",
             note="Exact adapter and primer bases are unavailable. Gel extraction retains 400–700-bp products for NovaSeq 6000 sequencing.")
    return wf


def final_rows():
    sc = Scene.duplex(list(FINAL_LIBRARY), label="PTA library")
    sc.labels("top")
    return sc.rows()


def protocol_graph():
    source="Defining source: https://doi.org/10.1073/pnas.2024176118"
    graph=ProtocolGraph("PTA",authority=source)
    genome=graph.input("genome","denatured single-cell genomic DNA",
                       MaterialKind.BIOLOGICAL,primary=True)
    primer=graph.input("random_primer","exonuclease-resistant random primer",
                       MaterialKind.OLIGO,inferred=True,
                       source=source+" (bases, length and modifications not disclosed)")
    phi29=graph.input("phi29","phi29 DNA polymerase",MaterialKind.CATALYST)
    terminators=graph.input("terminators","alpha-thio-ddATP/ddCTP/ddGTP/ddTTP",
                            MaterialKind.OLIGO)
    amplified_state=pta.amplification_product()
    amplified=graph.product("pta_product","terminated primary-template copies",
                            MaterialKind.MOLECULAR,state=amplified_state)
    graph.react("pta_amplification","Random-prime and displace primary-template copies",
                StepKind.EXTENSION,
                ((genome,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),
                 (terminators,InputRole.STRUCTURAL),(phi29,InputRole.CATALYST)),
                (amplified,))
    repair=graph.input("end_repair_mix","KAPA HyperPlus end-repair/dA-tailing mix",
                       MaterialKind.CATALYST)
    prepared_model=repair_and_dA_tail(Construct([
        seg("PTA-amplified genomic insert","X"*48,placeholder=True)],name="PTA product"))
    prepared_scene=dA_tailed_scene(prepared_model,label="ligation-ready PTA product")
    prepared=graph.product("prepared_product","end-repaired, dA-tailed PTA products",
                           MaterialKind.MOLECULAR,state=prepared_scene)
    graph.react("end_prepare","End-repair and dA-tail PTA products without fragmentation",
                StepKind.END_REPAIR,
                ((amplified,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(prepared,))
    adapters=graph.input("udi_adapters","IDT UDI adapters",MaterialKind.ADAPTER,
                         inferred=True,source=source+" (oligo sequences unavailable)")
    ligase=graph.input("ligase","KAPA HyperPlus DNA ligase",MaterialKind.CATALYST)
    ligated=graph.product("ligated_library","UDI-adapter-ligated PTA library",
                          MaterialKind.MOLECULAR,state=FINAL_LIBRARY,inferred=True,
                          source=source+" (adapter bases inferred from standard topology)")
    adapter_geometry=AdapterLigationGeometry.unavailable(
        "The paper identifies IDT UDI adapters but does not print their oligos or ligation-ready end geometry.")
    graph.react("adapter_ligation","Ligate IDT UDI adapters",StepKind.LIGATION,
                ((prepared,InputRole.TEMPLATE),(adapters,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(ligated,),geometry=adapter_geometry)
    forward=graph.input("pcr_forward","UDI library PCR forward primer",MaterialKind.OLIGO,
                        inferred=True,source=source+" (sequence unavailable)")
    reverse=graph.input("pcr_reverse","UDI library PCR reverse primer",MaterialKind.OLIGO,
                        inferred=True,source=source+" (sequence unavailable)")
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY,inferred=True,
                          source=source+" (canonical adapter bases inferred)")
    graph.react("library_pcr","Amplify the UDI library for 15 cycles",StepKind.PCR,
                ((ligated,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The defining paper does not print the IDT UDI PCR primers."))
    graph.sequence(library,SequencingPlan(
        library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),
        acquired_roles=("Read 1","Index 1 (i7)","Index 2 (i5)","Read 2"),
        read_lengths=READ_LENGTHS,
        read_uses={"Read 1":"include","Index 1 (i7)":"demultiplex",
                   "Index 2 (i5)":"demultiplex","Read 2":"include"},
        limitation="Exact UDI oligos and genomic read lengths are not reported."))
    graph.declare_reviewed_scope((
        "genome","random_primer","phi29","terminators","pta_amplification",
        "end_repair_mix","end_prepare","udi_adapters","ligase","adapter_ligation",
        "pcr_forward","pcr_reverse","pcr_polymerase","library_pcr","library"),source)
    return graph.require_publishable()
