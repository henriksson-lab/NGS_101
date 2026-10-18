"""Molecular construct model for scDNase-seq (Pico-Seq).

The paper's cited method identifies Illumina's Multiplexing Sample Prep Oligonucleotide
Kit (ref. 1005709); Illumina publishes that obsolete single-index oligo set.
"""
from __future__ import annotations

from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import closed_circle, role
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature
from endprep import EndPreparedDuplex, dA_tailed_scene

INDEX_NT = 6


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def digested_fragment_scene() -> Scene:
    """Representative short DNase-I-released genomic fragment."""
    return Scene.duplex([seg("DNase-accessible genomic fragment", "X" * 36,
                             placeholder=True)], label="released DNA")


def a_tailed_fragment_scene() -> Scene:
    """End-repaired fragment with one unpaired 3'-A at each end."""
    core = Construct([seg("end-repaired fragment", "X" * 32, placeholder=True)])
    return dA_tailed_scene(EndPreparedDuplex(core,five_prime_phosphate=True,three_prime_overhang="A",overhang_inferred=True), label="fragment")


def adapter_bottom() -> list[Segment]:
    return [
        seg("Read 1 arm", il.TRUSEQ_READ1[:-(len(il.STEM_COMPLEMENT) + 1)], "r1"),
        seg("stem complement", il.STEM_COMPLEMENT, "r1"),
        seg("3' T", "T", "r1"),
    ]


def adapter_top() -> list[Segment]:
    return [
        seg("stem", il.STEM, "r2"),
        seg("unpaired Read 2 arm", il.INDEX1_PRIMER[len(il.STEM):20], "r2"),
    ]


def adapter_scene() -> Scene:
    """Historical Multiplexing-kit fork; Scene enforces the 12-bp stem."""
    sc = Scene()
    sc.strand("bottom", adapter_bottom(), label="3'-T strand")
    sc.anneal("top", adapter_top(), to="bottom", pair=("stem", "stem complement"),
              label="5'-phosphorylated strand", mod5="p",
              unpaired=("unpaired Read 2 arm",))
    sc.mark("bottom", "3' T", "ligation overhang")
    return sc


def p5_primer() -> list[Segment]:
    return [seg("P5", il.P5, "p5"),
            seg("Read 1 continuation", il.TRUSEQ_READ1[4:], "r1")]


def indexed_p7_primer() -> list[Segment]:
    return [seg("P7", il.P7, "p7"),
            seg("i7", "N" * INDEX_NT, "cbc", placeholder=True,
                feature=feature("cell_i7","cell_barcode","whitelist",whitelist="published indexed primer set")),
            seg("indexed primer 3-prime end", il.TRUSEQ_READ2[:13], "r2")]


SEQ_PRIMERS = (
    sp.custom("Read 1", "Multiplexing-kit Read 1", il.TRUSEQ_READ1,
              'Illumina "Multiplexing Sample Preparation Oligonucleotide Kit" ref. 1005709'),
    sp.custom("Index 1 (i7)", "Multiplexing-kit index primer", il.INDEX1_PRIMER,
              'Illumina "Multiplexing Sample Preparation Oligonucleotide Kit" ref. 1005709'),
    sp.custom("Read 2", "Multiplexing-kit Read 2", il.TRUSEQ_READ2,
              'Illumina "Multiplexing Sample Preparation Oligonucleotide Kit" ref. 1005709'),
)


def final_library(insert_nt: int = 36) -> Construct:
    """PCR-completed single-index library; required primer sites are validated here."""
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("Read 1 arm", il.TRUSEQ_READ1[4:], "r1"),
        seg("DNase-accessible genomic insert", "X" * insert_nt, placeholder=True),
        seg("dA junction", "A",inferred=True),
        seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        seg("i7 reverse complement", "N" * INDEX_NT, "cbc", placeholder=True,
            feature=feature("cell_i7", "cell_barcode", "whitelist", whitelist="published indexed primer set")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="scDNase-seq library")
    problems = sp.verify(lib, SEQ_PRIMERS,
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems:
        raise ValueError("invalid scDNase-seq library: " + "; ".join(problems))
    return lib


def primer_landings():
    lib = final_library()
    return [(primer, sp.locate(lib, primer)) for primer in SEQ_PRIMERS]


def carrier_mixture_scene(target_label, carrier_label):
    """Independent target fragments and carrier/background, never one strand."""
    scene = Scene()
    scene.strand("target", [role(target_label)], label=target_label)
    scene.strand("carrier", [role(carrier_label)],
                 label="independent carrier/background")
    scene.note("carrier", "Added carrier is circular; its later background products remain unresolved.")
    scene.labels("target")
    scene.labels("carrier")
    return scene


def protocol_graph():
    source="Jin et al. 2015, doi:10.1038/nature15740, Methods DNase I digestion/Pico-Seq; kit oligos from Illumina Multiplexing reference1005709"
    g=ProtocolGraph("scDNase-seq",authority=source)
    cell=g.input("cell","one cell lysed in Triton/MgCl2 buffer",MaterialKind.BIOLOGICAL,primary=True)
    enzyme=g.input("nuclease","DNase I Roche04716728001,0.2–1U,37°C5min",MaterialKind.CATALYST)
    digested=g.product("digested_fragments","DNase-I-released genomic fragments",MaterialKind.MOLECULAR,state=digested_fragment_scene())
    g.react("digest","Digest one lysed cell with DNase I",StepKind.DIGESTION,
        ((cell,InputRole.TEMPLATE),(enzyme,InputRole.CATALYST)),(digested,))
    carrier=g.input("carrier","30ng circular plasmid carrier (identity unassigned)",MaterialKind.MOLECULAR,
        state=closed_circle("circular carrier DNA",(role("plasmid sequence unassigned"),),closure_label="covalently closed carrier"))
    pk=g.input("proteinase_k","1µl20mg/ml proteinase K",MaterialKind.CATALYST)
    stop=g.input("stop_buffer","80µl Stop Buffer,0.15%SDS/10mM EDTA;55°C1h",MaterialKind.CONDITION)
    stopped=g.product("stopped_dna","stopped digest with carrier DNA",MaterialKind.MOLECULAR,
        state=carrier_mixture_scene("DNase-released genomic fragments", "30ng circular carrier"))
    g.react("stop_digest","Stop DNase and digest protein; add circular carrier",StepKind.DEPROTEINIZATION,
        ((digested,InputRole.TEMPLATE),(carrier,InputRole.STRUCTURAL),(pk,InputRole.CATALYST),
         (stop,InputRole.CONDITION)),(stopped,))
    extraction=g.input("dna_extraction","phenol–chloroform;ethanol precipitation with20µg glycogen",MaterialKind.CONDITION)
    recovered=g.product("recovered_dna","purified genomic/carrier DNA mixture",MaterialKind.MOLECULAR,state=stopped.state)
    g.react("recover_dna","Recover DNA by phenol–chloroform and ethanol precipitation",StepKind.EXTRACTION,
        ((stopped,InputRole.TEMPLATE),(extraction,InputRole.CONDITION)),(recovered,),changes_molecule=False)
    repair=g.input("end_repair_mix","cited Illumina library end-preparation enzymes",MaterialKind.CATALYST)
    tailed=g.product("a_tailed_fragments","end-repaired dA-tailed genomic fragments (representative target subset)",
        MaterialKind.MOLECULAR,state=a_tailed_fragment_scene(),inferred=True,
        source=source+"; dA-tail inferred from the cited adapter compatibility")
    g.react("end_prepare","INFERRED — end-prepare genomic fragments for the cited Illumina adapter",StepKind.END_REPAIR,
        ((recovered,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(tailed,))
    adapters=g.input("adapters","obsolete Multiplexing kit adapter strands",MaterialKind.ADAPTER,state=adapter_scene())
    ligase=g.input("ligase","cited Illumina library DNA ligase",MaterialKind.CATALYST)
    ligated=g.product("ligated_fragments","adapter-ligated target fragments; carrier background unresolved",MaterialKind.MOLECULAR,
        state=carrier_mixture_scene("fork-adapted genomic fragments", "carrier background; products unresolved"))
    g.react("adapter_ligation","Ligate the Multiplexing adapter",StepKind.LIGATION,
        ((tailed,InputRole.TEMPLATE),(adapters,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),
        (ligated,),geometry=reported_adapter_end("Multiplexing adapter",kind="cohesive",polarity="3-prime",
            terminal_sequence="T",substrate_scenes=(a_tailed_fragment_scene(),),
            statement="Vendor obsolete Multiplexing oligos have a3′T strand for the source-qualified end-preparation reference."))
    p5=g.input("p5_primer","Multiplexing PCR Primer1.0",MaterialKind.OLIGO,state=tuple(p5_primer()))
    p7=g.input("indexed_p7_primer","Multiplexing indexed primer architecture (six-base index)",MaterialKind.OLIGO,state=tuple(indexed_p7_primer()))
    polymerase=g.input("pcr_polymerase","library PCR polymerase (formulation unassigned)",MaterialKind.CATALYST)
    amplified=g.product("indexed_amplicons","indexed target amplicons with unresolved carrier-derived background",MaterialKind.MOLECULAR,state=final_library())
    g.react("indexed_pcr","First index PCR:6 cycles",StepKind.PCR,
        ((ligated,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),(p7,InputRole.STRUCTURAL),
         (polymerase,InputRole.CATALYST)),(amplified,),
        geometry=PCRGeometry.unavailable(ligated.state,
            reason="The vendor index primer ends with13nt GTGACTGGAGTTC, not the full Read2 primer. The20nt adapter/Primer2.0 extension and exact assay primer mixture are not reconstructed; the completed library cannot stand in for this precursor."))
    selector=g.input("size_selection","2%E-gel,160–300bp",MaterialKind.SELECTION)
    selected=g.product("selected_amplicons","160–300bp first-PCR fraction",MaterialKind.MOLECULAR,state=amplified.state)
    rejected=g.product("rejected_sizes","first E-gel off-size DNA",MaterialKind.WASTE)
    g.react("size_select","Select160–300bp fragments on2%E-gel",StepKind.SELECTION,
        ((amplified,InputRole.TEMPLATE),(selector,InputRole.SELECTION)),(selected,rejected))
    g.discard(rejected,"E-gel fraction outside the160–300bp selected window is not recovered")
    final_p5=g.input("final_p5_primer","second-round P5 primer (exact physical sequence unassigned)",MaterialKind.OLIGO)
    final_p7=g.input("final_p7_primer","second-round P7 primer (exact physical sequence unassigned)",MaterialKind.OLIGO)
    final_product=g.product("final_amplicons","second-PCR amplified fragments",MaterialKind.MOLECULAR,state=final_library())
    g.react("final_pcr","Second PCR with P5/P7:22 cycles",StepKind.PCR,
        ((selected,InputRole.TEMPLATE),(final_p5,InputRole.STRUCTURAL),(final_p7,InputRole.STRUCTURAL),
         (polymerase,InputRole.CATALYST)),(final_product,),geometry=PCRGeometry.unavailable(selected.state,
            reason="Methods name P5 and P7 primers without assigning their exact physical sequences in this model. These are distinct from the first index-PCR kit-primer declarations."))
    final_selector=g.input("final_size_selection","second160–300bp E-gel selection",MaterialKind.SELECTION)
    library=g.product("library","twice E-gel-selected scDNase library",MaterialKind.LIBRARY,state=final_product.state)
    final_reject=g.product("final_rejected_sizes","second E-gel off-size DNA",MaterialKind.WASTE)
    g.react("final_size_select","Select160–300bp after second PCR",StepKind.SELECTION,
        ((final_product,InputRole.TEMPLATE),(final_selector,InputRole.SELECTION)),(library,final_reject))
    g.discard(final_reject,"second E-gel off-size fractions are not recovered")
    g.sequence(library,SequencingPlan(library=final_library(),primers=tuple(SEQ_PRIMERS),
        acquisition="unspecified",limitation="Methods name HiSeq2500; actual acquired read/index roles and cycle mask remain unassigned. Kit primer sites do not establish acquisition."))
    g.declare_reviewed_scope(("cell","nuclease","digest","carrier","proteinase_k","stop_buffer","stop_digest",
        "dna_extraction","recover_dna","end_repair_mix","end_prepare","adapters","ligase","adapter_ligation",
        "p5_primer","indexed_p7_primer","pcr_polymerase","indexed_pcr",
        "size_selection","size_select","final_p5_primer","final_p7_primer","final_pcr","final_size_selection","final_size_select","library"),source)
    return g.require_publishable()
