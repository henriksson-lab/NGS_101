"""Molecular construct model for scMNase-seq.

The inspected supplement specifies universal adapters through a cited method.
Their exact study assignment remains unchecked; adapter segments are inferred references.
"""
from __future__ import annotations

from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature
from endprep import EndPreparedDuplex, dA_tailed_scene

PCR_CYCLES = 23
SIZE_SELECTION_BP = (160, 300)


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def inferred(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    """An unverified adapter reference region, always dotted."""
    return Segment(name=name, top=top, tag=tag, inferred=True, **kw)


def digested_fragment_scene() -> Scene:
    """Representative nucleosome-protected fragment left by MNase."""
    return Scene.duplex([seg("MNase-protected genomic fragment", "X" * 44,
                             placeholder=True)], label="protected DNA")


def a_tailed_fragment_scene() -> Scene:
    """End-repaired fragment with one unpaired 3'-A at each end."""
    core = Construct([inferred("end-repaired fragment", "X" * 40, placeholder=True)])
    return dA_tailed_scene(EndPreparedDuplex(core,five_prime_phosphate=True,three_prime_overhang="A",overhang_inferred=True), label="fragment")


def adapter_bottom() -> list[Segment]:
    return [
        inferred("Read 1 arm", il.TRUSEQ_READ1[:-(len(il.STEM_COMPLEMENT) + 1)], "r1"),
        inferred("stem complement", il.STEM_COMPLEMENT, "r1"),
        inferred("3' T", "T", "r1"),
    ]


def adapter_top() -> list[Segment]:
    return [
        inferred("stem", il.STEM, "r2"),
        inferred("unpaired Read 2 arm", il.INDEX1_PRIMER[len(il.STEM):20], "r2"),
    ]


def adapter_scene() -> Scene:
    """Inferred TruSeq fork; Scene enforces the 12-bp stem."""
    sc = Scene()
    sc.strand("bottom", adapter_bottom(), label="3'-T strand")
    sc.anneal("top", adapter_top(), to="bottom", pair=("stem", "stem complement"),
              label="5'-phosphorylated strand", mod5="p",
              unpaired=("unpaired Read 2 arm",))
    sc.mark("bottom", "3' T", "ligation overhang")
    return sc


def p5_primer() -> list[Segment]:
    return [inferred("P5", il.P5, "p5"),
            inferred("Read 1 continuation", il.TRUSEQ_READ1[4:], "r1")]


def indexed_p7_primer() -> list[Segment]:
    return [inferred("P7", il.P7, "p7"),
            role("unassigned cell index",tag="cbc",molecular_feature=feature("cell_i7","cell_barcode","unknown"),inferred=True),
            inferred("Read 2", il.TRUSEQ_READ2, "r2")]


SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])


def final_library(insert_nt: int = 44) -> Construct:
    """PCR-completed single-index library; required primer sites are validated here."""
    lib = Construct([
        inferred("P5", il.P5, "p5"),
        inferred("Read 1 arm", il.TRUSEQ_READ1[4:], "r1"),
        seg("MNase-protected genomic insert", "X" * insert_nt, placeholder=True),
        inferred("dA junction", "A"),
        inferred("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        role("i7 reverse complement",tag="cbc",molecular_feature=feature("cell_i7","cell_barcode","unknown"),inferred=True),
        inferred("P7 reverse complement", il.P7_RC, "p7"),
    ], name="scMNase-seq library")
    for primer in SEQ_PRIMERS:
        if sp.locate(lib, primer) is None:
            raise ValueError(f"{primer.name} does not land on {lib.name}")
    return lib


def protocol_graph():
    source="Lai et al.2018,doi:10.1038/s41586-018-0567-3,Supplementary Methods pp2–3; cited universal-adapter method33 remains unassigned"
    inferred_source=source+"; conditional TruSeq reference, not verified study reagents"
    g=ProtocolGraph("scMNase-seq",authority=source)
    cell=g.input("cell","one cell lysed in Triton/CaCl2 buffer",MaterialKind.BIOLOGICAL,primary=True)
    enzyme=g.input("nuclease","MNase1:3000 in lysis buffer,37°C5min",MaterialKind.CATALYST)
    digested=g.product("digested_fragments","MNase-protected genomic fragments",MaterialKind.MOLECULAR,state=digested_fragment_scene())
    g.react("digest","Digest one lysed cell with micrococcal nuclease",StepKind.DIGESTION,
        ((cell,InputRole.TEMPLATE),(enzyme,InputRole.CATALYST)),(digested,))
    pk=g.input("proteinase_k","1µl20mg/ml proteinase K",MaterialKind.CATALYST)
    stop=g.input("stop_buffer","80µl0.15%SDS/10mM EGTA Stop Buffer;55°C1h",MaterialKind.CONDITION)
    stopped=g.product("stopped_dna","protein-depleted MNase digest",MaterialKind.MOLECULAR,state=digested.state)
    g.react("stop_digest","Stop MNase with EGTA/SDS and digest proteins",StepKind.DEPROTEINIZATION,
        ((digested,InputRole.TEMPLATE),(pk,InputRole.CATALYST),(stop,InputRole.CONDITION)),(stopped,))
    extraction=g.input("dna_extraction","phenol–chloroform;ethanol precipitation with glycogen",MaterialKind.CONDITION)
    recovered=g.product("recovered_dna","purified MNase-protected DNA",MaterialKind.MOLECULAR,state=stopped.state)
    g.react("recover_dna","Recover DNA before universal-adapter ligation",StepKind.EXTRACTION,
        ((stopped,InputRole.TEMPLATE),(extraction,InputRole.CONDITION)),(recovered,),changes_molecule=False)
    repair=g.input("end_repair_mix","unassigned end-preparation enzymes (conditional reference)",MaterialKind.CATALYST,inferred=True,source=inferred_source)
    tailed=g.product("a_tailed_fragments","conditional dA-tailed reference",MaterialKind.MOLECULAR,state=a_tailed_fragment_scene(),inferred=True,source=inferred_source)
    g.react("end_prepare","INFERRED — end-repair and dA-tail reference",StepKind.END_REPAIR,
        ((recovered,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(tailed,),source=inferred_source)
    adapters=g.input("adapters","universal adapters (TruSeq reference only)",MaterialKind.ADAPTER,state=adapter_scene(),inferred=True,source=inferred_source)
    ligase=g.input("ligase","DNA ligase for cited universal-adapter method",MaterialKind.CATALYST)
    ligated=g.product("ligated_fragments","universal-adapter-ligated protected DNA",MaterialKind.MOLECULAR,
        state=molecule("universal-adapter PCR input",(role("adapter-ligated protected DNA"),)))
    g.react("adapter_ligation","Ligate universal adapters; physical ends unresolved",StepKind.LIGATION,
        ((tailed,InputRole.TEMPLATE),(adapters,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(ligated,),
        geometry=AdapterLigationGeometry.unavailable("The defining supplement names universal adapters via method33; the inferred TruSeq reference does not establish the supplied adapters or their physical end geometry."))
    p5=g.input("p5_primer","conditional TruSeq P5-side primer reference",MaterialKind.OLIGO,state=tuple(p5_primer()),inferred=True,source=inferred_source)
    p7=g.input("indexed_p7_primer","conditional indexed P7 reference; index length unassigned",MaterialKind.OLIGO,state=tuple(indexed_p7_primer()),inferred=True,source=inferred_source)
    polymerase=g.input("pcr_polymerase","index PCR polymerase (formulation unassigned)",MaterialKind.CATALYST)
    amplified=g.product("indexed_amplicons","indexed amplified protected fragments (conditional arms)",MaterialKind.MOLECULAR,state=final_library(),inferred=True,source=inferred_source)
    g.react("indexed_pcr","Index-amplify adapter-ligated fragments:23 cycles",StepKind.PCR,
        ((ligated,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),(p7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(amplified,),
        geometry=PCRGeometry.unavailable(ligated.state,
            reason="The study-specific universal adapters/indexed primers and their precursor sites are unassigned. Conditional final TruSeq sites are not evidence for the preceding PCR template."))
    selector=g.input("size_selection","E-gel160–300bp fraction",MaterialKind.SELECTION)
    library=g.product("library","size-selected scMNase library (conditional arms)",MaterialKind.LIBRARY,state=amplified.state,inferred=True,source=inferred_source)
    rejected=g.product("rejected_sizes","E-gel off-size DNA",MaterialKind.WASTE)
    g.react("size_select","Select160–300bp library fragments by E-gel",StepKind.SELECTION,
        ((amplified,InputRole.TEMPLATE),(selector,InputRole.SELECTION)),(library,rejected))
    g.discard(rejected,"off-size E-gel fractions are not recovered")
    g.sequence(library,SequencingPlan(library=final_library(),acquired_roles=("Read 1","Read 2"),
        limitation="Supplementary Methods specify paired-end HiSeq2500. Run lengths/index acquisition and study adapters remain unassigned; exact reference sites do not establish an experimental run mask."))
    g.declare_reviewed_scope(("cell","nuclease","digest","proteinase_k","stop_buffer","stop_digest","dna_extraction","recover_dna",
        "end_repair_mix","end_prepare","adapters","ligase","adapter_ligation","p5_primer","indexed_p7_primer","pcr_polymerase","indexed_pcr","size_selection","size_select","library"),source)
    return g.require_publishable()
