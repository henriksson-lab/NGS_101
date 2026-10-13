"""Liu et al. 2019 TAPS chemistry."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from base_conversion import taps_path
from chemdraw import Construct, MolecularState, Row, Scene, Segment, Workflow, feature
from endprep import (ThreePrimeOverhangAdapter, cohesive_ligation_geometry,
                     dA_tailed_scene, ligation_product_scene, repair_and_dA_tail)
from pcr import PCRGeometry
from protocol_graph import (Disposition, Input, InputRole, Material, MaterialKind,
                            ProtocolGraph, SequencingPlan, SourceCoverage, SourceItem,
                            Step, StepKind, Terminal)

TITLE = "TAPS — bisulfite-free methylation sequencing"
NOTES = "01_taps.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41587-019-0041-2">Liu et al., <i>Nature Biotechnology</i> (2019)</a>.'
SUMMARY = "TET oxidation and pyridine–borane reduction make 5mC and 5hmC read as T, while unmodified C remains C."
CAVEAT = ""
INDEX6 = "GCCAAT"  # TruSeq Index 6, printed in the defining supplement
TRUSEQ_ADAPTER = ThreePrimeOverhangAdapter(
    "TAPS TruSeq adapter ligating stem", "GCTCTTCCGATCT", "GATCGGAAGAGC")


def library():
    lib = Construct([
        Segment("P5 before shared ACAC", il.P5[:-4], "p5"),
        Segment("Read 1 site", il.TRUSEQ_READ1, "r1"),
        Segment("TAPS-treated genomic insert", "X" * 36, placeholder=True),
        Segment("dA junction / Read 2 site start", "A"),
        Segment("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        Segment("i7 reverse complement", INDEX6, "cbc",
                feature=feature("sample_i7", "sample_index", "whitelist", whitelist="TruSeq index set")),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="TAPS TruSeq library")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems: raise ValueError("invalid TAPS library: " + "; ".join(problems))
    return lib, primers


FINAL_LIBRARY, SEQ_PRIMERS = library()
FINAL_CAPTION = "Single-index TruSeq library after TAPS treatment and PCR."
SEQUENCING_INTRO = "Paired reads interrogate the treated insert; i7 identifies the library."


def workflow():
    insert = Construct([Segment("TAPS genomic insert", "X" * 36, placeholder=True)])
    tailed = repair_and_dA_tail(insert)
    w = Workflow(MolecularState("fragmented genomic DNA", tuple(Scene.duplex(list(insert)).rows())))
    w.react("End-repair and add 3′ dA", dA_tailed_scene(tailed).rows(), name="dA-tailed DNA",
            note="The single-base insert overhang is the substrate for the printed TruSeq adapter.")
    geometry = cohesive_ligation_geometry(tailed, TRUSEQ_ADAPTER,
        statement="The TAPS supplement prints the complete TruSeq construct; its adapter stem leaves a 3′-T complementary to insert 3′-dA.")
    lig = ligation_product_scene(FINAL_LIBRARY, "TAPS-treated genomic insert", label="adapter-ligated molecule")
    w.ligate_adapter("Ligate TruSeq adapters", geometry, lig.rows(), name="adapter-ligated TAPS molecule",
                     note="The supplement prints the complete single-index construct.")
    w.react("Oxidize and reduce modified cytosines", [
            Row(chunks=[(taps_path("5mC").text(), "w1", False)]),
            Row(chunks=[(taps_path("5hmC").text(), "w1", False)]),
            Row(chunks=[(taps_path("C").text(), None, False)]),
        ], name="TAPS-converted library", note="TET oxidation makes 5caC; pyridine–borane makes DHU, which PCR copies as T.")
    w.react("PCR", [Row(chunks=[("DHU → T in amplified library; C remains C", None, False)])],
            name="amplified TAPS library", note="The conversion polarity is opposite to bisulfite sequencing.")
    return w


def protocol_graph() -> ProtocolGraph:
    source = "Liu et al. 2019, doi:10.1038/s41587-019-0041-2"
    graph = ProtocolGraph(TITLE, authority=source)
    insert = Construct([Segment("TAPS genomic insert", "X" * 36, placeholder=True)])
    tailed = repair_and_dA_tail(insert)
    graph.add_input(Material("genomic_dna", "fragmented genomic DNA",
                             MaterialKind.BIOLOGICAL,
                             state=Scene.duplex(list(insert)), primary=True, source=source))
    graph.add_input(Material("end_repair_enzymes", "end-repair and dA-tailing enzymes",
                             MaterialKind.CATALYST, source=source))
    graph.add_input(Material("datp", "dATP", MaterialKind.CONDITION, source=source))
    tailed_material = Material("da_tailed_dna", "end-repaired, 3-prime-dA-tailed DNA",
                               MaterialKind.MOLECULAR, state=dA_tailed_scene(tailed),
                               source=source)
    graph.add_step(Step(
        "end_repair_and_tailing", "End-repair and add 3-prime dA", StepKind.TAILING,
        (Input("genomic_dna", InputRole.TEMPLATE),
         Input("end_repair_enzymes", InputRole.CATALYST),
         Input("datp", InputRole.CONDITION)),
        (tailed_material.id,), source=source), (tailed_material,))

    geometry = cohesive_ligation_geometry(
        tailed, TRUSEQ_ADAPTER,
        statement="The TAPS supplement prints the complete TruSeq construct; its adapter stem leaves a 3-prime-T complementary to insert 3-prime-dA.")
    graph.add_input(Material("truseq_adapter", "3-prime-dT TruSeq adapter",
                             MaterialKind.ADAPTER, state=TRUSEQ_ADAPTER.scene(), source=source))
    graph.add_input(Material("dna_ligase", "DNA ligase", MaterialKind.CATALYST,
                             source=source))
    ligated_scene = ligation_product_scene(
        FINAL_LIBRARY, "TAPS-treated genomic insert", label="adapter-ligated molecule")
    ligated = Material("ligated_library", "adapter-ligated TAPS molecule",
                       MaterialKind.MOLECULAR, state=ligated_scene, source=source)
    graph.add_step(Step(
        "adapter_ligation", "Ligate TruSeq adapters", StepKind.LIGATION,
        (Input(tailed_material.id, InputRole.TEMPLATE),
         Input("truseq_adapter", InputRole.STRUCTURAL),
         Input("dna_ligase", InputRole.CATALYST)),
        (ligated.id,), geometry=geometry, source=source), (ligated,))

    graph.add_input(Material("tet", "TET oxidation enzymes", MaterialKind.CATALYST,
                             source=source))
    graph.add_input(Material("pyridine_borane", "pyridine-borane",
                             MaterialKind.CONDITION, source=source))
    converted = Material(
        "converted_library", "TAPS-converted adapter-ligated library",
        MaterialKind.MOLECULAR,
        state=(Row(chunks=[(taps_path("5mC").text(), "w1", False)]),
               Row(chunks=[(taps_path("5hmC").text(), "w1", False)]),
               Row(chunks=[(taps_path("C").text(), None, False)])), source=source)
    graph.add_step(Step(
        "convert_modified_c", "Oxidize and reduce modified cytosines",
        StepKind.CHEMICAL_CONVERSION,
        (Input(ligated.id, InputRole.TEMPLATE), Input("tet", InputRole.CATALYST),
         Input("pyridine_borane", InputRole.CONDITION)),
        (converted.id,), source=source), (converted,))

    graph.add_input(Material("library_forward_primer", "library PCR forward primer",
                             MaterialKind.OLIGO, source=source))
    graph.add_input(Material("library_reverse_primer", "indexed library PCR reverse primer",
                             MaterialKind.OLIGO, source=source))
    pcr_geometry = PCRGeometry.unavailable(
        FINAL_LIBRARY,
        reason="The current source transcription establishes the completed TruSeq construct but does not independently transcribe both library-PCR primer sequences.")
    final = Material("taps_library", "amplified TAPS sequencing library",
                     MaterialKind.LIBRARY, state=FINAL_LIBRARY, source=source)
    graph.add_step(Step(
        "library_pcr", "Amplify the converted library", StepKind.PCR,
        (Input(converted.id, InputRole.TEMPLATE),
         Input("library_forward_primer", InputRole.STRUCTURAL),
         Input("library_reverse_primer", InputRole.STRUCTURAL)),
        (final.id,), geometry=pcr_geometry, source=source), (final,))
    graph.mark_terminal(Terminal(final.id, Disposition.SEQUENCING), SequencingPlan(
        FINAL_LIBRARY, tuple(SEQ_PRIMERS), tuple(p.role for p in SEQ_PRIMERS)))

    graph.declare_reviewed_scope(dict((
            ("genomic_dna_input", "genomic_dna"),
            ("end_repair_enzyme_reagents", "end_repair_enzymes"),
            ("datp_reagent", "datp"),
            ("end_repair_and_tailing", "end_repair_and_tailing"),
            ("truseq_adapter_reagent", "truseq_adapter"),
            ("dna_ligase_reagent", "dna_ligase"),
            ("adapter_ligation", "adapter_ligation"),
            ("tet_reagent", "tet"),
            ("pyridine_borane_reagent", "pyridine_borane"),
            ("tet_oxidation", "convert_modified_c"),
            ("pyridine_borane_reduction", "convert_modified_c"),
            ("forward_pcr_primer_reagent", "library_forward_primer"),
            ("reverse_pcr_primer_reagent", "library_reverse_primer"),
            ("library_pcr", "library_pcr"),
            ("sequencing", "taps_library"))), source)
    return graph.require_publishable()
