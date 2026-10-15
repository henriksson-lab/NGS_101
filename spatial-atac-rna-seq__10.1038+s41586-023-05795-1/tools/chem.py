"""Protocol-owned spatial assay and reaction graph."""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
from protocol_families.spatial import (ident, SPATIAL_BEAD, UMI, CELL_A, CELL_B, HASH, SAMPLE,
    SECTION_HASH, SPACE_HASH, SPATIAL_A, SPATIAL_B, SPATIAL_A_ROUND, SPATIAL_B_ROUND,
    _finish_pcr, _spatial_dna_libraries, _spatial_rows, hdst_library, slide_seq_library,
    sci_space_library)
from adapter_ligation import AdapterLigationGeometry
from chemdraw import MolecularState, Scene, Workflow, duplex_rows
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import chain, molecule, role, role_library

from types import SimpleNamespace

TITLE = 'spatial-ATAC–RNA-seq'
NOTES = '01_spatial_atac_rna_seq.html'
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41586-023-05795-1">Zhang et al., Nature (2023)</a>.'
SUMMARY = 'Tn5 tags accessible DNA and in-situ poly(dT) RT captures RNA before orthogonal microfluidic barcode ligations assign both branches the same tissue pixel.'
CAVEAT = ''
FINAL_CAPTION = 'The exact representative panel is the ATAC branch: Read 1 enters accessible DNA, while Read 2 enters the B–A spatial chain from the N701 end. The separately recovered RNA branch retains the same spatial address plus a UMI and is Nextera-tagmented after pre-amplification.'
SEQUENCING_INTRO = 'The exact representative panel is the ATAC branch: Read 1 enters accessible DNA, while Read 2 enters the B–A spatial chain from the N701 end. The separately recovered RNA branch retains the same spatial address plus a UMI and is Nextera-tagmented after pre-amplification.'
READ_LENGTHS = None
READ_USES = None
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
oligos = lambda: ("<p>RT primer (5′-phosphate; one dT is internal biotin): 5'- CATCGGCGTACGACT<umi>NNNNNNNNNN</umi>TTTTTTTTTTTTTTTVN -3'</p>", "<p>Tn5ME-B: 5'- CATCGGCGTACGACTAGATGTGTATAAGAGACAG -3'</p>", "<p>Barcode A-1 (5′-phosphate): 5'- AGGCCAGAGCATTCG<cbc>AACGTGAT</cbc>GTGGCCGATGTTTCG -3'</p>", "<p>Ligation linker 1: 5'- AGTCGTACGCCGATGCGAAACATCGGCCAC -3'</p>", "<p>Barcode B-1: 5'- CAAGCGTTGGCTTCTCGCATCT<cbc>AACGTGAT</cbc>ATCCACGTGCTTGAG -3'</p>", "<p>Ligation linker 2: 5'- CGAATGCTCTGGCCTCTCAAGCACGTGGAT -3'</p>")
_template, FINAL_LIBRARY, _geometry, SEQ_PRIMERS = _spatial_dna_libraries("accessible DNA")
_definition = SimpleNamespace(title=TITLE, doi='10.1038/s41586-023-05795-1', citation='Zhang et al., Nature (2023)', endpoint='The exact representative panel is the ATAC branch: Read 1 enters accessible DNA, while Read 2 enters the B–A spatial chain from the N701 end. The separately recovered RNA branch retains the same spatial address plus a UMI and is Nextera-tagmented after pre-amplification.')


def _dbit(*, cuttag: bool):
    enzyme = "antibody-tethered pA–Tn5" if cuttag else "Tn5"
    dna_kind = "histone-mark DNA" if cuttag else "accessible DNA"
    start = Scene(); start.strand("targets", [role(dna_kind), role("poly(A) RNA")],
                                 label="fixed permeabilized tissue")
    wf = Workflow(MolecularState("Fixed tissue targets", tuple(start.rows())))
    wf.react(f"Tag {dna_kind} with {enzyme} and reverse-transcribe RNA",
             _spatial_rows(0, dna_kind), name="shared ligation handle on both branches",
             note="Tn5ME-B and the 5′-phosphorylated RT primer carry the same exact CATCGGCGTACGACT handle.")
    g1 = AdapterLigationGeometry.known(
        "splinted", SPATIAL_A_ROUND.adapter_scene(),
        statement="Linker 1 is exactly the reverse complement of the shared RT/Tn5ME-B handle followed by the barcode-A ligating arm.")
    wf.ligate_adapter("Ligate barcode A in the first microfluidic channel set", g1,
                      _spatial_rows(1, dna_kind), name="A-barcoded DNA and RNA branches",
                      note="The A-oligo 3′-OH is sealed to the shared 5′-phosphate; A-1 is shown from the 100-member plate.")
    g2 = AdapterLigationGeometry.known(
        "splinted", SPATIAL_B_ROUND.adapter_scene(),
        statement="Linker 2 is exactly the reverse complement of the barcode-A outgoing arm followed by the barcode-B ligating arm.")
    wf.ligate_adapter("Rotate the device and ligate barcode B", g2,
                      _spatial_rows(2, dna_kind), name="two-dimensional pixel address",
                      note="The B-oligo 3′-OH is sealed to the barcode-A 5′-phosphate; the ordered A–B pair identifies one channel intersection.")
    wf.react("Recover the biotin-selected RNA branch and amplify both libraries",
             _spatial_rows(2, dna_kind), name="separate matched DNA and RNA libraries",
             note="Both physical libraries retain the same exact B–A spatial chain; only the RNA branch carries the RT UMI and poly(dT).")
    _template, final, geometry, _primers = _spatial_dna_libraries(dna_kind)
    wf.pcr("Add the indexed Illumina shell to the DNA branch", geometry,
           duplex_rows(final, label="DNA library", unpaired=(f"{dna_kind} genomic DNA'",)),
           name="indexed DNA-branch library",
           note="N501 supplies P5/i5 at the Tn5ME-A end; N701 supplies P7/i7 at the spatial-barcode end. The independently recovered RNA branch is Nextera-tagmented after pre-amplification and is not collapsed into this fixed cut site.")
    return wf

workflow = lambda: _dbit(cuttag=False)

def _dbit_graph(d, final, primers, *, cuttag):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    dna_kind = "histone-mark DNA" if cuttag else "accessible DNA"
    tissue = graph.input("tissue", "fixed permeabilized tissue containing DNA and poly(A) RNA",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("fixed tissue targets", (role(dna_kind), role("poly(A) RNA"))),
                         primary=True)
    substrate = tissue
    reviewed = ["tissue"]
    if cuttag:
        antibody = graph.input("antibody", "histone-mark primary antibody",
                               MaterialKind.SELECTION)
        bound = graph.product("antibody_bound", "antibody-bound tissue chromatin",
                              MaterialKind.MOLECULAR, state=tissue.state)
        unbound = graph.product("unbound_chromatin", "unbound chromatin", MaterialKind.WASTE)
        graph.react("bind_antibody", "Bind antibody to the selected histone mark",
                    StepKind.BINDING,
                    ((tissue, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                    (bound, unbound))
        graph.discard(unbound, "does not retain the selected histone-mark complex")
        substrate = bound
        reviewed.extend(("antibody", "bind_antibody"))
    transposome = graph.input("transposome", "RT-handle-bearing loaded Tn5 adapter",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "antibody-tethered pA–Tn5" if cuttag else "Tn5 transposase",
                      MaterialKind.CATALYST)
    tagged = graph.product("tagged_tissue", "tissue with tagged DNA and retained poly(A) RNA",
                           MaterialKind.MOLECULAR,
                           state=molecule("shared ligation handle on both branches", (
                               role(dna_kind + " carrying the shared handle"), role("poly(A) RNA"))))
    waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tag " + dna_kind + " with the shared-handle transposome",
                StepKind.TAGMENTATION,
                ((substrate, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not retain a spatial-barcode-compatible genomic fragment")
    rt_primer = graph.input("rt_primer", "phosphorylated shared-handle oligo(dT) RT primer",
                            MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    paired = graph.product("paired_branches", "DNA and cDNA branches bearing the shared handle",
                           MaterialKind.MOLECULAR,
                           state=molecule("shared-handle DNA and RNA branches", (
                               role(dna_kind + " fragment"), UMI, role("transcript cDNA"))))
    graph.react("reverse_transcribe", "Reverse-transcribe RNA through the shared RT handle",
                StepKind.REVERSE_TRANSCRIPTION,
                ((tagged, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (paired,))
    a_adapter = graph.input("barcode_a", "barcode-A oligo pool and linker 1",
                            MaterialKind.ADAPTER)
    a_ligase = graph.input("a_ligase", "barcode-A DNA ligase", MaterialKind.CATALYST)
    a_product = graph.product("a_product", "A-barcoded DNA and cDNA branches",
                              MaterialKind.MOLECULAR, state=paired.state)
    g1 = AdapterLigationGeometry.known(
        "splinted", SPATIAL_A_ROUND.adapter_scene(),
        statement="Linker 1 is the reverse complement of the shared RT/Tn5ME-B handle followed by the barcode-A ligating arm.")
    graph.react("ligate_a", "Ligate barcode A in the first channel set", StepKind.LIGATION,
                ((paired, InputRole.TEMPLATE), (a_adapter, InputRole.STRUCTURAL),
                 (a_ligase, InputRole.CATALYST)), (a_product,), geometry=g1)
    b_adapter = graph.input("barcode_b", "barcode-B oligo pool and linker 2",
                            MaterialKind.ADAPTER)
    b_ligase = graph.input("b_ligase", "barcode-B DNA ligase", MaterialKind.CATALYST)
    b_product = graph.product("b_product", "two-dimensionally addressed DNA and cDNA branches",
                              MaterialKind.MOLECULAR, state=a_product.state)
    g2 = AdapterLigationGeometry.known(
        "splinted", SPATIAL_B_ROUND.adapter_scene(),
        statement="Linker 2 is the reverse complement of the barcode-A outgoing arm followed by the barcode-B ligating arm.")
    graph.react("ligate_b", "Rotate the device and ligate barcode B", StepKind.LIGATION,
                ((a_product, InputRole.TEMPLATE), (b_adapter, InputRole.STRUCTURAL),
                 (b_ligase, InputRole.CATALYST)), (b_product,), geometry=g2)
    template, _same_final, geometry, _same_primers = _spatial_dna_libraries(dna_kind)
    dna = graph.product("dna_branch", dna_kind + " spatial branch", MaterialKind.MOLECULAR,
                        state=template)
    rna = graph.product("rna_branch", "matched spatial RNA branch", MaterialKind.MOLECULAR,
                        state=molecule("spatial RNA branch", (SPATIAL_A, SPATIAL_B, UMI,
                            role("transcript cDNA"))))
    graph.react("split_branches", "Recover matched DNA and RNA branches",
                StepKind.SPLIT, ((b_product, InputRole.TEMPLATE),), (dna, rna),
                changes_molecule=False)
    reviewed.extend(("transposome", "tn5", "tagment", "rt_primer", "reverse_transcriptase",
                     "reverse_transcribe", "barcode_a", "a_ligase", "ligate_a",
                     "barcode_b", "b_ligase", "ligate_b", "split_branches"))
    _finish_pcr(graph, "dna", dna, final, primers, "DNA branch: " + d.endpoint,
                reviewed, label=d.title + " DNA", geometry=geometry)
    rna_transposome = graph.input("rna_transposome", "RNA-branch library transposome",
                                  MaterialKind.ADAPTER)
    rna_tn5 = graph.input("rna_tn5", "RNA-branch Tn5 transposase", MaterialKind.CATALYST)
    rna_tagged = graph.product("rna_tagged", "tagmented spatial cDNA branch",
                               MaterialKind.MOLECULAR, state=rna.state)
    rna_waste = graph.product("rna_tn5_waste", "non-library RNA-branch tagmentation products",
                              MaterialKind.WASTE)
    graph.react("tagment_rna", "Tagment the pre-amplified RNA branch",
                StepKind.TAGMENTATION,
                ((rna, InputRole.TEMPLATE), (rna_transposome, InputRole.STRUCTURAL),
                 (rna_tn5, InputRole.CATALYST)), (rna_tagged, rna_waste))
    graph.discard(rna_waste, "does not retain a spatially addressed transcript insert")
    reviewed.extend(("rna_transposome", "rna_tn5", "tagment_rna"))
    rna_final = role_library(d.title + " RNA-library roles",
                             (SPATIAL_A, SPATIAL_B, UMI, role("transcript insert")), index=SAMPLE)
    _finish_pcr(graph, "rna", rna_tagged, rna_final, (), "RNA branch: " + d.endpoint,
                reviewed, label=d.title + " RNA")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    return _dbit_graph(_definition, FINAL_LIBRARY, SEQ_PRIMERS, cuttag=False)
