"""Protocol-owned spatial assay and reaction graph."""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
from protocol_families.spatial import (UMI, SAMPLE, SPATIAL_A, SPATIAL_B,
    SPATIAL_A_ROUND, SPATIAL_B_ROUND, _finish_pcr, _spatial_dna_libraries,
    _spatial_branch_scene, _SPATIAL_A_FEATURE, _SPATIAL_B_FEATURE, SPATIAL_RT_HANDLE,
    SPATIAL_RT_POLYT, SPATIAL_RT_UMI, SPATIAL_TN5_A, SPATIAL_TN5_B, SPATIAL_TN5_REV,
    SPATIAL_LINKER_1, SPATIAL_LINKER_2, SPATIAL_A_PREFIX, SPATIAL_A1, SPATIAL_A_SUFFIX,
    SPATIAL_B_PREFIX, SPATIAL_B1, SPATIAL_B_SUFFIX, SPATIAL_N501, SPATIAL_N701)
import illumina as il
import nextera as nx
from batch_ngs import joined_scene, seg
from chemdraw import Construct, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import Disposition, Terminal
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
READ_LENGTHS = {'Read 1': 150, 'Read 2': 150}  # Methods: NovaSeq 6000 paired-end 150
READ_USES = None
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
# Zhang et al. 2023 Methods ("Preparation of transposome"; cDNA library construction)
# and Supplementary Table 4 print these as ordered.  The shared constants above are
# letter-identical to the 2023 Table 4/5/6 entries checked in the 2026-10-17 pass.
TSO_WRITTEN = "AAGCAGTGGTATCAACGCAGAGTGAATrGrG+G"  # Methods, cDNA library construction
PCR_PRIMER_2 = "AAGCAGTGGTATCAACGCAGAGT"             # Supplementary Table 4
_TSO_DNA = "AAGCAGTGGTATCAACGCAGAGTGAATGGG"

def oligos():
    p = lambda name, body: f"<p>{name}: 5'- {body} -3'</p>"
    return (
        p("RT primer (/5Phos/; /iBiodT/ is the first T after the UMI)",
          f"/5Phos/{SPATIAL_RT_HANDLE}<umi>{SPATIAL_RT_UMI}</umi>/iBiodT/{SPATIAL_RT_POLYT}"),
        p("Tn5ME-A", SPATIAL_TN5_A),
        p("Tn5MErev", "/5Phos/" + SPATIAL_TN5_REV),
        p("Tn5ME-B", "/5Phos/" + SPATIAL_TN5_B),
        p("Barcode A-1 (of A-1…A-100)", f"/5Phos/{SPATIAL_A_PREFIX}<cbc>{SPATIAL_A1}</cbc>{SPATIAL_A_SUFFIX}"),
        p("Ligation linker 1", SPATIAL_LINKER_1),
        p("Barcode B-1 (of B-1…B-100)", f"{SPATIAL_B_PREFIX}<cbc>{SPATIAL_B1}</cbc>{SPATIAL_B_SUFFIX}"),
        p("Ligation linker 2", SPATIAL_LINKER_2),
        p("Template switch primer", TSO_WRITTEN),
        p("PCR Primer 1 (cDNA pre-amplification)", SPATIAL_B_PREFIX),
        p("PCR Primer 2 (cDNA pre-amplification)", PCR_PRIMER_2),
    )
_template, FINAL_LIBRARY, _geometry, SEQ_PRIMERS = _spatial_dna_libraries("accessible DNA")
_definition = SimpleNamespace(title=TITLE, doi='10.1038/s41586-023-05795-1', citation='Zhang et al., Nature (2023)', endpoint='The exact representative panel is the ATAC branch: Read 1 enters accessible DNA, while Read 2 enters the B–A spatial chain from the N701 end. The separately recovered RNA branch retains the same spatial address plus a UMI and is Nextera-tagmented after pre-amplification.')


def rt_primer_scene() -> Scene:
    """Zhang et al. 2023 Supplementary Table 4 RT primer as ordered:
    /5Phos/ handle, N10, /iBiodT/, T15, VN (the biotin-dT is drawn as its base T)."""
    assert SPATIAL_RT_POLYT == "T" * 15 + "VN"
    sc = Scene()
    sc.strand("RT primer", [
        seg("RT handle", SPATIAL_RT_HANDLE, "r1"),
        seg("10-nt UMI", SPATIAL_RT_UMI, "umi", placeholder=True,
            feature=feature("spatial_rna_umi", "umi", "random")),
        seg("internal biotin-dT (/iBiodT/)", "T"),
        seg("oligo(dT)15", "T" * 15),
        seg("degenerate 3' anchor", "VN", placeholder=True),
    ], label="RT primer (/5Phos/, /iBiodT/)", mod5="p")
    return sc


def transposome_adapters() -> tuple[Scene, Scene]:
    """Zhang et al. 2023 Methods ("Preparation of transposome"): Tn5ME-A and
    /5Phos/Tn5ME-B, each annealed to /5Phos/Tn5MErev, loaded together on Tn5."""
    if not (SPATIAL_TN5_A == nx.S5 + nx.ME and SPATIAL_TN5_B == SPATIAL_RT_HANDLE + nx.ME
            and SPATIAL_TN5_REV == nx.ME_RC):
        raise ValueError("printed Tn5ME-A/B/MErev do not share the 19-bp mosaic end")
    out = []
    for name, outer, tag, mod5 in (("Tn5ME-A", nx.S5, "s5", ""),
                                   ("Tn5ME-B", SPATIAL_RT_HANDLE, "r1", "p")):
        sc = Scene()
        sc.strand(name, [seg(f"{name} outer arm", outer, tag),
                         seg(f"{name} mosaic end", nx.ME, "me")],
                  label=name + (" (/5Phos/)" if mod5 else ""), mod5=mod5)
        sc.anneal(f"Tn5MErev for {name}", [seg("Tn5MErev", SPATIAL_TN5_REV, "me")],
                  to=name, pair=("Tn5MErev", f"{name} mosaic end"),
                  label="Tn5MErev (/5Phos/)", mod5="p")
        out.append(sc)
    return tuple(out)


def barcode_adapter_scene(round_, linker: str, *, phosphorylated: bool) -> Scene:
    """One plate's barcode oligo pre-annealed to its ligation linker (Supplementary
    Tables 4-6); the barcode window is the featured 100-member plate placeholder."""
    splint_prior = revcomp(round_.previous_5_arm)
    if linker != splint_prior + revcomp(round_.barcode_3_arm):
        raise ValueError(f"{round_.name}: linker does not bridge its two arms")
    sc = Scene()
    sc.strand(f"{round_.name} linker",
              [seg(f"{round_.name} prior-molecule splint", splint_prior),
               seg(f"{round_.name} barcode splint", linker[len(splint_prior):])],
              label=f"ligation linker for {round_.name}")
    sc.anneal(f"{round_.name} oligo", round_.product_segments(),
              to=f"{round_.name} linker",
              pair=(f"{round_.name} ligating arm", f"{round_.name} barcode splint"),
              label=f"{round_.name} oligo" + (" (/5Phos/)" if phosphorylated else ""),
              mod5="p" if phosphorylated else "")
    return sc


def index_primers():
    """Representative N501 and N701 (of N701-N707), Supplementary Table 4."""
    i5 = feature("spatial_i5", "sample_index", "whitelist",
                 whitelist="N501 index from Zhang et al. 2023 Supplementary Table 4")
    i7 = feature("spatial_i7", "sample_index", "whitelist",
                 whitelist="N701 index (of N701-N707) from Zhang et al. 2023 Supplementary Table 4")
    n501 = Construct([seg("P5", il.P5, "p5"), seg("i5 (N501)", "TAGATCGC", "cbc", feature=i5),
                      seg("Tn5ME-A outer arm", nx.S5, "s5"), seg("mosaic end", nx.ME, "me")],
                     name="N501")
    n701 = Construct([seg("P7", il.P7, "p7"), seg("i7 (N701)", "TCGCCTTA", "cbc", feature=i7),
                      seg("S7", nx.S7, "s7"), seg("mosaic end", nx.ME, "me"),
                      seg("barcode B PCR arm", SPATIAL_B_PREFIX)], name="N701")
    if "".join(s.top for s in n501) != SPATIAL_N501 or "".join(s.top for s in n701) != SPATIAL_N701:
        raise ValueError("N501/N701 segments do not reassemble the printed primers")
    return n501, n701


def _rna_scene(stage: int, *, switched: bool) -> Scene:
    """Protocol-local RNA branch: the TSO-derived end exists only after on-bead switching.

    Zhang et al. 2023 Methods add the template switch primer to streptavidin-bound cDNA
    after both in situ ligations, reverse crosslinking and DNA/cDNA separation.
    """
    core = [seg("RT handle", SPATIAL_RT_HANDLE, "r1"),
            seg("10-nt UMI", "U" * 10, "umi", placeholder=True,
                feature=feature("spatial_rna_umi", "umi", "random")),
            seg("internal biotin-dT (/iBiodT/)", "T"),
            seg("poly(dT)VN", SPATIAL_RT_POLYT, placeholder=True),
            seg("transcript cDNA", "X" * 28, placeholder=True)]
    if switched:
        core.append(seg("TSO-derived end", revcomp(_TSO_DNA)))
    parts, junctions = core, []
    if stage >= 1:
        parts = [seg("barcode A outgoing arm", SPATIAL_A_PREFIX),
                 seg("barcode A", "A" * 8, "cbc", placeholder=True, feature=_SPATIAL_A_FEATURE),
                 seg("barcode A ligating arm", SPATIAL_A_SUFFIX), *parts]
        junctions.append(("barcode A ligating arm", "RT handle", "ligation A"))
    if stage >= 2:
        parts = [seg("barcode B PCR arm", SPATIAL_B_PREFIX),
                 seg("barcode B", "B" * 8, "cbc", placeholder=True, feature=_SPATIAL_B_FEATURE),
                 seg("barcode B ligating arm", SPATIAL_B_SUFFIX), *parts]
        junctions.insert(0, ("barcode B ligating arm", "barcode A outgoing arm", "ligation B"))
    return joined_scene(parts, tuple(junctions), label="RNA branch", duplex=False)


def _rows(stage: int, dna_kind: str, *, switched: bool = False) -> list:
    return [*_spatial_branch_scene(dna_kind, stage).rows(),
            *_rna_scene(stage, switched=switched).rows()]


def _dbit(*, cuttag: bool):
    dna_kind = "histone-mark DNA" if cuttag else "accessible DNA"
    tag_note = ("Primary antibody (4 °C overnight), guinea-pig anti-rabbit secondary antibody, then pA–Tn5 adaptor complex in 300 mM NaCl wash buffer; tagmentation is activated with 10 mM MgCl2 at 37 °C for 1 h and stopped with 40 mM EDTA."
                if cuttag else
                "Home-made Tn5 transposome (Tn5ME-A and Tn5ME-B, each annealed to /5Phos/Tn5MErev) tags fixed, permeabilized tissue at 37 °C for 30 min; 40 mM EDTA stops transposition.")
    start = Scene(); start.strand("targets", [role(dna_kind), role("poly(A) RNA")],
                                 label="fixed permeabilized tissue")
    wf = Workflow(MolecularState("Fixed tissue targets", tuple(start.rows())))
    wf.react(f"Tag {dna_kind}, then reverse-transcribe RNA in situ",
             _rows(0, dna_kind), name="shared ligation handle on both branches",
             note=tag_note + " In situ RT then uses Maxima H Minus with the /5Phos/ biotinylated RT primer (room temperature 30 min, 42 °C 90 min). Tn5ME-B and the RT primer carry the same CATCGGCGTACGACT handle.")
    g1 = AdapterLigationGeometry.known(
        "splinted", SPATIAL_A_ROUND.adapter_scene(),
        statement="Linker 1 is exactly the reverse complement of the shared RT/Tn5ME-B handle followed by the barcode-A ligating arm.")
    wf.ligate_adapter("Ligate barcode A in the first microfluidic channel set", g1,
                      _rows(1, dna_kind), name="A-barcoded DNA and RNA branches",
                      note="Each barcode A is pre-annealed with linker 1 and ligated with T4 DNA ligase at 37 °C for 30 min. The A-oligo 3′-OH is sealed to the shared 5′-phosphate; A-1 is shown from the 100-member plate (Supplementary Table 5).")
    g2 = AdapterLigationGeometry.known(
        "splinted", SPATIAL_B_ROUND.adapter_scene(),
        statement="Linker 2 is exactly the reverse complement of the barcode-A outgoing arm followed by the barcode-B ligating arm.")
    wf.ligate_adapter("Rotate the device and ligate barcode B", g2,
                      _rows(2, dna_kind), name="two-dimensional pixel address",
                      note="The B-oligo 3′-OH is sealed to the /5Phos/ barcode-A 5′ end; the ordered A–B pair identifies one channel intersection (Supplementary Table 6).")
    wf.react("Reverse-crosslink and separate cDNA from DNA with streptavidin",
             _rows(2, dna_kind), name="DNA in supernatant; biotinylated cDNA on beads",
             note="Proteinase K/SDS lysis at 58 °C for 2 h then 65 °C overnight; Zymo DCC-5 cleanup; Dynabeads MyOne Streptavidin C1 bind the biotinylated RT-primer cDNA, and the supernatant carries the tagmented DNA branch.")
    wf.react("Template-switch the bead-bound cDNA", _rows(2, dna_kind, switched=True),
             name="template-switched cDNA branch",
             note=f"Maxima H Minus with template switch primer {TSO_WRITTEN} on the beads (room temperature 30 min, 42 °C 90 min). PCR Primers 1 and 2 (KAPA HiFi, 5 cycles plus qPCR-determined cycles) then pre-amplify the cDNA between the barcode-B PCR arm and the TSO end; 0.8× AMPure. 1 ng is Nextera XT-tagmented and amplified 12 cycles with a P5 primer and an indexed P7 primer; 0.7× AMPure.")
    _template, final, geometry, _primers = _spatial_dna_libraries(dna_kind)
    pcr_note = ("NEBNext PCR from the purified supernatant: 58 °C 5 min, 72 °C 5 min, then 12 cycles; 1.3× AMPure."
                if cuttag else
                "NEBNext PCR from the purified supernatant: 72 °C 5 min, then 5 cycles plus qPCR-determined cycles; 1× AMPure.")
    wf.pcr("Add the indexed Illumina shell to the DNA branch", geometry,
           duplex_rows(final, label="DNA library", unpaired=(f"{dna_kind} genomic DNA'",)),
           name="indexed DNA-branch library",
           note=pcr_note + " N501 supplies P5/i5 at the Tn5ME-A end; N701 (one of N701–N707) supplies P7/i7 at the barcode-B PCR arm. The initial 72 °C step precedes cycling; its gap-fill product is not drawn separately.")
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
        secondary = graph.input("secondary_antibody", "guinea pig anti-rabbit IgG secondary antibody",
                                MaterialKind.SELECTION)
        bound = graph.product("antibody_bound", "antibody-bound tissue chromatin",
                              MaterialKind.MOLECULAR, state=tissue.state)
        unbound = graph.product("unbound_antibody", "unbound antibody removed by washing", MaterialKind.WASTE)
        graph.react("bind_antibody", "Bind primary and secondary antibodies to the selected histone mark",
                    StepKind.BINDING,
                    ((tissue, InputRole.TEMPLATE), (antibody, InputRole.SELECTION),
                     (secondary, InputRole.SELECTION)),
                    (bound, unbound))
        graph.discard(unbound, "washed away before pA–Tn5 binding")
        substrate = bound
        reviewed.extend(("antibody", "secondary_antibody", "bind_antibody"))
    transposome = graph.input("transposome", "Tn5ME-A / /5Phos/Tn5ME-B adapters annealed to /5Phos/Tn5MErev",
                              MaterialKind.ADAPTER, state=transposome_adapters())
    tn5 = graph.input("tn5", "Diagenode pA–Tn5 (C01070002)" if cuttag else "Diagenode unloaded Tn5 (C01070010)",
                      MaterialKind.CATALYST)
    tagged = graph.product("tagged_tissue", "tissue with tagged DNA and retained poly(A) RNA",
                           MaterialKind.MOLECULAR,
                           state=molecule("shared ligation handle on both branches", (
                               role(dna_kind + " carrying the shared handle"), role("poly(A) RNA"))))
    offtarget = graph.product("tn5_offproducts", "A–A or otherwise non-B-tagged fragments",
                              MaterialKind.MOLECULAR, inferred=True)
    graph.react("tagment", "Tag " + dna_kind + " with the shared-handle transposome",
                StepKind.TAGMENTATION,
                ((substrate, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, offtarget))
    graph.mark_terminal(Terminal("tn5_offproducts", Disposition.UNRESOLVED,
        "Candidate products inferred from mixed A/B transposome loading; not physically separated. They stay with the DNA branch and lack a ligatable B handle, so N501/N701 PCR is not expected to amplify them exponentially; fate not observed in the source."))
    rt_primer = graph.input("rt_primer", "/5Phos/ shared-handle N10 /iBiodT/ oligo(dT)VN RT primer",
                            MaterialKind.OLIGO, state=rt_primer_scene())
    rt = graph.input("reverse_transcriptase", "Maxima H Minus reverse transcriptase", MaterialKind.CATALYST)
    paired = graph.product("paired_branches", "DNA and cDNA branches bearing the shared handle",
                           MaterialKind.MOLECULAR,
                           state=molecule("shared-handle DNA and RNA branches", (
                               role(dna_kind + " fragment"), UMI, role("transcript cDNA"))))
    graph.react("reverse_transcribe", "Reverse-transcribe RNA through the shared RT handle",
                StepKind.REVERSE_TRANSCRIPTION,
                ((tagged, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (paired,))
    a_adapter = graph.input("barcode_a", "/5Phos/ barcode-A oligo pool annealed to linker 1",
                            MaterialKind.ADAPTER,
                            state=barcode_adapter_scene(SPATIAL_A_ROUND, SPATIAL_LINKER_1,
                                                        phosphorylated=True))
    a_ligase = graph.input("a_ligase", "T4 DNA ligase", MaterialKind.CATALYST)
    a_product = graph.product("a_product", "A-barcoded DNA and cDNA branches",
                              MaterialKind.MOLECULAR, state=paired.state)
    g1 = AdapterLigationGeometry.known(
        "splinted", SPATIAL_A_ROUND.adapter_scene(),
        statement="Linker 1 is the reverse complement of the shared RT/Tn5ME-B handle followed by the barcode-A ligating arm.")
    graph.react("ligate_a", "Ligate barcode A in the first channel set", StepKind.LIGATION,
                ((paired, InputRole.TEMPLATE), (a_adapter, InputRole.STRUCTURAL),
                 (a_ligase, InputRole.CATALYST)), (a_product,), geometry=g1)
    b_adapter = graph.input("barcode_b", "barcode-B oligo pool annealed to linker 2",
                            MaterialKind.ADAPTER,
                            state=barcode_adapter_scene(SPATIAL_B_ROUND, SPATIAL_LINKER_2,
                                                        phosphorylated=False))
    b_ligase = graph.input("b_ligase", "T4 DNA ligase", MaterialKind.CATALYST)
    b_product = graph.product("b_product", "two-dimensionally addressed DNA and cDNA branches",
                              MaterialKind.MOLECULAR, state=a_product.state)
    g2 = AdapterLigationGeometry.known(
        "splinted", SPATIAL_B_ROUND.adapter_scene(),
        statement="Linker 2 is the reverse complement of the barcode-A outgoing arm followed by the barcode-B ligating arm.")
    graph.react("ligate_b", "Rotate the device and ligate barcode B", StepKind.LIGATION,
                ((a_product, InputRole.TEMPLATE), (b_adapter, InputRole.STRUCTURAL),
                 (b_ligase, InputRole.CATALYST)), (b_product,), geometry=g2)
    proteinase = graph.input("proteinase_k", "proteinase K / SDS reverse-crosslinking mix",
                             MaterialKind.CATALYST)
    lysate = graph.product("lysate", "released barcoded DNA and cDNA",
                           MaterialKind.MOLECULAR, state=b_product.state)
    graph.react("reverse_crosslink", "Reverse-crosslink and digest the tissue region (58 °C 2 h, 65 °C overnight)",
                StepKind.DEPROTEINIZATION,
                ((b_product, InputRole.TEMPLATE), (proteinase, InputRole.CATALYST)),
                (lysate,), changes_molecule=False)
    beads = graph.input("streptavidin_beads", "Dynabeads MyOne Streptavidin C1", MaterialKind.SELECTION)
    template, _same_final, geometry, _same_primers = _spatial_dna_libraries(dna_kind)
    dna = graph.product("dna_branch", dna_kind + " spatial branch (supernatant)", MaterialKind.MOLECULAR,
                        state=template)
    rna = graph.product("rna_branch", "bead-bound biotinylated spatial cDNA", MaterialKind.MOLECULAR,
                        state=molecule("spatial RNA branch", (SPATIAL_A, SPATIAL_B, UMI,
                            role("transcript cDNA"))))
    graph.react("split_branches", "Capture biotinylated cDNA on streptavidin beads; keep DNA in the supernatant",
                StepKind.AFFINITY_CAPTURE,
                ((lysate, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (dna, rna),
                changes_molecule=False)
    reviewed.extend(("transposome", "tn5", "tagment", "rt_primer", "reverse_transcriptase",
                     "reverse_transcribe", "barcode_a", "a_ligase", "ligate_a",
                     "barcode_b", "b_ligase", "ligate_b", "proteinase_k", "reverse_crosslink",
                     "streptavidin_beads", "split_branches"))
    n501, n701 = index_primers()
    _finish_pcr(graph, "dna", dna, final, primers, "DNA branch: " + d.endpoint,
                reviewed, label=d.title + " DNA", geometry=geometry,
                forward_state=n501, reverse_state=n701)
    tso = graph.input("tso", "template switch primer " + TSO_WRITTEN, MaterialKind.OLIGO,
                      state=_TSO_DNA)
    tso_rt = graph.input("tso_rt", "Maxima H Minus reverse transcriptase (on beads)", MaterialKind.CATALYST)
    switched_state = molecule("template-switched spatial cDNA", (SPATIAL_A, SPATIAL_B, UMI,
                              role("transcript cDNA"), role("TSO-derived end")))
    switched = graph.product("rna_switched", "template-switched bead-bound cDNA",
                             MaterialKind.MOLECULAR, state=switched_state)
    graph.react("template_switch", "Template-switch the bead-bound cDNA", StepKind.TEMPLATE_SWITCHING,
                ((rna, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL), (tso_rt, InputRole.CATALYST)),
                (switched,))
    pre_f = graph.input("cdna_primer_1", "PCR Primer 1 " + SPATIAL_B_PREFIX, MaterialKind.OLIGO,
                        state=SPATIAL_B_PREFIX)
    pre_r = graph.input("cdna_primer_2", "PCR Primer 2 " + PCR_PRIMER_2, MaterialKind.OLIGO,
                        state=PCR_PRIMER_2)
    pre_pol = graph.input("cdna_polymerase", "KAPA HiFi HotStart ReadyMix", MaterialKind.CATALYST)
    amplified = graph.product("rna_amplified", "pre-amplified spatial cDNA", MaterialKind.MOLECULAR,
                              state=switched_state)
    graph.react("cdna_pcr", "Pre-amplify cDNA with PCR Primers 1 and 2", StepKind.PCR,
                ((switched, InputRole.TEMPLATE), (pre_f, InputRole.STRUCTURAL),
                 (pre_r, InputRole.STRUCTURAL), (pre_pol, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(switched_state,
                    reason="Primer 1 matches the barcode-B PCR arm and Primer 2 the TSO 5′ region, but the transcript insert is variable, so no single exact amplicon is drawn."))
    rna_transposome = graph.input("rna_transposome", "kit-loaded Nextera XT transposome adapters",
                                  MaterialKind.ADAPTER)
    rna_tn5 = graph.input("rna_tn5", "Nextera XT Amplicon Tagment Mix transposase", MaterialKind.CATALYST)
    rna_tagged = graph.product("rna_tagged", "tagmented spatial cDNA branch",
                               MaterialKind.MOLECULAR, state=rna.state)
    rna_other = graph.product("rna_tn5_offproducts", "internal cDNA fragments lacking the spatial end",
                              MaterialKind.MOLECULAR, inferred=True)
    graph.react("tagment_rna", "Tagment 1 ng pre-amplified cDNA (Nextera XT, 55 °C 5 min)",
                StepKind.TAGMENTATION,
                ((amplified, InputRole.TEMPLATE), (rna_transposome, InputRole.STRUCTURAL),
                 (rna_tn5, InputRole.CATALYST)),
                (rna_tagged, rna_other))
    graph.mark_terminal(Terminal("rna_tn5_offproducts", Disposition.UNRESOLVED,
        "Candidate products; not physically separated. Which P5 and indexed P7 primers the 12-cycle Nextera PCR used is not printed in the Methods, so their amplification fate is not modelled."))
    reviewed.extend(("tso", "tso_rt", "template_switch", "cdna_primer_1", "cdna_primer_2",
                     "cdna_polymerase", "cdna_pcr", "rna_transposome", "rna_tn5", "tagment_rna"))
    rna_final = role_library(d.title + " RNA-library roles",
                             (SPATIAL_A, SPATIAL_B, UMI, role("transcript insert")), index=SAMPLE)
    _finish_pcr(graph, "rna", rna_tagged, rna_final, (), "RNA branch: " + d.endpoint,
                reviewed, label=d.title + " RNA")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    return _dbit_graph(_definition, FINAL_LIBRARY, SEQ_PRIMERS, cuttag=False)
