"""Shared chemistry for nascent-RNA ends, cap selection and 4sU recoding.

The constructors encode the distinctions that matter experimentally: whether a molecule
was selected because it was capped, because a polymerase extended it in a run-on, or
because a metabolic label was enriched or chemically recoded.  Protocol pages only choose
a route; they cannot silently collapse these mechanisms into generic RNA-seq.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, Workflow, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import chain, molecule, role


@dataclass(frozen=True)
class NascentEndCapture:
    """A source of nascent RNA and the end whose coordinate is informative."""
    source: str
    selected_end: str
    cap_selected: bool = False

    def __post_init__(self) -> None:
        if self.selected_end not in {"5-prime", "3-prime", "both"}:
            raise ValueError("nascent-RNA capture must declare 5-prime, 3-prime or both")

    def RNA(self, name: str = "selected nascent RNA") -> Scene:
        parts = [role("5-prime cap" if self.cap_selected else "native 5-prime end"),
                 role("nascent RNA body"), role("polymerase-proximal 3-prime end")]
        sc = molecule(name, tuple(parts))
        sc.note("molecule", f"{self.source}; informative end: {self.selected_end}")
        return sc


@dataclass(frozen=True)
class FourSUReadout:
    """A validated 4sU route: enrichment and base recoding are independent choices."""
    chemistry: str
    output: str
    enrichment: bool

    def __post_init__(self) -> None:
        if self.output not in {"selected 4sU-RNA", "U-to-C mismatch", "native RNA ends"}:
            raise ValueError("unknown 4sU readout")
        if self.output == "U-to-C mismatch" and not self.chemistry:
            raise ValueError("base recoding requires a named chemistry")

    def labelled(self) -> Scene:
        return molecule("pulse-labelled RNA", (role("pre-existing RNA"),
                        Segment("incorporated 4sU", "4sU", placeholder=True, length_bp=1),
                        role("new RNA")))

    def product(self) -> Scene:
        if self.output == "U-to-C mismatch":
            mark = Segment("recoded 4sU read as C", "C", "me", placeholder=True,
                           length_bp=1)
        else:
            mark = Segment(self.output, "U", "me", placeholder=True, length_bp=1)
        sc = molecule(self.output, (role("RNA before label"), mark,
                                    role("RNA after label")))
        sc.note("molecule", self.chemistry)
        return sc


SAMPLE_INDEX = role("sample index", tag="cbc", length_bp=8,
                    molecular_feature=feature("sample_index", "sample_index", "unknown"))


def role_library(name: str, payload: tuple[Segment, ...], *, indexed: bool = True,
                 platform: str = "Illumina-compatible") -> Construct:
    parts = [role(f"{platform} left end")]
    if indexed:
        parts.append(SAMPLE_INDEX)
    parts.extend(payload)
    parts.append(role(f"{platform} right end"))
    return Construct(parts, name=name)


def _end_capture(key: str) -> Workflow:
    if key == "mnet-seq":
        cap = NascentEndCapture("MNase-solubilized Pol II CTD immunoprecipitate", "3-prime")
        return chain(molecule("native chromatin", (role("Pol II elongation complex"),
                                                    role("nascent RNA"))), (
            ("Digest native chromatin with MNase and immunoprecipitate a Pol II CTD state",
             cap.RNA(), "The protected RNA 3′ end is the Pol II active-site or a cotranscriptional-processing end."),
            ("Phosphorylate RNA 5′ ends, size-select and ligate the two RNA adapters",
             molecule("adapter-ligated nascent RNA", (role("5-prime RNA adapter"),
                 role("Pol II-associated RNA"), role("3-prime RNA adapter")),
                 junctions=(("5-prime RNA adapter", "Pol II-associated RNA", "RNA ligation"),
                            ("Pol II-associated RNA", "3-prime RNA adapter", "RNA ligation"))),
             "The original implementation used the TruSeq Small RNA workflow after end repair."),
            ("Reverse-transcribe, PCR-index and sequence from the RNA 3′ end",
             molecule("mNET-seq library", (role("left sequencing end"), SAMPLE_INDEX,
                                           role("nascent RNA cDNA"), role("right sequencing end"))),
             "Read orientation preserves the nucleotide at the selected RNA 3′ end."),
        ))
    if key == "start-seq":
        cap = NascentEndCapture("short capped RNA from isolated nuclei", "both", True)
        return chain(molecule("nuclear RNA", (role("short capped Start-RNA"),
                                              role("other nuclear RNA"))), (
            ("Size-select 25–65-nt short nuclear RNAs and retain capped Start-RNAs",
             cap.RNA(), "The capped 5′ end marks initiation; the native 3′ end marks the paused polymerase."),
            ("Prepare separate 5′-end and 3′-end adapter libraries",
             molecule("two directional libraries", (role("5-prime Start-RNA library"),
                                                     role("3-prime Start-RNA library"))),
             "Purpose-designed adapters keep both libraries compatible with Illumina sequencing."),
            ("Sequence the two ends of the same Start-RNA population",
             molecule("Start-seq endpoints", (role("transcription start nucleotide"),
                                               role("pause-site nucleotide"))),
             "The paired measurements distinguish initiation position from promoter-proximal pausing."),
        ))
    cap = NascentEndCapture("biotin-NTP nuclear run-on product", "5-prime", True)
    name = "GRO-cap" if key == "gro-cap" else "PRO-cap"
    runon = ("BrUTP run-on and anti-BrdU enrichment" if key == "gro-cap"
             else "single-biotin-NTP precision run-on and streptavidin enrichment")
    return chain(molecule("isolated nuclei", (role("engaged RNA polymerase"),
                                              role("nascent RNA"))), (
        (runon, molecule("enriched run-on RNA", (role("capped nascent 5-prime end"),
                                                  role("labelled run-on 3-prime end"))),
         "Only transcripts extended by an engaged polymerase acquire the affinity label."),
        ("Remove 5′ phosphates from uncapped RNA, decap with TAP and ligate the 5′ adapter",
         cap.RNA(f"{name} cap-selected RNA"),
         "Decapping converts the protected cap into the ligatable monophosphate; a no-TAP control measures background."),
        ("Install the other adapter, reverse-transcribe, index and amplify",
         molecule(f"{name} library", (role("left sequencing end"), SAMPLE_INDEX,
                                      role("cap-selected nascent cDNA"),
                                      role("right sequencing end"))),
         "The read begins at the nascent capped transcription start site."),
    ))


def _four_su(key: str) -> Workflow:
    cfg = {
        "timelapse-seq": FourSUReadout("oxidation plus 2,2,2-trifluoroethylamine nucleophilic substitution", "U-to-C mismatch", False),
        "tuc-seq": FourSUReadout("osmium-tetroxide-mediated 4sU-to-cytidine conversion", "U-to-C mismatch", False),
        "snu-seq": FourSUReadout("thiol biotinylation, streptavidin purification and bacterial poly(A) polymerase", "native RNA ends", True),
        "sf4su-seq": FourSUReadout("thiol biotinylation, streptavidin purification and ≈100-nt denaturing-gel selection", "native RNA ends", True),
    }[key]
    if key in {"timelapse-seq", "tuc-seq"}:
        return chain(molecule("living cells", (role("RNA synthesis during pulse"),)), (
            ("Pulse-label newly synthesized RNA with 4-thiouridine", cfg.labelled(),
             "Pre-existing RNA lacks 4sU; newly made RNA incorporates it."),
            ("Chemically recode 4sU before reverse transcription", cfg.product(),
             "The chemistry changes base-pairing so labelled uridines are read as cytidines."),
            ("Prepare an indexed RNA-seq library", molecule("recoding library", (
                role("left sequencing end"), SAMPLE_INDEX, role("cDNA with U-to-C marks"),
                role("right sequencing end"))),
             "New molecules are recognized from distributed T-to-C mismatches rather than a separate barcode."),
        ))
    extra = ("Polyadenylate the intact enriched RNA with bacterial poly(A) polymerase"
             if key == "snu-seq" else "Gel-select the abundant ≈100-nt enriched RNA fraction")
    result = ("3-prime-anchored intact 4sU-RNA" if key == "snu-seq"
              else "size-resolved 5-prime and 3-prime RNA ends")
    return chain(molecule("living cells", (role("RNA synthesis during pulse"),)), (
        ("Pulse-label RNA with 4sU, thiol-biotinylate and enrich on streptavidin", cfg.labelled(),
         "This route selects labelled molecules; it does not recode 4sU."),
        (extra, cfg.product(),
         "SNU-Seq and sf4sU-Seq deliberately create different endpoint information from the same labelled pool."),
        ("Build the directional sequencing library", molecule(result, (role("left sequencing end"),
            *(() if key == "snu-seq" else (SAMPLE_INDEX,)),
            role("4sU-enriched RNA-derived cDNA"), role("right sequencing end"))),
         "The final read coordinate reflects the selected intact or size-fractionated RNA end."),
    ))


def _cagescan() -> Workflow:
    return chain(molecule("total RNA", (role("capped transcript"),)), (
        ("Random-prime reverse transcription with the CAGEscan N6 primer",
         molecule("first-strand cDNA", (role("RT handle"), role("random-primed transcript cDNA"))),
         "The downstream mate is fixed by random priming rather than restriction cleavage."),
        ("Template-switch only at the capped RNA 5′ end",
         molecule("cap-linked cDNA", (role("template-switch 5-prime tag"),
                                      role("transcript cDNA"), role("random-primer end"))),
         "The cap-dependent template switch connects a promoter end to its transcript body."),
        ("PCR with paired-end Illumina primers",
         molecule("CAGEscan library", (role("P5 / cap-end sequencing handle"), SAMPLE_INDEX,
                                       role("5-prime transcript cDNA"),
                                       role("downstream transcript cDNA"),
                                       role("P7 / downstream sequencing handle"))),
         "Paired reads connect a transcription start site to a downstream exon or transcript annotation."),
    ))


@dataclass(frozen=True)
class Definition:
    title: str
    doi: str
    citation: str
    summary: str
    endpoint: str
    notes: str
    caveat: str = ""


DEFS = {
    "mnet-seq": Definition("mNET-seq — mammalian native elongating transcript sequencing", "10.1016/j.cell.2015.03.027", "Nojima et al., Cell (2015)", "MNase-solubilized native chromatin and CTD-state-specific Pol II immunoprecipitation expose the 3′ ends of polymerase-associated nascent RNA at nucleotide resolution.", "Read orientation reports the RNA 3′ end held at the Pol II active site or generated by cotranscriptional processing.", "01_mnet_seq.html"),
    "start-seq": Definition("Start-seq — start-site-associated RNA sequencing", "10.1126/science.1181421", "Nechaev et al., Science (2010)", "Short capped nuclear RNAs are sequenced from both ends so their 5′ ends report initiation and their 3′ ends report promoter-proximal pausing.", "Separate compatible libraries read the 5′ and 3′ ends of 25–65-nt capped Start-RNAs.", "01_start_seq.html"),
    "gro-cap": Definition("GRO-cap — capped global run-on sequencing", "10.1038/ng.3142", "Core et al., Nature Genetics (2014)", "BrUTP-labelled global run-on RNA is subjected to enzymatic cap selection, enriching the capped 5′ ends of nascent promoter and enhancer transcripts.", "The informative base is the first transcribed nucleotide immediately downstream of the ligated 5′ adapter.", "01_gro_cap.html"),
    "pro-cap": Definition("PRO-cap — capped precision run-on sequencing", "10.1126/science.1229386", "Kwak et al., Science (2013)", "Precision run-on marks engaged polymerases with one biotin-NTP while cap selection redirects the library to nascent transcription start sites.", "The cap-selected read begins at the initiating nucleotide; PRO-seq from the same pool reports the active-site 3′ end.", "01_pro_cap.html"),
    "timelapse-seq": Definition("TimeLapse-seq — nucleoside-recoding RNA dynamics", "10.1038/nmeth.4582", "Schofield et al., Nature Methods (2018)", "Oxidative nucleophilic substitution converts 4sU in newly synthesized RNA into a cytidine-like base, producing T-to-C mismatches in ordinary RNA-seq reads.", "Standard RNA-seq reads identify newly synthesized molecules from distributed T-to-C conversions.", "01_timelapse_seq.html"),
    "tuc-seq": Definition("TUC-seq — thiouridine-to-cytidine conversion sequencing", "10.1002/anie.201707465", "Riml et al., Angewandte Chemie (2017)", "Osmium-mediated chemistry converts incorporated 4sU into native cytidine, encoding newly transcribed RNA as U-to-C sequence changes.", "RNA-seq reads distinguish labelled molecules through U-to-C conversions created before reverse transcription.", "01_tuc_seq.html"),
    "snu-seq": Definition("SNU-Seq — single-nucleotide-resolution 4sU sequencing", "10.1093/nar/gkag652", "Gerlevik et al., Nucleic Acids Research (2026)", "Biotin-selected intact 4sU RNA is polyadenylated with bacterial poly(A) polymerase and captured by a 3′ RNA-seq workflow to preserve native endpoint information.", "The 3′-anchored read reports the native end of an intact pulse-labelled RNA molecule.", "01_snu_seq.html"),
    "sf4su-seq": Definition("sf4sU-Seq — size-fractionated 4sU sequencing", "10.1093/nar/gkag652", "Gerlevik et al., Nucleic Acids Research (2026)", "Biotin-enriched 4sU RNA is gel-selected around 100 nt and sequenced directionally from both ends to map initiation and promoter-proximal pausing.", "The two read ends report the native 5′ and 3′ coordinates of the selected short labelled RNA.", "01_sf4su_seq.html"),
    "cagescan": Definition("CAGEscan — paired-end cap-to-transcript sequencing", "10.1038/nmeth.1470", "Plessy et al., Nature Methods (2010)", "Cap-dependent template switching marks the transcript 5′ end while a random RT primer supplies a downstream mate, linking promoters to transcript bodies without tag cleavage.", "One mate enters at the cap-derived 5′ end and the other enters from the random-primer-defined downstream end.", "01_cagescan.html"),
}


CAGE_OLIGOS = (
    "5′-TAGTCGAACTGAAGGTCTCCAGCAXnrGrGrG-3′ — CAGEscan template-switch primer (Xn is a 3- or 6-base DNA barcode; terminal rGrGrG is RNA)",
    "5′-TAGTCGAACGAAGGTCTCCGAACCGCTCTTCCGATCTN6-3′ — CAGEscan RT primer",
    "5′-TAGTCGAACTGAAGGTCTCCAGC-3′ — CAGEscan semi-suppressive PCR forward primer",
    "5′-TGACGTCGTCTAGTCGAACTGAAGGTCTCCGAACC-3′ — CAGEscan semi-suppressive PCR reverse primer",
    "5′-AATGATACGGCGACCACCGAGATCTACACTAGTCGAACTGAAGG-3′ — CAGEscan forward PCR primer",
    "5′-CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT-3′ — CAGEscan reverse PCR primer",
    "5′-TAGTCGAACTGAAGGTCTCCAGCA-3′ — CAGEscan forward paired-end sequencing primer",
    "5′-CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT-3′ — CAGEscan reverse paired-end sequencing primer",
)


def _ligate_rna_adapter(graph: ProtocolGraph, prefix: str, template, label: str,
                        reviewed: list[str]):
    adapter = graph.input(prefix + "_adapter", label + " RNA adapter", MaterialKind.ADAPTER)
    ligase = graph.input(prefix + "_ligase", label + " RNA ligase", MaterialKind.CATALYST)
    product = graph.product(prefix + "_product", label + " adapter-ligated RNA",
                            MaterialKind.MOLECULAR,
                            state=molecule(label + " adapter-ligated RNA", (
                                role(label + " adapter"), role("selected RNA"))))
    graph.react(prefix + "_ligation", "Ligate the " + label + " RNA adapter",
                StepKind.LIGATION,
                ((template, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (product,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The defining source establishes this RNA-adapter role but not one immutable ligating-end sequence for the page model."))
    reviewed.extend((prefix + "_adapter", prefix + "_ligase", prefix + "_ligation"))
    return product


def _finish_cdna(graph: ProtocolGraph, prefix: str, template, final: Construct,
                 endpoint: str, reviewed: list[str], *, rt_label: str = "adapter RT primer"):
    rt_primer = graph.input(prefix + "_rt_primer", rt_label, MaterialKind.OLIGO)
    reverse_transcriptase = graph.input(prefix + "_reverse_transcriptase", "reverse transcriptase",
                                       MaterialKind.CATALYST)
    cdna = graph.product(prefix + "_cdna", "assay-derived first-strand cDNA",
                         MaterialKind.MOLECULAR,
                         state=molecule("assay-derived cDNA", (role("selected RNA cDNA"),)))
    graph.react(prefix + "_rt", "Reverse-transcribe the selected RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((template, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (cdna,))
    forward = graph.input(prefix + "_pcr_forward", "library PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input(prefix + "_pcr_reverse", "indexed library PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input(prefix + "_pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", "indexed sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr", "Index and amplify the cDNA library", StepKind.PCR,
                ((cdna, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The source establishes read direction and library roles but not one immutable complete primer geometry."))
    graph.sequence(library, SequencingPlan(limitation=endpoint))
    reviewed.extend((
        prefix + "_rt_primer", prefix + "_reverse_transcriptase", prefix + "_rt",
        prefix + "_pcr_forward", prefix + "_pcr_reverse", prefix + "_pcr_polymerase",
        prefix + "_pcr", prefix + "_library",
    ))
    return library


def _mnet_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    chromatin = graph.input(
        "native_chromatin", "native chromatin containing Pol II elongation complexes",
        MaterialKind.BIOLOGICAL,
        state=molecule("native chromatin", (role("Pol II elongation complex"),
                                             role("nascent RNA"))), primary=True)
    mnase = graph.input("mnase", "micrococcal nuclease", MaterialKind.CATALYST)
    solubilized = graph.product("solubilized_complexes", "MNase-solubilized Pol II complexes",
                                MaterialKind.MOLECULAR, state=chromatin.state)
    graph.react("mnase_digest", "Digest native chromatin with MNase",
                StepKind.DIGESTION,
                ((chromatin, InputRole.TEMPLATE), (mnase, InputRole.CATALYST)),
                (solubilized,))
    antibody = graph.input("ctd_antibody", "Pol II CTD-state-specific antibody",
                           MaterialKind.SELECTION)
    selected = graph.product(
        "polii_rna", "Pol II-associated nascent RNA complexes",
        MaterialKind.MOLECULAR,
        state=NascentEndCapture(
            "MNase-solubilized Pol II CTD immunoprecipitate", "3-prime").RNA())
    unbound = graph.product("unbound_complexes", "unbound chromatin complexes",
                            MaterialKind.WASTE)
    graph.react("immunoprecipitate", "Immunoprecipitate the selected Pol II CTD state",
                StepKind.AFFINITY_CAPTURE,
                ((solubilized, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                (selected, unbound))
    graph.discard(unbound, "does not carry the selected Pol II CTD state")
    extraction = graph.input("rna_extraction", "protein removal and RNA-extraction condition",
                             MaterialKind.CONDITION)
    rna = graph.product("selected_rna", "purified Pol II-associated nascent RNA",
                        MaterialKind.MOLECULAR, state=selected.state)
    graph.react("extract_rna", "Remove protein and recover nascent RNA",
                StepKind.DEPROTEINIZATION,
                ((selected, InputRole.TEMPLATE), (extraction, InputRole.CONDITION)), (rna,))
    kinase = graph.input("rna_kinase", "RNA 5-prime polynucleotide kinase",
                         MaterialKind.CATALYST)
    phosphorylated = graph.product("phosphorylated_rna", "5-prime-phosphorylated nascent RNA",
                                   MaterialKind.MOLECULAR, state=rna.state)
    graph.react("phosphorylate", "Phosphorylate RNA 5-prime ends",
                StepKind.PHOSPHORYLATION,
                ((rna, InputRole.TEMPLATE), (kinase, InputRole.CATALYST)),
                (phosphorylated,))
    size_rule = graph.input("size_selection", "nascent-RNA size-selection rule",
                            MaterialKind.SELECTION)
    sized = graph.product("sized_rna", "size-selected Pol II-associated RNA",
                          MaterialKind.MOLECULAR, state=phosphorylated.state)
    rejected = graph.product("rejected_rna", "RNA outside the selected size range",
                             MaterialKind.WASTE)
    graph.react("select_size", "Size-select Pol II-associated RNA", StepKind.SELECTION,
                ((phosphorylated, InputRole.TEMPLATE), (size_rule, InputRole.SELECTION)),
                (sized, rejected))
    graph.discard(rejected, "falls outside the library's nascent-RNA size window")
    reviewed = [
        "native_chromatin", "mnase", "mnase_digest", "ctd_antibody",
        "immunoprecipitate", "rna_extraction", "extract_rna", "rna_kinase",
        "phosphorylate", "size_selection", "select_size",
    ]
    ligated3 = _ligate_rna_adapter(graph, "three_prime", sized, "3-prime", reviewed)
    ligated5 = _ligate_rna_adapter(graph, "five_prime", ligated3, "5-prime", reviewed)
    _finish_cdna(graph, "library", ligated5, final, d.endpoint, reviewed)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _runon_cap_graph(key: str, d: Definition, final: Construct) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    nuclei = graph.input(
        "nuclei", "isolated nuclei containing engaged RNA polymerases",
        MaterialKind.BIOLOGICAL,
        state=molecule("isolated nuclei", (role("engaged RNA polymerase"),
                                            role("nascent RNA"))), primary=True)
    label_name = "BrUTP" if key == "gro-cap" else "single biotin-NTP"
    label = graph.input("labelled_ntp", label_name + " run-on nucleotide",
                        MaterialKind.MOLECULAR)
    polymerase = graph.input("engaged_polymerase", "engaged endogenous RNA polymerase",
                             MaterialKind.CATALYST)
    runon = graph.product(
        "runon_rna", "affinity-labelled nuclear run-on RNA",
        MaterialKind.MOLECULAR,
        state=molecule("labelled run-on RNA", (role("capped nascent 5-prime end"),
                                                role("labelled run-on 3-prime end"))))
    graph.react("run_on", "Extend engaged polymerases with " + label_name,
                StepKind.RUN_ON,
                ((nuclei, InputRole.TEMPLATE), (label, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (runon,))
    capture_label = ("anti-BrdU affinity reagent" if key == "gro-cap"
                     else "streptavidin affinity reagent")
    capture = graph.input("runon_capture", capture_label, MaterialKind.SELECTION)
    enriched = graph.product("enriched_runon", "affinity-enriched run-on RNA",
                             MaterialKind.MOLECULAR, state=runon.state)
    unlabelled = graph.product("unlabelled_rna", "RNA not extended in the run-on",
                               MaterialKind.WASTE)
    graph.react("enrich_runon", "Affinity-enrich labelled run-on RNA",
                StepKind.AFFINITY_CAPTURE,
                ((runon, InputRole.TEMPLATE), (capture, InputRole.SELECTION)),
                (enriched, unlabelled))
    graph.discard(unlabelled, "did not acquire the run-on affinity label")
    phosphatase = graph.input("phosphatase", "uncapped-RNA 5-prime phosphatase",
                              MaterialKind.CATALYST)
    dephosphorylated = graph.product("dephosphorylated_rna", "run-on RNA after uncapped-end blocking",
                                    MaterialKind.MOLECULAR, state=enriched.state)
    graph.react("dephosphorylate", "Remove 5-prime phosphates from uncapped RNA",
                StepKind.DEPHOSPHORYLATION,
                ((enriched, InputRole.TEMPLATE), (phosphatase, InputRole.CATALYST)),
                (dephosphorylated,))
    tap = graph.input("tap", "tobacco acid pyrophosphatase", MaterialKind.CATALYST)
    cap_selected = graph.product(
        "cap_selected_rna", "decapped nascent RNA with ligatable 5-prime monophosphate",
        MaterialKind.MOLECULAR,
        state=molecule(
            ("GRO-cap" if key == "gro-cap" else "PRO-cap") + " cap-selected RNA",
            (role("ligatable decapped 5-prime monophosphate"),
             role("nascent RNA body"), role("labelled run-on 3-prime end"))))
    graph.react("decap", "Decap protected RNA ends with TAP",
                StepKind.CHEMICAL_CONVERSION,
                ((dephosphorylated, InputRole.TEMPLATE), (tap, InputRole.CATALYST)),
                (cap_selected,))
    reviewed = [
        "nuclei", "labelled_ntp", "engaged_polymerase", "run_on", "runon_capture",
        "enrich_runon", "phosphatase", "dephosphorylate", "tap", "decap",
    ]
    ligated5 = _ligate_rna_adapter(graph, "five_prime", cap_selected, "5-prime", reviewed)
    ligated3 = _ligate_rna_adapter(graph, "three_prime", ligated5, "3-prime", reviewed)
    _finish_cdna(graph, "library", ligated3, final, d.endpoint, reviewed)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _start_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    nuclear_rna = graph.input(
        "nuclear_rna", "nuclear RNA containing short capped Start-RNAs",
        MaterialKind.BIOLOGICAL,
        state=molecule("nuclear RNA", (role("short capped Start-RNA"),
                                       role("other nuclear RNA"))), primary=True)
    size_rule = graph.input("size_selection", "25–65-nt RNA size-selection rule",
                            MaterialKind.SELECTION)
    short_rna = graph.product("short_rna", "25–65-nt short nuclear RNA",
                              MaterialKind.MOLECULAR, state=nuclear_rna.state)
    long_rna = graph.product("rejected_sizes", "nuclear RNA outside 25–65 nt",
                             MaterialKind.WASTE)
    graph.react("select_size", "Select 25–65-nt short nuclear RNAs",
                StepKind.SELECTION,
                ((nuclear_rna, InputRole.TEMPLATE), (size_rule, InputRole.SELECTION)),
                (short_rna, long_rna))
    graph.discard(long_rna, "falls outside the Start-RNA size interval")
    cap_reagent = graph.input("cap_selection", "capped-RNA retention chemistry",
                              MaterialKind.SELECTION)
    capped = graph.product(
        "capped_start_rna", "capped Start-RNA population", MaterialKind.MOLECULAR,
        state=NascentEndCapture("short capped RNA from isolated nuclei", "both", True).RNA())
    uncapped = graph.product("uncapped_short_rna", "uncapped short nuclear RNA",
                             MaterialKind.WASTE)
    graph.react("select_caps", "Retain capped Start-RNAs", StepKind.SELECTION,
                ((short_rna, InputRole.TEMPLATE), (cap_reagent, InputRole.SELECTION)),
                (capped, uncapped))
    graph.discard(uncapped, "lacks the protected capped end required by Start-seq")
    five_input = graph.product("five_end_input", "Start-RNA assigned to 5-prime library",
                               MaterialKind.MOLECULAR, state=capped.state)
    three_input = graph.product("three_end_input", "Start-RNA assigned to 3-prime library",
                                MaterialKind.MOLECULAR, state=capped.state)
    graph.react("split_end_libraries", "Split Start-RNA into 5-prime and 3-prime libraries",
                StepKind.SPLIT, ((capped, InputRole.TEMPLATE),),
                (five_input, three_input), changes_molecule=False)
    reviewed = [
        "nuclear_rna", "size_selection", "select_size", "cap_selection", "select_caps",
        "split_end_libraries",
    ]
    five3 = _ligate_rna_adapter(graph, "five_library_three", five_input,
                                "5-prime-library 3-prime", reviewed)
    five5 = _ligate_rna_adapter(graph, "five_library_five", five3,
                                "5-prime-library 5-prime", reviewed)
    _finish_cdna(graph, "five_library", five5, final,
                 "5-prime Start-seq library: " + d.endpoint, reviewed)
    three3 = _ligate_rna_adapter(graph, "three_library_three", three_input,
                                 "3-prime-library 3-prime", reviewed)
    three5 = _ligate_rna_adapter(graph, "three_library_five", three3,
                                 "3-prime-library 5-prime", reviewed)
    _finish_cdna(graph, "three_library", three5, final,
                 "3-prime Start-seq library: " + d.endpoint, reviewed)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _four_su_graph(key: str, d: Definition, final: Construct) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input(
        "cells", "living cells during the 4sU pulse",
        MaterialKind.BIOLOGICAL,
        state=molecule("living cells", (role("RNA synthesis during pulse"),)), primary=True)
    four_su = graph.input("four_su", "4-thiouridine pulse-label nucleotide",
                          MaterialKind.MOLECULAR)
    polymerase = graph.input("cellular_rna_polymerase", "cellular RNA polymerases",
                             MaterialKind.CATALYST)
    cfg = {
        "timelapse-seq": FourSUReadout(
            "oxidation plus 2,2,2-trifluoroethylamine nucleophilic substitution",
            "U-to-C mismatch", False),
        "tuc-seq": FourSUReadout(
            "osmium-tetroxide-mediated 4sU-to-cytidine conversion",
            "U-to-C mismatch", False),
        "snu-seq": FourSUReadout(
            "thiol biotinylation, streptavidin purification and bacterial poly(A) polymerase",
            "native RNA ends", True),
        "sf4su-seq": FourSUReadout(
            "thiol biotinylation, streptavidin purification and ≈100-nt denaturing-gel selection",
            "native RNA ends", True),
    }[key]
    labelled = graph.product("labelled_rna", "4sU pulse-labelled RNA",
                             MaterialKind.MOLECULAR, state=cfg.labelled())
    graph.react("pulse_label", "Pulse-label newly synthesized RNA with 4sU",
                StepKind.TRANSCRIPTION,
                ((cells, InputRole.TEMPLATE), (four_su, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (labelled,))
    reviewed = ["cells", "four_su", "cellular_rna_polymerase", "pulse_label"]

    if key in {"timelapse-seq", "tuc-seq"}:
        chemistry = graph.input("recoding_chemistry", cfg.chemistry,
                                MaterialKind.CONDITION)
        selected = graph.product("recoded_rna", "chemically recoded 4sU RNA",
                                 MaterialKind.MOLECULAR, state=cfg.product())
        graph.react("recode", "Chemically recode 4sU before reverse transcription",
                    StepKind.CHEMICAL_CONVERSION,
                    ((labelled, InputRole.TEMPLATE), (chemistry, InputRole.CONDITION)),
                    (selected,))
        reviewed.extend(("recoding_chemistry", "recode"))
    else:
        biotin = graph.input("thiol_biotin", "4sU thiol-biotinylation reagent",
                             MaterialKind.MOLECULAR)
        condition = graph.input("biotinylation_condition", "thiol-biotinylation condition",
                                MaterialKind.CONDITION)
        biotinylated = graph.product("biotinylated_rna", "biotinylated 4sU RNA",
                                    MaterialKind.MOLECULAR, state=labelled.state)
        graph.react("biotinylate", "Thiol-biotinylate incorporated 4sU",
                    StepKind.CHEMICAL_CONVERSION,
                    ((labelled, InputRole.TEMPLATE), (biotin, InputRole.STRUCTURAL),
                     (condition, InputRole.CONDITION)), (biotinylated,))
        streptavidin = graph.input("streptavidin", "streptavidin RNA-capture reagent",
                                  MaterialKind.SELECTION)
        enriched = graph.product("enriched_rna", "streptavidin-enriched 4sU RNA",
                                 MaterialKind.MOLECULAR, state=cfg.product())
        unlabelled = graph.product("unlabelled_rna", "unlabelled RNA background",
                                   MaterialKind.WASTE)
        graph.react("enrich", "Affinity-purify biotinylated 4sU RNA",
                    StepKind.AFFINITY_CAPTURE,
                    ((biotinylated, InputRole.TEMPLATE),
                     (streptavidin, InputRole.SELECTION)), (enriched, unlabelled))
        graph.discard(unlabelled, "does not carry the pulse-installed 4sU affinity label")
        reviewed.extend((
            "thiol_biotin", "biotinylation_condition", "biotinylate",
            "streptavidin", "enrich",
        ))
        if key == "snu-seq":
            pap = graph.input("poly_a_polymerase", "bacterial poly(A) polymerase",
                              MaterialKind.CATALYST)
            atp = graph.input("tailing_atp", "ATP for RNA polyadenylation",
                              MaterialKind.MOLECULAR)
            selected = graph.product("selected_rna", "polyadenylated intact 4sU RNA",
                                     MaterialKind.MOLECULAR, state=enriched.state)
            graph.react("polyadenylate", "Polyadenylate intact enriched RNA",
                        StepKind.TAILING,
                        ((enriched, InputRole.TEMPLATE), (atp, InputRole.STRUCTURAL),
                         (pap, InputRole.CATALYST)), (selected,))
            reviewed.extend(("poly_a_polymerase", "tailing_atp", "polyadenylate"))
        else:
            size_rule = graph.input("size_selection", "approximately 100-nt gel-selection rule",
                                    MaterialKind.SELECTION)
            selected = graph.product("selected_rna", "approximately 100-nt enriched 4sU RNA",
                                     MaterialKind.MOLECULAR, state=enriched.state)
            rejected = graph.product("rejected_sizes", "enriched RNA outside the selected size",
                                     MaterialKind.WASTE)
            graph.react("select_size", "Gel-select the approximately 100-nt RNA fraction",
                        StepKind.SELECTION,
                        ((enriched, InputRole.TEMPLATE), (size_rule, InputRole.SELECTION)),
                        (selected, rejected))
            graph.discard(rejected, "falls outside the size-resolved endpoint population")
            reviewed.extend(("size_selection", "select_size"))

    ligated = _ligate_rna_adapter(graph, "library_adapter", selected,
                                  "directional library", reviewed)
    _finish_cdna(graph, "library", ligated, final, d.endpoint, reviewed)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _cagescan_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    rna = graph.input(
        "total_rna", "total RNA containing capped transcripts",
        MaterialKind.BIOLOGICAL,
        state=molecule("total RNA", (role("capped transcript"),)), primary=True)
    n6 = graph.input("n6_rt_primer", "published CAGEscan N6 reverse-transcription primer",
                     MaterialKind.OLIGO)
    reverse_transcriptase = graph.input("reverse_transcriptase", "template-switching reverse transcriptase",
                                       MaterialKind.CATALYST)
    first = graph.product(
        "first_strand", "random-primed first-strand cDNA",
        MaterialKind.MOLECULAR,
        state=molecule("first-strand cDNA", (role("RT handle"),
                                              role("random-primed transcript cDNA"))))
    graph.react("reverse_transcription", "Random-prime reverse transcription",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna, InputRole.TEMPLATE), (n6, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (first,))
    tso = graph.input("cap_tso", "cap-dependent CAGEscan template-switch oligo",
                      MaterialKind.OLIGO)
    switched = graph.product(
        "cap_linked_cdna", "cap-linked CAGEscan cDNA",
        MaterialKind.MOLECULAR,
        state=molecule("cap-linked cDNA", (role("template-switch 5-prime tag"),
                                            role("transcript cDNA"), role("random-primer end"))))
    graph.react("template_switch", "Template-switch at the capped RNA 5-prime end",
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (switched,))
    forward = graph.input("pcr_forward", "published CAGEscan forward PCR primer",
                          MaterialKind.OLIGO)
    reverse = graph.input("pcr_reverse", "published CAGEscan reverse PCR primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("sequencing_library", "paired-end CAGEscan library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("library_pcr", "PCR with paired-end Illumina primers", StepKind.PCR,
                ((switched, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The historical PCR oligos are published, but this role-only construct does not assert modern instrument-primer landing geometry."))
    graph.sequence(library, SequencingPlan(limitation=d.endpoint))
    graph.declare_reviewed_scope((
        "total_rna", "n6_rt_primer", "reverse_transcriptase", "reverse_transcription",
        "cap_tso", "template_switch", "pcr_forward", "pcr_reverse", "pcr_polymerase",
        "library_pcr", "sequencing_library",
    ), source)
    return graph.require_publishable()


def get(key: str):
    d = DEFS[key]
    workflow_fn = (_cagescan if key == "cagescan" else
                   (lambda: _four_su(key)) if key in {"timelapse-seq", "tuc-seq", "snu-seq", "sf4su-seq"}
                   else (lambda: _end_capture(key)))
    payload = (role("cap-linked transcript cDNA"), role("downstream transcript cDNA")) if key == "cagescan" else (role("assay-derived cDNA"),)
    final = role_library(d.title + " final library", payload,
                         indexed=key != "snu-seq",
                         platform="platform-specific" if key == "snu-seq" else "Illumina-compatible")
    unavailable = ("The authoritative paper establishes the library roles and read direction, but does not expose one complete immutable sequencing-primer geometry; adapters are therefore shown by role."
                   if key != "cagescan" else
                   "The original CAGEscan supplement publishes the RT and PCR oligos, but the historical instrument-primer geometry is not asserted as a modern invariant.")
    spec = SimpleNamespace(
        TITLE=d.title, NOTES=d.notes,
        SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
        SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final, SEQ_PRIMERS=(),
        FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint, READ_LENGTHS=None,
        SEQUENCING_ENDING=d.endpoint, SEQUENCING_UNAVAILABLE=unavailable,
        workflow=workflow_fn, oligos=(lambda: CAGE_OLIGOS if key == "cagescan" else ()),
    )
    if key == "mnet-seq":
        spec.protocol_graph = lambda: _mnet_graph(d, final)
    elif key == "start-seq":
        spec.protocol_graph = lambda: _start_graph(d, final)
    elif key in {"gro-cap", "pro-cap"}:
        spec.protocol_graph = lambda: _runon_cap_graph(key, d, final)
    elif key in {"timelapse-seq", "tuc-seq", "snu-seq", "sf4su-seq"}:
        spec.protocol_graph = lambda: _four_su_graph(key, d, final)
    else:
        spec.protocol_graph = lambda: _cagescan_graph(d, final)
    return spec
