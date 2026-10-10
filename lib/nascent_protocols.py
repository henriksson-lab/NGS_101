"""Shared chemistry for nascent-RNA ends, cap selection and 4sU recoding.

The constructors encode the distinctions that matter experimentally: whether a molecule
was selected because it was capped, because a polymerase extended it in a run-on, or
because a metabolic label was enriched or chemically recoded.  Protocol pages only choose
a route; they cannot silently collapse these mechanisms into generic RNA-seq.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

from chemdraw import Construct, Scene, Segment, Workflow, feature
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


def role_library(name: str, payload: tuple[Segment, ...], *, indexed: bool = True) -> Construct:
    parts = [role("Illumina-compatible left end")]
    if indexed:
        parts.append(SAMPLE_INDEX)
    parts.extend(payload)
    parts.append(role("Illumina-compatible right end"))
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
            SAMPLE_INDEX, role("4sU-enriched RNA-derived cDNA"), role("right sequencing end"))),
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
    "5′-TAGTCGAACGAAGGTCTCCGAACCGCTCTTCCGATCTN6-3′ — CAGEscan RT primer",
    "5′-AATGATACGGCGACCACCGAGATCTACACTAGTCGAACTGAAGG-3′ — CAGEscan forward PCR primer",
    "5′-CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT-3′ — CAGEscan reverse PCR primer",
)


def get(key: str):
    d = DEFS[key]
    workflow_fn = (_cagescan if key == "cagescan" else
                   (lambda: _four_su(key)) if key in {"timelapse-seq", "tuc-seq", "snu-seq", "sf4su-seq"}
                   else (lambda: _end_capture(key)))
    payload = (role("cap-linked transcript cDNA"), role("downstream transcript cDNA")) if key == "cagescan" else (role("assay-derived cDNA"),)
    final = role_library(d.title + " final library", payload)
    unavailable = ("The authoritative paper establishes the library roles and read direction, but does not expose one complete immutable sequencing-primer geometry; adapters are therefore shown by role."
                   if key != "cagescan" else
                   "The original CAGEscan supplement publishes the RT and PCR oligos, but the historical instrument-primer geometry is not asserted as a modern invariant.")
    return SimpleNamespace(
        TITLE=d.title, NOTES=d.notes,
        SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
        SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final, SEQ_PRIMERS=(),
        FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint, READ_LENGTHS=None,
        SEQUENCING_ENDING=d.endpoint, SEQUENCING_UNAVAILABLE=unavailable,
        workflow=workflow_fn, oligos=(lambda: CAGE_OLIGOS if key == "cagescan" else ()),
    )
