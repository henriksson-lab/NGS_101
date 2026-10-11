"""Reusable models for R-loop capture and differential cytosine readout protocols."""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import seqprimers as sp
from base_conversion import ace_path, oxbs_path, tab_path
from batch_ngs import seg, truseq_library
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow


def _part(name: str, symbol: str = "X", n: int = 24, **kw) -> Segment:
    return Segment(name, symbol * n, placeholder=True, **kw)


@dataclass(frozen=True)
class RLoopCapture:
    """A three-stranded R-loop plus a capture route and sequenced component."""
    probe: str
    sequenced_component: str
    in_vivo: bool = False

    def __post_init__(self) -> None:
        if self.sequenced_component not in {"DNA", "RNA"}:
            raise ValueError("an R-loop capture must sequence either its DNA or RNA component")
        if not self.probe.strip():
            raise ValueError("an R-loop capture requires an affinity probe")

    def rloop(self) -> Scene:
        sc = Scene()
        sc.strand("template DNA", [_part("DNA in RNA:DNA hybrid")], label="R-loop")
        sc.anneal("hybrid RNA", [_part("RNA in RNA:DNA hybrid", "x")],
                  to="template DNA", pair=("RNA in RNA:DNA hybrid", "DNA in RNA:DNA hybrid"),
                  label="R-loop")
        sc.strand("displaced DNA", [_part("displaced single-stranded DNA", "Y")],
                  label="displaced DNA")
        sc.note("hybrid RNA", f"captured by {self.probe}")
        return sc

    def captured_component(self) -> Scene:
        sc = Scene()
        sc.strand("selected molecule", [_part(f"R-loop {self.sequenced_component} component")],
                  label=f"selected {self.sequenced_component}")
        sc.note("selected molecule", f"{self.probe}; " +
                ("captured in fixed cells" if self.in_vivo else "captured after gentle extraction"))
        return sc


def _workflow(initial: Scene, steps: tuple[tuple[str, Scene, str], ...]) -> Workflow:
    w = Workflow(MolecularState("Starting material", tuple(initial.rows())))
    for action, scene, note in steps:
        w.react(action, scene.rows(), name=action, note=note)
    return w


def _role_library(name: str, payload: tuple[str, ...]) -> Construct:
    return Construct([_part("Illumina-compatible left adapter", "L", 18),
                      *[_part(x, "X" if i % 2 == 0 else "Y", 26)
                        for i, x in enumerate(payload)],
                      _part("Illumina-compatible right adapter", "R", 18)], name=name)


def _conversion_scene(states: tuple[tuple[str, str], ...], label: str) -> Scene:
    sc = Scene()
    for i, (origin, state) in enumerate(states):
        sc.strand(origin, [Segment(state, state, "me" if "m" in state else None,
                                   placeholder=True, length_bp=1)], label=f"{origin} → {state}", at=0)
    return sc


def _conversion_workflow(key: str) -> Workflow:
    initial = _conversion_scene((("original C", "C"), ("original 5mC", "5mC"),
                                 ("original 5hmC", "5hmC")), "genomic cytosines")
    if key == "oxbs-seq":
        tagged = _role_library_scene("adapter-ligated oxBS precursor", ("unconverted genomic insert",))
        oxidized = _conversion_scene((("original C", "C"), ("original 5mC", "5mC"),
                                      ("original 5hmC", "5fC")), "after oxidation")
        readout = _conversion_scene(tuple((b, oxbs_path(b).states[-1])
                                          for b in ("C", "5mC", "5hmC")), "PCR readout")
        return _workflow(initial, (
            ("Sonicate, end-repair and add one 3-prime dA to each end",
             _role_library_scene("dA-tailed oxBS insert", ("unconverted genomic insert",), ligated=False),
             "The dA-tailed insert is the substrate for the following adapter ligation."),
            ("Ligate TruSeq adapters", tagged,
             "The Nature Protocols workflow installs the library adapters before oxidation and bisulfite."),
            ("Oxidize 5hmC selectively with potassium perruthenate", oxidized,
             "5hmC becomes 5fC; 5mC remains protected."),
            ("Bisulfite-convert, then amplify with uracil-tolerant polymerase", readout,
             "C and oxidized 5hmC read as T; only 5mC reads as C."),
        ))
    if key == "tab-seq":
        protected = _conversion_scene((("original C", "C"), ("original 5mC", "5mC"),
                                       ("original 5hmC", "5gmC")), "after β-GT")
        oxidized = _conversion_scene((("original C", "C"), ("original 5mC", "5caC"),
                                      ("original 5hmC", "5gmC")), "after mTet1")
        readout = _conversion_scene(tuple((b, tab_path(b).states[-1])
                                          for b in ("C", "5mC", "5hmC")), "PCR readout")
        return _workflow(initial, (
            ("Sonicate genomic DNA and add conversion-control molecules", initial,
             "C, 5mC and 5hmC spike-ins measure conversion, oxidation and protection."),
            ("Glucosylate 5hmC with T4 β-glucosyltransferase", protected,
             "5gmC is protected from the following TET oxidation."),
            ("Oxidize 5mC to 5caC with excess mTet1", oxidized,
             "The protected 5hmC branch remains 5gmC."),
            ("Prepare the Illumina WGBS library, bisulfite-convert and amplify", readout,
             "Only original 5hmC reads as C; C and original 5mC read as T."),
        ))
    protected = _conversion_scene((("original C", "C"), ("original 5mC", "5mC"),
                                   ("original 5hmC", "5gmC")), "after β-GT")
    readout = _conversion_scene(tuple((b, ace_path(b).states[-1])
                                      for b in ("C", "5mC", "5hmC")), "PCR readout")
    return _workflow(initial, (
        ("Shear genomic DNA and add C/5mC/5hmC controls", initial,
         "The defining whole-genome experiment used approximately 300-bp input fragments."),
        ("Glucosylate 5hmC with T4 β-glucosyltransferase", protected,
         "Glucosylation strengthens the resistance of 5hmC to APOBEC3A."),
        ("Denature, then deaminate C and 5mC with APOBEC3A", readout,
         "C and 5mC read as T; protected 5hmC remains C without bisulfite."),
        ("Install Accel Methyl-NGS adapters and index by PCR",
         _role_library_scene("ACE-seq library", ("ACE-treated insert",), ligated=False),
         "Adapter construction follows deamination; the commercial oligo sequences are not public."),
    ))


def _hybrid_workflow(key: str) -> Workflow:
    if key == "drip-seq":
        cap = RLoopCapture("S9.6 antibody", "DNA")
        return _workflow(cap.rloop(), (
            ("Gently extract nucleic acid and restriction-digest genomic DNA", cap.rloop(),
             "The five-enzyme cocktail preserves RNA:DNA hybrids while reducing fragment size."),
            ("Immunoprecipitate RNA:DNA hybrids with S9.6", cap.captured_component(),
             "An RNase H-treated aliquot is the hybrid-specificity control."),
            ("Sonicate captured DNA, end-repair and add one 3-prime dA to each end",
             _role_library_scene("dA-tailed captured DNA", ("captured DNA fragment",), ligated=False),
             "The dA-tailed captured fragment is shown before adapter ligation."),
            ("Ligate adapters to the dA-tailed captured DNA",
             _role_library_scene("non-strand-specific DRIP-seq library", ("captured DNA fragment",)),
             "The DNA component is sequenced; the original implementation used a standard dsDNA library."),
        ))
    if key == "dripc-seq":
        cap = RLoopCapture("S9.6 antibody", "RNA")
        return _workflow(cap.rloop(), (
            ("Gently extract, restriction-digest and immunoprecipitate with S9.6", cap.rloop(),
             "The same hybrid capture precedes a different molecular readout."),
            ("Digest DNA and reverse-transcribe the retained R-loop RNA", cap.captured_component(),
             "Poly(dT) and random hexamers copy the RNA component into cDNA."),
            ("Make the second strand with dUTP, sonicate and ligate indexed adapters",
             _role_library_scene("dUTP-marked DRIPc-seq library", ("R-loop RNA-derived cDNA", "dUTP strand")),
             "UNG removes the dUTP-containing strand before PCR, preserving RNA strand identity."),
        ))
    cap = RLoopCapture("V5-tagged catalytically inactive RNASEH1 D210N", "DNA", in_vivo=True)
    primed = Scene(); primed.strand("copied molecule", [
        _part("N9-primer handle"), _part("N9 random primer", "N", 9),
        _part("R-loop DNA-derived copy")], label="unidirectional R-ChIP copy")
    primed.note("copied molecule", "the 5′ primer overhang blocks adapter ligation at this end")
    return _workflow(cap.rloop(), (
        ("Crosslink cells, shear chromatin and immunoprecipitate V5–RNASEH1", cap.captured_component(),
         "The D210N mutant binds R-loops; the binding-defective WKKD mutant is the control."),
        ("Denature and copy captured DNA with the tailed N9 primer", primed,
         "Random priming converts the DNA component to dsDNA and creates a blocking 5′ overhang."),
        ("Add one 3-prime dA to the unobstructed end", primed,
         "The tailed N9-primer end remains blocked; only the opposite end becomes ligation-ready."),
        ("Ligate one Illumina adapter and PCR-amplify",
         _rchip_library_scene(),
         "Only the unobstructed end accepts the adapter, so read orientation is retained."),
    ))


def _role_library_scene(name: str, payload: tuple[str, ...], *, ligated: bool = True) -> Scene:
    lib = _role_library(name, payload)
    sc = Scene.duplex(list(lib), label=name)
    if ligated:
        sc.junction("top", "Illumina-compatible left adapter", payload[0], "adapter ligation")
        sc.junction("top", payload[-1], "Illumina-compatible right adapter", "adapter ligation")
    return sc


def _rchip_library() -> Construct:
    return Construct([
        Segment("tailed N9 handle", "CAAGCAGAAGACGGCATACGAG", "p7"),
        _part("N9 random-primer remnant", "N", 9),
        _part("R-loop DNA-derived insert"),
        _part("one ligated Illumina adapter", "R", 18),
    ], name="strand-specific R-ChIP library")


def _rchip_library_scene() -> Scene:
    sc = Scene.duplex(list(_rchip_library()), label="strand-specific R-ChIP library")
    sc.junction("top", "R-loop DNA-derived insert", "one ligated Illumina adapter",
                "adapter ligation")
    return sc


@dataclass(frozen=True)
class Definition:
    title: str
    doi: str
    source: str
    summary: str
    payload: tuple[str, ...]
    sequencing: str
    notes: str
    truseq: bool = False
    caveat: str = ""


DEFS = {
    "drip-seq": Definition("DRIP-seq — DNA:RNA immunoprecipitation sequencing",
        "10.1016/j.molcel.2012.01.017", "Ginno et al., Molecular Cell (2012)",
        "Gentle extraction and restriction digestion preserve R-loops for S9.6 immunoprecipitation; a conventional dsDNA library sequences the captured DNA component.",
        ("captured R-loop DNA",), "Single-end Illumina reads report captured restriction fragments without strand identity.", "01_drip_seq.html",
        caveat="The defining paper names the library reagents but does not expose complete adapter or sequencing-primer sequences; adapter regions are shown by role."),
    "dripc-seq": Definition("DRIPc-seq — strand-specific sequencing of R-loop RNA",
        "10.1016/j.molcel.2016.05.032", "Sanz et al., Molecular Cell (2016)",
        "After S9.6 capture, DNase releases the R-loop RNA for reverse transcription; dUTP second-strand synthesis and UNG make the final library strand-specific.",
        ("R-loop RNA-derived cDNA",), "Single-end reads derive from the retained cDNA strand and therefore identify the RNA strand in the original hybrid.", "01_dripc_seq.html",
        caveat="The paper and detailed protocol specify indexed Illumina construction but do not print complete adapter sequences; adapter regions are shown by role."),
    "r-chip": Definition("R-ChIP — in-vivo R-loop capture with inactive RNASEH1",
        "10.1016/j.molcel.2017.10.008", "Chen et al., Molecular Cell (2017)",
        "V5-tagged RNASEH1 D210N binds R-loops in fixed cells; a tailed N9 primer and one-ended adapter ligation produce a strand-specific DNA-component library.",
        ("tailed N9 end", "R-loop DNA-derived insert"), "Single-end Illumina reads preserve the orientation of the DNA strand captured from each R-loop.", "01_r_chip.html",
        caveat="The defining paper prints the tailed N9 primer, but identifies the opposite adapter only as an Illumina adapter; that adapter is shown by role."),
    "oxbs-seq": Definition("oxBS-seq — oxidative bisulfite sequencing",
        "10.1126/science.1220671", "Booth et al., Science (2012)",
        "Potassium perruthenate converts 5hmC to 5fC before bisulfite, so 5hmC and unmodified C read as T while 5mC alone reads as C.",
        ("oxidative-bisulfite-converted insert",), "TruSeq reads interrogate the converted insert; an i7 read identifies the library.", "01_oxbs_seq.html", truseq=True),
    "tab-seq": Definition("TAB-seq — TET-assisted bisulfite sequencing",
        "10.1016/j.cell.2012.04.027", "Yu et al., Cell (2012)",
        "β-GT protects 5hmC as 5gmC, then mTet1 oxidizes 5mC to 5caC; after bisulfite, only original 5hmC reads as C.",
        ("TAB-treated insert",), "TruSeq reads interrogate the converted insert; an i7 read identifies the library.", "01_tab_seq.html", truseq=True),
    "ace-seq": Definition("ACE-seq — APOBEC-coupled epigenetic sequencing",
        "10.1038/nbt.4204", "Schutsky et al., Nature Biotechnology (2018)",
        "β-GT protects 5hmC before APOBEC3A deaminates C and 5mC, giving a bisulfite-free readout in which only original 5hmC remains C.",
        ("ACE-treated insert",), "The defining experiment used single-end NextSeq sequencing of Accel Methyl-NGS libraries.", "01_ace_seq.html",
        caveat="Accel Methyl-NGS adapter oligos are proprietary and are shown by role; no exact sequencing-primer sequence is asserted."),
}


def get(key: str):
    d = DEFS[key]
    if d.truseq:
        final, primers = truseq_library([seg(d.payload[0], "X" * 34, placeholder=True)],
                                        d.title + " final library", dual_index=False)
    else:
        final, primers = ((_rchip_library() if key == "r-chip" else
                           _role_library(d.title + " final library", d.payload)), ())
    return SimpleNamespace(
        TITLE=d.title, NOTES=d.notes,
        SOURCE=(f'Defining source: <a href="https://doi.org/{d.doi}">{d.source}</a>.' +
                (' Detailed chemistry was checked against the authors’ later protocol.'
                 if key in {"drip-seq", "dripc-seq", "oxbs-seq", "tab-seq", "r-chip"} else '')),
        SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final, SEQ_PRIMERS=primers,
        FINAL_CAPTION=d.sequencing, SEQUENCING_INTRO=d.sequencing, READ_LENGTHS=None,
        SEQUENCING_ENDING=d.sequencing,
        SEQUENCING_UNAVAILABLE=("The authoritative source does not expose the exact sequencing-primer sequence for this library implementation."
                                if not primers else ""),
        workflow=lambda: (_conversion_workflow(key) if key in {"oxbs-seq", "tab-seq", "ace-seq"}
                          else _hybrid_workflow(key)),
        oligos=(lambda: ('5′-/invddt/CAAGCAGAAGACGGCATACGAGNNNNNNNNN-3′ — tailed N9 random primer',)
                if key == "r-chip" else ()),
    )
