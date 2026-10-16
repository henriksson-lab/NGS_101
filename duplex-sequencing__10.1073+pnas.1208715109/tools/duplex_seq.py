"""Schmitt2012 Duplex Sequencing, using the March2013-corrected original SI."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceCoverage, SourceItem, StepKind)

TITLE = "Duplex Sequencing — complementary-strand consensus"
NOTES = "01_duplex-sequencing.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1073/pnas.1208715109">Schmitt et al., '
          '<i>PNAS</i> (2012)</a>; original Supporting Information corrected March 27, 2013.')
SUMMARY = "Complementary random duplex tags identify both strands of each original DNA molecule. Strand families are paired by reciprocal end tags to form duplex consensuses after sequencing."
CAVEAT = ("This page uses the original Supporting Information explicitly labelled ‘Corrected March 27, 2013’, recovered in a public USPTO exhibit of the article and SI. The uncorrected SI has not been compared, so these oligos are attributed to that corrected source version. It supplies dA-tailed adapters, dT-tailed fragments and a six-position multiplexing index. The later 2014 redesign has different tail polarity and is a separate chemistry. The corrected SI explicitly prints the five-base CAGTA read spacer; it also follows from copying the template’s ACTG and adding adapter dA. The adapter-copy product contains blunt and dA-tailed molecules; a separate dATP-only treatment completes A-tailing. The physical pre-PCR Y-shaped molecule is shown separately from the duplex copied during PCR.")
TAG_ALPHA = feature("duplex_tag_alpha", "umi", "random", group="duplex_tag", part="alpha")
TAG_BETA = feature("duplex_tag_beta", "umi", "random", group="duplex_tag", part="beta")
TAG_ADAPTER = feature("duplex_tag_adapter", "umi", "random")
SAMPLE_INDEX = feature("sample_index", "sample_index", "unknown")

# Exact factorizations checked against the corrected original SI in selftest.py.
PRIMER_STRAND = il.P5[:-4] + il.TRUSEQ_READ1
TEMPLATE_FIXED = "ACTG"
TEMPLATE_STRAND = TEMPLATE_FIXED + "N" * 12 + revcomp(il.TRUSEQ_READ2)
PCR_FORWARD = il.P5[:21]
PCR_REVERSE = il.P7 + "N" * 6 + il.TRUSEQ_READ2[:23]
COPIED_FIXED = revcomp(TEMPLATE_FIXED)
READ_SPACER = COPIED_FIXED + "A"
# Exact shared overlap of the two source oligos; no complete-arm complement assumed.
CORE = il.TRUSEQ_READ1[-13:]
if CORE != il.TRUSEQ_READ2[-13:]:
    raise ValueError("source adapter oligos do not share the declared terminal overlap")


def tag(name, semantic):
    return Segment(name, "N" * 12, "umi", placeholder=True, feature=semantic)


def roles(name):
    return Construct([Segment(name, "[" + name + "]", placeholder=True)], name=name)


def adapter_scene(*, copied=False, a_tailed=False):
    if a_tailed and not copied:
        raise ValueError("adapter A-tailing follows copying the random template")
    top = [Segment("primer Y arm", PRIMER_STRAND[:-13], "p5"), Segment("paired core", CORE)]
    if copied:
        top += [tag("copied duplex tag", TAG_ADAPTER), Segment("copied ACTG", COPIED_FIXED)]
    if a_tailed:
        top += [Segment("adapter 3-prime dA", "A")]
    bottom = [Segment("template ACTG", TEMPLATE_FIXED), tag("template duplex tag", TAG_ADAPTER),
              Segment("core complement", revcomp(CORE)),
              Segment("template Y arm", revcomp(il.TRUSEQ_READ2)[13:], "r2")]
    scene = Scene()
    scene.strand("primer", top, label="primer-derived strand")
    scene.anneal("template", bottom, to="primer", pair=("core complement", "paired core"),
                 unpaired=("template Y arm",), label="source template strand", mod5="/5phos/")
    if a_tailed:
        scene.overhang("primer", "adapter 3-prime dA", polarity="3-prime", owner="adapter",
                       text="unpaired 3′ dA; added after copying")
    scene.labels("primer")
    return scene


def ligation_geometry():
    substrate = Scene()
    substrate.strand("top", [Segment("fragment top", "X" * 24, placeholder=True),
                              Segment("fragment 3-prime dT", "T")], label="end-repaired fragment")
    substrate.anneal("bottom", [Segment("fragment bottom", "X" * 24, placeholder=True)],
                     to="top", pair=("fragment bottom", "fragment top"), label="fragment complement")
    substrate.overhang("top", "fragment 3-prime dT", polarity="3-prime", owner="substrate",
                       text="unpaired fragment 3′ dT")
    return AdapterLigationGeometry.known("cohesive", adapter_scene(copied=True, a_tailed=True), substrate,
        statement="Corrected original SI ligates copied, separately dA-tailed Y adapters to end-repaired, dT-tailed DNA with T4 DNA ligase. Both one-base tails are unpaired in the ligation-ready substrates.")


def central_segments():
    return [tag("duplex tag α", TAG_ALPHA), Segment("left fixed read spacer", READ_SPACER),
            Segment("genomic insert", "X" * 38, placeholder=True),
            Segment("right fixed junction", revcomp(READ_SPACER)), tag("duplex tag β", TAG_BETA)]


def ligated_y_scene():
    """The two original DNA strands retain different non-complementary Y arms."""
    middle = central_segments()
    top = [Segment("left primer Y arm", PRIMER_STRAND[:-13]),
           Segment("left paired core", CORE), *middle,
           Segment("right paired core", revcomp(CORE)),
           Segment("right template Y arm", revcomp(il.TRUSEQ_READ2)[13:])]
    bottom = [Segment("opposite primer Y arm", PRIMER_STRAND[:-13]),
              Segment("opposite left paired core", CORE), *complement_segments(middle),
              Segment("opposite right paired core", revcomp(CORE)),
              Segment("opposite template Y arm", revcomp(il.TRUSEQ_READ2)[13:])]
    scene = Scene()
    scene.strand("top", top, label="one original strand and adapters")
    scene.anneal("bottom", bottom, to="top", pair=("genomic insert'", "genomic insert"),
                 unpaired=("opposite primer Y arm", "opposite template Y arm"),
                 label="complementary original strand and adapters")
    scene.labels("top")
    return scene


# The original top strand provides the template for a reverse-primer first copy.
# The copied duplex core acquires the P5-complement site used by the forward primer.
PCR_TEMPLATE = Construct([Segment("source primer-derived left arm", PRIMER_STRAND),
                          *central_segments(), Segment("source template-derived right arm", revcomp(il.TRUSEQ_READ2))],
                         name="Duplex core after the first PCR copy; original Y molecule shown separately")
PCR_GEOMETRY = PCRGeometry.exact(PCR_TEMPLATE,
    PCRPrimer("source 21-base forward primer", PCR_FORWARD, anneal_nt=21),
    PCRPrimer("source six-index reverse primer", PCR_REVERSE, anneal_nt=23),
    statement="The source reverse primer has a 23-base site on the original template-derived right arm. Its first extension copies the left primer arm, providing the inward-facing 21-base forward-primer site. This panel shows that copied duplex core; the reverse primer's P7/index overhang is shown separately. It does not make the pre-PCR Y arms complementary.")
FINAL_LIBRARY = Construct([*PCR_TEMPLATE.segments,
                           Segment("six-position sample index; reverse-complement orientation", "N" * 6, "i7",
                                   placeholder=True, feature=SAMPLE_INDEX),
                           Segment("P7 complement", revcomp(il.P7), "p7")],
                          name="Corrected-original Duplex Sequencing amplified library")
FINAL_CAPTION = "Each tag is followed in its sequencing read by the derived five-base CAGTA spacer, then genomic DNA. The right top-strand junction TACTG is the reverse complement of CAGTA. α/β identify the complementary tag pair; the six-position sample index is a separate feature."
SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
SEQUENCING_INTRO = "Corrected SI describes HiSeq 2000 sequencing without supplying a cycle count here. Both genomic reads begin with a 12-base tag followed by five fixed bases at cycles 13–17. Processing removes these and the next four genomic bases; that filtering is separate from physical library preparation."
SHOW_READ_CONTENT = True


def oligos():
    return [f"Adapter primer strand: 5′-{PRIMER_STRAND}-3′",
            f"Adapter template strand: 5′-/5phos/{TEMPLATE_STRAND}-3′",
            f"PCR forward: 5′-{PCR_FORWARD}-3′",
            f"PCR reverse family (source XXXXXX = six-position fixed multiplexing index): 5′-{il.P7}[XXXXXX]{il.TRUSEQ_READ2[:23]}-3′"]


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    dna = g.input("dna", "Sheared genomic DNA", MaterialKind.MOLECULAR, primary=True)
    end_repair = g.input("end_repair_kit", "NEB DNA End Repair kit", MaterialKind.CATALYST)
    repaired = g.product("repaired", "End-repaired DNA", MaterialKind.MOLECULAR)
    g.react("end_repair", "Repair sheared DNA ends", StepKind.END_REPAIR,
            ((dna, InputRole.TEMPLATE), (end_repair, InputRole.CATALYST)), (repaired,))
    size_beads = g.input("size_beads", "AMPure XP beads", MaterialKind.SELECTION)
    selected = g.product("selected", "Size-selected fragments", MaterialKind.MOLECULAR)
    size_waste = g.product("size_waste", "Fragments excluded by source size selection", MaterialKind.WASTE)
    g.react("size_select", "Select the source fragment-size interval", StepKind.SELECTION,
            ((repaired, InputRole.TEMPLATE), (size_beads, InputRole.SELECTION)), (selected, size_waste))
    g.discard(size_waste, "Source bead selection physically excludes large and small fragments; this is not consensus filtering.")
    klenow = g.input("klenow", "Klenow exo− DNA polymerase", MaterialKind.CATALYST)
    dttp = g.input("dttp", "dTTP-only tailing condition", MaterialKind.CONDITION)
    t_tailed = g.product("t_tailed", "dT-tailed genomic duplexes", MaterialKind.MOLECULAR)
    g.react("tail_insert", "Add fragment 3′ dT tails", StepKind.TAILING,
            ((selected, InputRole.TEMPLATE), (klenow, InputRole.CATALYST), (dttp, InputRole.CONDITION)), (t_tailed,))
    adapter_primer = g.input("adapter_primer", "Source adapter primer strand", MaterialKind.OLIGO,
                            state=Construct([Segment("adapter primer strand", PRIMER_STRAND)]))
    adapter_template = g.input("adapter_template", "Source 5′-phosphorylated adapter template", MaterialKind.OLIGO,
        state=Construct([Segment("ACTG", TEMPLATE_FIXED), tag("random adapter tag", TAG_ADAPTER),
                         Segment("template read2 arm", revcomp(il.TRUSEQ_READ2))]))
    annealed = g.product("annealed_adapter", "Annealed source Y-adapter oligos", MaterialKind.MOLECULAR,
                         state=adapter_scene().rows())
    g.react("anneal_adapter", "Anneal the two source adapter oligos", StepKind.HYBRIDIZATION,
            ((adapter_template, InputRole.TEMPLATE), (adapter_primer, InputRole.STRUCTURAL)), (annealed,))
    dntps = g.input("dntps", "dNTPs for adapter copying", MaterialKind.CONDITION)
    copied = g.product("copied_adapter", "Copied duplex-tag adapter; blunt/dA mixture", MaterialKind.MOLECULAR,
                       state=adapter_scene(copied=True).rows())
    g.react("copy_adapter", "Copy the template's random tag and ACTG end", StepKind.EXTENSION,
            ((annealed, InputRole.TEMPLATE), (klenow, InputRole.CATALYST), (dntps, InputRole.CONDITION)), (copied,))
    datp = g.input("datp", "dATP-only adapter-tail condition", MaterialKind.CONDITION)
    ready_adapter = g.product("tag_adapter", "dA-tailed Y adapter", MaterialKind.ADAPTER,
                              state=adapter_scene(copied=True, a_tailed=True).rows())
    g.react("tail_adapter", "A-tail residual blunt adapter molecules", StepKind.TAILING,
            ((copied, InputRole.TEMPLATE), (klenow, InputRole.CATALYST), (datp, InputRole.CONDITION)), (ready_adapter,))
    ligase = g.input("ligase", "T4 DNA ligase", MaterialKind.CATALYST)
    tagged = g.product("tagged", "Original duplex with ligated Y adapters", MaterialKind.MOLECULAR,
                       state=ligated_y_scene().rows())
    g.react("ligate_tags", "Ligate dA adapters to dT-tailed fragment ends", StepKind.LIGATION,
            ((t_tailed, InputRole.TEMPLATE), (ready_adapter, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)),
            (tagged,), geometry=ligation_geometry())
    forward = g.input("pcr_forward", "Source 21-base PCR forward primer", MaterialKind.OLIGO,
                       state=Construct([Segment("PCR forward", PCR_FORWARD)]))
    reverse = g.input("pcr_reverse", "Source six-index PCR reverse primer", MaterialKind.OLIGO,
        state=Construct([Segment("P7", il.P7), Segment("six-position multiplexing index", "N" * 6,
                         placeholder=True, feature=SAMPLE_INDEX), Segment("source 23-base read2 arm", il.TRUSEQ_READ2[:23])]))
    kapa = g.input("kapa", "KAPA HiFi PCR polymerase", MaterialKind.CATALYST)
    library = g.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    g.react("library_pcr", "PCR-copy the Y arms and amplify strand families", StepKind.PCR,
            ((tagged, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
             (reverse, InputRole.STRUCTURAL), (kapa, InputRole.CATALYST)),
            (library,), geometry=PCR_GEOMETRY)
    g.sequence(library, SequencingPlan(library=FINAL_LIBRARY, primers=SEQ_PRIMERS))
    g.declare_reviewed_scope(("dna", "end_repair_kit", "end_repair", "size_beads", "size_select",
        "klenow", "dttp", "tail_insert", "adapter_primer", "adapter_template", "anneal_adapter",
        "dntps", "copy_adapter", "datp", "tail_adapter", "ligase", "ligate_tags", "pcr_forward",
        "pcr_reverse", "kapa", "library_pcr", "library"), SOURCE)
    g.add_source_item(SourceItem("consensus", "Tag/spacer filtering, SSCS and DCS computation", SOURCE))
    g.cover_source_item(SourceCoverage("consensus", "excluded",
        reason="Computational processing removes read bases and rejects discordant observations; it does not physically remove library molecules."))
    return g.require_publishable()
