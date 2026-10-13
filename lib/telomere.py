"""Shared structures for chromosome-end tagging and telomere capture workflows."""
from __future__ import annotations

from dataclasses import dataclass

from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, complement_segments, revcomp
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role

HUMAN_G_REPEAT = "TTAGGG"
HUMAN_C_REPEAT = revcomp(HUMAN_G_REPEAT)


def phases(sequence: str) -> tuple[str, ...]:
    """Every cyclic phase of a repeat, in deterministic order."""
    sequence = sequence.upper()
    if not sequence or set(sequence) - set("ACGT"):
        raise ValueError("repeat must be non-empty DNA")
    return tuple(sequence[i:] + sequence[:i] for i in range(len(sequence)))


HUMAN_G_PHASES = phases(HUMAN_G_REPEAT)
HUMAN_C_PHASES = phases(HUMAN_C_REPEAT)


@dataclass(frozen=True)
class TelomereEnd:
    """A duplex telomere ending in a native G-rich 3-prime overhang."""

    repeat: str = HUMAN_G_REPEAT
    duplex_copies: int = 3
    overhang_copies: int = 2
    subtelomere_bases: int = 16

    def __post_init__(self) -> None:
        if self.duplex_copies < 1 or self.overhang_copies < 1:
            raise ValueError("telomere drawing needs duplex repeat and overhang")
        if self.subtelomere_bases < 1:
            raise ValueError("telomere drawing needs a subtelomeric anchor")

    @property
    def duplex(self) -> str:
        return self.repeat * self.duplex_copies

    @property
    def overhang(self) -> str:
        return self.repeat * self.overhang_copies

    def scene(self) -> Scene:
        top = [Segment("subtelomere", "X" * self.subtelomere_bases, placeholder=True),
               Segment("duplex telomere", self.duplex, "r1"),
               Segment("G-rich 3-prime overhang", self.overhang, "r2")]
        bottom = complement_segments(top[:2])
        sc = Scene(); sc.strand("G-rich strand", top, label="chromosome")
        sc.anneal("C-rich strand", bottom, to="G-rich strand",
                  pair=("duplex telomere'", "duplex telomere"), label="chromosome")
        sc.mark("G-rich strand", "G-rich 3-prime overhang", "native terminal overhang")
        return sc


@dataclass(frozen=True)
class PhasedTelomereAdapter:
    """Adapter family whose terminal arms cover every repeat phase."""

    core: str
    arm_repeat: str = HUMAN_C_REPEAT
    arm_copies: int = 3
    phosphorylated_5: bool = True
    core_segments: tuple[Segment, ...] = ()

    def __post_init__(self) -> None:
        if not self.core or set(self.core.upper()) - set("ACGTNX"):
            raise ValueError("adapter core must be DNA or explicit N/X placeholders")
        if not self.arm_repeat or set(self.arm_repeat.upper()) - set("ACGT"):
            raise ValueError("telomere arm repeat must be real DNA")
        if self.arm_copies < 1:
            raise ValueError("adapter needs at least one repeat copy")
        if self.core_segments and "".join(s.top for s in self.core_segments) != self.core:
            raise ValueError("structured adapter-core segments must exactly rebuild the core")

    @property
    def arms(self) -> tuple[str, ...]:
        return tuple((p + self.arm_repeat * self.arm_copies)[:
                     len(self.arm_repeat) * self.arm_copies]
                     for p in phases(self.arm_repeat))

    def annealed_scene(self, end: TelomereEnd | None = None, phase: int = 0) -> Scene:
        end = end or TelomereEnd()
        arm = self.arms[phase]
        top = [Segment("subtelomere", "X" * end.subtelomere_bases, placeholder=True),
               Segment("duplex telomere", end.duplex, "r1"),
               Segment("G-rich overhang", end.overhang, "r2")]
        bottom = complement_segments(top[:2])
        sc = Scene(); sc.strand("G-rich strand", top, label="chromosome")
        sc.anneal("C-rich strand", bottom, to="G-rich strand",
                  pair=("duplex telomere'", "duplex telomere"), label="chromosome")
        core = list(self.core_segments) or [
            Segment("adapter core", self.core, "r3",
                    placeholder="N" in self.core or "X" in self.core)
        ]
        adapter = [*core, Segment("telomere-complement arm", arm, "r2")]
        sc.anneal("terminal adapter", adapter, to="G-rich strand",
                  pair=("telomere-complement arm", "G-rich overhang"),
                  shift=len(end.overhang) - len(arm), label="adapter",
                  unpaired=tuple(s.name for s in core),
                  mod5="p" if self.phosphorylated_5 else "")
        sc.ligation_end("terminal adapter", core[0].name, kind="cohesive",
                        polarity="5-prime", owner="adapter",
                        text="adapter 5-prime phosphate at native C-rich 3-prime end")
        sc.footer("** ligase seals the adapter 5′ end to the native C-rich 3′ end",
                  "terminal adapter", core[0].name)
        return sc


@dataclass(frozen=True)
class AffinityCapture:
    """A tagged molecule selected on a solid phase and released by a named cut."""

    affinity: str
    release_enzyme: str | None = None

    def steps(self) -> tuple[str, ...]:
        out = (f"bind {self.affinity}-tagged molecules", "wash away untagged DNA")
        return out + ((f"release with {self.release_enzyme}",) if self.release_enzyme else ())


def nanopore_telomere_graph(kind: str, title: str, authority: str) -> ProtocolGraph:
    """Typed provenance for the four terminal-telomere nanopore workflows."""
    if kind not in {"telotag", "telo_seq", "dtm", "profiling"}:
        raise ValueError(f"unknown telomere workflow {kind!r}")
    graph = ProtocolGraph(title, authority=authority)
    dna = graph.input("genomic_dna", "high-molecular-weight genomic DNA with native telomeres",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("chromosome ends", (role("subtelomere"), role("native telomere terminus"))),
                      primary=True)
    reviewed = ["genomic_dna"]
    current = dna

    if kind == "profiling":
        nuclease = graph.input("internal_nuclease", "internal genomic restriction enzyme",
                               MaterialKind.CATALYST)
        digested = graph.product("internally_digested", "restriction-shortened telomere-bearing DNA",
                                 MaterialKind.MOLECULAR, state=current.state)
        graph.react("digest_internal", "Digest genomic DNA away from native telomere ends",
                    StepKind.DIGESTION,
                    ((current, InputRole.TEMPLATE), (nuclease, InputRole.CATALYST)),
                    (digested,))
        reviewed.extend(("internal_nuclease", "digest_internal")); current = digested

    if kind == "telotag":
        tdt = graph.input("tdt", "terminal deoxynucleotidyl transferase",
                          MaterialKind.CATALYST)
        tailed = graph.product("polya_tailed", "poly(A)-tailed native chromosome ends",
                               MaterialKind.MOLECULAR, state=current.state)
        graph.react("polya_tail", "Add poly(A) at the native 3-prime chromosome end",
                    StepKind.TAILING,
                    ((current, InputRole.TEMPLATE), (tdt, InputRole.CATALYST)), (tailed,))
        reviewed.extend(("tdt", "polya_tail")); current = tailed

    tag = graph.input("terminal_tag", "phased telomere-capture or TeloTag oligo",
                      MaterialKind.OLIGO,
                      inferred=kind in {"dtm", "profiling"},
                      source=authority + (" (some ordered bases unavailable)" if kind in {"dtm", "profiling"} else ""))
    hybrid = graph.product("tag_hybrid", "capture oligo annealed at the native telomere end",
                           MaterialKind.MOLECULAR, state=current.state)
    graph.react("anneal_tag", "Anneal the terminal capture oligo",
                StepKind.HYBRIDIZATION,
                ((current, InputRole.TEMPLATE), (tag, InputRole.STRUCTURAL)), (hybrid,))
    reviewed.extend(("terminal_tag", "anneal_tag")); current = hybrid

    if kind != "telotag":
        terminal_ligase = graph.input("terminal_ligase", "terminal-tag DNA ligase",
                                      MaterialKind.CATALYST)
        tagged = graph.product("terminal_tagged", "covalently terminal-tagged telomere DNA",
                               MaterialKind.MOLECULAR, state=current.state)
        graph.react("ligate_terminal_tag", "Seal the tag to the native telomere terminus",
                    StepKind.LIGATION,
                    ((current, InputRole.TEMPLATE), (tag, InputRole.STRUCTURAL),
                     (terminal_ligase, InputRole.CATALYST)), (tagged,),
                    geometry=AdapterLigationGeometry.unavailable(
                        "The terminal telomere pairing is drawn from the phased oligo family; it is not a universal adapter-end geometry."))
        reviewed.extend(("terminal_ligase", "ligate_terminal_tag")); current = tagged

    if kind in {"telotag", "dtm"}:
        polymerase = graph.input("gap_fill_polymerase", "Sulfolobus DNA polymerase IV",
                                 MaterialKind.CATALYST)
        filled = graph.product("gap_filled", "gap-filled terminal-tagged telomere DNA",
                               MaterialKind.MOLECULAR, state=current.state)
        graph.react("gap_fill", "Fill the terminal tag–chromosome gap",
                    StepKind.EXTENSION,
                    ((current, InputRole.TEMPLATE), (polymerase, InputRole.CATALYST)),
                    (filled,))
        reviewed.extend(("gap_fill_polymerase", "gap_fill")); current = filled
        if kind == "telotag":
            nick_ligase = graph.input("nick_ligase", "T4 DNA ligase", MaterialKind.CATALYST)
            sealed = graph.product("sealed", "nick-sealed TeloTag chromosome end",
                                   MaterialKind.MOLECULAR, state=current.state)
            graph.react("seal_nick", "Seal the filled TeloTag–chromosome nick",
                        StepKind.NICK_SEALING,
                        ((current, InputRole.TEMPLATE), (nick_ligase, InputRole.CATALYST)),
                        (sealed,))
            reviewed.extend(("nick_ligase", "seal_nick")); current = sealed

    if kind == "telo_seq":
        ecov = graph.input("ecorv", "EcoRV", MaterialKind.CATALYST)
        cut = graph.product("ecorv_cut", "telorette-tagged DNA with distal blunt EcoRV ends",
                            MaterialKind.MOLECULAR, state=current.state)
        graph.react("digest_ecorv", "Cut distal genomic DNA with EcoRV", StepKind.DIGESTION,
                    ((current, InputRole.TEMPLATE), (ecov, InputRole.CATALYST)), (cut,))
        klenow = graph.input("klenow", "Klenow exo-minus dA-tailing enzyme",
                             MaterialKind.CATALYST)
        tailed = graph.product("distal_tailed", "dA-tailed internal EcoRV ends",
                               MaterialKind.MOLECULAR, state=cut.state)
        graph.react("tail_distal", "dA-tail only the new distal EcoRV ends",
                    StepKind.TAILING,
                    ((cut, InputRole.TEMPLATE), (klenow, InputRole.CATALYST)), (tailed,))
        reviewed.extend(("ecorv", "digest_ecorv", "klenow", "tail_distal")); current = tailed

    if kind == "dtm":
        shear = graph.input("fragmentation_condition", "mechanical fragmentation condition",
                            MaterialKind.CONDITION)
        fragmented = graph.product("fragmented", "mechanically fragmented tagged telomere DNA",
                                   MaterialKind.MOLECULAR, state=current.state)
        graph.react("fragment", "Mechanically fragment the tagged telomere DNA",
                    StepKind.FRAGMENTATION,
                    ((current, InputRole.TEMPLATE), (shear, InputRole.CONDITION)),
                    (fragmented,))
        reviewed.extend(("fragmentation_condition", "fragment")); current = fragmented

    if kind == "profiling":
        beads = graph.input("streptavidin", "streptavidin capture beads",
                            MaterialKind.SELECTION)
        captured = graph.product("captured", "biotin-TeloTag-enriched telomere DNA",
                                 MaterialKind.MOLECULAR, state=current.state)
        rejected = graph.product("uncaptured", "untagged genomic fragments", MaterialKind.WASTE)
        graph.react("capture", "Select biotinylated TeloTag molecules",
                    StepKind.AFFINITY_CAPTURE,
                    ((current, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                    (captured, rejected))
        graph.discard(rejected, "does not carry the biotinylated terminal tag")
        release = graph.input("release_nuclease", "adapter-site release nuclease",
                              MaterialKind.CATALYST)
        released = graph.product("released", "restriction-released tagged telomere molecules",
                                 MaterialKind.MOLECULAR, state=captured.state)
        graph.react("release", "Release captured molecules at the adapter restriction site",
                    StepKind.DIGESTION,
                    ((captured, InputRole.TEMPLATE), (release, InputRole.CATALYST)),
                    (released,))
        reviewed.extend(("streptavidin", "capture", "release_nuclease", "release")); current = released

    if kind == "telo_seq":
        splint = graph.input("s1_splint", "S1 motor-adapter recruitment splint",
                             MaterialKind.OLIGO)
        recruited = graph.product("motor_recruited", "S1-recruited telomere library substrate",
                                  MaterialKind.MOLECULAR, state=current.state)
        graph.react("recruit_motor", "Hybridize S1 to recruit the nanopore adapter",
                    StepKind.HYBRIDIZATION,
                    ((current, InputRole.TEMPLATE), (splint, InputRole.STRUCTURAL)),
                    (recruited,))
        reviewed.extend(("s1_splint", "recruit_motor")); current = recruited

    motor = graph.input("motor_adapter", "motor-loaded Oxford Nanopore sequencing adapter",
                        MaterialKind.ADAPTER, inferred=True,
                        source=authority + " (proprietary terminal bases unavailable)")
    motor_ligase = graph.input("motor_ligase", "Oxford Nanopore adapter ligase",
                               MaterialKind.CATALYST)
    final_state = Construct([
        Segment("motor-loaded ONT adapter", "[motor adapter]", placeholder=True, inferred=True),
        Segment("terminal-tagged telomere fragment", "[tagged telomere DNA]", placeholder=True),
    ], name=title + " nanopore library")
    library = graph.product("library", final_state.name, MaterialKind.LIBRARY,
                            state=final_state)
    graph.react("ligate_motor", "Attach the motor-loaded nanopore sequencing adapter",
                StepKind.LIGATION,
                ((current, InputRole.TEMPLATE), (motor, InputRole.STRUCTURAL),
                 (motor_ligase, InputRole.CATALYST)), (library,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The cited workflow does not disclose the proprietary motor-adapter terminal bases."))
    graph.sequence(library, SequencingPlan(
        acquired_roles=("Sequencing output",),
        limitation="Oxford Nanopore motor entry is platform-driven; no sequencing primer binds the library."))
    reviewed.extend(("motor_adapter", "motor_ligase", "ligate_motor", "library"))
    graph.declare_reviewed_scope(reviewed, authority)
    return graph.require_publishable()
