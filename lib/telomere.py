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
        # A phase is a window into one infinitely repeated motif.  Prefixing a
        # phase to a fresh repeat would duplicate bases at the join (for example
        # CCTAAC + CCCTAA -> CCTAACCCCTAA).  Slice one long repeat instead.
        width = len(self.arm_repeat) * self.arm_copies
        repeated = self.arm_repeat * (self.arm_copies + 1)
        return tuple(repeated[offset:offset + width]
                     for offset in range(len(self.arm_repeat)))

    def annealed_scene(self, end: TelomereEnd | None = None, phase: int = 0) -> Scene:
        """Draw one phased arm annealed to the native G-rich 3-prime overhang.

        The arm pairs antiparallel with the overhang, so the adapter's 3-prime end
        (the arm) abuts the native C-rich strand's 5-prime end at the start of the
        overhang; that nick is the ligatable junction.  The adapter core, with its
        optional 5-prime phosphate, stays single-stranded beyond the chromosome
        terminus.  The default drawn overhang is as long as one arm.
        """
        arm = self.arms[phase]
        end = end or TelomereEnd(overhang_copies=self.arm_copies)
        if len(end.overhang) < len(arm):
            raise ValueError("the drawn G-rich overhang must be at least as long as the "
                             "adapter arm, or the arm would overlap the native C-rich strand")
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
                  shift=0, label="adapter",
                  unpaired=tuple(s.name for s in core),
                  mod5="p" if self.phosphorylated_5 else "")
        sc.ligation_end("terminal adapter", "telomere-complement arm", kind="splinted",
                        polarity="3-prime", owner="adapter",
                        text="adapter 3-prime end abuts the native C-rich strand 5-prime end")
        sc.footer("** ligase can seal the adapter 3′ end to the native C-rich 5′ end",
                  "terminal adapter", "telomere-complement arm")
        return sc


@dataclass(frozen=True)
class AffinityCapture:
    """A tagged molecule selected on a solid phase and released by a named cut."""

    affinity: str
    release_enzyme: str | None = None

    def steps(self) -> tuple[str, ...]:
        out = (f"bind {self.affinity}-tagged molecules", "wash away untagged DNA")
        return out + ((f"release with {self.release_enzyme}",) if self.release_enzyme else ())


def start_telomere_graph(graph: ProtocolGraph):
    return graph.input("genomic_dna", "high-molecular-weight genomic DNA with native telomeres",
                       MaterialKind.BIOLOGICAL,
                       state=molecule("chromosome ends", (role("subtelomere"),
                                                         role("native telomere terminus"))),
                       primary=True)


def anneal_terminal_tag(graph: ProtocolGraph, current, *, label: str):
    tag = graph.input("terminal_tag", label, MaterialKind.OLIGO, inferred=False)
    hybrid = graph.product("tag_hybrid", "capture oligo annealed at the native telomere end",
                           MaterialKind.MOLECULAR, state=current.state)
    graph.react("anneal_tag", "Anneal the terminal capture oligo", StepKind.HYBRIDIZATION,
                ((current, InputRole.TEMPLATE), (tag, InputRole.STRUCTURAL)), (hybrid,))
    return tag, hybrid


def fill_terminal_gap(graph: ProtocolGraph, current, *, polymerase: str,
                      source: str | None = None):
    """Gap-fill a terminal-tagged telomere; the caller names the source's polymerase."""
    if not polymerase.strip():
        raise ValueError("fill_terminal_gap needs the source-named polymerase")
    polymerase = graph.input("gap_fill_polymerase", polymerase,
                             MaterialKind.CATALYST, source=source)
    filled = graph.product("gap_filled", "gap-filled terminal-tagged telomere DNA",
                           MaterialKind.MOLECULAR, state=current.state)
    graph.react("gap_fill", "Fill the terminal tag–chromosome gap", StepKind.EXTENSION,
                ((current, InputRole.TEMPLATE), (polymerase, InputRole.CATALYST)), (filled,))
    return filled


def attach_nanopore_motor(graph: ProtocolGraph, current, *, title: str,
                          authority: str, reviewed: tuple[str, ...]):
    motor = graph.input("motor_adapter", "motor-loaded Oxford Nanopore sequencing adapter",
                        MaterialKind.ADAPTER, inferred=True,
                        source=authority + " (the complete kit- and revision-specific motor-bearing assembly is unavailable)")
    ligase = graph.input("motor_ligase", "Oxford Nanopore adapter ligase",
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
                 (ligase, InputRole.CATALYST)), (library,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The cited workflow does not disclose the complete kit- and revision-specific motor-bearing adapter assembly."))
    graph.sequence(library, SequencingPlan(
        acquired_roles=("Sequencing output",),
        limitation="Oxford Nanopore motor entry is platform-driven; no sequencing primer binds the library."))
    graph.declare_reviewed_scope((*reviewed, "motor_adapter", "motor_ligase", "ligate_motor", "library"), authority)
    return graph.require_publishable()
