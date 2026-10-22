"""Independent, deliberately narrow checks of declared molecular chemistry.

This reads a ProtocolGraph after construction.  It does not decide whether its
declarations faithfully transcribe a paper, and a locally established claim never
means that the whole reaction or protocol has been proved.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from chemdraw import Construct, Scene, bases_pair
from circular import ClosedCircle, RCAProduct, rolling_circle
from cohesive_assembly import CircularCohesiveAssembly, CohesiveFragment
from end_combinations import EndCombinationGeometry, LoadedEndPool, validate_groups
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, MolecularRelation,
                            ProtocolGraph, StepKind)
from restriction import Digest, FilledDigest
from topology import HairpinTopology
import seqprimers


class Status(str, Enum):
    ESTABLISHED = "established"
    CONDITIONAL = "conditional"
    UNRESOLVED = "unresolved"
    CONTRADICTORY = "contradictory"


@dataclass(frozen=True)
class Claim:
    target: str
    question: str
    status: Status
    reason: str


@dataclass(frozen=True)
class ValidationWarning:
    target: str
    message: str


def warnings(graph: ProtocolGraph) -> tuple[ValidationWarning, ...]:
    """Explicit gaps and failed or uncheckable claims for the audit and page."""
    found = [ValidationWarning(step.id, step.reason) for step in graph.steps.values()
             if step.kind is StepKind.MAGIC]
    for claim in check_graph(graph):
        if claim.status is Status.CONTRADICTORY:
            found.append(ValidationWarning(
                claim.target, f"Contradictory {claim.question}: {claim.reason}"))
        elif (claim.status is Status.UNRESOLVED and claim.target in graph.steps
              and graph.steps[claim.target].relations):
            found.append(ValidationWarning(
                claim.target, f"Unresolved {claim.question}: {claim.reason}"))
    return tuple(found)


def _chemical(material):
    return material.chemical_state if material.chemical_state is not None else material.state


def _uncertain_construct(value: Construct) -> bool:
    return any(seg.placeholder or seg.inferred or set(seg.top.upper()) - set("ACGTU")
               for seg in value)


def _uncertain_scene(a, b, columns: set[int]) -> bool:
    for strand in (a, b):
        pos = strand.col
        for seg in strand.drawn():
            if (seg.placeholder or seg.inferred or set(seg.top.upper()) - set("ACGTU")):
                if columns.intersection(range(pos, pos + len(seg))):
                    return True
            pos += len(seg)
    return False


def _uncertain_construct_window(value: Construct, start: int, end: int) -> bool:
    pos = 0
    for seg in value:
        if pos < end and start < pos + len(seg):
            if seg.placeholder or seg.inferred or set(seg.top.upper()) - set("ACGTU"):
                return True
        pos += len(seg)
    return False


def _pairing(target: str, scene: Scene) -> list[Claim]:
    claims = []
    for strand in scene.strands.values():
        if strand.partner is None:
            continue
        partner = scene.strands[strand.partner]
        mine, theirs = strand.cells(), partner.cells()
        columns = set(mine) & set(theirs)
        checked = {col for col in columns
                   if mine[col][1] not in strand.unpaired
                   and theirs[col][1] not in partner.unpaired}
        bad = [col for col in sorted(checked)
               if not bases_pair(mine[col][0], theirs[col][0])]
        question = f"pairing {strand.name} ↔ {partner.name}"
        if bad:
            claims.append(Claim(target, question, Status.CONTRADICTORY,
                                f"incompatible bases at drawn column {bad[0]}"))
        elif not checked:
            claims.append(Claim(target, question, Status.UNRESOLVED,
                                "no non-exempt paired columns"))
        elif _uncertain_scene(strand, partner, checked):
            claims.append(Claim(target, question, Status.CONDITIONAL,
                                f"{len(checked)} compatible columns include unknown or inferred bases"))
        else:
            claims.append(Claim(target, question, Status.ESTABLISHED,
                                f"{len(checked)} exact Watson–Crick columns"))
    return claims


def _pcr(target: str, geometry: PCRGeometry) -> Claim:
    question = "primer binding within declared PCR geometry"
    if geometry.availability != "exact":
        return Claim(target, question, Status.UNRESOLVED,
                     f"{geometry.availability} primer geometry: {geometry.statement}")
    if (geometry.forward_hit is None or geometry.reverse_hit is None
            or geometry.forward_hit.strand != "bottom"
            or geometry.reverse_hit.strand != "top"
            or geometry.forward_hit.end > geometry.reverse_hit.start):
        return Claim(target, question, Status.CONTRADICTORY,
                     "claimed primer sites do not face inward")
    uncertain = False
    for primer, claimed in ((geometry.forward, geometry.forward_hit),
                            (geometry.reverse, geometry.reverse_hit)):
        if primer is None or claimed is None:
            return Claim(target, question, Status.CONTRADICTORY,
                         "exact geometry lacks a primer or landing site")
        query = seqprimers.SeqPrimer("PCR", primer.name, primer.sequence)
        hits = seqprimers.locate_all(geometry.template, query,
                                    min_anneal=primer.anneal_nt)
        if claimed not in hits:
            return Claim(target, question, Status.CONTRADICTORY,
                         f"{primer.name} claimed site is absent from the declared template")
        uncertain |= (bool(set(primer.sequence.upper()) - set("ACGT"))
                      or _uncertain_construct_window(geometry.template,
                                                     claimed.start, claimed.end))
    return Claim(target, question,
                 Status.CONDITIONAL if uncertain else Status.ESTABLISHED,
                 "both inward primer sites recomputed; template contains unknown or inferred bases"
                 if uncertain else "both inward primer sites recomputed from exact bases")


def _cut(target: str, digest: Digest) -> Claim:
    site = digest.enzyme.site
    top_cut, bottom_cut = digest.enzyme.top_cut, digest.enzyme.bottom_cut
    overhang = site[min(top_cut, bottom_cut):max(top_cut, bottom_cut)]
    if (digest.left_top + digest.right_top != digest.top
            or digest.left_bottom_drawn + digest.right_bottom_drawn != digest.bottom_drawn
            or digest.enzyme.overhang != overhang):
        return Claim(target, "restriction cut follows declared coordinates",
                     Status.CONTRADICTORY, "cut fragments or overhang disagree with the input site")
    return Claim(target, "restriction cut follows declared coordinates",
                 Status.ESTABLISHED,
                 "both strand paths and the cohesive end follow the declared site cuts")


def _transition(graph: ProtocolGraph, step,
                relation: MolecularRelation | None) -> Claim:
    target = step.id
    question = "declared product follows from molecular inputs"
    if step.kind is StepKind.MAGIC:
        return Claim(target, question, Status.UNRESOLVED,
                     f"unmodelled transition: {step.reason}")
    if relation is MolecularRelation.COHESIVE_ASSEMBLY:
        question = relation.value
        geometry = step.geometry
        if not isinstance(geometry, CircularCohesiveAssembly):
            return Claim(target, question, Status.UNRESOLVED,
                         "explicit circular cohesive-fragment geometry is required")
        ports = [port.material for port in step.inputs
                 if port.role in {InputRole.TEMPLATE, InputRole.STRUCTURAL}]
        declared = dict(geometry.fragments)
        if len(ports) != len(set(ports)) or set(ports) != set(declared):
            return Claim(target, question, Status.CONTRADICTORY,
                         "assembly geometry does not account for the molecular input IDs")
        fragments = []
        for identifier, declared_fragment in geometry.fragments:
            actual = _chemical(graph.materials[identifier])
            if not isinstance(actual, CohesiveFragment):
                return Claim(target, question, Status.UNRESOLVED,
                             f"input {identifier} has no declared cohesive fragment state")
            actual_core = "".join(seg.top for seg in actual.core).upper()
            declared_core = "".join(seg.top for seg in declared_fragment.core).upper()
            if (actual.left, actual.right, actual.left_polarity, actual.right_polarity) != (
                    declared_fragment.left, declared_fragment.right,
                    declared_fragment.left_polarity, declared_fragment.right_polarity):
                return Claim(target, question, Status.CONTRADICTORY,
                             f"geometry changes the known cohesive ends of {identifier}")
            if actual_core != declared_core:
                uncertain_core = any(seg.placeholder or seg.inferred or set(seg.top.upper()) - set("ACGT")
                                     for fragment in (actual, declared_fragment) for seg in fragment.core)
                return Claim(target, question,
                             Status.UNRESOLVED if uncertain_core else Status.CONTRADICTORY,
                             f"geometry and input core differ for {identifier}; "
                             + ("unspecified cores do not establish sequence identity"
                                if uncertain_core else "known core DNA is changed"))
            fragments.append(actual)
        if any((a.right, a.right_polarity) != (b.left, b.left_polarity)
               for a, b in zip(fragments, fragments[1:] + fragments[:1])):
            return Claim(target, question, Status.CONTRADICTORY,
                         "consecutive fragments have incompatible cohesive ends")
        products = [_chemical(graph.materials[out]) for out in step.outputs
                    if graph.materials[out].kind is not MaterialKind.WASTE]
        if len(products) != 1 or not isinstance(products[0], Construct):
            return Claim(target, question, Status.UNRESOLVED,
                         "one molecular product strand is required")
        expected = "".join(fragment.left + "".join(seg.top for seg in fragment.core)
                           for fragment in fragments)
        observed = products[0].top().upper()
        expected = expected.upper()
        uncertain = (_uncertain_construct(products[0]) or any(
            seg.placeholder or seg.inferred or set(seg.top.upper()) - set("ACGT")
            for fragment in fragments for seg in fragment.core))
        if len(observed) != len(expected) or observed not in expected + expected:
            return Claim(target, question,
                         Status.UNRESOLVED if uncertain else Status.CONTRADICTORY,
                         "unspecified sequence prevents an exact product comparison"
                         if uncertain else
                         "product changes fragment order, core sequence or fusion multiplicity")
        return Claim(target, question,
                     Status.CONDITIONAL if uncertain else Status.ESTABLISHED,
                     "declared ends and ordered cores agree with one fusion per seam; "
                     + ("cores include unspecified or inferred sequence; " if uncertain else "")
                     + "digestion, reaction efficiency and competing products are not established")
    if relation is MolecularRelation.END_COMBINATIONS:
        question = "all loaded end combinations match their two-primer library fate"
        pools = [graph.materials[port.material].chemical_state for port in step.inputs
                 if port.role is InputRole.STRUCTURAL and isinstance(
                     graph.materials[port.material].chemical_state, LoadedEndPool)]
        if len(pools) != 1 or not isinstance(step.geometry, EndCombinationGeometry):
            return Claim(target, question, Status.UNRESOLVED,
                         "one loaded-end pool and end-combination geometry are required")
        pcr_step = graph.steps.get(step.geometry.pcr_step)
        if (pcr_step is None or pcr_step.kind is not StepKind.PCR
                or not isinstance(pcr_step.geometry, PCRGeometry)):
            return Claim(target, question, Status.UNRESOLVED,
                         "the declared downstream PCR step or its geometry is unavailable")
        pcr_primers = (pcr_step.geometry.forward, pcr_step.geometry.reverse)
        if pcr_step.geometry.availability != "exact" or any(
                primer is None for primer in pcr_primers):
            return Claim(target, question, Status.UNRESOLVED,
                         "the downstream PCR primer sites are not exact")
        if {(primer.sequence, primer.anneal_nt) for primer in pcr_primers} != {
                (primer.sequence, primer.anneal_nt) for primer in step.geometry.primers}:
            return Claim(target, question, Status.CONTRADICTORY,
                         "end-class primers differ from the downstream PCR primers")
        retained = [group.output for group in step.geometry.groups
                    if group.output in graph.materials
                    and graph.materials[group.output].kind is not MaterialKind.WASTE]
        frontier = list(retained)
        seen: set[str] = set()
        reaches_pcr = False
        while frontier:
            material = frontier.pop()
            if material in seen:
                continue
            seen.add(material)
            for consumer in graph.consumers_of(material):
                if not any(port.material == material and port.role is InputRole.TEMPLATE
                           for port in consumer.inputs):
                    continue
                if consumer.id == pcr_step.id:
                    reaches_pcr = True
                frontier.extend(consumer.outputs)
        if not reaches_pcr:
            return Claim(target, question, Status.CONTRADICTORY,
                         "retained end classes do not feed the declared PCR step")
        mismatch = validate_groups(
            pools[0], step.geometry,
            {output: graph.materials[output].kind is MaterialKind.WASTE
             for output in step.outputs})
        if mismatch:
            return Claim(target, question, Status.CONTRADICTORY, mismatch)
        return Claim(target, question, Status.ESTABLISHED,
                     "all end-identity classes are accounted for by exact primer-entry "
                     "footprints; this does not establish reaction yield or PCR efficiency")
    templates = [_chemical(graph.materials[p.material]) for p in step.inputs
                 if p.role is InputRole.TEMPLATE]
    products = [_chemical(graph.materials[id]) for id in step.outputs
                if graph.materials[id].kind.value != "discarded material"]
    if relation is MolecularRelation.CIRCLE_CLOSURE:
        before = next((x for x in templates if isinstance(x, Construct)), None)
        after = next((x for x in products if isinstance(x, ClosedCircle)), None)
        if before is not None and after is not None:
            if before.top() != after.sequence:
                return Claim(target, question, Status.CONTRADICTORY,
                             "closed circle sequence differs from its linear substrate")
            return Claim(target, "circle preserves linear strand order",
                         Status.CONDITIONAL if _uncertain_construct(before)
                         else Status.ESTABLISHED,
                         "closed sequence equals substrate; end-group chemistry is outside this claim")
    if relation is MolecularRelation.RESTRICTION_FILL:
        digest = next((x for x in templates if isinstance(x, Digest)), None)
        filled = next((x for x in products if isinstance(x, FilledDigest)), None)
        if digest is not None and filled is not None:
            if filled.digest != digest:
                return Claim(target, question, Status.CONTRADICTORY,
                             "filled ends refer to a different restriction substrate")
            try:
                expected = digest.fill_in(filled.biotin_base)
            except ValueError as exc:
                return Claim(target, question, Status.CONTRADICTORY, str(exc))
            if (filled.left_new_top != expected.left_new_top
                    or filled.right_new_bottom != expected.right_new_bottom):
                return Claim(target, question, Status.CONTRADICTORY,
                             "filled bases disagree with the restriction overhang")
            return Claim(target, "recessed ends copy the restriction overhang",
                         Status.ESTABLISHED,
                         "both fill-in tracts derive from the declared cut coordinates")
    if relation is MolecularRelation.ROLLING_CIRCLE_COPY:
        before = next((x for x in templates if isinstance(x, ClosedCircle)), None)
        after = next((x for x in products if isinstance(x, RCAProduct)), None)
        if before is not None and after is not None:
            try:
                expected = rolling_circle(before, after.primer, after.copies)
            except ValueError as exc:
                return Claim(target, question, Status.CONTRADICTORY, str(exc))
            if (after.template.sequence != before.sequence
                    or after.sequence != expected.sequence
                    or after.binding_start != expected.binding_start):
                return Claim(target, question, Status.CONTRADICTORY,
                             "RCA repeat or primer landing disagrees with input circle")
            return Claim(target, "RCA repeats copy the input circle",
                         Status.CONDITIONAL if _uncertain_construct(before.linear)
                         else Status.ESTABLISHED,
                         "repeat sequence and primer landing derived from the input circle")
    if relation is MolecularRelation.CONTACT_JUNCTION:
        filled = next((x for x in templates if isinstance(x, FilledDigest)), None)
        product = next((x for x in products if isinstance(x, Construct)), None)
        if filled is not None and product is not None:
            count = product.top().count(filled.junction)
            if count != 1:
                return Claim(target, question, Status.CONTRADICTORY,
                             f"expected one derived restriction junction, found {count}")
            return Claim(target, "contact product contains the derived junction",
                         Status.ESTABLISHED,
                         "one exact fill-in junction occurs in the declared product")
    if relation in {MolecularRelation.HAIRPIN_STEM, MolecularRelation.STRAND_ORDER}:
        before = next((x for x in templates if isinstance(x, Construct)), None)
        fold = next((x for x in products if isinstance(x, HairpinTopology)), None)
        if relation is MolecularRelation.HAIRPIN_STEM and before is not None and fold is not None:
            if len(before) != len(fold.sequence):
                return Claim(target, question, Status.CONTRADICTORY,
                             "folded strand length differs from its input")
            pos = 0
            for seg in before:
                if not seg.placeholder:
                    for base, folded in zip(seg.top.upper(), fold.sequence[pos:pos + len(seg)]):
                        if base in "ACGT" and folded in "ACGT" and base != folded:
                            return Claim(target, question, Status.CONTRADICTORY,
                                         "folded strand changes a known input base")
                pos += len(seg)
            first = fold.sequence[fold.first_stem_start:fold.first_stem_start + fold.stem_length]
            second = fold.sequence[fold.second_stem_start:fold.second_stem_start + fold.stem_length]
            if any(not bases_pair(a, b) for a, b in zip(first, reversed(second))):
                return Claim(target, question, Status.CONTRADICTORY,
                             "declared stem cannot base-pair antiparallel")
            uncertain = (_uncertain_construct(before) or
                         set(first + second) - set("ACGT"))
            return Claim(target, "declared intramolecular stem can pair",
                         Status.CONDITIONAL if uncertain else Status.ESTABLISHED,
                         "stem pairing is possible given unknown bases; fold preference is unproved"
                         if uncertain else "exact stem bases pair; fold preference is unproved")
        after = next((x for x in products if isinstance(x, Construct)), None)
        if relation is MolecularRelation.STRAND_ORDER and before is not None and after is not None:
            if before.top() != after.top():
                return Claim(target, question, Status.CONTRADICTORY,
                             "a sequence-preserving operation changes the declared strand")
            return Claim(target, "covalent strand order is preserved",
                         Status.CONDITIONAL if (_uncertain_construct(before)
                                                or _uncertain_construct(after))
                         else Status.ESTABLISHED,
                         "fold preference or cross-link formation is not established")
    return Claim(target, relation.value if relation is not None else question,
                 Status.UNRESOLVED,
                 "declared molecular states are insufficient for this relation"
                 if relation is not None else
                 "no molecular input-to-output relation is declared")


def check_graph(graph: ProtocolGraph) -> tuple[Claim, ...]:
    """Check only claims supported by graph declarations; never infer missing states."""
    claims = []
    for material in graph.materials.values():
        state = _chemical(material)
        if isinstance(state, Digest):
            claims.append(_cut(material.id, state))
        if isinstance(material.chemical_state, RCAProduct) and isinstance(material.state, Construct):
            claims.append(Claim(
                material.id, "rendered concatemer agrees with chemical state",
                (Status.CONTRADICTORY if material.state.top() != state.sequence
                 else Status.CONDITIONAL if _uncertain_construct(material.state)
                 else Status.ESTABLISHED),
                "drawn strand differs from derived RCA repeats"
                if material.state.top() != state.sequence else
                "drawn strand equals the derived repeats"))
        if isinstance(material.state, Scene):
            claims.extend(_pairing(material.id, material.state))
        if isinstance(state, Scene) and state is not material.state:
            claims.extend(_pairing(material.id, state))
    for step in graph.steps.values():
        if isinstance(step.geometry, PCRGeometry):
            claims.append(_pcr(step.id, step.geometry))
        scenes = getattr(step.geometry, "scenes", ())
        for i, scene in enumerate(scenes, 1):
            claims.extend(_pairing(f"{step.id}:geometry:{i}", scene))
        claims.extend(_transition(graph, step, relation)
                      for relation in (step.relations or (None,)))
    return tuple(claims)
