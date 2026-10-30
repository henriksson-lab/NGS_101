"""Independent, deliberately narrow checks of declared molecular chemistry.

This reads a ProtocolGraph after construction.  It does not decide whether its
declarations faithfully transcribe a paper, and a locally established claim never
means that the whole reaction or protocol has been proved.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from chemdraw import Construct, Scene, bases_pair, revcomp
from circular import ClosedCircle, RCAProduct, rolling_circle
from assembly_geometry import StrandRetention, BridgedAssembly, OverlapExtension, ActivatedEndTransfer, DuplexEndState, LinearFivePrimeAssembly, CircularOverlapReference, DuplexPartition, strands as cohesive_strands
from cohesive_assembly import CircularCohesiveAssembly, CohesiveFragment
from end_combinations import EndCombinationGeometry, LoadedEndPool, validate_groups
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, MolecularRelation,
                            ProtocolGraph, StepKind)
from restriction import Digest, FilledDigest
from topology import HairpinTopology
from rna_lariat import RNALariatTopology
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


def _duplex_reference_construct(value):
    """Read an explicitly supplied full duplex view as one reference axis."""
    if isinstance(value,Construct):
        return value
    if isinstance(value,Scene) and set(value.strands)=={'top','bottom'}:
        top,bottom=value.strands['top'],value.strands['bottom']
        a=''.join(seg.top for seg in top.segs);b=''.join(seg.top for seg in bottom.segs)
        if (not top.rev and bottom.rev and top.col==bottom.col and len(a)==len(b)
                and bottom.partner=='top' and not top.unpaired and not bottom.unpaired
                and a.upper()==revcomp(b).upper()):
            return Construct(tuple(top.segs),'explicit duplex reference')
    return None


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
        uncertain |= (primer.inferred or bool(set(primer.sequence.upper()) - set("ACGT"))
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


def _circular_overlap_reference(graph, step) -> Claim:
    """Check intended reference arithmetic, independently of unmodelled chemistry."""
    target=step.id
    geom=step.geometry
    if not isinstance(geom,CircularOverlapReference):
        return Claim(target,"conditional circular reference",Status.UNRESOLVED,'explicit duplex-map overlap reference required')
    ports={p.material for p in step.inputs if p.role in {InputRole.TEMPLATE,InputRole.STRUCTURAL} and graph.materials[p.material].kind is MaterialKind.MOLECULAR}
    if ports!={i for i,_,_,_ in geom.parts}:
        return Claim(target,"conditional circular reference",Status.CONTRADICTORY,'reference order changes molecular contributors')
    maps=[]
    for identifier,declared,left,right in geom.parts:
        actual=_chemical(graph.materials[identifier])
        if not isinstance(actual,Construct):
            return Claim(target,"conditional circular reference",Status.UNRESOLVED,'an actual duplex reference map is required')
        if actual.top()!=declared.top():
            return Claim(target,"conditional circular reference",Status.CONTRADICTORY,'actual duplex map differs from intended reference')
        maps.append((actual,left,right))
    if any(ra!=lb or a.top()[-ra:]!=b.top()[:lb] for (a,_,ra),(b,lb,_) in zip(maps,maps[1:]+maps[:1])):
        return Claim(target,"conditional circular reference",Status.CONTRADICTORY,'reference overlaps differ')
    expected=''.join(c.top()[:-r] for c,_,r in maps)
    products=[_chemical(graph.materials[o]) for o in step.outputs if isinstance(_chemical(graph.materials[o]),Construct)]
    if len(products)!=1 or products[0].top()!=expected:
        return Claim(target,"conditional circular reference",Status.CONTRADICTORY,'intended reference loses, changes or duplicates an overlap')
    return Claim(target,"conditional circular reference",Status.CONDITIONAL,'Only intended circular reference order is checked; resection, recombination, repair, closure and biological recovery remain unestablished')

def _transition(graph: ProtocolGraph, step,
                relation: MolecularRelation | None) -> Claim:
    target = step.id
    question = "declared product follows from molecular inputs"
    if step.kind is StepKind.MAGIC:
        return Claim(target, question, Status.UNRESOLVED,
                     f"unmodelled transition: {step.reason}")
    if relation is MolecularRelation.LINEAR_COHESIVE_ASSEMBLY:
        geom=step.geometry
        if not isinstance(geom,LinearFivePrimeAssembly):
            return Claim(target,relation.value,Status.UNRESOLVED,'explicit linear fragment order required')
        ports={p.material for p in step.inputs if p.role in {InputRole.TEMPLATE,InputRole.STRUCTURAL}}
        if ports!={i for i,_ in geom.fragments}:
            return Claim(target,relation.value,Status.CONTRADICTORY,'linear geometry changes contributors')
        fragments=[]
        for identifier,declared in geom.fragments:
            actual=_chemical(graph.materials[identifier])
            if not isinstance(actual,CohesiveFragment):
                return Claim(target,relation.value,Status.UNRESOLVED,'explicit actual cohesive fragments required')
            if cohesive_strands(actual)!=cohesive_strands(declared):
                return Claim(target,relation.value,Status.CONTRADICTORY,'linear assembly changes an input strand')
            fragments.append(actual)
        if any(a.right!=b.left or (a.right_polarity,b.left_polarity)!=('5-prime','5-prime') for a,b in zip(fragments,fragments[1:])):
            return Claim(target,relation.value,Status.CONTRADICTORY,'linear seams do not pair')
        products=[_chemical(graph.materials[o]) for o in step.outputs if isinstance(_chemical(graph.materials[o]),CohesiveFragment)]
        expected=''.join(('' if i==0 else f.left)+''.join(s.top for s in f.core) for i,f in enumerate(fragments))
        if len(products)!=1 or ''.join(s.top for s in products[0].core)!=expected or (products[0].left,products[0].right,products[0].left_polarity,products[0].right_polarity)!=(fragments[0].left,fragments[-1].right,'5-prime','5-prime'):
            return Claim(target,relation.value,Status.CONTRADICTORY,'linear product loses external ends or changes seam multiplicity')
        return Claim(target,relation.value,Status.CONDITIONAL if products[0].inferred or any(s.inferred or s.placeholder for s in products[0].core) else Status.ESTABLISHED,'Ordered cores and internal seams agree; external ends are retained, not joined. End-group declarations are preconditions, not experimental evidence')
    if relation is MolecularRelation.ACTIVATED_END_TRANSFER:
        geom=step.geometry
        if not isinstance(geom,ActivatedEndTransfer):
            return Claim(target,relation.value,Status.UNRESOLVED,'explicit activated end chemistry required')
        ports={p.material for p in step.inputs if p.role in {InputRole.TEMPLATE,InputRole.STRUCTURAL}}
        if ports!={geom.vector_id,geom.insert_id}:
            return Claim(target,relation.value,Status.CONTRADICTORY,'transfer changes molecular contributors')
        for identifier,declared in ((geom.vector_id,geom.vector),(geom.insert_id,geom.insert)):
            actual=_chemical(graph.materials[identifier])
            if not isinstance(actual,DuplexEndState) or (actual.top5,actual.top3,actual.bottom5,actual.bottom3)!=(declared.top5,declared.top3,declared.bottom5,declared.bottom3) or cohesive_strands(actual.fragment)!=cohesive_strands(declared.fragment):
                return Claim(target,relation.value,Status.CONTRADICTORY,'actual substrate end groups or strands differ from transfer geometry')
        if (geom.vector.top3,geom.vector.bottom3)!=('phosphotyrosyl','phosphotyrosyl') or (geom.insert.top5,geom.insert.bottom5)!=('OH','OH'):
            return Claim(target,relation.value,Status.CONTRADICTORY,'activated donor/hydroxyl acceptor requirement is not met')
        a,b=geom.insert.fragment,geom.vector.fragment
        if (a.right,a.right_polarity)!=(b.left,b.left_polarity) or (b.right,b.right_polarity)!=(a.left,a.left_polarity):
            return Claim(target,relation.value,Status.CONTRADICTORY,'declared donor and insert ends do not pair')
        expected=a.left+''.join(s.top for s in a.core)+b.left+''.join(s.top for s in b.core)
        products=[_chemical(graph.materials[o]) for o in step.outputs if isinstance(_chemical(graph.materials[o]),Construct)]
        if len(products)!=1 or products[0].top()!=expected:
            return Claim(target,relation.value,Status.CONTRADICTORY,'transfer reference loses or changes a contributor')
        return Claim(target,relation.value,Status.CONDITIONAL if _uncertain_construct(products[0]) else Status.ESTABLISHED,'Activated donor and hydroxyl acceptor end declarations and paired reference order agree. Remaining strand closures and biological recovery are not established')
    if relation is MolecularRelation.PCR_PRODUCT:
        geom=step.geometry
        if not isinstance(geom,PCRGeometry) or geom.availability!='exact':
            return Claim(target,relation.value,Status.UNRESOLVED,'exact primer-binding geometry required')
        templates=[_chemical(graph.materials[p.material]) for p in step.inputs if p.role is InputRole.TEMPLATE]
        oligos=[_chemical(graph.materials[p.material]) for p in step.inputs if p.role is InputRole.STRUCTURAL]
        if geom.amplicon_id and geom.amplicon_id not in step.outputs:
            return Claim(target,relation.value,Status.CONTRADICTORY,'declared amplicon is not an output of this PCR')
        output_ids=(geom.amplicon_id,) if geom.amplicon_id else step.outputs
        actual=[_chemical(graph.materials[o]) for o in output_ids if isinstance(_chemical(graph.materials[o]),Construct)]
        if len(templates)!=1 or not isinstance(templates[0],Construct) or templates[0].top()!=geom.template.top():
            return Claim(target,relation.value,Status.CONTRADICTORY,'PCR template differs from binding geometry')
        if not all(isinstance(o,Construct) for o in oligos) or sorted(o.top() for o in oligos)!=sorted((geom.forward.sequence,geom.reverse.sequence)):
            return Claim(target,relation.value,Status.CONTRADICTORY,'input primer polymers differ from geometry')
        expected=geom.forward.sequence+templates[0].top()[geom.forward_hit.end:geom.reverse_hit.start]+revcomp(geom.reverse.sequence)
        if len(actual)!=1 or actual[0].top()!=expected:
            return Claim(target,relation.value,Status.CONTRADICTORY,'product changes primer tails or retained template interval')
        return Claim(target,relation.value,Status.CONDITIONAL if _uncertain_construct(actual[0]) else Status.ESTABLISHED,'Reference primer polymers and selected template interval agree; terminal additions, errors and yield are outside scope')
    if relation is MolecularRelation.OVERLAP_EXTENSION:
        geom=step.geometry
        if not isinstance(geom,OverlapExtension):
            return Claim(target,relation.value,Status.UNRESOLVED,'explicit terminal overlap required')
        declared=((geom.left_id,geom.left),(geom.right_id,geom.right))
        ports={p.material for p in step.inputs if p.role in {InputRole.TEMPLATE,InputRole.STRUCTURAL}}
        if ports!={geom.left_id,geom.right_id}:
            # A separately drawn active pair retains its molecular parents; it is
            # not a new covalent template or an excuse to skip their provenance.
            if len(ports)!=1:
                return Claim(target,relation.value,Status.CONTRADICTORY,'overlap geometry changes contributing inputs')
            active_id=next(iter(ports));active=_chemical(graph.materials[active_id])
            producer=graph.producer_of(active_id)
            if (not isinstance(active,OverlapExtension) or producer is None
                    or producer.kind is not StepKind.HYBRIDIZATION
                    or {p.material for p in producer.inputs if p.role in {InputRole.TEMPLATE,InputRole.STRUCTURAL}}!={geom.left_id,geom.right_id}
                    or active.left.top()!=geom.left.top() or active.right.top()!=geom.right.top()
                    or active.overlap_nt!=geom.overlap_nt):
                return Claim(target,relation.value,Status.CONTRADICTORY,'active pair does not retain the declared contributing fragment geometry')
        for identifier,con in declared:
            actual=_chemical(graph.materials[identifier])
            if not isinstance(actual,Construct) or actual.top()!=con.top():
                return Claim(target,relation.value,Status.CONTRADICTORY,'fragment differs from declared overlap')
        products=[_chemical(graph.materials[o]) for o in step.outputs if isinstance(_chemical(graph.materials[o]),Construct)]
        n=geom.overlap_nt
        if geom.left.top()[-n:]!=geom.right.top()[:n] or len(products)!=1 or products[0].top()!=geom.left.top()+geom.right.top()[n:]:
            return Claim(target,relation.value,Status.CONTRADICTORY,'product duplicates, loses or changes the shared interval')
        return Claim(target,relation.value,Status.CONDITIONAL if _uncertain_construct(products[0]) else Status.ESTABLISHED,'Two reference fragments share the declared overlap once; successful extension and competing products are not established')
    if relation is MolecularRelation.DUPLEX_PARTITION:
        geom=step.geometry
        if not isinstance(geom,DuplexPartition):
            return Claim(target,relation.value,Status.UNRESOLVED,'explicit reference partition required')
        if geom.input_id not in {p.material for p in step.inputs} or set(step.outputs)!={i for i,_,_ in geom.intervals}:
            return Claim(target,relation.value,Status.CONTRADICTORY,'partition changes input or output accounting')
        parent=_duplex_reference_construct(_chemical(graph.materials[geom.input_id]))
        if not isinstance(parent,Construct):
            return Claim(target,relation.value,Status.UNRESOLVED,'explicit parent reference required')
        if parent.top()!=geom.parent.top():
            return Claim(target,relation.value,Status.CONTRADICTORY,'partition parent differs')
        for identifier,start,end in geom.intervals:
            actual=_duplex_reference_construct(_chemical(graph.materials[identifier]))
            if not isinstance(actual,Construct):
                return Claim(target,relation.value,Status.UNRESOLVED,'each output needs a reference map')
            if actual.top()!=parent.top()[start:end]:
                return Claim(target,relation.value,Status.CONTRADICTORY,'partition output changes or duplicates parent bases')
        return Claim(target,relation.value,Status.CONDITIONAL,'Reference intervals account for each parent base once; physical cut offsets, end chemistry, recognition and yield remain unassigned')
    if relation is MolecularRelation.STRAND_RETENTION:
        geometry=step.geometry
        if not isinstance(geometry,StrandRetention):
            return Claim(target,relation.value,Status.UNRESOLVED,'explicit retained strand intervals required')
        if geometry.input_id not in {p.material for p in step.inputs} or geometry.output_id not in step.outputs:
            return Claim(target,relation.value,Status.CONTRADICTORY,'interval geometry names different reaction materials')
        before=_chemical(graph.materials[geometry.input_id]);after=_chemical(graph.materials[geometry.output_id])
        if not isinstance(before,Scene) or not isinstance(after,CohesiveFragment):
            return Claim(target,relation.value,Status.UNRESOLVED,'parent duplex Scene and cohesive product required')
        actual=tuple(''.join(s.top for s in strand.segs) for strand in before.strands.values())
        if actual != (geometry.top,geometry.bottom) or cohesive_strands(after) != cohesive_strands(geometry.fragment):
            return Claim(target,relation.value,Status.CONTRADICTORY,'physical parent/product strands differ from retained intervals')
        kept=(actual[0][slice(*geometry.top_span)],actual[1][slice(*geometry.bottom_span)])
        uncertain = after.inferred or any(seg.inferred or seg.placeholder for strand in before.strands.values() for seg in strand.segs)
        return Claim(target,relation.value,Status.CONDITIONAL if uncertain else Status.ESTABLISHED,
                     'Explicit physical strands match retained intervals; enzyme specificity and the experimental distribution of ends are not established') if kept==cohesive_strands(after) else Claim(target,relation.value,Status.CONTRADICTORY,'retained bases differ')
    if relation is MolecularRelation.BRIDGED_ASSEMBLY:
        geometry=step.geometry
        if not isinstance(geometry,BridgedAssembly):
            return Claim(target,relation.value,Status.UNRESOLVED,'explicit bridge geometry required')
        for identifier,con in geometry.fragments:
            if identifier not in graph.materials:
                return Claim(target,relation.value,Status.CONTRADICTORY,'geometry names an unknown reaction material')
            actual=_chemical(graph.materials[identifier])
            if identifier not in {p.material for p in step.inputs} or not isinstance(actual,Construct) or actual.top()!=con.top():
                return Claim(target,relation.value,Status.CONTRADICTORY,'fragment input differs from bridge geometry')
        for identifier,sequence in geometry.bridges:
            if identifier not in graph.materials:
                return Claim(target,relation.value,Status.CONTRADICTORY,'geometry names an unknown reaction material')
            actual=_chemical(graph.materials[identifier])
            if identifier not in {p.material for p in step.inputs} or not isinstance(actual,Construct) or actual.top()!=sequence:
                return Claim(target,relation.value,Status.CONTRADICTORY,'bridge input differs from declared geometry')
        products=[_chemical(graph.materials[o]) for o in step.outputs if isinstance(_chemical(graph.materials[o]),Construct)]
        expected=geometry.product('expected').top()
        if len(products)!=1 or products[0].top()!=expected:
            return Claim(target,relation.value,Status.CONTRADICTORY,'product incorporates bridge DNA or changes fragment order')
        return Claim(target,relation.value,Status.CONDITIONAL if _uncertain_construct(products[0]) else Status.ESTABLISHED,
                     'Bridge arms complement adjacent ends; circular product concatenates only the declared fragments. Ligation chemistry, yield and other products remain outside this relation')
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
        products = [_duplex_reference_construct(_chemical(graph.materials[out])) for out in step.outputs
                    if graph.materials[out].kind is not MaterialKind.WASTE]
        if len(products) != 1 or products[0] is None:
            return Claim(target, question, Status.UNRESOLVED,
                         "one molecular product strand is required")
        expected = "".join(fragment.left + "".join(seg.top for seg in fragment.core)
                           for fragment in fragments)
        observed = products[0].top().upper()
        expected = expected.upper()
        uncertain = (_uncertain_construct(products[0]) or any(fragment.inferred for fragment in fragments) or any(
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
        if isinstance(state, RNALariatTopology):
            try:
                state.validate_connectivity()
            except ValueError as error:
                claims.append(Claim(material.id, "declared RNA lariat backbone connectivity",
                                    Status.CONTRADICTORY, str(error)))
            else:
                claims.append(Claim(
                    material.id, "declared RNA lariat backbone connectivity", Status.ESTABLISHED,
                    "supplied RNA retains its ordinary backbone and one specified 2′–5′ bond; "
                    "this does not establish formation, enzyme behavior or source identity"))
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
        if isinstance(step.geometry,CircularOverlapReference):
            claims.append(_circular_overlap_reference(graph,step))
        if isinstance(step.geometry, PCRGeometry):
            claims.append(_pcr(step.id, step.geometry))
        scenes = getattr(step.geometry, "scenes", ())
        for i, scene in enumerate(scenes, 1):
            claims.extend(_pairing(f"{step.id}:geometry:{i}", scene))
        claims.extend(_transition(graph, step, relation)
                      for relation in (step.relations or (None,)))
    return tuple(claims)
