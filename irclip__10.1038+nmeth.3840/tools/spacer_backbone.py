"""Local ordered DNA/C18 backbone and explicitly conditional copying operation.

A linker occupies a covalent-backbone node, never a base or an empty string. These
objects cannot be passed to ordinary DNA sequence/circle routines by flattening.
The traversal engine supplies no biological rule about polymerase/linker behavior:
the caller must provide the stop hypothesis and its evidence boundary explicitly.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
from typing import TypeAlias
from chemdraw import MolecularFeature, Scene, Segment, iupac_complement


def revcomp(sequence: str) -> str:
    """Reverse-complement nucleotide classes, never schematic placeholders."""
    return iupac_complement(sequence)[::-1]

@dataclass(frozen=True)
class DNA:
    name: str
    sequence: str
    feature: MolecularFeature | None = None
    def __post_init__(self):
        if not self.name or not self.sequence or set(self.sequence)-set('ACGTRYSWKMBDHVN'):
            raise ValueError('named DNA blocks require nonempty IUPAC nucleotide sequences')

@dataclass(frozen=True)
class VariableDNA:
    """RNA-derived DNA of unassigned length/sequence, retaining its molecular origin."""
    name: str
    origin: str
    reverse_complemented: bool
    def __post_init__(self):
        if not self.name or not self.origin: raise ValueError('variable DNA requires its origin')

@dataclass(frozen=True)
class C18:
    name: str
    code: str = '/iSp18/'
    def __post_init__(self):
        if not self.name or self.code != '/iSp18/':
            raise ValueError('this local construct supports only source-written C18 linkers')

Node: TypeAlias = DNA | VariableDNA | C18

@dataclass(frozen=True)
class Backbone:
    nodes: tuple[Node, ...]
    circular: bool = False
    five_prime: str | None = 'phosphate'
    three_prime: str | None = 'OH'
    inferred: bool = False
    def __post_init__(self):
        if not isinstance(self.nodes,tuple):
            raise TypeError('ordered backbone nodes must be an immutable tuple')
        if not self.nodes or len({n.name for n in self.nodes}) != len(self.nodes):
            raise ValueError('backbone needs nonempty uniquely named nodes')
        if not all(isinstance(n, (DNA, VariableDNA, C18)) for n in self.nodes):
            raise TypeError('unsupported backbone node')
        if isinstance(self.nodes[0], C18) or isinstance(self.nodes[-1], C18):
            raise ValueError('source constructs have nucleotide blocks at their free/closure ends')
        if self.circular and (self.five_prime is not None or self.three_prime is not None):
            raise ValueError('a closed backbone has no free 5-prime/3-prime termini')
        if not self.circular and (self.five_prime not in {'OH', 'phosphate'} or self.three_prime != 'OH'):
            raise ValueError('linear source backbone requires declared free ends')
    def close(self) -> Backbone:
        if self.circular or (self.five_prime,self.three_prime) != ('phosphate','OH'):
            raise ValueError('closure requires linear 5-prime phosphate and 3-prime OH')
        # Preserve the ordered node objects, including every non-nucleotide linker.
        return replace(self,circular=True,five_prime=None,three_prime=None)
    def written(self) -> str:
        body=''.join(n.sequence if isinstance(n,DNA) else n.code if isinstance(n,C18)
                     else '['+n.origin+(' RC' if n.reverse_complemented else ' sense')+']'
                     for n in self.nodes)
        return ('circle('+body+')' if self.circular else
                ('/5phos/' if self.five_prime=='phosphate' else '/5OH/')+body)
    def exact_sequence(self) -> str:
        if any(not isinstance(n,DNA) for n in self.nodes):
            raise ValueError('linkers and unassigned RNA copies cannot be flattened into DNA bases')
        return ''.join(n.sequence for n in self.nodes)
    def display_segments(self):
        """Annotation/rendering segments, retaining non-nucleotide C18 placeholders.

        A Construct of these segments can validate feature definitions. It is not
        a conversion of this chemical backbone into an ordinary DNA sequence.
        """
        parts=[]
        for n in self.nodes:
            if isinstance(n,DNA):
                parts.append(Segment(n.name,n.sequence,feature=n.feature,
                                     placeholder=bool(set(n.sequence)-set("ACGT")),inferred=self.inferred))
            elif isinstance(n,C18):
                parts.append(Segment(n.name,'[C18 linker]',placeholder=True,inferred=self.inferred,
                                     note='Non-nucleotide covalent linker; no base-pair length assigned.'))
            else:
                parts.append(Segment(n.name,'['+n.origin+(' RC' if n.reverse_complemented else ' sense')+']',
                                     placeholder=True,inferred=self.inferred))
        return tuple(parts)

    def rows(self):
        sc=Scene()
        sc.strand('backbone',self.display_segments(),label='circle opened only for display' if self.circular else 'linear backbone',
                  mod5='' if self.circular else self.five_prime,mod3='' if self.circular else self.three_prime)
        sc.note('backbone', 'C18 boxes are chemical links, not bases. '+
                ('Covalently closed: last node joins first; both C18 links remain.' if self.circular
                 else 'Nucleotide and non-nucleotide nodes retain their original order.'))
        return sc.rows()

@dataclass(frozen=True)
class Landing:
    node_index: int
    start: int
    end: int
    primer: str


def bind_suffix(template: Backbone, primer: str, *, block_name: str) -> Landing:
    """Longest exact 3-prime primer suffix on one specified DNA block.

    This is a sequence relation, not a claim about affinity, kinetics or an enzyme's
    minimum annealing length. A repeated longest site is rejected as ambiguous.
    """
    if not primer or set(primer)-set('ACGTRYSWKMBDHVN'):
        raise ValueError('primer must be a nonempty IUPAC DNA string')
    candidates=[(i,n) for i,n in enumerate(template.nodes) if n.name==block_name]
    if len(candidates)!=1 or not isinstance(candidates[0][1],DNA):
        raise ValueError('landing requires one named nucleotide block')
    index,node=candidates[0]
    for width in range(min(len(primer),len(node.sequence)),0,-1):
        target=revcomp(primer[-width:])
        hits=[i for i in range(len(node.sequence)-width+1) if node.sequence[i:i+width]==target]
        if len(hits)>1: raise ValueError('longest exact landing is ambiguous')
        if hits: return Landing(index,hits[0],hits[0]+width,primer)
    raise ValueError('primer has no exact terminal suffix on the specified block')

@dataclass(frozen=True)
class StopHypothesis:
    linker_codes: frozenset[str]
    statement: str
    def __post_init__(self):
        if not isinstance(self.linker_codes,frozenset):
            raise TypeError('stop codes must be an immutable explicit set')
        if not self.statement.strip(): raise ValueError('copying requires an explicit hypothesis/evidence boundary')

@dataclass(frozen=True)
class PathPiece:
    node_index: int
    template_start: int | None
    template_end: int | None

@dataclass(frozen=True)
class Extension:
    product: Backbone
    copied_path: tuple[PathPiece,...]
    stopped_at: int
    closure_crossings: int
    hypothesis: StopHypothesis


def extend_to_linker(template: Backbone, landing: Landing, *, hypothesis: StopHypothesis) -> Extension:
    """Copy 3′→5′ along the template, wrapping the real circle, stopping before a linker.

    Traversal starts immediately before the primer's 3′ landing base. Copying a
    variable DNA block reverses its origin orientation without inventing its bases
    or length. An unqualified linker or a complete turn without a linker raises;
    neither is silently skipped or treated as rolling-circle amplification.
    """
    if not template.circular: raise ValueError('this operation requires an explicitly closed backbone')
    if not 0<=landing.node_index<len(template.nodes): raise ValueError('landing outside backbone')
    node=template.nodes[landing.node_index]
    if not isinstance(node,DNA) or not 0<=landing.start<landing.end<=len(node.sequence):
        raise ValueError('landing interval must lie inside a nucleotide block')
    if node.sequence[landing.start:landing.end]!=revcomp(landing.primer[-(landing.end-landing.start):]):
        raise ValueError('primer suffix does not pair with supplied landing interval')
    copied: list[Node]=[DNA('extension primer',landing.primer)]
    trace=[];index=landing.node_index;end=landing.start;wraps=0
    # At most one complete circuit is searched. No unbounded no-stop circle copying.
    for _ in range(len(template.nodes)+1):
        n=template.nodes[index]
        if isinstance(n,C18):
            if n.code not in hypothesis.linker_codes:
                raise ValueError('linker traversal is unqualified; no automatic bypass')
            return Extension(Backbone(tuple(copied),five_prime='OH',inferred=True),tuple(trace),index,wraps,hypothesis)
        if isinstance(n,VariableDNA):
            copied.append(VariableDNA('copy '+n.name,n.origin,not n.reverse_complemented))
            trace.append(PathPiece(index,None,None))
        elif end:
            copied.append(DNA('copy '+n.name,revcomp(n.sequence[:end]),n.feature))
            trace.append(PathPiece(index,0,end))
        index-=1
        if index<0:
            index=len(template.nodes)-1;wraps+=1
        if index==landing.node_index:
            raise ValueError('no qualified linker before a full circular turn; no RCA assumed')
        next_node=template.nodes[index]
        end=len(next_node.sequence) if isinstance(next_node,DNA) else 0
    raise AssertionError('bounded backbone traversal failed to terminate')
