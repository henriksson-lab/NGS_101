"""Molecular relations for caller-declared DNA examples.

No enzymes, recognition sites, overlap lengths, assembly orders or repair mechanisms
are selected here. Callers own those choices and their source qualifications.
"""
from dataclasses import dataclass, replace

from chemdraw import Construct, Scene, Segment, revcomp
from cohesive_assembly import CohesiveFragment


def strands(fragment: CohesiveFragment) -> tuple[str, str]:
    """Physical 5′→3′ strands of an explicitly supplied cohesive fragment."""
    core = ''.join(s.top for s in fragment.core)
    top = (fragment.left if fragment.left_polarity == '5-prime' else '') + core + (fragment.right if fragment.right_polarity == '3-prime' else '')
    bottom = (revcomp(fragment.right) if fragment.right_polarity == '5-prime' else '') + revcomp(core) + (revcomp(fragment.left) if fragment.left_polarity == '3-prime' else '')
    return top, bottom


@dataclass(frozen=True)
class StrandRetention:
    """Retain two explicit intervals; this does not predict an enzyme's cuts."""
    input_id: str
    output_id: str
    top: str
    bottom: str
    top_span: tuple[int, int]
    bottom_span: tuple[int, int]
    fragment: CohesiveFragment
    statement: str

    def __post_init__(self):
        if not self.statement.strip():
            raise ValueError('strand retention needs caller scope')
        for seq, (start, end) in ((self.top, self.top_span), (self.bottom, self.bottom_span)):
            if not 0 <= start < end <= len(seq):
                raise ValueError('retained interval outside physical strand')
        retained = (self.top[slice(*self.top_span)], self.bottom[slice(*self.bottom_span)])
        if retained != strands(self.fragment):
            raise ValueError('retained strands do not make the declared cohesive fragment')

    def removed(self):
        return tuple(part for seq, (start, end) in
                     ((self.top, self.top_span), (self.bottom, self.bottom_span))
                     for part in (seq[:start], seq[end:]) if part)

    def rows(self):
        return self.fragment.rows()


@dataclass(frozen=True)
class BridgedAssembly:
    """Declared circular order of blunt duplex parts and nonincorporated splints."""
    fragments: tuple[tuple[str, Construct], ...]
    bridges: tuple[tuple[str, str], ...]  # physical 5′→3′ splints, one per seam
    arm_nt: int
    five_prime_phosphorylated: bool
    three_prime_hydroxyl: bool
    statement: str
    availability = 'known'
    ready_label = 'Declared bridge-supported ends before sealing'
    action_label = 'Pair the caller-declared splints across each seam'

    def __post_init__(self):
        if not self.five_prime_phosphorylated or not self.three_prime_hydroxyl:
            raise ValueError('ligation requires caller-declared 5-prime phosphate and 3-prime hydroxyl ends')
        if len(self.fragments) < 2 or len(self.bridges) != len(self.fragments):
            raise ValueError('one bridge is required for each circular seam')
        if self.arm_nt <= 0 or not self.statement.strip():
            raise ValueError('bridge arms and scope must be explicit')
        if len({i for i,_ in self.fragments}) != len(self.fragments) or len({i for i,_ in self.bridges}) != len(self.bridges):
            raise ValueError('fragment and bridge IDs must be unique')
        for (_,a), (_,b), (_,bridge) in zip(self.fragments,self.fragments[1:]+self.fragments[:1],self.bridges):
            if min(len(a),len(b)) < self.arm_nt:
                raise ValueError('bridge arm exceeds its fragment')
            expected = revcomp(a.top()[-self.arm_nt:] + b.top()[:self.arm_nt])
            if bridge != expected:
                raise ValueError('bridge does not complement the declared adjacent ends')

    def product(self, name):
        return Construct(tuple(replace(s,name=i+': '+s.name) for i,c in self.fragments for s in c),name)

    def rows(self):
        result=[]
        for (ai,a),(bi,b),(bridge_id,bridge) in zip(self.fragments,self.fragments[1:]+self.fragments[:1],self.bridges):
            left=a.top()[-self.arm_nt:];right=b.top()[:self.arm_nt]
            sc=Scene()
            sc.strand('left',[Segment('left end',left,inferred=True)],label=ai)
            sc.anneal('bridge',[Segment('right arm',bridge[:self.arm_nt],inferred=True),Segment('left arm',bridge[self.arm_nt:],inferred=True)],to='left',pair=('left arm','left end'),label=bridge_id)
            sc.anneal('right',[Segment('right end',right,inferred=True)],to='bridge',pair=('right end','right arm'),label=bi)
            sc.footer('Unsealed seam; bridge is a template, not product DNA.','left')
            result.extend(sc.rows())
        return result


def dna_window(con: Construct, start: int, end: int, *, prefix: str) -> tuple[Segment, ...]:
    """Slice explicit DNA coordinates while retaining component annotations."""
    if not 0 <= start <= end <= len(con) or any(s.placeholder for s in con):
        raise ValueError('coordinate slicing requires explicit, bounded DNA')
    result=[];pos=0
    for seg in con:
        lo=max(start-pos,0);hi=min(end-pos,len(seg));pos+=len(seg)
        if lo < hi:
            result.append(replace(seg,name=prefix+seg.name,top=seg.top[lo:hi],
                                  bottom=seg.bottom[lo:hi] if seg.bottom is not None else None,
                                  length_bp=hi-lo if seg.length_bp is not None else None))
    return tuple(result)


def pcr_product(geometry, forward: Construct, reverse: Construct, *, name: str) -> Construct:
    """Reference amplicon from explicit primer polymers and inward binding sites.

    Untemplated terminal additions, error rates and PCR yield are not predicted.
    """
    from chemdraw import complement_segments
    if geometry.availability!='exact' or not isinstance(geometry.template,Construct):
        raise ValueError('amplicon derivation needs exact template and binding sites')
    if forward.top()!=geometry.forward.sequence or reverse.top()!=geometry.reverse.sequence:
        raise ValueError('primer polymers differ from the binding geometry')
    middle=dna_window(geometry.template,geometry.forward_hit.end,geometry.reverse_hit.start,prefix='template: ')
    return Construct(tuple(replace(s,name='forward: '+s.name) for s in forward)+middle+
                     tuple(replace(s,name='reverse copy: '+s.name) for s in complement_segments(tuple(reverse))),name)


@dataclass(frozen=True)
class OverlapExtension:
    """Two supplied linear references share a specified terminal overlap."""
    left_id: str
    right_id: str
    left: Construct
    right: Construct
    overlap_nt: int
    statement: str

    def __post_init__(self):
        n=self.overlap_nt
        if not self.statement.strip() or not 0<n<min(len(self.left),len(self.right)):
            raise ValueError('overlap and scope must be explicit and bounded')
        if self.left.top()[-n:]!=self.right.top()[:n] or set(self.left.top()+self.right.top())-set('ACGT'):
            raise ValueError('terminal overlap requires equal explicit DNA')
        if any(s.feature is not None for s in dna_window(self.left,len(self.left)-n,len(self.left),prefix='')+dna_window(self.right,0,n,prefix='')):
            raise ValueError('overlapping identifier features require an explicit annotation decision')

    def product(self,name):
        return Construct(tuple(replace(s,name='left: '+s.name) for s in self.left)+
                         dna_window(self.right,self.overlap_nt,len(self.right),prefix='right: '),name)

    def rows(self):
        n=self.overlap_nt
        left_body=dna_window(self.left,0,len(self.left)-n,prefix='left ')
        overlap=dna_window(self.left,len(self.left)-n,len(self.left),prefix='overlap ')
        right_body=dna_window(self.right,n,len(self.right),prefix='right ')
        from chemdraw import complement_segments
        top=left_body+overlap
        # Bottom is physical 5′→3′: right-body complement followed by overlap complement.
        bottom=tuple(complement_segments(overlap+right_body))
        sc=Scene();sc.strand('left',top,label='left fragment active strand')
        sc.anneal('right',bottom,to='left',pair=(overlap[0].name+"'",overlap[0].name),label='right fragment active strand')
        sc.footer('Paired 3-prime ends can extend on the opposite template; one overlap is retained.','left')
        return sc.rows()


@dataclass(frozen=True)
class NucleotideBoundary(StrandRetention):
    """Caller-declared 3′ retention boundary and excluded-base suffix relation.

    The boundary base is data; no enzyme or nucleotide identity is selected here.
    This checks composition, not an enzyme's stopping distribution.
    """
    boundary_base: str

    def __post_init__(self):
        super().__post_init__()
        if self.boundary_base not in 'ACGT' or len(self.boundary_base)!=1:
            raise ValueError('one explicit boundary base required')
        for seq,span,retained in zip((self.top,self.bottom),(self.top_span,self.bottom_span),strands(self.fragment)):
            if span[0]!=0 or not retained.endswith(self.boundary_base) or self.boundary_base in seq[span[1]:]:
                raise ValueError('declared 3-prime boundary or excluded-base suffix is inconsistent')


@dataclass(frozen=True)
class DuplexEndState:
    """Explicit physical strand end groups on a caller-supplied duplex."""
    fragment: CohesiveFragment
    top5: str
    top3: str
    bottom5: str
    bottom3: str

    def __post_init__(self):
        allowed={'OH','phosphate','phosphotyrosyl','blocked','unassigned'}
        if any(group not in allowed for group in (self.top5,self.top3,self.bottom5,self.bottom3)):
            raise ValueError('end groups must be explicit or unassigned')

    def rows(self):
        scene=self.fragment.scene()
        scene.strands['top'].mod5=self.top5;scene.strands['top'].mod3=self.top3
        scene.strands['bottom'].mod5=self.bottom5;scene.strands['bottom'].mod3=self.bottom3
        return scene.rows()


@dataclass(frozen=True)
class ActivatedEndTransfer:
    """Declared 3′-phosphotyrosyl donors and incoming 5′-OH acceptors.

    This relation neither chooses an enzyme nor asserts that all circular seams
    become sealed. The caller must model remaining strand-end fates separately.
    """
    vector_id: str
    insert_id: str
    vector: DuplexEndState
    insert: DuplexEndState
    statement: str

    def __post_init__(self):
        if not self.statement.strip() or self.vector_id==self.insert_id:
            raise ValueError('distinct contributors and a scope statement required')
        if (self.vector.top3,self.vector.bottom3)!=('phosphotyrosyl','phosphotyrosyl'):
            raise ValueError('both donor 3-prime ends must be explicitly activated')
        if (self.insert.top5,self.insert.bottom5)!=('OH','OH'):
            raise ValueError('incoming physical 5-prime ends must be hydroxyl acceptors')
        from cohesive_assembly import CircularCohesiveAssembly
        CircularCohesiveAssembly(((self.insert_id,self.insert.fragment),(self.vector_id,self.vector.fragment)),self.statement)

    def reference(self,name):
        from cohesive_assembly import CircularCohesiveAssembly
        return CircularCohesiveAssembly(((self.insert_id,self.insert.fragment),(self.vector_id,self.vector.fragment)),self.statement).product(name)


@dataclass(frozen=True)
class LinearFivePrimeAssembly:
    """Join a supplied linear order, retaining both external 5′ overhangs."""
    fragments: tuple[tuple[str, CohesiveFragment], ...]
    five_prime_phosphorylated: bool
    three_prime_hydroxyl: bool
    statement: str
    availability = 'known'
    ready_label = 'Declared linear cohesive seams before sealing'
    action_label = 'Pair adjacent supplied cohesive ends'

    def __post_init__(self):
        if len(self.fragments)<2 or not self.statement.strip():
            raise ValueError('a linear assembly needs fragments and scope')
        if len({i for i,_ in self.fragments})!=len(self.fragments):
            raise ValueError('linear contributor IDs must be unique')
        if not self.five_prime_phosphorylated or not self.three_prime_hydroxyl:
            raise ValueError('supplied joining ends must be phosphorylated/hydroxyl terminated')
        if any(f.left_polarity!='5-prime' or f.right_polarity!='5-prime' for _,f in self.fragments):
            raise ValueError('this geometry explicitly supports 5-prime cohesive ends')
        if any(a.right!=b.left for (_,a),(_,b) in zip(self.fragments,self.fragments[1:])):
            raise ValueError('adjacent linear cohesive ends differ')

    def product(self,name):
        core=[]
        for index,(identifier,fragment) in enumerate(self.fragments):
            if index:core.append(Segment(identifier+' retained seam',fragment.left,inferred=fragment.inferred))
            core.extend(replace(s,name=identifier+': '+s.name) for s in fragment.core)
        return CohesiveFragment(name,self.fragments[0][1].left,tuple(core),self.fragments[-1][1].right,'5-prime','5-prime',inferred=any(f.inferred for _,f in self.fragments))

    def rows(self):
        # Each input duplex is individually paired by Scene. The constructor checks
        # the adjacent port identities; no circular last-to-first seam is implied.
        result=[]
        for _,fragment in self.fragments:result.extend(fragment.rows())
        return result


@dataclass(frozen=True)
class CircularOverlapReference:
    """A conditional circular reference from explicitly overlapping duplex maps.

    This does not select resection, strand exchange, repair or a biochemical path.
    """
    parts: tuple[tuple[str, Construct, int, int], ...]
    statement: str

    def __post_init__(self):
        if len(self.parts)<2 or not self.statement.strip() or len({i for i,_,_,_ in self.parts})!=len(self.parts):
            raise ValueError('distinct overlapping maps and scope required')
        for _,con,left,right in self.parts:
            if not 0<left<len(con)-right or right<=0 or set(con.top())-set('ACGT'):
                raise ValueError('overlap windows require explicit bounded DNA and an internal body')
            if any(s.feature is not None for s in dna_window(con,0,left,prefix='')+dna_window(con,len(con)-right,len(con),prefix='')):
                raise ValueError('overlapping identifier features require explicit annotation decisions')
        for (_,a,_,right),(_,b,left,_) in zip(self.parts,self.parts[1:]+self.parts[:1]):
            if right!=left or a.top()[-right:]!=b.top()[:left]:
                raise ValueError('declared circular reference overlaps differ')

    def reference(self,name):
        return Construct(tuple(replace(s,inferred=True) for identifier,con,_,right in self.parts for s in dna_window(con,0,len(con)-right,prefix=identifier+': ')),name)


@dataclass(frozen=True)
class DuplexPartition:
    """Partition an explicit duplex reference at supplied top-axis boundaries.

    This compares base accounting only. It assigns no recognition sites, physical
    end offsets, end groups, nuclease specificity or product frequencies.
    """
    input_id: str
    parent: Construct
    intervals: tuple[tuple[str, int, int], ...]
    statement: str

    def __post_init__(self):
        if not self.input_id or not self.statement.strip() or not self.intervals:
            raise ValueError('a partition requires a parent, intervals and scope')
        if len({i for i,_,_ in self.intervals}) != len(self.intervals):
            raise ValueError('partition output IDs must be distinct')
        if set(self.parent.top())-set('ACGT'):
            raise ValueError('partition requires an explicit DNA reference')
        cursor=0
        for _,start,end in self.intervals:
            if start!=cursor or not start<end<=len(self.parent):
                raise ValueError('partition intervals must cover the parent once in order')
            cursor=end
        if cursor!=len(self.parent):
            raise ValueError('partition loses the parent suffix')

    def products(self):
        return tuple((identifier,Construct(dna_window(self.parent,start,end,prefix=identifier+': '),identifier))
                     for identifier,start,end in self.intervals)
