"""Portable molecular drawing/timing declarations; no rendering dependencies.

Coordinates describe a shot, not atomic geometry. Sequences are always stored 5'–3'.
Exact extension claims are checked independently of presentation choreography.
"""
from __future__ import annotations
from dataclasses import asdict, dataclass, field, replace
import math
import re

SCHEMA = 'ngs101.animation/4'
DEFAULT_SPEED = 1.5
KINDS = {'extend', 'anneal', 'separate', 'seal', 'cleave', 'fold', 'tail', 'bind', 'release', 'rca', 'transcribe', 'resect', 'convert', 'capture', 'transfer', 'highlight', 'displace', 'pair', 'unpair', 'strand_transfer', 'resect3', 'follow', 'peel', 'introduce'}
BASES = set('ACGTUNRYSWKMBDHV')


def finite(value, name):
    if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value):
        raise ValueError(f'{name} must be a finite number')
    return value


def complement(seq):
    return seq.upper().replace('U', 'T').translate(str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN'))


@dataclass
class Region:
    start: int
    end: int
    label: str
    color: str = 'adapter'
    side: int = -1
    # Optional explicit callout position for densely packed molecular anatomy.
    label_position: tuple[float, float] | None = None
    label_offset: tuple[float, float] | None = None


@dataclass
class Strand:
    id: str
    sequence: str
    start: tuple[float, float]
    end: tuple[float, float]
    color: str = 'dna'
    label: str = ''
    polymer: str = 'DNA'
    five_prime: str = 'OH'
    three_prime: str = 'OH'
    initially_visible: int | None = None
    schematic: bool = False
    # Explicit piecewise route: (nucleotide-boundary index, x, y), stored 5'→3'.
    path: list[tuple[float, float, float]] = field(default_factory=list)
    circular: bool = False
    nicks: list[int] = field(default_factory=list)
    marks: dict[str, str] = field(default_factory=dict)
    mark_label_offsets: dict[str,tuple[float,float]] = field(default_factory=dict)
    # Kept intervals beyond a growing gap; coordinates are half-open.
    preserve: list[tuple[int, int]] = field(default_factory=list)
    # Region colors/labels do not infer chemical function.
    bands: list[tuple[int, int, str]] = field(default_factory=list)
    label_position: tuple[float, float] | None = None
    # Native mixed oligos must declare each physical sugar interval explicitly.
    sugars: list[tuple[int, int, str]] = field(default_factory=list)
    # Retained interval -> source molecule id and its 5'→3' starting coordinate.
    lineage: list[tuple[int, int, str, int]] = field(default_factory=list)

    regions: list[Region] = field(default_factory=list)
    pieces: list[tuple[int, int, float, float]] = field(default_factory=list)
    continued_5prime: bool = False
    continued_3prime: bool = False
    nucleotide_ids: list[str] = field(default_factory=list)
    sealed_marks: list[int] = field(default_factory=list)

    def validate(self):
        if any(not isinstance(i,int) or not 0<=i<len(self.sequence) for i in self.sealed_marks):raise ValueError("junction mark leaves strand")
        if self.nucleotide_ids and len(self.nucleotide_ids)!=len(self.sequence):raise ValueError('one nucleotide identity per base required')
        for key,offset in self.mark_label_offsets.items():
            if key not in self.marks or len(offset)!=2:raise ValueError('mark label offset needs an existing mark and two coordinates')
            for v in offset:finite(v,'mark label offset')
        for a,b,dx,dy in self.pieces:
            if not 0<=a<b<=len(self.sequence):raise ValueError('initial piece leaves strand')
            finite(dx,'piece displacement');finite(dy,'piece displacement')
        for r in self.regions:
            if not 0 <= r.start < r.end <= len(self.sequence): raise ValueError("region leaves strand")
            if r.side not in {-1, 1}: raise ValueError("region side must be -1 or 1")
        if not re.fullmatch(r'[a-z][a-z0-9_]*', self.id):
            raise ValueError('invalid strand id')
        if not self.sequence or set(self.sequence.upper()) - BASES:
            raise ValueError(f'{self.id}: use an IUPAC sequence, with schematic=True for unknown classes')
        if self.polymer not in {'DNA', 'RNA', 'mixed'}:
            raise ValueError('polymer must be DNA, RNA or mixed')
        if self.sequence != self.sequence.upper():
            raise ValueError('sequence must be uppercase; annotations belong in marks')
        if self.polymer == 'RNA' and 'T' in self.sequence:
            raise ValueError('RNA should use U rather than T')
        if self.polymer == 'mixed':
            cursor = 0
            for a, b, sugar in self.sugars:
                if a != cursor or not a < b <= len(self.sequence) or sugar not in {'DNA', 'RNA', 'LNA'}:
                    raise ValueError('mixed strand requires a contiguous physical sugar map')
                cursor = b
            if cursor != len(self.sequence): raise ValueError('mixed sugar map must cover strand')

        for p in (self.start, self.end):
            if len(p) != 2:
                raise ValueError('positions need x/y')
            for v in p: finite(v, 'position')
        n = len(self.sequence) if self.initially_visible is None else self.initially_visible
        if not isinstance(n, int) or isinstance(n, bool) or not 0 <= n <= len(self.sequence):
            raise ValueError('visible length outside strand')
        for a,b in self.preserve:
            if not 0 <= a < b <= len(self.sequence): raise ValueError('invalid preserved interval')
        if self.path:
            if self.path[0][0] != 0 or self.path[-1][0] != len(self.sequence):
                raise ValueError('route must cover whole nucleotide-boundary interval')
            if any(b[0] <= a[0] for a,b in zip(self.path,self.path[1:])):
                raise ValueError('route coordinates must increase in strand order')
            for knot in self.path:
                if len(knot)!=3: raise ValueError('route knot needs boundary/x/y')
                for v in knot: finite(v, 'route')
        if any(not isinstance(i,int) or not 0<=i<len(self.sequence) for i in self.nicks):
            raise ValueError('invalid nick coordinate')
        if 0 in self.nicks and not (self.path or self.circular):
            raise ValueError('closure nick needs routed or circular geometry')
        for i in self.marks:
            if not str(i).isdigit() or not 0<=int(i)<len(self.sequence): raise ValueError('invalid base mark')
        if self.circular and (not self.path or math.dist(self.path[0][1:],self.path[-1][1:])>1e-6):
            raise ValueError('closed circle needs a geometrically closed ordered route')
        if self.circular and (self.five_prime!='none' or self.three_prime!='none' or 0 in self.nicks):
            raise ValueError('closed circle has no free termini or closure nick')
        return self


@dataclass
class Junction:
    bond: int
    acceptor_3prime: str = 'OH'
    donor_5prime: str = 'P'
    mode: str = 'nick'
    source: str = ''

    def validate(self):
        if self.acceptor_3prime!='OH' or self.donor_5prime not in {'P','App'}:
            raise ValueError('ligation requires a free 3primeOH and 5prime phosphate or adenylated phosphate')
        if self.mode not in {'nick','splinted','single-strand','blunt'} or not self.source:
            raise ValueError('ligation requires explicit source-scoped end geometry')


@dataclass
class Pairing:
    a: str
    a_start: int
    b: str
    b_start: int
    length: int
    b_step: int = -1
    exact: bool = True
    source: str = ''
    id: str = ''
    initially_bound: bool = True

    @property
    def key(self):
        return self.id or f'{self.a}:{self.a_start}:{self.b}:{self.b_start}:{self.length}'


@dataclass
class Decoration:
    kind: str
    x: float
    y: float
    radius: float
    label: str = ''
    color: str = 'muted'


@dataclass
class Event:
    kind: str
    strand: str
    enzyme: str = ''
    template: str = ''
    # Complementary template bases in their stored 5'→3' coordinates; stepping may be -1.
    template_start: int = 0
    template_step: int = -1
    primer_length: int = 0
    primer_binding_length: int | None = None
    target: tuple[float, float] | None = None
    bond: int | None = None
    fraction_start: float = 0
    fraction_end: float = 1
    exact: bool = False
    note: str = ''
    copy_end: int | None = None
    junction: Junction | None = None
    # New route for a fold; closure is a separate seal event, never inferred.
    route: list[tuple[float, float, float]] = field(default_factory=list)
    positions: list[int] = field(default_factory=list)
    replacement: str = ''
    amount: int = 0
    other: str = ''
    span: tuple[int, int] | None = None
    complex: str = ''
    persistent: bool = False
    enzyme_prebound: bool = False
    enzyme_class: str = ''
    enzyme_label: str | None = None
    actor_offset: tuple[float, float] = (0., 0.)
    enzyme_label_offset: tuple[float, float] = (0., -61.)
    pairing: str = ""
    molecular_operation: dict = field(default_factory=dict)


@dataclass
class Shot:
    id: str
    title: str
    strands: list[Strand]
    events: list[Event]
    narration: str
    cue: str = ''
    duration: float | None = None
    narration_seconds: float | None = None
    establish: float = 1.8
    action: float = 5.5
    hold: float = 2.5
    scope: str = 'Illustrative geometry'
    source: str = ''
    pairings: list[Pairing] = field(default_factory=list)
    decorations: list[Decoration] = field(default_factory=list)
    inventory: dict[str, str] = field(default_factory=dict)
    camera: tuple[float,float,float] = (800.,500.,1.)
    transition: float = 0.
    product_pool: bool = False
    chemical_pairs: list[tuple[str,str]] = field(default_factory=list)
    chemistry_generated: bool = False
    initial_nucleotides: list[str] = field(default_factory=list)
    initial_atom_records: list[dict] = field(default_factory=list)
    initial_backbone: list[tuple[str,str]] = field(default_factory=list)
    initial_released: list[str] = field(default_factory=list)

    def timing(self, words_per_minute=125, *, estimate_narration=True):
        finite(words_per_minute, 'speaking rate')
        if words_per_minute <= 0: raise ValueError('speaking rate must be positive')
        for name in ('establish', 'action', 'hold'):
            finite(getattr(self, name), name)
            if getattr(self, name) < 0: raise ValueError('negative shot timing')
        if self.action <= 0: raise ValueError('action interval must be positive')
        speech = self.narration_seconds
        if speech is None:
            speech = len(self.narration.split()) * 60 / words_per_minute if estimate_narration else 0
        finite(speech, 'narration duration')
        if speech < 0: raise ValueError('negative narration duration')
        minimum = self.establish + self.action + self.hold
        total = max(minimum, speech + 1)
        if self.duration is not None:
            finite(self.duration, 'shot duration')
            if self.duration < minimum or self.duration < speech:
                raise ValueError(f'{self.id}: duration would truncate narration or minimum holds')
            total = self.duration
        # Longer narration extends the final hold, not the polymerase motion.
        return self.establish, self.action, total

    def validate(self):
        if not self.id or not self.title or not self.narration or not self.scope or not self.source:
            raise ValueError('shot needs id, title, narration, scope and source')
        self.timing(estimate_narration=False)
        finite(self.transition,'transition')
        if not 0<=self.transition<=self.establish:raise ValueError('transition must fit before chemistry starts')
        if len(self.camera)!=3 or self.camera[2]<=0:raise ValueError('camera needs center x/y and positive zoom')
        for v in self.camera:finite(v,'camera')
        by_id = {s.id: s.validate() for s in self.strands}
        if self.chemistry_generated:
            ids=[i for strand in self.strands for i in strand.nucleotide_ids]
            if any(not strand.nucleotide_ids for strand in self.strands):raise ValueError('simulated view needs physical nucleotide identities')
            for strand in self.strands:
                n=len(strand.sequence)if strand.initially_visible is None else strand.initially_visible
                present=list(strand.nucleotide_ids[:n])+[id for a,b in strand.preserve for id in strand.nucleotide_ids[a:b]]
                if set(present)-set(self.initial_nucleotides):raise ValueError('visible molecule is absent from initial simulated state')
            if len(set(ids))!=len(ids):raise ValueError('one physical nucleotide cannot appear in two drawn molecules')
        if len(by_id) != len(self.strands): raise ValueError('duplicate strand ids')
        for strand in self.strands:
            for a,b,origin,start in strand.lineage:
                if origin not in self.inventory or not 0<=a<b<=len(strand.sequence) or start<0:
                    raise ValueError('invalid molecular identity correspondence')
                if strand.sequence[a:b]!=self.inventory[origin][start:start+b-a]:
                    raise ValueError('retained interval differs from original molecule')
        if len({p.key for p in self.pairings})!=len(self.pairings):raise ValueError('duplicate pairing identity')
        for pairing in self.pairings:
            if pairing.a not in by_id or pairing.b not in by_id or pairing.b_step!=-1 or pairing.length<=0:
                raise ValueError('pairing requires existing antiparallel strands and positive length')
            a,b = by_id[pairing.a],by_id[pairing.b]
            ia=list(range(pairing.a_start,pairing.a_start+pairing.length))
            ib=[pairing.b_start-i for i in range(pairing.length)]
            if min(ia)<0 or max(ia)>=len(a.sequence) or min(ib)<0 or max(ib)>=len(b.sequence):
                raise ValueError('pairing interval leaves molecule')
            if pairing.exact:
                if a.schematic or b.schematic or set(''.join(a.sequence[i]for i in ia)+''.join(b.sequence[i]for i in ib))-set('ACGTU'): raise ValueError('unknown classes cannot establish exact pairing')
                if any(a.sequence[i].replace('U','T') != complement(b.sequence[j]) for i,j in zip(ia,ib)):
                    raise ValueError('declared antiparallel pairing is not complementary')
        for d in self.decorations:
            if d.kind not in {'bead','droplet','cell','nucleosome'}: raise ValueError('unknown decoration')
            for x in (d.x,d.y,d.radius): finite(x,'decoration geometry')
            if d.radius<=0: raise ValueError('decoration radius must be positive')
        affected=set()
        for e in self.events:
            if e.kind not in KINDS or e.strand not in by_id:
                raise ValueError('unknown event/strand')
            for v in (e.fraction_start, e.fraction_end): finite(v, 'event time')
            if not 0 <= e.fraction_start < e.fraction_end <= 1:
                raise ValueError('event times must be ordered within action interval')
            s = by_id[e.strand]
            if e.enzyme_class not in {'','polymerase','ligase','nuclease','transposase'}:raise ValueError('unknown enzyme class')
            if e.kind in {'follow','peel'} and (e.template not in by_id or e.amount<=0):raise ValueError('binding layout requires a template and positive nucleotide spacing')
            if e.kind=='introduce':
                if not self.chemistry_generated or e.molecular_operation.get('kind')!='reagent':raise ValueError('intact reagent introduction requires a generated reagent operation')
                if e.molecular_operation.get('pairing_added') or e.molecular_operation.get('backbone_removed'):raise ValueError('reagent introduction cannot hide reactions')
            if self.chemistry_generated and e.kind in {'extend','tail','seal','cleave','strand_transfer','rca','transcribe','resect','resect3','convert','pair','unpair'} and not e.molecular_operation:
                raise ValueError('chemical animation event lacks a generated reaction')
            if e.kind in {'pair','unpair'} and e.pairing not in {p.key for p in self.pairings}:
                raise ValueError('binding transition needs an explicit pairing identity')
            if e.kind in {'extend','rca','transcribe'}:
                end=len(s.sequence) if e.copy_end is None else e.copy_end
                if not 0 <= e.primer_length < end <= len(s.sequence):
                    raise ValueError('invalid copy interval')
                if e.kind!='transcribe' and (e.primer_length==0 or s.three_prime!='OH'):
                    raise ValueError('extension requires an existing primer with extendable 3primeOH')
                if e.kind=='transcribe' and (s.polymer!='RNA' or e.primer_length!=0):
                    raise ValueError('transcription creates RNA de novo')
                if s.initially_visible!=e.primer_length:
                    predecessors=[x for x in self.events if x.strand==s.id and x.kind=='extend' and x.copy_end==e.primer_length and x.fraction_end<=e.fraction_start]
                    if not predecessors:raise ValueError('extension must start from actual initial primer or an earlier generated extension')
                if e.template not in by_id or e.template_step!=-1:
                    raise ValueError('copying requires explicit antiparallel template coordinates')
                template=by_id[e.template]
                if e.kind=='rca' and not template.circular:
                    raise ValueError('RCA requires a covalently closed template')
                indices=[e.template_start+i*e.template_step for i in range(end-e.primer_length)]
                if e.kind=='rca': indices=[i%len(template.sequence) for i in indices]
                if any(i<0 or i>=len(template.sequence) for i in indices): raise ValueError('copy mapping leaves template')
                binding=e.primer_length if e.primer_binding_length is None else e.primer_binding_length
                if not isinstance(binding,int)or not 0<=binding<=e.primer_length or (e.kind!='transcribe'and not binding):
                    raise ValueError('primer binding length must fit the existing primer')
                if e.exact:
                    begin=e.primer_length-binding
                    primer_indices=[e.template_start+(i-e.primer_length)*e.template_step for i in range(begin,e.primer_length)]
                    if e.kind=='rca':primer_indices=[i%len(template.sequence)for i in primer_indices]
                    if any(i<0 or i>=len(template.sequence)for i in primer_indices):raise ValueError('primer mapping leaves template')
                    template_bases=''.join(template.sequence[i]for i in primer_indices+indices)
                    if s.schematic or template.schematic or set(s.sequence[begin:end]+template_bases)-set('ACGTU'):
                        raise ValueError('unknown classes cannot establish exact copying')
                    copied=complement(template_bases)
                    if s.polymer=='RNA':copied=copied.replace('T','U')
                    elif s.polymer=='mixed':
                        copied=''.join(base.replace('T','U')if any(lo<=i<hi and sugar=='RNA'for lo,hi,sugar in s.sugars)else base for i,base in enumerate(copied,begin))
                    if s.sequence[begin:end]!=copied:raise ValueError('primer or new daughter bases do not complement template mapping')
            if e.kind=='strand_transfer':
                if e.other not in by_id or e.bond is None or not 0<e.bond<len(by_id[e.other].sequence):
                    raise ValueError('strand transfer requires an explicit internal target bond')
                if s.three_prime!='OH':raise ValueError('strand transfer needs adapter 3primeOH')
                if not e.note:raise ValueError('strand transfer requires source-scoped chemistry')
            if e.molecular_operation and e.kind in {'extend','tail','rca','transcribe'}:
                added=e.molecular_operation['added']
                begin=e.primer_length if e.kind!='tail' or e.primer_length else s.initially_visible
                end=e.copy_end or len(s.sequence)
                if ''.join(a['base']for a in added)!=s.sequence[begin:end]:
                    raise ValueError('drawn synthesis differs from generated reaction product')
                if s.nucleotide_ids and [a['id']for a in added]!=s.nucleotide_ids[begin:end]:
                    raise ValueError('drawn nucleotide identities differ from generated reaction')
                if e.kind in {'extend','rca','transcribe'} and by_id[e.template].nucleotide_ids:
                    expected=[]
                    for j in range(len(added)):
                        i=e.template_start-j
                        if e.kind=='rca':i%=len(by_id[e.template].nucleotide_ids)
                        expected.append(by_id[e.template].nucleotide_ids[i])
                    if [a['template']for a in added]!=expected:raise ValueError('drawn copy follows a different template than simulated')
            if self.chemistry_generated and e.kind=='cleave':
                removed=e.molecular_operation['backbone_removed']
                bond=[s.nucleotide_ids[e.bond-1],s.nucleotide_ids[e.bond]]
                if bond not in removed and tuple(bond)not in removed:raise ValueError('drawn cut is not a simulated backbone cleavage')
            if self.chemistry_generated and e.kind=='strand_transfer':
                added=e.molecular_operation['backbone_added']
                bond=[s.nucleotide_ids[-1],by_id[e.other].nucleotide_ids[e.bond]]
                if bond not in added and tuple(bond)not in added:raise ValueError('drawn transfer junction is not simulated')
            if e.kind in {'seal','cleave'}:
                if e.bond is None or not 0<=e.bond<len(s.sequence): raise ValueError('seal/cleave needs a nucleotide bond coordinate')
                if e.kind=='seal':
                    if e.bond not in s.nicks: raise ValueError('ligation seals a preexisting nick, not an intact bond or gap')
                    if e.junction is None: raise ValueError('ligation needs explicit two-end chemistry')
                    e.junction.validate()
                    if e.junction.bond!=e.bond: raise ValueError('junction coordinate mismatch')
                    if e.bond==0 and (s.five_prime!=e.junction.donor_5prime or s.three_prime!=e.junction.acceptor_3prime):
                        raise ValueError('closure requires actual declared terminal modifications')
                    if e.bond==0:
                        route=e.route or s.path
                        if not route or math.dist(route[0][1:],route[-1][1:])>1e-6:
                            raise ValueError('circle closure route must actually meet')
                    if (s.id,e.bond) in affected: raise ValueError('cannot seal one nick twice')
                    affected.add((s.id,e.bond))
                    if s.initially_visible is not None and s.initially_visible<(e.bond or len(s.sequence)):
                        fills=[x for x in self.events if x.strand==s.id and x.kind=='extend' and x.fraction_end<=e.fraction_start and (x.copy_end or len(s.sequence))>=(e.bond or len(s.sequence))]
                        if not fills: raise ValueError('cannot seal a gap before extension')
            if e.kind in {'highlight','displace'}:
                if e.span is None or not 0 <= e.span[0] < e.span[1] <= len(s.sequence): raise ValueError('event needs valid molecular span')
            if e.complex not in {'', 'homodimer'}: raise ValueError('unknown enzyme complex')
            if len(e.enzyme_label_offset)!=2:raise ValueError('enzyme label offset needs x/y')
            for v in e.enzyme_label_offset:finite(v,'enzyme label offset')
            if e.kind in {'anneal','separate','capture','transfer','displace'}:
                if e.target is None or len(e.target)!=2: raise ValueError('movement needs target displacement')
                for x in e.target: finite(x,'movement')
            if e.kind=='fold':
                if not e.route: raise ValueError('fold needs explicit routed strand geometry')
                replace(s,path=e.route).validate()
            if e.kind in {'resect','resect3'} and (not isinstance(e.amount,int) or not 0<e.amount<len(s.sequence)):
                raise ValueError('resection requires explicit removed 5prime interval')
            if e.kind=='convert':
                if not e.positions or e.replacement not in BASES: raise ValueError('conversion requires explicit sites and chemical product')
                for i in e.positions:
                    if not 0<=i<len(s.sequence): raise ValueError('conversion leaves strand')
                    if s.sequence[i]!='C' or s.marks.get(str(i))=='5mC': raise ValueError('selected C conversion cannot modify protected or non-C base')
            if e.kind=='tail' and (s.three_prime!='OH' or s.initially_visible is None):
                raise ValueError('tailing needs an explicit extendable initial strand')
        return self


@dataclass
class Movie:
    title: str
    shots: list[Shot]
    protocol: str = ''
    provenance: dict = field(default_factory=dict)
    schema: str = SCHEMA

    def validate(self):
        if self.schema != SCHEMA: raise ValueError('unsupported animation schema')
        if not self.title or not self.shots: raise ValueError('movie needs title and shots')
        if len({s.id for s in self.shots}) != len(self.shots): raise ValueError('duplicate shot ids')
        for s in self.shots: s.validate()
        return self

    def to_dict(self):
        self.validate()
        result=asdict(self)
        for shot in result["shots"]:shot["initial_released"]=sorted(shot["initial_released"])
        return result

    @classmethod
    def from_dict(cls, d):
        if not isinstance(d,dict)or d.get('schema')!=SCHEMA:raise ValueError('missing or unsupported animation schema')
        d = dict(d)
        shots = []
        for raw in d.pop('shots'):
            raw = dict(raw)
            strands = [Strand(**dict(s, regions=[Region(**r) for r in s.get("regions", [])])) for s in raw.pop('strands')]
            events = []
            for e in raw.pop('events'):
                e=dict(e)
                if e.get('junction') is not None:e['junction']=Junction(**e['junction'])
                events.append(Event(**e))
            raw['pairings']=[Pairing(**x)for x in raw.get('pairings',[])]
            raw['decorations']=[Decoration(**x)for x in raw.get('decorations',[])]
            shots.append(Shot(strands=strands, events=events, **raw))
        return cls(shots=shots, **d).validate()

    def duration(self, wpm=125, speed=DEFAULT_SPEED):
        finite(speed, 'speed')
        if speed <= 0: raise ValueError('speed must be positive')
        return sum(s.timing(wpm)[2] for s in self.shots) / speed

    def at(self, seconds, wpm=125, speed=DEFAULT_SPEED):
        finite(seconds, 'timestamp')
        self.duration(wpm, speed)
        remaining = max(0, seconds * speed)
        for shot in self.shots:
            if remaining < shot.timing(wpm)[2]: return shot, remaining
            remaining -= shot.timing(wpm)[2]
        last = self.shots[-1]
        return last, last.timing(wpm)[2]


def progress(shot, event, seconds, wpm=125):
    establish, action, _ = shot.timing(wpm)
    p = (seconds - establish) / action
    p = max(0., min(1., (p - event.fraction_start) / (event.fraction_end - event.fraction_start)))
    # Ease without frame-to-frame state: same t produces identical geometry.
    return p * p * (3 - 2 * p)
