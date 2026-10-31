"""Substrate-driven nucleotide/connectivity simulator, independent of drawing/timing.

Operations consume existing molecular state; no method accepts an expected product.
This is a selected-reaction simulator, not a predictor of reaction rates or yields.
Enzyme capabilities, substrate choices and unresolved branches are supplied locally.
"""
from __future__ import annotations
from dataclasses import dataclass,asdict,replace
from copy import deepcopy
from functools import wraps
from .model import complement


def atomic(method):
    """Failed preconditions never leave a partly changed molecular state."""
    @wraps(method)
    def wrapped(self,*args,**kw):
        before=self.snapshot();serial=self._serial;history=list(self.history)
        try:return method(self,*args,**kw)
        except Exception:
            for name in ('atoms','next','previous','paired','five','three','names','released','bound','ntp','ntmod','associations'):
                setattr(self,name,deepcopy(getattr(before,name)))
            self._serial=serial;self.history=history
            raise
    return wrapped


@dataclass(frozen=True)
class Ligase:
    name: str
    acceptor_sugars: tuple[str,...]
    donor_sugars: tuple[str,...]
    donor_ends: tuple[str,...]
    required_cofactors: tuple[str,...]
    splint_sugars: tuple[str,...]


@dataclass(frozen=True)
class Nucleotide:
    id: str
    base: str
    sugar: str
    origin: str
    template: str = ''
    features: tuple[str,...] = ()
    symbol: str = ''
    complemented: bool = False
    modification: str = ''


@dataclass(frozen=True)
class Polymerase:
    name: str
    template_sugars: tuple[str,...]
    product_sugar: str
    strand_displacement: bool = False
    displacement_sugars: tuple[str,...] = ("DNA",)
    terminal_bases: tuple[str,...] = ()
    de_novo: bool = False
    primer_sugars: tuple[str,...] = ('DNA',)

    def __post_init__(self):
        if self.product_sugar not in {'DNA','RNA'}:raise ValueError('polymerase product sugar required')
        if not self.template_sugars or set(self.template_sugars)-{'DNA','RNA','LNA'}:raise ValueError('explicit template sugar capabilities required')
        if not self.primer_sugars or set(self.primer_sugars)-{'DNA','RNA','LNA'}:raise ValueError('explicit primer sugar capabilities required')


@dataclass(frozen=True)
class Exonuclease:
    name: str
    direction: str
    target_sugars: tuple[str,...]
    required_end: tuple[str,...]
    substrate: str
    retained_residues: int
    maximum_3prime_overhang: int = 0

    def __post_init__(self):
        if self.direction not in {'5to3','3to5'}or self.substrate not in {'single','duplex','either'}:raise ValueError('explicit nuclease direction/substrate required')
        if self.retained_residues<0:raise ValueError('invalid terminal residue policy')


@dataclass
class Snapshot:
    atoms: dict
    next: dict
    previous: dict
    paired: dict
    five: dict
    three: dict
    names: dict
    released: set
    bound: dict
    ntp: dict
    ntmod: dict
    associations: dict


@dataclass
class Operation:
    kind: str
    before: Snapshot
    after: Snapshot
    added: tuple[str,...]
    source: str
    details: dict
    id: int = 0

    def to_dict(self):
        def edges(mapping):return {tuple(sorted((a,b)))for a,b in mapping.items()}
        return {'operation_id':self.id,'kind':self.kind,'source':self.source,'details':self.details,
                'added':[asdict(self.after.atoms[i])for i in self.added],
                'updated':[asdict(n)for i,n in self.after.atoms.items()if i in self.before.atoms and n!=self.before.atoms[i]],
                'backbone_added':sorted(set(self.after.next.items())-set(self.before.next.items())),
                'backbone_removed':sorted(set(self.before.next.items())-set(self.after.next.items())),
                'pairing_added':sorted(edges(self.after.paired)-edges(self.before.paired)),
                'pairing_removed':sorted(edges(self.before.paired)-edges(self.after.paired)),
                'released':sorted(self.after.released-self.before.released),
                'five_added':{i:v for i,v in self.after.five.items()if self.before.five.get(i)!=v},
                'five_removed':sorted(set(self.before.five)-set(self.after.five)),
                'three_added':{i:v for i,v in self.after.three.items()if self.before.three.get(i)!=v},
                'three_removed':sorted(set(self.before.three)-set(self.after.three)),
                'protein_added':{i:sorted(v)for i,v in self.after.bound.items()if self.before.bound.get(i)!=v},
                'protein_removed':sorted(set(self.before.bound)-set(self.after.bound)),
                'nucleotide_pools':{i:sorted(v)for i,v in self.after.ntp.items()},
                'input_names':{i:v for i,v in self.after.names.items()if i not in self.before.names}}


class Simulator:
    def __init__(self):
        self.atoms={};self.next={};self.previous={};self.paired={}
        self.five={};self.three={};self.names={};self.released=set();self.history=[];self._serial=0;self.bound={};self.ntp={};self.ntmod={};self.associations={}

    def snapshot(self):
        return Snapshot(dict(self.atoms),dict(self.next),dict(self.previous),dict(self.paired),dict(self.five),dict(self.three),dict(self.names),set(self.released),deepcopy(self.bound),deepcopy(self.ntp),deepcopy(self.ntmod),deepcopy(self.associations))

    def _record(self,kind,before,source,**details):
        if not source:raise ValueError('reaction needs explicit source or illustrative scope')
        op=Operation(kind,before,self.snapshot(),tuple(i for i in self.atoms if i not in before.atoms),source,details,len(self.history))
        self.history.append(op);self.validate();return op

    def _new(self,base,sugar,origin,**kw):
        if base not in 'ACGTUNRYSWKMBDHV' or sugar not in {'DNA','RNA','LNA'}:raise ValueError('invalid nucleotide')
        if sugar=='RNA' and base=='T':raise ValueError('RNA nucleotide must use U')
        self._serial+=1;id=f'n{self._serial}'
        self.atoms[id]=Nucleotide(id,base,sugar,origin,**kw);return id

    def _link(self,a,b):
        if a in self.next or b in self.previous:raise ValueError('backbone already occupied')
        self.next[a]=b;self.previous[b]=a;self.three.pop(a,None);self.five.pop(b,None)

    def _break(self,a,b):
        if self.next.get(a)!=b or self.previous.get(b)!=a:raise ValueError('backbone bond does not exist')
        del self.next[a];del self.previous[b];self.three[a]='OH';self.five[b]='P'

    @atomic
    def add(self,name,sequence,*,sugar='DNA',five='OH',three='OH',features=(),sugars=(),modifications=None,symbols=None,source):
        if name in self.names or not sequence:raise ValueError('new nonempty reagent identity required')
        before=self.snapshot();ids=[]
        for j,base in enumerate(sequence):
            typ=next((kind for lo,hi,kind in sugars if lo<=j<hi),sugar)
            tags=tuple(tag for lo,hi,tag in features if lo<=j<hi)
            symbol=(symbols or {}).get(j,f'{name}:{j}')if base not in 'ACGTU'else''
            ids.append(self._new(base,typ,name,features=tags,symbol=symbol,modification=(modifications or {}).get(j,'')))
        for a,b in zip(ids,ids[1:]):self._link(a,b)
        self.five[ids[0]]=five;self.three[ids[-1]]=three;self.names[name]=ids[0]
        self._record('reagent',before,source,name=name,ids=ids);return tuple(ids)

    def resolve(self,ref):
        id=self.names.get(ref,ref)
        if id not in self.atoms or id in self.released:raise ValueError('missing/released molecular substrate')
        return id

    def chain(self,ref):
        id=self.resolve(ref);seen=set();head=id
        while head in self.previous:
            if head in seen:
                head=id;break
            seen.add(head);head=self.previous[head]
        out=[];cursor=head
        while cursor not in out:
            out.append(cursor)
            if cursor not in self.next:break
            cursor=self.next[cursor]
        return tuple(out)

    def sequence(self,ref):return ''.join(self.atoms[i].base for i in self.chain(ref))
    def tail(self,ref):
        ids=self.chain(ref)
        if ids[-1]in self.next:raise ValueError('closed molecule has no free3prime end')
        return ids[-1]

    def _compatible(self,a,b,allow_ambiguous=False):
        a,b=self.atoms[a],self.atoms[b]
        if a.base in 'ACGTU'and b.base in 'ACGTU':return a.base.replace('U','T')==complement(b.base)
        if a.symbol and a.symbol==b.symbol:return a.complemented!=b.complemented
        return allow_ambiguous

    def _pair(self,a,b,allow_ambiguous=False):
        if a==b or a in self.paired or b in self.paired:raise ValueError('base already has a duplex partner')
        if not self._compatible(a,b,allow_ambiguous):raise ValueError('noncomplementary or unresolved base pair')
        self.paired[a]=b;self.paired[b]=a

    @atomic
    def unpair(self,ids,*,source):
        before=self.snapshot()
        for i in ids:
            if i in self.paired:
                other=self.paired.pop(i);del self.paired[other]
        return self._record('unpair',before,source,ids=list(ids))

    @atomic
    def anneal(self,primer,template,*,binding,template_site=None,allow_ambiguous=False,expand=False,source):
        """Bind a primer3prime suffix; locate its unique actual template site if omitted."""
        before=self.snapshot();p=self.chain(primer);t=self.chain(template)
        if not 0<binding<=len(p):raise ValueError('invalid physical primer-binding interval')
        circular=self.next.get(t[-1])==t[0]
        if template_site is None:
            hits=[j for j in range(0 if circular else binding-1,len(t))if all(self._compatible(p[len(p)-binding+k],t[(j-k)%len(t)],allow_ambiguous)for k in range(binding))]
            if len(hits)!=1:raise ValueError(f'primer has {len(hits)} physical binding sites; select explicitly')
            template_site=hits[0]
        if not (0 if circular else binding-1)<=template_site<len(t):raise ValueError('binding leaves physical template')
        if expand:
            while binding<len(p) and template_site+1<len(t) and self._compatible(p[len(p)-binding-1],t[template_site+1],allow_ambiguous):
                binding+=1;template_site+=1
        for k in range(binding):self._pair(p[len(p)-binding+k],t[(template_site-k)%len(t)],allow_ambiguous)
        return self._record('anneal',before,source,primer=list(p),template=list(t),binding=binding,template_site=template_site)

    @atomic
    def bind_spans(self,a,b,*,source,allow_ambiguous=False):
        before=self.snapshot();a=tuple(a);b=tuple(b)
        if not a or len(a)!=len(b):raise ValueError('matching nonempty binding spans required')
        if any(self.next.get(x)!=y for x,y in zip(a,a[1:]))or any(self.previous.get(x)!=y for x,y in zip(b,b[1:])):raise ValueError('binding spans must be covalently contiguous and antiparallel')
        for x,y in zip(a,b):self._pair(x,y,allow_ambiguous)
        return self._record('anneal',before,source)

    @atomic
    def duplex(self,a,b,*,allow_ambiguous=False,source):
        before=self.snapshot();a=self.chain(a);b=self.chain(b)
        if len(a)!=len(b):raise ValueError('duplex lengths differ')
        for x,y in zip(a,reversed(b)):self._pair(x,y,allow_ambiguous)
        return self._record('duplex',before,source)

    @atomic
    def extend(self,primer,enzyme,*,max_bases=None,source):
        """Walk actual template predecessors. No expected product or output length."""
        before=self.snapshot();tail=self.tail(primer)
        if self.atoms[tail].sugar not in enzyme.primer_sugars:raise ValueError('polymerase does not accept this primer sugar')
        if self.three.get(tail)!='OH':raise ValueError('primer3prime is not extendable')
        if tail not in self.paired:raise ValueError('extension needs an annealed primer3prime')
        template=self.paired[tail];copied=[];displaced=[];incorporations=[];visited=set();stop='physical_template_end'
        if max_bases is not None and (not isinstance(max_bases,int)or max_bases<0):raise ValueError('invalid synthesis budget')
        while template in self.previous:
            if max_bases is not None and len(copied)>=max_bases:stop='synthesis_budget';break
            next_template=self.previous[template]
            if next_template in visited and max_bases is None:raise ValueError('circular copying needs an explicit synthesis budget')
            visited.add(next_template)
            nt=self.atoms[next_template]
            if nt.sugar not in enzyme.template_sugars:raise ValueError('polymerase cannot read this template sugar')
            if any(next_template in ids for ids in self.bound.values()):stop='protein_occupied';break
            occupied=self.paired.get(next_template)
            if occupied is not None:
                if not enzyme.strand_displacement or self.atoms[occupied].sugar not in enzyme.displacement_sugars:stop='occupied_template';break
                del self.paired[occupied];del self.paired[next_template];displaced.append(occupied)
            base=complement(nt.base)
            if enzyme.product_sugar=='RNA':base=base.replace('T','U')
            self._require_nt(enzyme.product_sugar,base)
            new=self._new(base,enzyme.product_sugar,enzyme.name,template=next_template,modification=self.ntmod.get(enzyme.product_sugar,{}).get(base,''),features=nt.features,symbol=nt.symbol,complemented=not nt.complemented)
            self._link(tail,new);self.three[new]='OH';self._pair(new,next_template)
            copied.append(new);incorporations.append({'nucleotide':new,'template':next_template,'displaced':occupied or ''});tail=new;template=next_template
        return self._record('extend',before,source,enzyme=asdict(enzyme),copied=copied,displaced=displaced,incorporations=incorporations,stop=stop)

    def extension_sites(self,enzyme,substrates):
        """Discover accessible 3′ ends throughout the supplied molecular complexes.

        Substrates scope the reaction vessel, not which strands count as primers.
        Traversal follows both covalent and duplex bonds, so a template's recessed
        3′ end has exactly the same eligibility as a supplied oligo's 3′ end.
        """
        members=set();pending=[self.resolve(ref)for ref in substrates]
        while pending:
            id=pending.pop()
            if id in members:continue
            members.add(id)
            pending.extend(i for i in (self.next.get(id),self.previous.get(id),self.paired.get(id))if i is not None and i not in members)
        sites=[]
        for id in sorted(members):
            if id in self.next or self.three.get(id)!='OH' or self.atoms[id].sugar not in enzyme.primer_sugars:continue
            mate=self.paired.get(id);template=self.previous.get(mate)
            if template is None or self.atoms[template].sugar not in enzyme.template_sugars:continue
            if any(template in ids for ids in self.bound.values()):continue
            occupied=self.paired.get(template)
            if occupied is not None and (not enzyme.strand_displacement or self.atoms[occupied].sugar not in enzyme.displacement_sugars):continue
            sites.append(id)
        return sites

    @atomic
    def extend_available(self,enzyme,substrates,*,source):
        """Extend every eligible end, then rescan until the vessel is quiescent.

        Unbounded strand-displacing copying of circles fails via extend()'s
        finite-budget gate. RCA uses its explicitly budgeted operation.
        The returned operations retain independent physical strand identities.
        """
        substrates=tuple(substrates);operations=[]
        while sites:=self.extension_sites(enzyme,substrates):
            for site in sites:
                # Another extension can alter access: re-evaluate actual state.
                if site not in self.extension_sites(enzyme,substrates):continue
                op=self.extend(site,enzyme,source=source)
                if not op.added:raise ValueError('eligible extension made no progress')
                op.details['starting_3prime']=site
                operations.append(op)
        return operations

    @atomic
    def rolling_circle(self,primer,enzyme,*,max_bases,source):
        tail=self.tail(primer)
        if tail not in self.paired:raise ValueError('RCA requires an annealed primer3prime')
        template=self.chain(self.paired[tail])
        if self.next.get(template[-1])!=template[0]:raise ValueError('RCA requires a physically covalently closed template')
        if not enzyme.strand_displacement:raise ValueError('continued RCA requires explicit strand displacement capability')
        if not isinstance(max_bases,int)or max_bases<1:raise ValueError('RCA needs a finite illustrative synthesis budget')
        op=self.extend(primer,enzyme,max_bases=max_bases,source=source)
        op.kind='rca';op.details['closed_template']=list(template)
        return op

    @atomic
    def terminal_add(self,primer,base,count,enzyme,*,source):
        before=self.snapshot();tail=self.tail(primer)
        if self.three.get(tail)!='OH' or base not in enzyme.terminal_bases:raise ValueError('terminal addition unsupported for this end/enzyme')
        if not isinstance(count,int)or count<1:raise ValueError('terminal addition needs positive illustrative extent')
        for _ in range(count):
            self._require_nt(enzyme.product_sugar,base)
            new=self._new(base,enzyme.product_sugar,enzyme.name);self._link(tail,new);self.three[new]='OH';tail=new
        return self._record('terminal_add',before,source,count=count,base=base,enzyme=asdict(enzyme))

    @atomic
    def cleave(self,left,*,source):
        before=self.snapshot();left=self.resolve(left)
        if left not in self.next:raise ValueError('cannot cleave a missing backbone bond')
        self._break(left,self.next[left]);return self._record('cleave',before,source)

    @atomic
    def ligate(self,acceptor,donor,enzyme,*,mode,cofactors,source):
        before=self.snapshot();acceptor=self.resolve(acceptor);donor=self.resolve(donor)
        if acceptor in self.next or donor in self.previous:raise ValueError('ligation requires two actual free termini')
        if self.three.get(acceptor)!='OH' or self.five.get(donor)not in {'P','App'}:raise ValueError('ligation requires3OH and5P/App')
        if self.atoms[acceptor].sugar not in enzyme.acceptor_sugars or self.atoms[donor].sugar not in enzyme.donor_sugars or self.five[donor]not in enzyme.donor_ends:
            raise ValueError('ligase substrate/end specificity does not match')
        if set(enzyme.required_cofactors)-set(cofactors):raise ValueError('missing ligase cofactor')
        if mode not in {'nick','splinted','single-strand','blunt'}:raise ValueError('explicit ligation geometry required')
        if mode in {'nick','splinted'}:
            a,b=self.paired.get(acceptor),self.paired.get(donor)
            if a is not None and self.atoms[a].sugar not in enzyme.splint_sugars:raise ValueError('ligase splint sugar not supported')
            if a is None or b is None or self.previous.get(a)!=b:raise ValueError('a gap or unpaired ends cannot be sealed as a nick')
        if mode=='blunt':
            a,b=self.paired.get(acceptor),self.paired.get(donor)
            if a is None or b is None or a in self.previous or b in self.next:raise ValueError('blunt ligation requires two fully paired blunt termini')
        self._link(acceptor,donor);return self._record('ligate',before,source,mode=mode,enzyme=asdict(enzyme),cofactors=list(cofactors))

    @atomic
    def strand_transfer(self,adapter3,target5,*,source):
        """Atomic transesterification, distinct from hydrolysis plus ligase."""
        before=self.snapshot();a=self.resolve(adapter3);t=self.resolve(target5)
        if a in self.next or self.three.get(a)!='OH':raise ValueError('adapter needs a free3primeOH')
        if t not in self.previous:raise ValueError('strand transfer requires an intact internal target bond')
        upstream=self.previous[t];self._break(upstream,t);self._link(a,t)
        return self._record('strand_transfer',before,source,adapter3=a,target5=t,freed3=upstream)

    @atomic
    def tn5(self,adapter1,adapter2,target1,target2,*,stagger,complex_id,source):
        """Validate target duplex correspondence before either backbone changes."""
        a,b=self.resolve(target1),self.resolve(target2);chain=self.chain(a)
        mate=self.paired.get(b)
        if a not in self.previous or b not in self.previous or mate not in chain:raise ValueError('Tn5 needs two internal opposite target bonds')
        if chain.index(mate)-chain.index(a)+1!=stagger:raise ValueError('Tn5 target stagger differs from declared geometry')
        ends=[self.tail(adapter1),self.tail(adapter2)]
        if ends[0]==ends[1]or any(e in self.next or self.three.get(e)!='OH'for e in ends):raise ValueError('Tn5 requires two distinct adapter3OH ends')
        before=self.snapshot()
        if complex_id in self.bound:raise ValueError('complex identity already bound')
        self.bound[complex_id]=set(chain[chain.index(a):chain.index(mate)+1])|{self.paired[i]for i in chain[chain.index(a):chain.index(mate)+1]}
        for e,t in zip(ends,[a,b]):
            upstream=self.previous[t];self._break(upstream,t);self._link(e,t)
        return self._record('tn5',before,source,adapter3=ends,target5=[a,b],stagger=stagger)

    @atomic
    def deadenylate(self,ref,*,source):
        before=self.snapshot();head=self.chain(ref)[0]
        if self.five.get(head)!='App':raise ValueError('deadenylase requires a5adenylated substrate')
        self.five[head]='P'
        return self._record('deadenylate',before,source,released_small_molecule='AMP')

    @atomic
    def resect5(self,ref,count=None,*,required_five,required_sugars=('DNA',),require_duplex=False,source):
        before=self.snapshot();ids=self.chain(ref)
        if ids[0]in self.previous or self.five.get(ids[0])not in required_five:raise ValueError('exonuclease lacks an accessible compatible5prime end')
        if count is not None and(not isinstance(count,int)or not 0<count<=len(ids)):raise ValueError('invalid resection budget')
        removed=[];stop='substrate_end'
        for i in ids:
            if count is not None and len(removed)>=count:stop='illustrative_budget';break
            if self.atoms[i].sugar not in required_sugars:stop='incompatible_sugar';break
            if require_duplex and i not in self.paired:stop='unpaired_substrate';break
            if i in self.paired:
                other=self.paired.pop(i);del self.paired[other]
            if i in self.next:self._break(i,self.next[i])
            self.five.pop(i,None);self.three.pop(i,None);self.released.add(i);removed.append(i)
        remaining=[i for i in ids if i not in self.released]
        for name,id in list(self.names.items()):
            if id in removed:
                if remaining:self.names[name]=remaining[0]
                else:del self.names[name]
        return self._record('resect5',before,source,removed=removed,stop=stop)

    @atomic
    def convert(self,ids,*,from_base,to_base,protected=(),source):
        before=self.snapshot()
        for id in ids:
            nt=self.atoms[self.resolve(id)]
            if nt.base!=from_base or nt.modification in protected:raise ValueError('conversion targets an incompatible/protected nucleotide')
            self.atoms[id]=replace(nt,base=to_base)
        return self._record('convert',before,source,from_base=from_base,to_base=to_base)

    @atomic
    def transcribe(self,template,start,enzyme,*,promoter,source,max_bases=None):
        """Explicit initiation coordinate; promoter bases may themselves be transcribed."""
        if not enzyme.de_novo or enzyme.product_sugar!='RNA':raise ValueError('transcription needs a de-novo RNA polymerase')
        t=self.chain(template)
        if not 0<=start<len(t):raise ValueError('initiation outside physical template')
        if not promoter or any(i not in self.paired for i in promoter):raise ValueError('transcription requires an explicitly paired promoter')
        before=self.snapshot();current=t[start];out=[];seen=set()
        while True:
            if current in seen and max_bases is None:raise ValueError('circular transcription needs synthesis budget')
            seen.add(current);nt=self.atoms[current]
            if nt.sugar not in enzyme.template_sugars:raise ValueError('RNA polymerase cannot read template sugar')
            self._require_nt('RNA',complement(nt.base).replace('T','U'))
            new=self._new(complement(nt.base).replace('T','U'),'RNA',enzyme.name,template=current,features=nt.features,symbol=nt.symbol,complemented=not nt.complemented)
            if out:self._link(out[-1],new)
            out.append(new)
            if current not in self.previous or(max_bases is not None and len(out)>=max_bases):break
            current=self.previous[current]
        self.five[out[0]]='ppp';self.three[out[-1]]='OH'
        op=self._record('transcribe',before,source,enzyme=asdict(enzyme),output=out,initiation=t[start]);return tuple(out),op

    def _require_nt(self,sugar,base):
        if sugar not in self.ntp:raise ValueError('nucleotide substrates not supplied')
        classes={'N':set('ACGT'),'R':set('AG'),'Y':set('CT'),'S':set('CG'),'W':set('AT'),'K':set('GT'),'M':set('AC'),'B':set('CGT'),'D':set('AGT'),'H':set('ACT'),'V':set('ACG')}
        possible=classes.get(base,{base})
        if sugar=='RNA':possible={b.replace('T','U')for b in possible}
        if possible-set(self.ntp[sugar]):raise ValueError('required nucleotide substrate unavailable')

    @atomic
    def supply_nt(self,sugar,bases,*,modifications=None,source):
        before=self.snapshot()
        if sugar not in {'DNA','RNA'}or set(bases)-set('ACGTU'):raise ValueError('invalid nucleotide pool')
        self.ntp.setdefault(sugar,set()).update(bases)
        self.ntmod[sugar]=dict(modifications or {})
        return self._record('nucleotide_supply',before,source,sugar=sugar,bases=list(bases),modifications=dict(modifications or {}))

    @atomic
    def release_protein(self,complex_id,*,source):
        before=self.snapshot()
        if complex_id not in self.bound:raise ValueError('missing bound protein complex')
        del self.bound[complex_id]
        return self._record('protein_removal',before,source,complex_id=complex_id)

    @atomic
    def digest(self,ref,enzyme,*,max_bases=None,source):
        before=self.snapshot();ids=list(self.chain(ref));removed=[];stop='substrate_end'
        if self.next.get(ids[-1])==ids[0]:
            return self._record('digest',before,source,enzyme=asdict(enzyme),removed=[],stop='no_free_end')
        entry=ids[0]if enzyme.direction=='5to3'else ids[-1]
        end=self.five.get(entry)if enzyme.direction=='5to3'else self.three.get(entry)
        if end not in enzyme.required_end:return self._record('digest',before,source,enzyme=asdict(enzyme),removed=[],stop='incompatible_end')
        while len(ids)>max(1,enzyme.retained_residues):
            if max_bases is not None and len(removed)>=max_bases:stop='illustrative_budget';break
            i=ids[0]if enzyme.direction=='5to3'else ids[-1]
            if self.atoms[i].sugar not in enzyme.target_sugars:stop='incompatible_sugar';break
            paired=i in self.paired
            if enzyme.substrate=='single' and paired:stop='paired_substrate';break
            if enzyme.substrate=='duplex' and not paired:
                overhang=0
                for nt in reversed(ids)if enzyme.direction=='3to5'else ids:
                    if nt in self.paired:break
                    overhang+=1
                if overhang>enzyme.maximum_3prime_overhang or overhang==len(ids):stop='inaccessible_duplex_end';break
            old_five=self.five.get(i,'P')
            if i in self.paired:
                mate=self.paired.pop(i);del self.paired[mate]
            if i in self.next:self._break(i,self.next[i])
            if i in self.previous:self._break(self.previous[i],i)
            self.five[i]=old_five;self.three[i]='OH';self.released.add(i);removed.append(i)
            if enzyme.direction=='5to3':ids.pop(0)
            else:ids.pop()
        if len(ids)==enzyme.retained_residues and ids:stop='terminal_residue_limit'
        elif len(ids)==1:stop='no_remaining_backbone'
        for name,id in list(self.names.items()):
            if id in removed:
                if ids:self.names[name]=ids[0]
                else:del self.names[name]
        return self._record('digest',before,source,enzyme=asdict(enzyme),removed=removed,stop=stop)

    @atomic
    def rnase_h1(self,bonds,*,source):
        before=self.snapshot()
        for left in bonds:
            left=self.resolve(left);right=self.next.get(left)
            first=self.previous.get(left);last=self.next.get(right)
            window=[first,left,right,last]
            if any(i is None or self.atoms[i].sugar!='RNA'for i in window):raise ValueError('RNaseH1 needs two consecutive RNA nucleotides on each side of the cut')
            mates=[self.paired.get(i)for i in window]
            if any(i is None or self.atoms[i].sugar!='DNA'for i in mates):raise ValueError('RNaseH1 requires an RNA-DNA hybrid')
            if any(self.previous.get(a)!=b for a,b in zip(mates,mates[1:])):raise ValueError('RNaseH1 needs a continuous antiparallel DNA partner')
            self._break(left,right)
        return self._record('rnase_h1',before,source,cleaved=list(bonds))

    @atomic
    def restrict(self,top,site,recognition,top_cut,bottom_cut,*,source):
        before=self.snapshot();ids=self.chain(top)
        if not 0<=site<=len(ids)-len(recognition) or ''.join(self.atoms[i].base for i in ids[site:site+len(recognition)])!=recognition:
            raise ValueError('restriction recognition site absent from substrate')
        if not 0<site+top_cut<len(ids) or not 0<site+bottom_cut<len(ids):raise ValueError('restriction cuts outside substrate')
        top_right=ids[site+top_cut]
        # The bottom5prime side lies immediately left of its top-coordinate
        # boundary. Use original duplex partners, not display column guesses.
        bottom_right=self.paired.get(ids[site+bottom_cut-1])
        if bottom_right is None or top_right not in self.previous or bottom_right not in self.previous:raise ValueError('restriction needs an intact paired duplex')
        self._break(self.previous[top_right],top_right);self._break(self.previous[bottom_right],bottom_right)
        return self._record('restriction',before,source,recognition=recognition,site=site,top_cut=top_cut,bottom_cut=bottom_cut)

    @atomic
    def capture(self,name,ref,*,required_modification,source):
        before=self.snapshot();ids=self.chain(ref)
        ligands=[i for i in ids if self.atoms[i].modification==required_modification]
        if not ligands:raise ValueError('capture ligand absent from actual molecule')
        self.associations[name]=tuple(ligands)
        return self._record('capture',before,source,association=name,ligands=ligands)

    def validate(self):
        for a,b in self.next.items():
            if self.previous.get(b)!=a or a not in self.atoms or b not in self.atoms:raise ValueError('invalid covalent connectivity')
        for b,a in self.previous.items():
            if self.next.get(a)!=b:raise ValueError('asymmetric backbone')
        for a,b in self.paired.items():
            if self.paired.get(b)!=a or a in self.released or b in self.released:raise ValueError('invalid duplex connectivity')
        if any(i in self.next for i in self.three)or any(i in self.previous for i in self.five):raise ValueError('terminal chemistry on internal nucleotide')
        return self

    def trace(self):return [op.to_dict()for op in self.history]
