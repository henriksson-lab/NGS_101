"""Pure timestamp evaluation. No Cairo, mutable scene history or frame integration."""
from __future__ import annotations
from dataclasses import dataclass
import math
from .model import progress, Pairing


def route_point(strand, boundary, route=None, offset=(0.,0.)):
    path=route if route else strand.path
    if path:
        t=max(0,min(len(strand.sequence),boundary))
        for a,b in zip(path,path[1:]):
            if a[0]<=t<=b[0]:
                f=(t-a[0])/(b[0]-a[0]);return (a[1]+f*(b[1]-a[1])+offset[0],a[2]+f*(b[2]-a[2])+offset[1])
    f=boundary/len(strand.sequence)
    return (strand.start[0]+f*(strand.end[0]-strand.start[0])+offset[0],strand.start[1]+f*(strand.end[1]-strand.start[1])+offset[1])


@dataclass
class State:
    strand: object
    visible: list
    route: list
    offset: tuple
    nicks: dict
    closed: bool
    removed: float
    converted: dict
    synthesized: tuple | None
    displacements: list
    highlights: list
    molecular: object = None

    def point(self,boundary):return route_point(self.strand,boundary,self.route,self.offset)

    def segment_point(self,boundary,member):
        x,y=self.point(boundary)
        for a,b,dx,dy in self.displacements:
            if a <= member < b: x+=dx;y+=dy
        return x,y

    def fraction(self,index):
        if self.molecular is not None:
            id=self.strand.nucleotide_ids[index]
            return float(id in self.molecular.atoms and id not in self.molecular.released)
        return min(1.,sum(max(0.,min(index+1,b)-max(index,a))for a,b in self.visible))


@dataclass
class Actor:
    label: str
    position: tuple
    progress: float
    kind: str
    note: str
    complex: str = ""
    persistent: bool = False
    label_offset: tuple = (0., -61.)
    enzyme_class: str = ''


def evaluate(shot,seconds,wpm=125):
    states={};actors=[]
    molecular=None
    if shot.chemistry_generated:
        from .molecular_frame import molecular_frame
        molecular=molecular_frame(shot,seconds,wpm)
    for s in shot.strands:
        n=len(s.sequence);initial=n if s.initially_visible is None else s.initially_visible
        state=State(s,[(0.,float(initial)),*s.preserve],list(s.path),(0.,0.),{i:1. for i in s.nicks},s.circular,0.,{},None,list(s.pieces),[])
        for e in shot.events:
            if e.strand!=s.id:continue
            p=progress(shot,e,seconds,wpm)
            if seconds < shot.establish+shot.action*e.fraction_start:continue
            if e.kind in {'extend','tail','rca','transcribe'}:
                end=e.copy_end or n;begin=e.primer_length if e.kind in {'extend','tail','rca'} else initial
                visible=begin+(end-begin)*p
                state.visible=[(0.,visible),*s.preserve];state.synthesized=(initial,visible)
            elif e.kind=='fold':
                # Sample original/final routes at common source nucleotide coordinates.
                indices=sorted({0,n,*[k[0]for k in state.route],*[k[0]for k in e.route]})
                state.route=[(i,*tuple(a+(b-a)*p for a,b in zip(route_point(s,i,state.route),route_point(s,i,e.route))))for i in indices]
            elif e.kind in {'anneal','separate','capture','transfer'}:
                state.offset=(state.offset[0]+e.target[0]*p,state.offset[1]+e.target[1]*p)
            elif e.kind=='highlight':
                state.highlights.append((*e.span, .35+.65*math.sin(math.pi*p)**2))
            elif e.kind=='displace':
                state.displacements.append((*e.span,e.target[0]*p,e.target[1]*p))
            elif e.kind=='seal':
                state.nicks[e.bond]=1-p
                if e.bond==0 and p>=1:state.closed=True
            elif e.kind=='cleave':state.nicks[e.bond]=p
            elif e.kind=='resect':
                state.removed=e.amount*p;state.visible=[(state.removed,float(n))]
            elif e.kind=='resect3':
                state.removed=e.amount*p;state.visible=[(0.,n-state.removed)]
            elif e.kind=='convert':
                for j,i in enumerate(e.positions):
                    if p>=(j+1)/(len(e.positions)+1):state.converted[i]=e.replacement
        if molecular is not None:
            state.molecular=molecular
            active=[i for i,id in enumerate(s.nucleotide_ids)if id in molecular.atoms and id not in molecular.released]
            spans=[]
            for i in active:
                if spans and spans[-1][1]==i:spans[-1]=(spans[-1][0],i+1)
                else:spans.append((i,i+1))
            state.visible=spans or [(0,0)]
            state.nicks={i:1. for i in range(1,len(s.nucleotide_ids))if i in active and i-1 in active and (s.nucleotide_ids[i-1],s.nucleotide_ids[i])not in molecular.backbone}
            state.closed=bool(active)and(s.nucleotide_ids[-1],s.nucleotide_ids[0])in molecular.backbone
            state.converted={i:molecular.atoms[id]['base']for i,id in enumerate(s.nucleotide_ids)if id in molecular.atoms and molecular.atoms[id]['base']!=s.sequence[i]}
        states[s.id]=state
    from .layout import following_path,generated_rca_path
    for e in shot.events:
        if e.fraction_start>0 and seconds<shot.establish+shot.action*e.fraction_start:continue
        if e.kind=='follow':
            state=states[e.strand]
            state.route=[(i,x-state.offset[0],y-state.offset[1])for i,x,y in following_path(states[e.template],e.template_start+1,len(state.strand.sequence),e.amount)]
        elif e.kind=='peel':
            state=states[e.strand];template=states[e.template]
            if state.molecular is None:raise ValueError('peeling needs explicit molecular pairing state')
            ids=set(template.strand.nucleotide_ids);mates={a:b for a,b in state.molecular.duplex for a,b in [(a,b),(b,a)]}
            bound=[j for j,id in enumerate(state.strand.nucleotide_ids)if mates.get(id)in ids]
            lo=min(bound)if bound else len(state.strand.sequence)
            if lo:
                x,y=state.point(lo);last=state.point(len(state.strand.sequence));dx,dy=e.target;length=max(.001,math.hypot(dx,dy))
                points=[(0,x+lo*e.amount*dx/length,y+lo*e.amount*dy/length),(lo,x,y)]
                if lo<len(state.strand.sequence):points.append((len(state.strand.sequence),*last))
                state.route=[(i,x-state.offset[0],y-state.offset[1])for i,x,y in points]
    # RCA's 3prime end stays at the polymerase. Older copied intervals peel away
    # into the displaced tail, preserving order instead of growing disconnected
    # far from an enzyme on the template.
    for e in shot.events:
        if e.kind!='rca':continue
        out=states[e.strand];template=states[e.template];visible=out.visible[0][1]
        if shot.chemistry_generated:
            out.route=generated_rca_path(out,template,e)
            continue
        paired=min(20.,visible);peel=visible-paired;n=len(out.strand.sequence)
        def paired_point(i):
            site=(e.template_start+(i-e.primer_length)*e.template_step)%len(template.strand.sequence)
            x,y=template.point(site)
            center=(sum(k[1]for k in template.route)/len(template.route),sum(k[2]for k in template.route)/len(template.route))
            dx,dy=x-center[0],y-center[1];r=max(1,math.hypot(dx,dy))
            return x+dx/r*13,y+dy/r*13
        anchor=paired_point(peel);route=[]
        if peel>0:
            center=(sum(k[1]for k in template.route)/len(template.route),sum(k[2]for k in template.route)/len(template.route))
            angle=math.atan2(anchor[1]-center[1],anchor[0]-center[0])
            outlet_angle=.65
            delta=(angle-outlet_angle+math.pi)%(math.tau)-math.pi
            radius=max(math.dist(anchor,center),155)
            outlet=(center[0]+radius*math.cos(outlet_angle),center[1]+radius*math.sin(outlet_angle))
            vertices=[(outlet[0]+550*peel/n*(1-j/20),outlet[1]+45*math.sin((1-j/20)*peel/n*math.pi*7))for j in range(21)]
            vertices += [(center[0]+radius*math.cos(outlet_angle+delta*j/18),center[1]+radius*math.sin(outlet_angle+delta*j/18))for j in range(1,19)]
            vertices.append(anchor)
            length=sum(math.dist(a,b)for a,b in zip(vertices,vertices[1:]));walk=0
            route=[(0,*vertices[0])]
            for a,b in zip(vertices,vertices[1:]):
                walk+=math.dist(a,b)
                index=peel*walk/length
                if index>route[-1][0]:route.append((index,*b))
        for j in range(25):
            i=peel+paired*j/24
            if route and i<=route[-1][0]:continue
            route.append((i,*paired_point(i)))
        if visible<n:route.append((n,*route[-1][1:]))
        out.route=route
    actor_indices={}
    for e in shot.events:
        if not e.enzyme or (not e.enzyme_prebound and seconds < shot.establish+shot.action*e.fraction_start):continue
        p=progress(shot,e,seconds,wpm);s=states[e.strand]
        if e.kind=='rca':
            pos=s.point(s.visible[0][1])
        elif e.kind in {'extend','tail','transcribe'}:pos=s.point(s.visible[0][1])
        elif e.kind=='resect':pos=s.point(s.removed)
        elif e.kind=='resect3':pos=s.point(len(s.strand.sequence)-s.removed)
        elif e.kind in {'transfer','release'}:pos=s.point(len(s.strand.sequence))
        elif e.bond is not None:pos=s.point(e.bond)
        elif e.kind=='convert'and e.positions:
            pos=s.point(e.positions[min(len(e.positions)-1,int(p*len(e.positions)))]+.5)
        else:pos=s.point(len(s.strand.sequence)/2)
        pos=(pos[0]+e.actor_offset[0],pos[1]+e.actor_offset[1])
        enzyme_class=e.enzyme_class or ('polymerase'if e.kind in {'extend','tail','rca','transcribe'}else'ligase'if e.kind=='seal'else'')
        actor=Actor(e.enzyme if e.enzyme_label is None else e.enzyme_label,pos,p,e.kind,e.note,e.complex,e.persistent,e.enzyme_label_offset,enzyme_class)
        key=(e.strand,e.enzyme)
        if key in actor_indices:actors[actor_indices[key]]=actor
        else:actor_indices[key]=len(actors);actors.append(actor)
    return states,actors


@dataclass
class BindingState:
    pairing: object
    strength: float


def evaluate_pairings(shot,seconds,wpm=125,*,states=None):
    """Declared binding state is independent of strand distance or camera framing."""
    if shot.chemistry_generated:
        if states is None:states,_=evaluate(shot,seconds,wpm)
        locations={id:(s.id,j)for s in shot.strands for j,id in enumerate(s.nucleotide_ids) if states[s.id].fraction(j)>.95}
        molecular=next(iter(states.values())).molecular
        return [BindingState(Pairing(locations[a][0],locations[a][1],locations[b][0],locations[b][1],1,exact=False,id=a+':'+b),1.)
                for a,b in molecular.duplex if a in locations and b in locations]
    bindings=[]
    for pairing in shot.pairings:
        strength=1. if pairing.initially_bound else 0.
        for event in shot.events:
            if event.kind not in {'pair','unpair'} or event.pairing!=pairing.key:continue
            if seconds < shot.establish+shot.action*event.fraction_start:continue
            p=progress(shot,event,seconds,wpm)
            strength=p if event.kind=='pair' else 1-p
        bindings.append(BindingState(pairing,strength))
    return bindings


def evaluate_joins(shot,seconds,wpm=125):
    """Adapter 3prime to target 5prime covalent junctions, separate from docking."""
    return [(e,progress(shot,e,seconds,wpm)) for e in shot.events
            if e.kind=='strand_transfer' and seconds>=shot.establish+shot.action*e.fraction_start]


def presentation(movie,shot,seconds,wpm=125):
    """Interpolate conserved nucleotide geometry during an incoming establish hold.

    Molecular state comes exclusively from evaluate; this moves coordinates only.
    Paths use half-base knots so both letter centers and covalent joins match.
    """
    from dataclasses import replace
    states,actors=evaluate(shot,seconds,wpm);camera=shot.camera
    index=next(i for i,s in enumerate(movie.shots)if s is shot)
    if not index or not shot.transition or seconds>=shot.transition:return states,actors,camera
    previous=movie.shots[index-1]
    old,_=evaluate(previous,previous.timing(wpm)[2],wpm)
    before=next(iter(old.values())).molecular;after=next(iter(states.values())).molecular
    if before is None or before!=after:raise ValueError('reframing requires an unchanged molecular state')
    locations={id:(st,j)for st in old.values()for j,id in enumerate(st.strand.nucleotide_ids)if st.fraction(j)}
    t=max(0,seconds/shot.transition);t=t*t*(3-2*t)
    mix=lambda a,b:tuple(x+(y-x)*t for x,y in zip(a,b))
    # Capture actor attachment before replacing paths.
    anchors=[]
    for actor in actors:
        choices=[(math.dist(actor.position,st.point(j)),st.strand.id,j)for st in states.values()for j in range(len(st.strand.sequence)+1)]
        anchors.append(min(choices))
    geometries={};rigid={}
    for key,st in states.items():
        ids=st.strand.nucleotide_ids;path=[]
        knots=[]
        for k in range(2*len(ids)+1):
            q=k/2;j=min(int(q),len(ids)-1);id=ids[j]
            dest=st.segment_point(q,j);origin=dest
            if id in locations:
                before,at=locations[id];origin=before.segment_point(at+q-j,at)
            knots.append((q,origin,dest))
        # Factor out rigid rotation before deforming the path. A circle must
        # rotate and translate, never collapse halfway through a 180-degree turn.
        count=len(knots);oc=tuple(sum(k[1][i]for k in knots)/count for i in (0,1));dc=tuple(sum(k[2][i]for k in knots)/count for i in (0,1))
        dot=sum((o[0]-oc[0])*(d[0]-dc[0])+(o[1]-oc[1])*(d[1]-dc[1])for _,o,d in knots)
        cross=sum((o[0]-oc[0])*(d[1]-dc[1])-(o[1]-oc[1])*(d[0]-dc[0])for _,o,d in knots)
        angle=math.atan2(cross,dot);center=mix(oc,dc)
        def rotate(v,angle):
            co,si=math.cos(angle),math.sin(angle)
            return co*v[0]-si*v[1],si*v[0]+co*v[1]
        residual=max(math.dist(rotate((o[0]-oc[0],o[1]-oc[1]),angle),(d[0]-dc[0],d[1]-dc[1]))for _,o,d in knots)
        geometries[key]=knots
        if residual<1e-6:rigid[key]=(oc,dc,angle)
    owners={id:key for key,st in states.items()for id in st.strand.nucleotide_ids}
    frame=next(iter(states.values())).molecular
    mates={a:b for a,b in frame.duplex for a,b in [(a,b),(b,a)]}
    for key,st in states.items():
        path=[]
        for q,origin,dest in geometries[key]:
            j=min(int(q),len(st.strand.sequence)-1)
            template=owners.get(mates.get(st.strand.nucleotide_ids[j]))
            transform=rigid.get(key)or rigid.get(template)
            if transform:
                oc,dc,angle=transform;center=mix(oc,dc)
                source=(origin[0]-oc[0],origin[1]-oc[1]);target=rotate((dest[0]-dc[0],dest[1]-dc[1]),-angle)
                point=rotate(mix(source,target),angle*t)
                point=(point[0]+center[0],point[1]+center[1])
            else:point=mix(origin,dest)
            path.append((q,*point))
        st.route=path;st.offset=(0.,0.);st.displacements=[]
        if key in old and st.strand.label_position:
            prior=old[key].strand.label_position
            if prior:st.strand=replace(st.strand,label_position=mix(prior,st.strand.label_position))
    for actor,(_,key,q) in zip(actors,anchors):actor.position=states[key].point(q)
    camera=mix(previous.camera,camera)
    return states,actors,camera
