"""Pure timestamp evaluation. No Cairo, mutable scene history or frame integration."""
from __future__ import annotations
from dataclasses import dataclass
import math
from .model import progress


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

    def point(self,boundary):return route_point(self.strand,boundary,self.route,self.offset)

    def segment_point(self,boundary,member):
        x,y=self.point(boundary)
        for a,b,dx,dy in self.displacements:
            if a <= member < b: x+=dx;y+=dy
        return x,y

    def fraction(self,index):
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


def evaluate(shot,seconds,wpm=125):
    states={};actors=[]
    for s in shot.strands:
        n=len(s.sequence);initial=n if s.initially_visible is None else s.initially_visible
        state=State(s,[(0.,float(initial)),*s.preserve],list(s.path),(0.,0.),{i:1. for i in s.nicks},s.circular,0.,{},None,list(s.pieces),[])
        for e in shot.events:
            if e.strand!=s.id:continue
            p=progress(shot,e,seconds,wpm)
            if seconds < shot.establish+shot.action*e.fraction_start:continue
            if e.kind in {'extend','tail','rca','transcribe'}:
                end=e.copy_end or n;begin=e.primer_length if e.kind=='tail'and e.primer_length else initial
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
            elif e.kind=='convert':
                for j,i in enumerate(e.positions):
                    if p>=(j+1)/(len(e.positions)+1):state.converted[i]=e.replacement
        states[s.id]=state
    # RCA's 3prime end stays at the polymerase. Older copied intervals peel away
    # into the displaced tail, preserving order instead of growing disconnected
    # far from an enzyme on the template.
    for e in shot.events:
        if e.kind!='rca':continue
        out=states[e.strand];template=states[e.template];visible=out.visible[0][1]
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
        if not e.enzyme or (not e.persistent and seconds < shot.establish+shot.action*e.fraction_start):continue
        p=progress(shot,e,seconds,wpm);s=states[e.strand]
        if e.kind=='rca':
            template=states[e.template];count=max(0,s.visible[0][1]-e.primer_length)
            site=(e.template_start+e.template_step*count)%len(template.strand.sequence)
            pos=template.point(site)
        elif e.kind in {'extend','tail','transcribe'}:pos=s.point(s.visible[0][1])
        elif e.kind=='resect':pos=s.point(s.removed)
        elif e.kind=='transfer':pos=s.point(len(s.strand.sequence))
        elif e.bond is not None:pos=s.point(e.bond)
        elif e.kind=='convert'and e.positions:
            pos=s.point(e.positions[min(len(e.positions)-1,int(p*len(e.positions)))]+.5)
        else:pos=s.point(len(s.strand.sequence)/2)
        pos=(pos[0]+e.actor_offset[0],pos[1]+e.actor_offset[1])
        actor=Actor(e.enzyme,pos,p,e.kind,e.note,e.complex,e.persistent)
        key=(e.strand,e.enzyme)
        if key in actor_indices:actors[actor_indices[key]]=actor
        else:actor_indices[key]=len(actors);actors.append(actor)
    return states,actors
