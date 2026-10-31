"""Native Cairo paths/text; molecular state comes from the pure evaluator."""
from __future__ import annotations
from dataclasses import dataclass
import math
from .evaluate import evaluate
from .model import DEFAULT_SPEED

PALETTE={'background':'#101b2c','text':'#eef4fc','muted':'#9cacc2','dna':'#58d9cc',
         'rna':'#f8b869','primer':'#9d9bff','new':'#63c9ff','enzyme':'#6689e8',
         'ligase':'#ec9cba','adapter':'#c9a1ef','bio':'#eed17d','waste':'#6d7d94'}


@dataclass
class Style:
    width:int=1920
    height:int=1080
    bases:bool=False
    directions:bool=False
    evidence:bool=False
    coordinates:bool=False
    reagents:bool=False
    font:str='DejaVu Sans'


def rgb(key):
    value=PALETTE.get(key,key).lstrip('#')
    if len(value)!=6:raise ValueError(f'invalid animation color {key!r}')
    return tuple(int(value[i:i+2],16)/255 for i in (0,2,4))


class Renderer:
    def __init__(self,style=None):
        try:import cairo
        except ImportError as exc:raise RuntimeError('Animation rendering needs pycairo; see animation/README.md')from exc
        self.cairo=cairo;self.style=style or Style()
        for v in (self.style.width,self.style.height):
            if not isinstance(v,int)or isinstance(v,bool)or v<64:raise ValueError('frame dimensions must be integers >=64')
        if abs(self.style.width/self.style.height-16/9)>.01:
            raise ValueError('use 16:9 dimensions to preserve presentation layout')
        self.overflows=[]

    def text(self,c,text,x,y,size=24,color='text',align='left',max_width=None):
        c.set_source_rgb(*rgb(color));c.select_font_face(self.style.font);c.set_font_size(size)
        ext=c.text_extents(text)
        if max_width and ext[2]>max_width:
            desired=size*max_width/ext[2]
            size=max(14,desired);c.set_font_size(size)
            while len(text)>3 and c.text_extents(text)[2]>max_width:text=text[:-2]
            if desired<14:text=text.rstrip()+'…'
            ext=c.text_extents(text)
        if align=='center':x-=ext[2]/2+ext[0]
        elif align=='right':x-=ext[2]+ext[0]
        clip=c.clip_extents()
        if (x+ext[0]<max(clip[0],30)-1 or x+ext[0]+ext[2]>min(clip[2],1570)+1
                or y+ext[1]<max(clip[1],15)-1 or y+ext[1]+ext[3]>min(clip[3],885)+1):
            self.overflows.append(text)
        c.move_to(x,y);c.show_text(text)

    def line(self,c,points,color,width=5,alpha=1,rounded=True):
        if len(points)<2:return
        c.set_source_rgba(*rgb(color),alpha);c.set_line_width(width)
        c.set_line_cap(self.cairo.LINE_CAP_ROUND if rounded else self.cairo.LINE_CAP_BUTT);c.set_line_join(self.cairo.LINE_JOIN_ROUND)
        c.move_to(*points[0])
        for pt in points[1:]:c.line_to(*pt)
        c.stroke()

    def enzyme(self,c,actor):
        p=actor.progress
        # A brief docking/undocking motion, not measured kinetics or Brownian diffusion.
        dock=p if actor.kind=='bind'else 1-p if actor.kind=='release'else min(1,p/.13,(1-p)/.13)
        x,y=actor.position
        if actor.complex != 'homodimer':y-=12+70*(1-max(0,1 if actor.persistent else dock))
        color='ligase'if actor.kind=='seal'else'enzyme'
        if actor.complex == 'homodimer':
            alpha=1-p if actor.kind=='release' else 1
            if alpha<=0:return
            if actor.kind=='release':y-=100*p
            # Two identical coupled subunits, one complex, with a central DNA cleft.
            self.line(c,[(x-36,y-30),(x+36,y+30)],'enzyme',20,alpha)
            for dx,dy in [(-36,-30),(36,30)]:
                c.set_source_rgba(*rgb('enzyme'),.94*alpha)
                c.new_path();c.arc(x+dx,y+dy,32,0,math.tau);c.fill()
                self.line(c,[(x+dx-12,y+dy+12),(x+dx+12,y+dy-12)],'background',7)
            self.text(c,actor.label,x,max(280,y-64),22,align='center',max_width=390)
            return
        c.save();c.translate(x,y);c.rotate(-.12)
        # A soft, asymmetric clamp rather than a protein structural claim.
        c.set_source_rgba(*rgb(color),.96)
        c.new_path();c.move_to(-42,8);c.curve_to(-60,-33,-20,-51,17,-36)
        c.curve_to(58,-43,64,18,27,32);c.curve_to(9,49,-29,38,-42,8);c.fill()
        c.set_source_rgba(*rgb('background'),.55);c.set_line_width(3)
        c.new_path();c.move_to(-12,21);c.line_to(15,-15);c.set_line_width(8);c.set_line_cap(self.cairo.LINE_CAP_ROUND);c.stroke();c.restore()
        # Keep long enzyme names on screen at either end of a strand.
        label_x=max(210,min(1390,x))
        self.text(c,actor.label,label_x,max(280,y-61),22,align='center',max_width=390)

    def draw(self,surface,movie,seconds,wpm=125,speed=DEFAULT_SPEED):
        self.overflows=[]
        c=self.cairo.Context(surface);c.scale(self.style.width/1600,self.style.height/900)
        c.set_source_rgb(*rgb('background'));c.paint()
        shot,local=movie.at(seconds,wpm,speed);states,actors=evaluate(shot,local,wpm)
        self.text(c,movie.title,88,76,22,'muted',max_width=1400)
        self.text(c,shot.title,88,135,42,max_width=1400)
        self.text(c,shot.scope,88,180,20,'muted',max_width=1400)
        c.save();c.rectangle(40,215,1520,575);c.clip()
        cx,cy,zoom=shot.camera;c.translate(800,500);c.scale(zoom,zoom);c.translate(-cx,-cy)
        for d in shot.decorations:
            c.set_source_rgba(*rgb(d.color),.10 if d.kind=='droplet'else .28)
            c.new_path();c.arc(d.x,d.y,d.radius,0,math.tau);c.fill_preserve()
            c.set_source_rgba(*rgb(d.color),.65);c.set_line_width(2);c.stroke()
            if d.label:self.text(c,d.label,d.x,d.y+d.radius+36,22,d.color,align='center')
        # Binding bars only where the explicitly declared molecules are visible.
        for pairing in shot.pairings:
            a,b=states[pairing.a],states[pairing.b]
            for j in range(pairing.length):
                ia=pairing.a_start+j;ib=pairing.b_start-j
                if a.fraction(ia)<.9 or b.fraction(ib)<.9:continue
                pa=a.point(ia+.5);pb=b.point(ib+.5)
                # Avoid drawing bars during long-range approach; they appear at docking.
                if math.dist(pa,pb)<75:self.line(c,[pa,pb],'muted',1.4,.35)
        for st in states.values():
            s=st.strand;n=len(s.sequence)
            # Route is linear within each interval. Subsample long backbones but keep
            # exact boundaries, color bands, growing frontier and nick sites.
            stride=max(1,n//150)
            bounds={0,n,*range(0,n,stride),*[int(i)for i in st.nicks]}
            for a,b in st.visible:bounds.update((max(0,a),min(n,b)))
            for a,b,_ in s.bands:bounds.update((a,b))
            for r in s.regions:bounds.update((r.start,r.end))
            for a,b,_,_ in st.displacements:bounds.update((a,b))
            bounds.update(k[0]for k in st.route)
            ordered=sorted(x for x in bounds if 0<=x<=n)
            for a,b in zip(ordered,ordered[1:]):
                mid=(a+b)/2
                if not any(v0<=mid<v1 for v0,v1 in st.visible):continue
                color=s.color
                for v0,v1,tag in s.bands:
                    if v0<=mid<v1:color=tag
                if st.synthesized and st.synthesized[0]<=mid<st.synthesized[1]:color='rna' if s.polymer=='RNA' else 'new'
                for r in s.regions:
                    if r.start<=mid<r.end:color=r.color
                left,right=st.segment_point(a,mid),st.segment_point(b,mid)
                for v0,v1,intensity in st.highlights:
                    if v0<=mid<v1:self.line(c,[left,right],color,20,.18*intensity,rounded=False)
                # Bonds are gaps BETWEEN retained bases, not deleted base glyphs.
                gap_a=st.nicks.get(a,0);gap_b=st.nicks.get(b if b<n else 0,0)
                length=max(1,math.dist(left,right));ga=min(.35,8*gap_a/length);gb=min(.35,8*gap_b/length)
                aa=tuple(x+(y-x)*ga for x,y in zip(left,right));bb=tuple(y+(x-y)*gb for x,y in zip(left,right))
                self.line(c,[aa,bb],color,6)
            for r in s.regions:
                if st.fraction(r.start)<.9:continue
                end=min(r.end,max(b for a,b in st.visible))
                if end<=r.start:continue
                pa=st.segment_point(r.start,(r.start+end)/2);pb=st.segment_point(end,(r.start+end)/2)
                mid=((pa[0]+pb[0])/2,(pa[1]+pb[1])/2)
                tx,ty=r.label_position or (mid[0],mid[1]+r.side*65)
                self.line(c,[mid,(tx,ty-r.side*17)],r.color,1.3,.65)
                self.text(c,r.label,tx,ty,20,r.color,align='center',max_width=340)
            # Chemical marks are always shown; these distinguish a converted base from
            # a copied base even when letter-level detail is disabled.
            for raw,label in s.marks.items():
                i=int(raw)
                if st.fraction(i)<.9:continue
                x,y=st.point(i+.5);color='bio'if label=='biotin'else'ligase'
                c.set_source_rgb(*rgb(color));c.new_path();c.arc(x,y,7,0,math.tau);c.fill()
                self.text(c,label,x,y-25,17,color,align='center')
            for i,base in st.converted.items():
                x,y=st.point(i+.5);c.set_source_rgb(*rgb('new'));c.new_path();c.arc(x,y,8,0,math.tau);c.fill()
                self.text(c,base,x,y-23,21,'new',align='center')
            if st.removed:
                # Cleavage products persist as an unexpanded product population.
                for j in range(min(14,max(1,int(st.removed)))):
                    x=s.start[0]+12*j;y=s.start[1]+100+14*math.sin(j)
                    c.set_source_rgba(*rgb('waste'),.65);c.new_path();c.arc(x,y,3,0,math.tau);c.fill()
            if self.style.bases and n<=90 and not s.schematic:
                for i,base in enumerate(s.sequence):
                    if st.fraction(i)<.9:continue
                    pos=st.point(i+.5)
                    self.text(c,st.converted.get(i,base),pos[0],pos[1]+28,17,s.color,align='center')
            if s.label:
                if s.label_position:
                    x,y=s.label_position;align='center'
                elif st.route:
                    xs=[k[1]for k in st.route];ys=[k[2]for k in st.route]
                    x=(min(xs)+max(xs))/2+st.offset[0];y=max(ys)+62+st.offset[1];align='center'
                else:
                    x=min(s.start[0],s.end[0])+st.offset[0]-30
                    y=(s.start[1]+s.end[1])/2+st.offset[1]+8;align='right'
                self.text(c,s.label,x,y,22,s.color,align=align,max_width=380)
            if self.style.directions and not st.closed:
                start=st.point(st.removed);end=st.point(max(b for a,b in st.visible))
                dy=-25 if s.polymer=='RNA'else 52
                self.text(c,'5′',start[0]-20,start[1]+dy,18,s.color)
                self.text(c,'3′',end[0]+12,end[1]+dy,18,s.color)

        for actor in actors:self.enzyme(c,actor)
        c.restore()
        if self.style.coordinates:
            for j,s in enumerate(shot.strands):
                detail=f'{s.id}: {len(s.sequence)} nt'if not s.schematic else f'{s.id}: schematic'
                self.text(c,detail,88,242+21*j,15,'muted',max_width=300)
        if self.style.reagents:
            notes=list(dict.fromkeys(e.note for e in shot.events if e.note))
            for j,note in enumerate(notes[:2]):self.text(c,note,800,730+30*j,19,'muted',align='center',max_width=1360)
        index=movie.shots.index(shot)
        self.text(c,f'{index+1:02d} / {len(movie.shots):02d}',1510,840,19,'muted',align='right')
        if self.style.evidence:self.text(c,shot.source,88,818,14,'muted',max_width=1400)
        self.line(c,[(88,870),(1510,870)],'muted',2,.25)
        self.line(c,[(88,870),(88+1422*min(1,max(0,seconds/movie.duration(wpm,speed))),870)],'dna',3)
        surface.flush();return surface

    def image(self,movie,seconds,wpm=125,speed=DEFAULT_SPEED):
        surface=self.cairo.ImageSurface(self.cairo.FORMAT_RGB24,self.style.width,self.style.height)
        return self.draw(surface,movie,seconds,wpm,speed)

    def png(self,path,movie,seconds,wpm=125,speed=DEFAULT_SPEED):self.image(movie,seconds,wpm,speed).write_to_png(str(path))

    def svg(self,path,movie,seconds,wpm=125,speed=DEFAULT_SPEED):
        surface=self.cairo.SVGSurface(str(path),self.style.width,self.style.height)
        self.draw(surface,movie,seconds,wpm,speed);surface.finish()
