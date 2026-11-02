"""Native Cairo paths/text; molecular state comes from the pure evaluator."""
from __future__ import annotations
from dataclasses import dataclass
import math
from .evaluate import evaluate, evaluate_pairings, evaluate_joins, route_point, presentation
from .model import DEFAULT_SPEED

PALETTE={'background':'#101b2c','text':'#eef4fc','muted':'#9cacc2','dna':'#58d9cc',
         'rna':'#f8b869','primer':'#9d9bff','new':'#63c9ff','enzyme':'#6689e8',
         'ligase':'#ec9cba','adapter':'#c9a1ef','bio':'#eed17d','waste':'#6d7d94'}


@dataclass
class Style:
    width:int=1920
    height:int=1080
    bases:bool=True
    directions:bool=True
    evidence:bool=False
    coordinates:bool=False
    reagents:bool=False
    font:str='DejaVu Sans'
    timeline:bool=False


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
        self.in_molecule_layer=False

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
        if self.in_molecule_layer and (x+ext[0]<clip[0] or x+ext[0]+ext[2]>clip[2] or y+ext[1]<clip[1] or y+ext[1]+ext[3]>clip[3]):return
        if (x+ext[0]<max(clip[0],30)-1 or x+ext[0]+ext[2]>min(clip[2],1570)+1
                or y+ext[1]<max(clip[1],15)-1 or y+ext[1]+ext[3]>min(clip[3],885)+1):
            self.overflows.append(text)
        c.move_to(x,y);c.show_text(text)
        if self.in_molecule_layer and text:
            self.molecule_text_boxes.append((x+ext[0],y+ext[1],x+ext[0]+ext[2],y+ext[1]+ext[3]))

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
        if not actor.persistent and actor.kind not in {'bind','release'} and p>=1:return
        dock=p if actor.kind=='bind'else 1-p if actor.kind=='release'else min(1,p/.13,(1-p)/.13)
        x,y=actor.position
        if actor.complex != 'homodimer':y-=12+70*(1-max(0,1 if actor.persistent else dock))
        color='ligase'if actor.enzyme_class=='ligase'else'enzyme'
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
            self.text(c,actor.label,x,y-145,22,align='center',max_width=390)
            return
        c.save();c.translate(x,y);c.rotate(-.12)
        if actor.enzyme_class in {'polymerase','ligase'}:c.scale(.5,.5)
        # A soft, asymmetric clamp rather than a protein structural claim.
        c.set_source_rgba(*rgb(color),.96)
        c.new_path();c.move_to(-42,8);c.curve_to(-60,-33,-20,-51,17,-36)
        c.curve_to(58,-43,64,18,27,32);c.curve_to(9,49,-29,38,-42,8);c.fill()
        c.set_source_rgba(*rgb('background'),.55);c.set_line_width(3)
        c.new_path();c.move_to(-12,21);c.line_to(15,-15);c.set_line_width(8);c.set_line_cap(self.cairo.LINE_CAP_ROUND);c.stroke();c.restore()
        # Keep long enzyme names on screen at either end of a strand.
        label_x=x+actor.label_offset[0]
        clip=c.clip_extents()
        if clip[0]<=x<=clip[2]:
            c.save();c.select_font_face(self.style.font);c.set_font_size(22)
            half=min(390,c.text_extents(actor.label)[2])/2+12;c.restore()
            label_x=max(clip[0]+half,min(clip[2]-half,label_x))
        self.text(c,actor.label,label_x,y+actor.label_offset[1],22,align='center',max_width=390)

    def draw(self,surface,movie,seconds,wpm=125,speed=DEFAULT_SPEED):
        self.overflows=[]
        self.in_molecule_layer=False
        self.molecule_text_boxes=[]
        c=self.cairo.Context(surface);c.scale(self.style.width/1600,self.style.height/900)
        c.set_source_rgb(*rgb('background'));c.paint()
        shot,local=movie.at(seconds,wpm,speed);states,actors,camera=presentation(movie,shot,local,wpm)
        self.text(c,movie.title,88,76,22,'muted',max_width=1400)
        self.text(c,shot.title,88,135,42,max_width=1400)
        self.text(c,shot.scope,88,180,20,'muted',max_width=1400)
        self.in_molecule_layer=True
        c.save();c.rectangle(40,215,1520,575);c.clip()
        cx,cy,zoom=camera;c.translate(800,500);c.scale(zoom,zoom);c.translate(-cx,-cy)
        for d in shot.decorations:
            c.set_source_rgba(*rgb(d.color),.10 if d.kind=='droplet'else .28)
            c.new_path();c.arc(d.x,d.y,d.radius,0,math.tau);c.fill_preserve()
            c.set_source_rgba(*rgb(d.color),.65);c.set_line_width(2);c.stroke()
            if d.label:self.text(c,d.label,d.x,d.y+d.radius+36,22,d.color,align='center')
        # Pairing is an explicit molecular relation; no distance threshold can
        # create/destroy it while molecules move or the camera reframes.
        for bound in evaluate_pairings(shot,local,wpm,states=states):
            pairing=bound.pairing
            if bound.strength<=0:continue
            a,b=states[pairing.a],states[pairing.b]
            for j in range(pairing.length):
                ia=pairing.a_start+j;ib=pairing.b_start-j
                if a.fraction(ia)<.95 or b.fraction(ib)<.95:continue
                pa=a.segment_point(ia+.5,ia);pb=b.segment_point(ib+.5,ib)
                dx,dy=pb[0]-pa[0],pb[1]-pa[1];length=max(1,math.hypot(dx,dy))
                # Letter boxes reserve eight pixels on either strand.
                inset=min(8,length*.3)
                self.line(c,[(pa[0]+dx/length*inset,pa[1]+dy/length*inset),
                             (pb[0]-dx/length*inset,pb[1]-dy/length*inset)],'muted',1.1,.6*bound.strength)
        joined_adapters={e.strand for e,p in evaluate_joins(shot,local,wpm) if p>.5}
        locations={id:(st,j)for st in states.values()for j,id in enumerate(st.strand.nucleotide_ids)if st.fraction(j)>.95}
        end_labels=[]
        for st in states.values():
            s=st.strand;n=len(s.sequence)
            # A planned output is not an extant molecule. Captions, termini
            # and regions appear only once at least one nucleotide exists.
            if not any(st.fraction(i)>0 for i in range(n)):continue
            # Draw individual nucleotide cells, with a directional box like
            # the site's strand polygons. Unknown bases remain literal N/IUPAC.
            for i,base in enumerate(s.sequence):
                fraction=st.fraction(i)
                if fraction<=0:continue
                center=st.segment_point(i+.5,i)
                left=st.segment_point(i,i);right=st.segment_point(i+1,i)
                cell=math.dist(left,right)
                if cell<.1:continue
                color=s.color
                if st.synthesized and st.synthesized[0]<=i<st.synthesized[1]:color='rna' if s.polymer=='RNA' else 'new'
                for a,b,tag in s.bands:
                    if a<=i<b:color=tag
                for region in s.regions:
                    if region.start<=i<region.end:color=region.color
                glow=0
                for a,b,intensity in st.highlights:
                    if a<=i<b:glow=max(glow,intensity)
                c.save();c.translate(*center)
                angle=math.atan2(right[1]-left[1],right[0]-left[0])
                # Keep letters readable on both antiparallel strands.
                if angle>math.pi/2:angle-=math.pi
                if angle< -math.pi/2:angle+=math.pi
                c.rotate(angle)
                c.set_source_rgba(*rgb(color),(.13+.2*glow)*fraction)
                c.rectangle(-cell/2,-10,cell,19);c.fill()
                if self.style.bases:
                    c.select_font_face('DejaVu Sans Mono');c.set_font_size(min(15,cell*1.5))
                    letter=st.converted.get(i,base)
                    ext=c.text_extents(letter);c.set_source_rgba(*rgb(color),fraction)
                    c.move_to(-ext[2]/2-ext[0],5);c.show_text(letter)
                c.restore()
            for bond in s.sealed_marks:
                a,b=(bond-1)%n,bond
                if st.molecular is None:continue
                if (s.nucleotide_ids[a],s.nucleotide_ids[b]) not in st.molecular.backbone:continue
                left=st.segment_point(a+.5,a);right=st.segment_point(b+.5,b)
                self.line(c,[(left[0],left[1]+13),(right[0],right[1]+13)],'ligase',3)
                for j in (a,b):
                    x,y=st.segment_point(j+.5,j)
                    c.set_source_rgb(*rgb('ligase'));c.set_line_width(2)
                    c.rectangle(x-7,y-11,14,23);c.stroke()
            for nick,strength in st.nicks.items():
                if strength<.5 or not 0<nick<n:continue
                xx,yy=st.segment_point(nick,nick-.1)
                c.set_source_rgb(*rgb('background'));c.rectangle(xx-1,yy-11,2,22);c.fill()
            # Arrowheads mark physical 3prime termini, including both sides
            # of a cleaved backbone. They move with the retained fragments.
            if self.style.directions and not st.closed:
                for a,b in st.visible:
                    cuts=[a,*sorted(i for i,v in st.nicks.items() if a<i<b and v>.5),b]
                    for lo,hi in zip(cuts,cuts[1:]):
                        if hi-lo<.2 or (hi==len(s.sequence) and (s.continued_3prime or s.id in joined_adapters)):continue
                        end=st.segment_point(hi,hi-.1);near=st.segment_point(hi-.2,hi-.1)
                        dx,dy=end[0]-near[0],end[1]-near[1];length=max(.01,math.hypot(dx,dy));dx/=length;dy/=length
                        internal=hi<n and hi in st.nicks
                        tip=end if internal else (end[0]+dx*8,end[1]+dy*8)
                        base=(end[0]-dx*6,end[1]-dy*6) if internal else end
                        c.set_source_rgba(*rgb(s.color),.85);c.new_path();c.move_to(*tip)
                        c.line_to(base[0]-dy*8,base[1]+dx*8);c.line_to(base[0]+dy*8,base[1]-dx*8);c.close_path();c.fill()
            for r in s.regions:
                if st.fraction(r.start)<.9:continue
                end=min(r.end,max(b for a,b in st.visible))
                if end<=r.start:continue
                pa=st.segment_point(r.start,(r.start+end)/2);pb=st.segment_point(end,(r.start+end)/2)
                mid=st.segment_point((r.start+end)/2,(r.start+end)/2)
                if r.label_offset is not None:
                    dx,dy=r.label_offset
                elif r.label_position is not None:
                    original=route_point(s,(r.start+r.end)/2)
                    dx,dy=r.label_position[0]-original[0],r.label_position[1]-original[1]
                else:dx,dy=0,r.side*46
                tx,ty=mid[0]+dx,mid[1]+dy
                self.line(c,[mid,(tx,ty-r.side*17)],r.color,1.3,.65)
                self.text(c,r.label,tx,ty,20,r.color,align='center',max_width=340)
            # Chemical marks are always shown; these distinguish a converted base from
            # a copied base even when letter-level detail is disabled.
            for raw,label in s.marks.items():
                i=int(raw)
                if st.fraction(i)<.9:continue
                x,y=st.point(i+.5);color='bio'if label=='biotin'else'ligase'
                c.set_source_rgb(*rgb(color));c.new_path();c.arc(x,y,7,0,math.tau);c.fill()
                dx,dy=s.mark_label_offsets.get(raw,(0,-25))
                if raw in s.mark_label_offsets:self.line(c,[(x,y),(x+dx,y+dy+8)],color,1,.6)
                self.text(c,label,x+dx,y+dy,17,color,align='center')
            for i,base in st.converted.items():
                x,y=st.point(i+.5);c.set_source_rgb(*rgb('new'));c.new_path();c.arc(x,y,8,0,math.tau);c.fill()
                self.text(c,base,x,y-23,21,'new',align='center')
            if st.removed:
                # Released nucleotide IDs remain in molecular state; visual
                # particles leave the viewport instead of being erased.
                event=next((e for e in shot.events if e.strand==s.id and e.kind in {'resect','resect3'}),None)
                if event is not None:
                    op=event.molecular_operation;removed=op.get('details',{}).get('removed',[])
                    count=int(st.removed)
                    for j in range(count):
                        index=j if event.kind=='resect'else n-1-j
                        if removed and s.nucleotide_ids:
                            id=removed[j];index=s.nucleotide_ids.index(id)
                        phase=(local-shot.establish)/shot.action
                        fraction=(j+1)/max(1,event.amount)
                        timing=.5-math.sin(math.asin(max(-1,min(1,1-2*fraction)))/3)
                        release=event.fraction_start+(event.fraction_end-event.fraction_start)*timing
                        age=max(0,(phase-release)*shot.action)
                        xx,yy=st.segment_point(index+.5,index)
                        xx+=((j%9)-4)*24*age;yy+=65*age+30*age*age
                        if not 40<xx<1560 or not 225<yy<785:continue
                        c.set_source_rgba(*rgb('waste'),.25);c.new_path();c.arc(xx,yy,10,0,math.tau);c.fill()
                        self.text(c,s.sequence[index],xx,yy+5,14,'muted',align='center')
            if s.label:
                if s.label_position:
                    x,y=s.label_position;x+=st.offset[0];y+=st.offset[1];align='center'
                elif st.route:
                    xs=[k[1]for k in st.route];ys=[k[2]for k in st.route]
                    x=(min(xs)+max(xs))/2+st.offset[0];y=max(ys)+62+st.offset[1];align='center'
                else:
                    x=min(s.start[0],s.end[0])+st.offset[0]-30
                    y=(s.start[1]+s.end[1])/2+st.offset[1]+8;align='right'
                self.text(c,s.label,x,y,22,s.color,align=align,max_width=380)
            if self.style.directions and not st.closed:
                last=max(b for a,b in st.visible)
                start=st.segment_point(st.removed,st.removed);end=st.segment_point(last,max(st.removed,last-.1))
                dy=-25 if s.polymer=='RNA'else 52
                if st.molecular is not None:
                    mates={a:b for a,b in st.molecular.duplex for a,b in [(a,b),(b,a)]}
                    for j,id in enumerate(s.nucleotide_ids):
                        mate=mates.get(id)
                        if mate not in locations:continue
                        other,k=locations[mate]
                        delta=other.segment_point(k+.5,k)[1]-st.segment_point(j+.5,j)[1]
                        if abs(delta)>1:
                            dy=-30 if delta>0 else 52
                            break
                direction=1 if end[0]>=start[0] else -1
                if not s.continued_5prime:end_labels.append(('5′',start[0]-direction*16,start[1]+dy,s.color,dy))
                if not s.continued_3prime and s.id not in joined_adapters:end_labels.append(('3′',end[0]+direction*16,end[1]+dy,s.color,dy))

        for actor in actors:self.enzyme(c,actor)
        # Closely opposed free ends can share a nick coordinate. Keep their
        # labels distinct without moving a molecule or changing its bonding.
        for label,x,y,color,side in end_labels:
            c.select_font_face(self.style.font);c.set_font_size(18)
            box=c.text_extents(label)
            def collides():
                left,right=x-box[2]/2,x+box[2]/2
                top,bottom=y+box[1],y+box[1]+box[3]
                return any(left<b[2]+3 and right>b[0]-3 and top<b[3]+3 and bottom>b[1]-3 for b in self.molecule_text_boxes)
            while collides():
                y+=24 if side>0 else -24
            self.text(c,label,x,y,18,color,align='center')
        c.restore()
        self.in_molecule_layer=False
        if shot.product_pool:
            frame=next(iter(states.values())).molecular
            released=sorted(frame.released)if frame else[]
            if released:
                bases=' '.join(frame.atoms[id]['base']for id in released[:56])
                if len(released)>56:bases+=' …'
                self.text(c,f'Released nucleotides · {len(released)}',88,827,18,'muted',max_width=1400)
                self.text(c,bases,88,853,17,'waste',max_width=1400)
        if self.style.coordinates:
            for j,s in enumerate(shot.strands):
                detail=f'{s.id}: {len(s.sequence)} nt'if not s.schematic else f'{s.id}: schematic'
                self.text(c,detail,88,242+21*j,15,'muted',max_width=300)
        if self.style.reagents:
            notes=list(dict.fromkeys(e.note for e in shot.events if e.note))
            for j,note in enumerate(notes[:2]):self.text(c,note,800,730+30*j,19,'muted',align='center',max_width=1360)
        if self.style.timeline:
            index=movie.shots.index(shot)
            self.text(c,f'{index+1:02d} / {len(movie.shots):02d}',1510,840,19,'muted',align='right')
        if self.style.evidence:self.text(c,shot.source,88,885 if shot.product_pool else 818,14,'muted',max_width=1400)
        if self.style.timeline:
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
