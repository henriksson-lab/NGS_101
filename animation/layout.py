"""Nucleotide-indexed geometry; these functions never alter molecular bonds."""
import math


def offset_point(template,boundary,distance):
    n=len(template.strand.sequence)
    def point(q):return template.point(q%n if template.closed else max(0,min(n,q)))
    x,y=point(boundary);a=point(boundary-.01);b=point(boundary+.01)
    dx,dy=b[0]-a[0],b[1]-a[1];length=max(.0001,math.hypot(dx,dy))
    normal=(-dy/length,dx/length)
    return (x+distance*normal[0],y+distance*normal[1]),normal


def following_path(template,start,length,distance,*,offset=0):
    return [(offset+j/2,*offset_point(template,start-j/2,distance)[0])for j in range(2*length+1)]


def generated_rca_path(output,template,event):
    """Pair only actual bound nucleotides; displace only the unpaired prefix."""
    visible=int(output.visible[0][1]);n=len(output.strand.sequence)
    ids=set(template.strand.nucleotide_ids)
    mates={a:b for a,b in output.molecular.duplex for a,b in [(a,b),(b,a)]}
    bound=[j for j,id in enumerate(output.strand.nucleotide_ids[:visible])if mates.get(id)in ids]
    if not bound:raise ValueError('RCA layout requires the actual paired growing end')
    lo=min(bound)
    if bound!=list(range(lo,visible)):raise ValueError('RCA layout needs a contiguous paired suffix')
    start=event.template_start+event.primer_length+1
    route=[]
    if lo:
        # Each base keeps its template phase. Older displaced material spreads
        # out gradually instead of rotating a whole tail with the growing end.
        spacing=event.amount or 72
        route=[(j/2,*offset_point(template,start-j/2,24+(lo-j/2)*spacing/len(template.strand.sequence))[0])for j in range(2*lo)]
    route+=following_path(template,start-lo,visible-lo,24,offset=lo)
    if visible<n:route.append((n,*route[-1][1:]))
    return route
