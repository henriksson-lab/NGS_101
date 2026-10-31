"""Explicit selected Tn5 products: retained adapter strands, target 3prime ends, repair.

Canonical sequences are supplied by the protocol adapter, never selected here.
The target sequence and cut positions are an illustrative molecular example.
"""
from dataclasses import replace
from .model import Strand, Region, Pairing, Shot, Event, complement

MECHANISM_SOURCE=('Adey et al.2010, doi:10.1186/gb-2010-11-12-r119 Background/Fig1; '
                  'Kia et al.2017, doi:10.1186/s12896-016-0326-1 Background/FigS3; '
                  'target 5prime transfer and9bp stagger; no simulated kinetics.')


def reaction(target, left, right, me, cut_left, cut_right, *, cell=10, x0=100, y=480,derived=None):
    """Two insertion events; productive central pair retained, flanks exit viewport."""
    rc=lambda s:complement(s[::-1])
    n=len(target);a=cut_left;b=cut_right;gap=9
    if not 0<a<a+gap<b<b+gap<n:raise ValueError('target needs ordered internal9bp-staggered cut pairs')
    if not left.endswith(me) or not right.endswith(me):raise ValueError('each transferred end must terminate in declared ME')
    x=lambda i:x0+cell*i
    top=Strand('target_top',target,(x(0),y),(x(n),y),'dna','')
    if derived is not None:top.nucleotide_ids=list(derived['target_top_ids'])
    bottom=Strand('target_bottom',rc(target),(x(n),y+24),(x(0),y+24),'dna','')
    if derived is not None:bottom.nucleotide_ids=list(derived['target_bottom_ids'])
    events=[];pairings=[];strands=[top,bottom];landed_short_routes={}
    # Explicit binding regions ensure only the retained central overlap stays
    # paired during selection; flanking binding relations are unpaired before exit.
    intervals=[(0,a,'left_flank'),(a,a+gap,'left_bridge'),(a+gap,b,'central'),(b,b+gap,'right_bridge'),(b+gap,n,'right_flank')]
    for lo,hi,id in intervals:
        pairings.append(Pairing('target_top',lo,'target_bottom',n-lo-1,hi-lo,id=id,source=MECHANISM_SOURCE))
    for position,adapter,name,upper in [(a,left,'left',True),(b,right,'right',False)]:
        # Two separate transferred oligos per homodimer; ME complements are
        # separate nontransferred strands, not a second covalent adapter insertion.
        for suffix,seq,start,end,yy,central in [
            ('top',adapter,position-len(adapter),position,y,name=='left'),
            ('bottom',adapter,position+gap+len(adapter),position+gap,y+24,name=='right')]:
            id=f'{name}_{suffix}_adapter';short_id=id+'_nontransfer'
            offset=-900 if suffix=='top' else 900
            end_index=len(seq);is_top=suffix=='top'
            bend=-145 if is_top else 145
            st=Strand(id,seq,(x(start),yy+offset+bend),(x(end),yy+offset),'adapter','',
                      regions=[Region(0,len(seq)-len(me),('s5' if name=='left'else 's7'),'primer',-1 if is_top else 1,label_offset=(0,-42 if is_top else 45)),
                               Region(len(seq)-len(me),len(seq),'ME','adapter',-1 if is_top else 1,label_offset=(0,-42 if is_top else 45))])
            short_y=yy+24 if is_top else yy-24
            short=Strand(short_id,rc(me),(x(end),short_y+offset),
                         (x(end-len(me) if is_top else end+len(me)),short_y+offset+bend*len(me)/len(seq)),'waste','')
            if derived is not None:
                st.nucleotide_ids=list(derived['reagents'][id]);short.nucleotide_ids=list(derived['reagents'][short_id])
            strands.extend([st,short])
            events.extend([Event('transfer',id,'Tn5 homodimer'if central else '',target=(0,-offset),fraction_end=.35,
                                 complex='homodimer'if central else '',persistent=central,actor_offset=(45 if is_top else -45,12 if is_top else -12)),
                           Event('transfer',short_id,target=(0,-offset),fraction_end=.35)])
            pairings.append(Pairing(id,len(seq)-len(me),short_id,len(me)-1,len(me),id=id+'_me',source=MECHANISM_SOURCE))
            target_id='target_top' if is_top else 'target_bottom'
            bond=position if is_top else n-position-gap
            events.append(Event('strand_transfer',id,other=target_id,bond=bond,fraction_start=.35,fraction_end=.55,note=MECHANISM_SOURCE))
            if central:
                # Once the flanks have left, straighten the retained end without
                # changing any base identity or inventing a filled strand.
                events.append(Event('fold',id,route=[(0,x(start),yy+offset),(len(seq),x(end),yy+offset)],fraction_start=.85))
                short_final=[(0,x(end),short_y),(len(me),x(end-len(me) if is_top else end+len(me)),short_y)]
                landed_short_routes[short_id]=short_final
                events.append(Event('fold',short_id,route=[(i,xx,sy+offset)for i,xx,sy in short_final],fraction_start=.85))
                events.append(Event('release',id,'Tn5 homodimer',complex='homodimer',actor_offset=(45 if is_top else -45,12 if is_top else -12),fraction_start=.6,fraction_end=.7))
            else:
                dx=-2300 if name=='left' else 2300
                events.extend([Event('separate',id,target=(dx,0),fraction_start=.7),Event('separate',short_id,target=(dx,0),fraction_start=.7)])
    for id,bonds,intervals in [('target_top',[a,b],[(0,a,-2300),(b,n,2300)]),
                              ('target_bottom',[n-b-gap,n-a-gap],[(0,n-b-gap,2300),(n-a-gap,n,-2300)])]:
        events.extend(Event('cleave',id,bond=bond,fraction_start=.35,fraction_end=.55)for bond in bonds)
        events.extend(Event('displace',id,span=(lo,hi),target=(dx,0),fraction_start=.7)for lo,hi,dx in intervals)
    events.extend([Event('unpair','target_top',pairing='left_bridge',fraction_start=.6,fraction_end=.7),Event('unpair','target_top',pairing='right_bridge',fraction_start=.6,fraction_end=.7)])
    tag=Shot('tagment','Tn5 transfers adapter strands at staggered cuts',strands,events,
        'Each loaded Tn5 homodimer joins two adapter three-prime ends to new target five-prime ends, nine bases apart. After transposase removal, we follow one central fragment while the flanks move out of view. The recessed target three-prime ends are still unextended.',
        action=9,hold=3,scope='Illustrative cDNA target · canonical supplied adapters · retained strand identities',source=MECHANISM_SOURCE,pairings=pairings)
    # Full products describe future copies; only the actual transferred adapter and
    # original target interval exist initially. No completed duplex is substituted.
    if derived is None:raise ValueError('Tn5 rendering requires substrate-driven reaction results')
    insert=derived['insert']
    plus=derived['repaired_top'];minus=derived['repaired_bottom']
    top_initial=len(derived['tag_top_ids']);bottom_initial=len(derived['tag_bottom_ids'])
    plus_start=x(a-len(left));plus_end=x(b+gap+len(right))
    t=Strand('selected_top',plus,(plus_start,y),(plus_end,y),'dna','',initially_visible=top_initial,
             regions=[Region(0,len(left)-len(me),'s5','primer',-1,label_offset=(0,-45)),Region(len(left)-len(me),len(left),'ME','adapter',-1,label_offset=(0,-45)),
                      Region(len(left),len(left)+len(insert),'cDNA insert','dna',-1,label_offset=(0,-45)),Region(len(left)+len(insert),len(plus),'Copied right adapter','new',-1,label_offset=(0,-45))],
             lineage=[(0,len(left),'left_top_adapter',0),(len(left),top_initial,'target_top',a)])
    u=Strand('selected_bottom',minus,(plus_end,y+24),(plus_start,y+24),'dna','',initially_visible=bottom_initial,
             regions=[Region(0,len(right),'Transferred s7 + ME','adapter',1,label_offset=(0,70)),Region(bottom_initial,len(minus),'Copy 9 bases + left adapter','new',1,label_offset=(0,70))],
             lineage=[(0,len(right),'right_bottom_adapter',0),(len(right),bottom_initial,'target_bottom',n-b-gap)])
    inventory={s.id:s.sequence for s in strands}
    t.nucleotide_ids=list(derived['repaired_top_ids']);u.nucleotide_ids=list(derived['repaired_bottom_ids'])
    retained=[]
    for id in ['left_top_adapter_nontransfer','right_bottom_adapter_nontransfer']:
        source=next(s for s in strands if s.id==id)
        route=landed_short_routes[id]
        retained.append(replace(source,start=tuple(route[0][1:]),end=tuple(route[-1][1:]),path=route))
    bindings=[Pairing(t.id,0,u.id,len(minus)-1,len(plus),id='duplex',source=MECHANISM_SOURCE)]
    # Short ME complements are lost during repair. This selected illustration
    # leaves displacement versus digestion unassigned; no ligation is invented.
    bindings.extend([Pairing(t.id,len(left)-len(me),retained[0].id,len(me)-1,len(me),id='short_left',source=MECHANISM_SOURCE),
                     Pairing(u.id,len(right)-len(me),retained[1].id,len(me)-1,len(me),id='short_right',source=MECHANISM_SOURCE)])
    repair_events=[Event('unpair',t.id,pairing='short_left',fraction_start=.32,fraction_end=.65),Event('unpair',u.id,pairing='short_right',fraction_start=.32,fraction_end=.65),
                   Event('separate',retained[0].id,target=(-1700,-700),fraction_start=.65,fraction_end=1),Event('separate',retained[1].id,target=(1700,700),fraction_start=.65,fraction_end=1),
]
    for strand,template,prefix,initial in [(t,u,'top',top_initial),(u,t,'bottom',bottom_initial)]:
        ids=derived['repaired_bottom_ids' if prefix=='top' else 'repaired_top_ids']
        gap_op=derived[prefix+'_gap'];copy_op=derived[prefix+'_adapter_copy']
        after_gap=initial+len(gap_op.added)
        for op,start,end,p,limit in [(gap_op,.05,.32,initial,after_gap),(copy_op,.65,1,after_gap,len(strand.sequence))]:
            first_template=op.after.atoms[op.added[0]].template
            repair_events.append(Event('extend',strand.id,'DNA polymerase',template=template.id,template_start=ids.index(first_template),primer_length=p,primer_binding_length=1,copy_end=limit,exact=True,fraction_start=start,fraction_end=end,molecular_operation=op.to_dict(),enzyme_label_offset=(0,-120)))
    # Every eventual copied base has a declared complement; visibility still gates
    # bonds until each nucleotide has actually been synthesized.
    repair=Shot('fill','Extend the original recessed 3′ ends', [t,u,*retained],repair_events,
        'Each original target three-prime end extends across the nine-base stagger, then stops at the bound mosaic-end complement. After that short oligo loses binding, extension continues through the attached adapter. This selected path leaves the mechanism of oligo release unassigned. The original target and transferred adapter strands are retained.',
        action=8,hold=3,scope='9 nt stagger · original 3′ ends extend',source=MECHANISM_SOURCE,pairings=bindings,inventory=inventory)
    return tag,repair,replace(t,initially_visible=None),replace(u,initially_visible=None),insert
