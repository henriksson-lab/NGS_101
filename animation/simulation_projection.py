"""Bind visual tracks/events to generated nucleotide and bond changes."""

def extensions(shot,before,after,operations,*,fraction_start=.55):
    """Project every generated extension, irrespective of input/primer role."""
    from .model import Event
    from dataclasses import replace
    shot.strands=[replace(s)for s in shot.strands]
    for strand in shot.strands:
        old_length=len(strand.sequence);ids=[];seen=set();id=strand.nucleotide_ids[0]
        while id is not None:
            if id in seen:raise ValueError('linear extension projection cannot flatten a covalent circle')
            seen.add(id)
            ids.append(id);id=after.next.get(id)
        step=((strand.end[0]-strand.start[0])/old_length,(strand.end[1]-strand.start[1])/old_length)
        strand.nucleotide_ids=ids;strand.sequence=''.join(after.atoms[id].base for id in ids)
        sugars=[]
        for j,id in enumerate(ids):
            kind=after.atoms[id].sugar
            if sugars and sugars[-1][2]==kind:sugars[-1]=(sugars[-1][0],j+1,kind)
            else:sugars.append((j,j+1,kind))
        strand.polymer=sugars[0][2]if len(sugars)==1 and sugars[0][2]in {'DNA','RNA'}else'mixed'
        strand.sugars=sugars if strand.polymer=='mixed'else[]
        strand.end=(strand.start[0]+step[0]*len(ids),strand.start[1]+step[1]*len(ids))
        strand.initially_visible=sum(id in before.atoms for id in ids)
    locations={id:(strand,j)for strand in shot.strands for j,id in enumerate(strand.nucleotide_ids)}
    shot.events=[e for e in shot.events if e.kind!='extend']
    for op in operations:
        strand,j=locations[op.details['starting_3prime']]
        template,k=locations[op.after.atoms[op.added[0]].template]
        exact=not(set(strand.sequence+template.sequence)-set('ACGTU'))
        shot.events.append(Event('extend',strand.id,op.details['enzyme']['name'],template=template.id,template_start=k,primer_length=j+1,primer_binding_length=1,copy_end=j+1+len(op.added),exact=exact,fraction_start=fraction_start,molecular_operation=op.to_dict(),persistent=True))


def binding_spans(shot):
    """Derive validation spans from the actual binding/melting deltas."""
    from dataclasses import replace
    locations={id:(s.id,j)for s in shot.strands for j,id in enumerate(s.nucleotide_ids)}
    for j,p in enumerate(shot.pairings):
        e=next((e for e in shot.events if e.pairing==p.key and e.kind in {'pair','unpair'}),None)
        if e is None:continue
        pairs=e.molecular_operation['pairing_added'if e.kind=='pair'else'pairing_removed']
        columns=[]
        for a,b in pairs:
            if a not in locations or b not in locations:continue
            u,v=locations[a],locations[b]
            if u[0]==p.b and v[0]==p.a:u,v=v,u
            if u[0]==p.a and v[0]==p.b:columns.append((u[1],v[1]))
        columns.sort()
        if not columns or any(b!=(a[0]+1,a[1]-1)for a,b in zip(columns,columns[1:])):
            raise ValueError('binding projection requires one contiguous antiparallel span')
        shot.pairings[j]=replace(p,a_start=columns[0][0],b_start=columns[0][1],length=len(columns),exact=False)

def configure(shots,r):
    snapshots={'anchor':'anchor_before','rt':'rt_before','terminal':'terminal_before','switch':'switch_before','firstcopy':'firstcopy_before','pcr':'pcr_before','tagment':'tag_before','fill':'repair_before','index':'index_before','index_copy':'index_copy_before','index_duplex':'index_duplex_before','library':'final_state'}
    bindings={'anchor':{'rt_anchor':'anchor_bind'},'switch':{'tso_site':'switch_bind'},'firstcopy':{'first_primer':'firstcopy_bind'},'pcr':{'minus_site':'preamp_r2_bind','plus_site':'preamp_f2_bind'},'index':{'index_f':'index_f1_bind','index_r':'index_r1_bind'},'index_copy':{'next_f':'index_f2_bind','next_r':'index_r2_bind'},'index_duplex':{'third_f':'index_f3_anneal','third_r':'index_r3_anneal'}}
    melts={'firstcopy':{'hybrid':'rna_melt','tso_site':'tso_melt'},'pcr':{'parents':'preamp_melt'},'tagment':{'left_bridge':'left_bridge_melt','right_bridge':'right_bridge_melt'},'fill':{'short_left':'bottom_me_release','short_right':'top_me_release'},'index':{'index_parents':'index_melt'},'index_copy':{'old_f':'index_f1_melt','old_r':'index_r1_melt'},'index_duplex':{'old_f2':'index_f2_melt','old_r2':'index_r2_melt'}}
    for shot in shots:
        state=r[snapshots[shot.id]]
        shot.chemical_pairs=[(a,b)for a,b in state.paired.items()if a<b]
        from dataclasses import asdict
        shot.initial_atom_records=[asdict(a)for a in state.atoms.values()]
        shot.initial_backbone=list(state.next.items());shot.initial_released=list(state.released)
        shot.initial_nucleotides=list(state.atoms);shot.chemistry_generated=True
        tracks={s.id:s for s in shot.strands}
        for e in shot.events:
            if e.kind=='pair':e.molecular_operation=r[bindings[shot.id][e.pairing]].to_dict()
            elif e.kind=='unpair':e.molecular_operation=r[melts[shot.id][e.pairing]].to_dict()
            elif e.kind=='cleave':
                ids=tracks[e.strand].nucleotide_ids;bond=(ids[e.bond-1],ids[e.bond])
                op=next(op for op in [r['transfer_left'],r['transfer_right']]if bond in set(op.before.next.items())-set(op.after.next.items()))
                e.molecular_operation=op.to_dict()
            elif e.kind=='strand_transfer':e.molecular_operation=r['transfer_left'if e.strand.startswith('left')else'transfer_right'].to_dict()
            elif e.kind=='release'and shot.id=='tagment':e.molecular_operation=r['release_left'if e.strand.startswith('left')else'release_right'].to_dict()
