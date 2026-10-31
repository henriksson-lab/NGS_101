"""Project simulator results into existing mechanism layouts.

These functions never repair a molecule to match a drawing: an output mismatch
raises, so the source recipe or the old animation must be investigated.
"""
from .model import Event,Pairing,Strand,Shot,Region
from .other_simulations import mip as run_mip,smallseq as run_smallseq


def sequence(state,ids):return ''.join(state.atoms[i].base for i in ids)


def track(strand,state,ids):
    ids=list(ids);value=sequence(state,ids)
    if strand.sequence!=value:raise ValueError(f'{strand.id}: source simulation disagrees with legacy drawing')
    strand.nucleotide_ids=ids


def generated(shot,state):
    from dataclasses import asdict
    shot.initial_atom_records=[asdict(a)for a in state.atoms.values()]
    shot.initial_backbone=list(state.next.items());shot.initial_released=list(state.released)
    shot.chemistry_generated=True;shot.initial_nucleotides=list(state.atoms)
    shot.chemical_pairs=[(a,b)for a,b in state.paired.items()if a<b]


def evidence(movie,sim):
    movie.provenance['molecular_simulation']=sim.trace()
    movie.provenance['chemistry_authority']='Substrate-driven nucleotide/backbone/duplex reactions; drawing consumes generated products.'
    return movie.validate()


def mip(movie,m):
    r=run_mip(m);s=r['simulator'];fold,gap,close=movie.shots
    l=len(m.LIG_ARM);e=len(m.EXT_ARM);n=len(r['probe']);full=len(r['filled']);gap_nt=full-n
    x0=110;cell=8.8;yy=480;right=x0+len(r['target'])*cell;ext=x0+e*cell;lig=ext+gap_nt*cell
    route=[(0,lig,yy),(l,right,yy),(l+16,right,300),(n-e-16,x0,300),(n-e,x0,yy),(n,ext,yy),(full,lig,yy)]
    fold.strands[0].path=[(0,x0,360),(n,right,360)];fold.events[0].route=route[:-1]
    for shot in [fold,gap,close]:
        shot.strands[0].label_position=(780,255)
        shot.strands[1].start=(right,yy+24);shot.strands[1].end=(x0,yy+24)
    gap.strands[0].path=route;close.strands[0].path=route
    for shot,probe_ids,state in [(fold,r['probe'],r['initial']),(gap,r['filled'],r['before']),(close,r['filled'],r['filled_state'])]:
        track(shot.strands[0],s.snapshot(),probe_ids);track(shot.strands[1],s.snapshot(),r['target']);generated(shot,state)
    for p,key in zip(fold.pairings,['mip_lig','mip_ext']):p.id=key;p.initially_bound=False
    fold.events.extend([Event('pair','probe',pairing='mip_lig',fraction_start=.8,molecular_operation=r['bind5'].to_dict()),Event('pair','probe',pairing='mip_ext',fraction_start=.8,molecular_operation=r['bind3'].to_dict())])
    gap.events[0].molecular_operation=r['extend'].to_dict();close.events[0].molecular_operation=r['seal'].to_dict()
    from dataclasses import replace
    circle=replace(close.strands[0],circular=True,five_prime='none',three_prime='none',nicks=[])
    target=replace(close.strands[1],label='Linear target fragment')
    free=Strand('unused_probe',sequence(s.snapshot(),r['unused']),(-1680,710),(-320,710),'primer','',five_prime='P',nucleotide_ids=list(r['unused']),
                regions=[Region(0,2,'5′ probe dinucleotide · retained','primer',1,label_offset=(190,50))])
    exo=Shot('exonucleases','Digest from accessible 3′ ends; retain the closed circle',[circle,target,free],
        [Event('transfer','unused_probe',target=(1800,0),fraction_end=.2),Event('resect3','target','Exonuclease III',amount=len(r['target_digest'].details['removed']),molecular_operation=r['target_digest'].to_dict()),
         Event('resect3','unused_probe','Exonuclease I',amount=len(r['unused_digest'].details['removed']),fraction_start=.2,molecular_operation=r['unused_digest'].to_dict())],
        'Exonuclease I removes single-stranded DNA from its free three-prime end, retaining the source-specified five-prime dinucleotide. Exonuclease III acts from an accessible duplex three-prime end. The closed probe has no entry end for these exonuclease activities. Released nucleotides remain products, rather than being erased.',
        scope='Selected accessible molecules · no quantitative clearance or universal circle immunity',source=movie.provenance['authority'],action=8,hold=3)
    exo.inventory.update(close.inventory)
    generated(exo,r['exo_before']);movie.shots.append(exo);movie.provenance['steps_used']['exonucleases']=['digest_linear']
    return evidence(movie,s)



def smallseq(movie,m):
    r=run_smallseq(m);s=r['simulator'];mask,lig3,digest,lig5=movie.shots
    a=lig3.strands[0];track(a,r['lig3'].after,r['first']);generated(lig3,r['lig3_before'])
    next(e for e in lig3.events if e.kind=='seal').molecular_operation=r['lig3'].to_dict()
    track(digest.strands[0],r['digest_before'],r['free']);track(digest.strands[1],r['digest_before'],r['rtp']);generated(digest,r['digest_before'])
    digest.events[0].amount=len(r['digest'].details['removed']);digest.events[0].molecular_operation=r['digest'].to_dict()
    track(lig5.strands[0],s.snapshot(),r['product']);generated(lig5,r['lig5_before'])
    next(e for e in lig5.events if e.kind=='seal').molecular_operation=r['lig5'].to_dict()
    # Masking is an inferred binding view, not an enzymatic product. It remains
    # explicitly conditional until its own binding state projection is migrated.
    return evidence(movie,s)
