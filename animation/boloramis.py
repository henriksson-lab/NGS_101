"""Source-qualified BOLORAMIS padlock chemistry and its presentation projection."""
import math
from dataclasses import replace
from .model import Strand,Region,Shot,Event,Pairing,Junction
from .simulator import Simulator,Ligase,Polymerase
from .simulation_adapters import generated
from .protocols import load,qualify

SLUG='boloramis__10.1093+nar+gkab120'
SOURCE='BOLORAMIS2021 doi:10.1093/nar/gkab120; selected ACTN1 A1 source probe; conditional matched RNA window'
X=460;Y=480;PITCH=12


def simulate(m,*,budget=None,rna_displacement=True,phosphate='P',probe_sequence=None,rna_sequence=None,primer_sequence=None):
    s=Simulator()
    lengths=[len(m.LONG),4,len(m.ANCHOR),len(m.CODE),len(m.SHORT)]
    features=[];offset=0
    for size,tag in zip(lengths,['long','copied_code','anchor','code','short']):
        features.append((offset,offset+size,'probe.'+tag));offset+=size
    probe=s.add('probe',m.probe().top()if probe_sequence is None else probe_sequence,five=phosphate,features=features,source=SOURCE+(' synthetic perturbation'if probe_sequence is not None else''))
    rna=s.add('rna',m.ConditionalRNAMatch().sequence if rna_sequence is None else rna_sequence,sugar='RNA',five='unassigned',three='unassigned',source=SOURCE+(' synthetic target perturbation'if rna_sequence is not None else' inferred conditional window, not an isolated full transcript'))
    primer_source=SOURCE+' primer co-hybridized before ligation'+(' synthetic perturbation'if primer_sequence is not None else'')
    primer_sequence=m.RCA_PRIMER if primer_sequence is None else primer_sequence
    primer=s.add('primer',primer_sequence,features=[(0,len(primer_sequence),'rca.primer')],source=primer_source)
    initial=s.snapshot();long=s.bind_spans(probe[:lengths[0]],tuple(reversed(rna[:lengths[0]])),source=SOURCE)
    prime=s.anneal('primer','probe',binding=len(primer),source=SOURCE)
    fold_before=s.snapshot();short=s.anneal('probe','rna',binding=lengths[-1],template_site=len(rna)-1,source=SOURCE)
    nick=s.snapshot();ligase=Ligase('SplintR',('DNA',),('DNA',),('P',),('ATP',),('RNA',))
    seal=s.ligate(probe[-1],probe[0],ligase,mode='splinted',cofactors=('ATP',),source=SOURCE)
    s.supply_nt('DNA','ACGT',modifications={'T':'dT or aminoallyl-dU: occupancy unassigned'},source=SOURCE+' post-ligation source RCA mix; displayed DNA backbone notation')
    closed=s.snapshot()
    profile=Polymerase('Phi29',('DNA',),'DNA',strand_displacement=True,displacement_sugars=('DNA','RNA')if rna_displacement else('DNA',))
    if s.extension_sites(profile,['probe','primer','rna'])!=[s.tail('primer')]:raise ValueError('RCA trajectory must account for every eligible end')
    rca=s.rolling_circle('primer',profile,max_bases=len(probe)*3 if budget is None else budget,source=SOURCE+' selected successful RNA-hybrid passage, supported mechanistically by doi:10.1093/nar/gkad1101; applicability to NxGen37C inferred')
    return dict(simulator=s,initial=initial,long=long,prime=prime,fold_before=fold_before,short=short,nick=nick,seal=seal,closed=closed,rca=rca,after=s.snapshot(),probe=probe,rna=rna,primer=primer,product=s.chain('primer'),features=features)


def folded_route(n):
    # Rounded joins keep the loop tangent to both RNA-binding arms.
    points=[(X+216,Y)]
    points.extend((X+216+30*math.cos(math.pi/2-math.pi*j/40),Y-30+30*math.sin(math.pi/2-math.pi*j/40))for j in range(1,21))
    points.extend((X+66+180*math.cos(-math.pi*j/80),Y-30+100*math.sin(-math.pi*j/80))for j in range(1,81))
    points.extend((X-84+30*math.cos(math.pi-math.pi*j/40),Y-30+30*math.sin(math.pi-math.pi*j/40))for j in range(1,21))
    total=sum(math.dist(a,b)for a,b in zip(points,points[1:]));walk=0
    route=[(0,X,Y),(18,*points[0])]
    for a,b in zip(points,points[1:]):
        walk+=math.dist(a,b);route.append((18+35*walk/total,*b))
    route[-1]=(n-7,*points[-1])
    route.append((n,X,Y));return route


def probe_regions():
    return [Region(0,18,'5′ arm','adapter',1,label_offset=(160,80)),Region(18,22,'Code copy','bio',-1,label_offset=(100,-95)),Region(22,45,'RCA site','primer',-1,label_offset=(0,-55)),Region(45,53,'ACTN1 code','bio',-1,label_offset=(-15,-120)),Region(53,60,'3′ arm','adapter',1,label_offset=(-115,110))]


def tracks(r,state,*,folded=False,product=False):
    from .generated_views import strand
    probe=strand('probe',state,r['probe'],X,Y,cell=PITCH,regions=probe_regions())
    probe.sealed_marks=[0]
    if folded:probe.path=folded_route(len(r['probe']));probe.start=tuple(probe.path[0][1:]);probe.end=tuple(probe.path[-1][1:])
    rna=strand('rna',state,r['rna'],X+216,Y+24,cell=PITCH,direction=-1,regions=[Region(0,len(r['rna']),'Conditional RNA window','rna',1,label_offset=(0,75))])
    rna.continued_5prime=rna.continued_3prime=True
    ids=r['product']if product else r['primer']
    dna=strand('dna',state,ids,X+45*PITCH,Y+24,cell=PITCH,direction=-1,initial=len(r['primer']),regions=[Region(0,len(r['primer']),'RCA primer','primer',1,label_offset=(0,65))])
    dna.color='new'if product else'primer'
    return [probe,rna,dna]


def binding(a,ai,b,bi,length,id,bound=True):return Pairing(a,ai,b,bi,length,id=id,source=SOURCE,initially_bound=bound,exact=True)


def follow(r):
    site=list(r['probe']).index(r['prime'].after.paired[r['primer'][0]])
    return Event('follow','dna',template='probe',template_start=site,amount=24)


def final_view(r,event):
    """Lay out the same retained molecules; no melting or product replacement."""
    from types import SimpleNamespace
    from .generated_views import strand
    from .evaluate import route_point
    from .layout import offset_point,following_path
    state=r['after'];ids=list(r['product'])
    bound=[j for j,id in enumerate(ids)if state.paired.get(id)in set(r['probe'])];lo=min(bound)
    template=tracks(r,state,folded=True)[0];template.circular=True;template.five_prime=template.three_prime='none';template.regions=[]
    geo=lambda:SimpleNamespace(strand=template,closed=True,point=lambda q:route_point(template,q))
    start=event.template_start+event.primer_length+1
    _,normal=offset_point(geo(),start-lo,24);angle=math.pi-math.atan2(normal[1],normal[0]);co,si=math.cos(angle),math.sin(angle)
    def transform(x,y):
        dx,dy=x-(X+66),y-(Y-59)
        return 1180+co*dx-si*dy,420+si*dx+co*dy
    template.path=[(i,*transform(x,y))for i,x,y in template.path]
    template.start=tuple(template.path[0][1:]);template.end=tuple(template.path[-1][1:]);template.label='Closed DNA template';template.label_position=(1180,245)
    paired=following_path(geo(),start-lo,len(ids)-lo,24,offset=lo);anchor=tuple(paired[0][1:])
    vertices=[(100,650),(900,650)]
    vertices.extend((900+35*math.cos(math.pi/2-math.pi*j/20),615+35*math.sin(math.pi/2-math.pi*j/20))for j in range(1,21))
    vertices.append((100,580))
    vertices.extend((100+35*math.cos(math.pi/2+math.pi*j/20),545+35*math.sin(math.pi/2+math.pi*j/20))for j in range(1,21))
    vertices.extend([(900,510),anchor])
    length=sum(math.dist(a,b)for a,b in zip(vertices,vertices[1:]));path=[(0,*vertices[0])];walk=0
    for a,b in zip(vertices,vertices[1:]):
        walk+=math.dist(a,b);path.append((lo*walk/length,*b))
    path[-1]=(lo,*anchor);path.extend(paired[1:])
    product=strand('dna',state,ids,100,650,label='One continuous RCA product')
    product.path=path;product.end=tuple(path[-1][1:]);product.label_position=(490,765)
    colors={'probe.long':'adapter','probe.short':'adapter','probe.anchor':'primer','probe.code':'bio','probe.copied_code':'bio','rca.primer':'primer'}
    bands=[]
    for j,id in enumerate(ids):
        color=colors[state.atoms[id].features[0]]
        if bands and bands[-1][2]==color:bands[-1]=(bands[-1][0],j+1,color)
        else:bands.append((j,j+1,color))
    product.bands=bands
    rna=strand('rna',state,r['rna'],1400,685,cell=10,direction=-1,label='RNA window retained')
    rna.continued_5prime=rna.continued_3prime=True;rna.label_position=(1270,740)
    shot=Shot('product','One circle produces a continuous tandem product',[template,product,rna],
        [Event('bind','dna','Phi29',bond=len(ids),persistent=True,enzyme_prebound=True,enzyme_class='polymerase',enzyme_label='Phi29 polymerase',enzyme_label_offset=(-60,-70))],
        'This is the same DNA product, unfolded for inspection while its growing end remains paired to the original circle. Target-code regions repeat along one covalent strand. The displaced RNA window is retained separately. The original primer remains at the five-prime end; newly synthesized bases account for the rest.',action=9,hold=4,source=SOURCE,scope='Same molecules, rearranged for inspection · gold marks target-code sequence')
    generated(shot,state);return shot


def movie():
    m=load(SLUG,'chem.py');r=simulate(m);primer=len(r['primer']);n=len(r['probe'])
    initial=tracks(r,r['initial'])
    initial[1]=replace(initial[1],start=(initial[1].start[0]-1800,initial[1].start[1]),end=(initial[1].end[0]-1800,initial[1].end[1]))
    initial[2]=replace(initial[2],start=(initial[2].start[0]+1800,initial[2].start[1]),end=(initial[2].end[0]+1800,initial[2].end[1]))
    enter=Shot('bind','The probe and RCA primer find their binding sites',initial,
        [Event('anneal','rna',target=(1800,0),fraction_end=.4),Event('pair','probe',pairing='long',fraction_start=.4,fraction_end=.5,molecular_operation=r['long'].to_dict()),Event('anneal','dna',target=(-1800,0),fraction_start=.3,fraction_end=.65),Event('pair','dna',pairing='primer',fraction_start=.65,fraction_end=.75,molecular_operation=r['prime'].to_dict())],
        'The source supplies a phosphorylated DNA padlock and its RCA primer together during hybridization. The long arm binds this conditional matching RNA window, while the RCA primer binds the probe backbone. The short arm remains available to close the gap around the target.',
        action=8,source=SOURCE,scope='Source ACTN1 probe · inferred matching RNA window',pairings=[binding('probe',0,'rna',17,18,'long',False),binding('dna',0,'probe',44,primer,'primer',False)])
    generated(enter,r['initial'])
    f=follow(r)
    destination=folded_route(n)
    # The unbound short arm approaches above the existing hybrid instead of
    # passing straight through the long arm during a linear interpolation.
    lifted=[(i,x,y-50*min(1,max(0,(i-45)/8)))for i,x,y in destination]
    fold=Shot('fold','The second arm brings the probe ends together',tracks(r,r['fold_before']),
        [Event('fold','probe',route=lifted,fraction_end=.65),Event('fold','probe',route=destination,fraction_start=.65,fraction_end=.8),f,Event('pair','probe',pairing='short',fraction_start=.8,fraction_end=.95,molecular_operation=r['short'].to_dict())],
        'The same probe bends around the RNA. Its short arm binds beside the long arm, placing the three-prime hydroxyl next to the five-prime phosphate. Folding and hybridization do not seal that nick. The RCA primer remains attached throughout.',action=7,source=SOURCE,scope='One continuous probe · primer retained',pairings=[binding('probe',53,'rna',24,7,'short',False)])
    generated(fold,r['fold_before'])
    nicktracks=tracks(r,r['nick'],folded=True);nicktracks[0].nicks=[0]
    close=Shot('close','SplintR seals the RNA-splinted nick',nicktracks,
        [follow(r),Event('seal','probe','SplintR',bond=0,junction=Junction(0,'OH','P','splinted',SOURCE),fraction_start=.15,fraction_end=.8,enzyme_label='SplintR ligase',enzyme_label_offset=(-150,30),molecular_operation=r['seal'].to_dict())],
        'SplintR joins the adjacent DNA three-prime hydroxyl and five-prime phosphate on the RNA splint. Only this covalent bond makes the padlock a closed DNA circle. The RNA and the original RCA primer remain separate molecules; neither is converted into the circle.',action=5,source=SOURCE,scope='SplintR ligase · pink marks the sealed DNA junction · RNA retained')
    generated(close,r['nick'])
    close.camera=(X+66,Y-50,1.55)
    copytracks=tracks(r,r['after'],folded=True,product=True)
    # Initial circularity follows the actual sealed backbone, not folding.
    copytracks[0].circular=True;copytracks[0].five_prime=copytracks[0].three_prime='none'
    copytracks[0].bands=[(x.start,x.end,x.color)for x in copytracks[0].regions];copytracks[0].regions=[]
    copytracks[0].label='Closed padlock';copytracks[0].label_position=(X+66,Y-45)
    first_template=r['rca'].after.atoms[r['rca'].added[0]].template
    event=Event('rca','dna','Phi29',template='probe',template_start=list(r['probe']).index(first_template),primer_length=primer,primer_binding_length=primer,exact=True,amount=72,persistent=True,enzyme_label='',molecular_operation=r['rca'].to_dict())
    copy=Shot('rca','Phi29 copies the circle and displaces earlier strands',copytracks,[event,Event('peel','rna',template='probe',amount=PITCH,target=(1,1))],
        'Phi29 extends the retained primer around the covalently closed DNA template. Where it encounters the RNA hybrid, this selected successful trajectory displaces RNA while preserving its backbone. On returning to DNA already copied, it displaces the older daughter strand. The growing product stays continuous. Three template traversals after the primer illustrate extension, not measured yield.',action=24,hold=4,source=SOURCE,scope='Selected successful RCA · RNA passage inferred · DNA backbone notation')
    generated(copy,r['closed'])
    copy.camera=(X+66,Y-50,.85)
    copytracks[1].regions[0].label_offset=(90,75)
    # Move RNA aside only after incorporation has removed its final hybrid bond.
    collisions=[j for j,x in enumerate(r['rca'].details['incorporations'])if x['displaced']in set(r['rna'])]
    fraction=(max(collisions)+1)/len(r['rca'].added);low,high=0.,1.
    for _ in range(40):
        mid=(low+high)/2
        if mid*mid*(3-2*mid)<fraction:low=mid
        else:high=mid
    copy.events.append(Event('separate','rna',target=(-350,0),fraction_start=high+.005,fraction_end=min(.95,high+.2)))
    movie=qualify(m,SLUG,'BOLORAMIS',[enter,fold,close,copy,final_view(r,event)],{'bind':['hybridize'],'fold':['hybridize'],'close':['circle'],'rca':['rca'],'product':['rca']},limits=['RNA is the conditional25nt perfect-match window derived from source probe arms, not a verified occupied transcript or an isolated full RNA molecule.','RCA primer is co-hybridized and retained; no fresh primer or RNase is invented between ligation and RCA.','Selected successful Phi29 RNA-hybrid passage is mechanistically supported by doi:10.1093/nar/gkad1101, with applicability to the source NxGen37C reaction inferred; stalled branches and yields are not predicted.','Canonical DNA backbone lettering: each newly incorporated T-equivalent has unresolved dT versus aminoallyl-dU occupancy. No modification sites or crosslinks are invented.','Three template traversals after the23nt primer are an illustrative synthesis budget.'],extra={'molecular_simulation':r['simulator'].trace(),'presentation_version':'substrate-padlock-v1','probe_sequence':m.probe().top(),'rna_window':m.ConditionalRNAMatch().sequence,'source_primer':m.RCA_PRIMER})
    for shot in movie.shots[1:]:shot.transition=shot.establish
    return movie
