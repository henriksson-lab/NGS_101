"""Base-resolved SMART-seq2 teaching example, assembled from protocol-owned oligos.

A short synthetic transcript and source-era example indexes are explicitly not
sample data. The 2014 long TSO variant retains its source qualification. Reaction
products come from nucleotide connectivity, with the endpoint declaration checked
independently after simulation.
"""
from dataclasses import replace
from .model import Strand, Region, Pairing, Shot, Event, complement

CELL=10
DY=24
X0=100
Y=456
rc=lambda s:complement(s[::-1])


def shifted(s,dx=0,dy=0):
    return replace(s,start=(s.start[0]+dx,s.start[1]+dy),end=(s.end[0]+dx,s.end[1]+dy),
                   path=[(i,x+dx,y+dy)for i,x,y in s.path])


def bp(a,a_start,b,b_start,length,id,source,bound=True,exact=True):
    return Pairing(a,a_start,b,b_start,length,id=id,source=source,initially_bound=bound,exact=exact)


def bind(id,strand,start,end):
    return Event('pair',strand,pairing=id,fraction_start=start,fraction_end=end)


def movie():
    from .protocols import load,qualify
    from .tagmentation import reaction
    slug='smart-seq2__10.1038+nmeth.2639';m=load(slug,'smartseq2.py');source='Protocol SMART-seq2 oligos and2014 TSO variant; '+m.DOI
    h=m.ISPCR;H=len(h);dt=m.OLIGO_DT[:-2]+'GA';P=len(dt)
    body=('ACGUGCUAGCUACGACUGGUACGUACCGUAUGCUGACGUAGCUACGUA')[:48]
    from .smartseq2_simulation import simulate
    derived=simulate(m,body)
    rseq=''.join(derived['simulator'].atoms[i].base for i in derived['rna_ids'])
    tso_seq=''.join(s.top for s in m.drawn_tso());T=len(tso_seq)
    first=derived['first'];F=len(first);tail=derived['tail'];U=len(tail);full=derived['full'];N=len(full)
    x=lambda i:X0+CELL*i
    rna=Strand('rna',rseq,(x(F),Y-DY),(x(len(dt)-32),Y-DY),'rna','',polymer='RNA',
               nucleotide_ids=list(derived['rna_ids']),regions=[Region(0,len(body),'RNA','rna',-1,label_offset=(0,-45)),Region(len(body)+2,len(rseq),'poly(A)','bio',-1,label_offset=(0,-45))])
    regions=[Region(0,H,'ISPCR handle','primer',1,label_offset=(0,48)),Region(H,H+len(m.DT_LINKER),'AC','adapter',1,label_offset=(0,82)),Region(H+len(m.DT_LINKER),P,'oligo(dT)30VN','bio',1,label_offset=(0,48))]
    first_dna=Strand('cdna',first,(x(0),Y),(x(F),Y),'dna','',initially_visible=P,regions=regions,nucleotide_ids=list(derived['first_ids']))
    anchor=Shot('anchor','The RT primer binds poly(A)',[rna,shifted(first_dna,dx=-1800)],
        [Event('anneal','cdna',target=(1800,0),fraction_end=.55),bind('rt_anchor','cdna',.55,.7),Event('highlight','cdna',span=(0,P),fraction_start=.55)],
        'The reverse transcription primer enters and binds through oligo dT. Its amplification handle and AC linker remain unpaired. RNA and DNA are antiparallel.',
        action=5,pairings=[bp('cdna',P-32,'rna',len(rseq)-1,32,'rt_anchor',source,False)],source=source,scope='Published primer · synthetic transcript · GA is one VN member')
    rt=Shot('rt','Reverse transcriptase extends the bound primer',[rna,first_dna],
        [Event('extend','cdna','SuperScript II',template='rna',template_start=len(body)-1,primer_length=P,primer_binding_length=32,exact=True,persistent=True)],
        'SuperScript II extends the existing primer. The RNA remains paired with its growing complementary DNA, and the original five-prime handle stays attached.',
        action=6,pairings=[bp('cdna',P-32,'rna',len(rseq)-1,len(rseq),'hybrid',source)],source=source,scope='First-strand synthesis · RNA retained')
    tail_dna=replace(first_dna,sequence=tail,end=(x(U),Y),initially_visible=F,nucleotide_ids=list(derived['tail_ids']),regions=regions+[Region(F,U,'CCC','new',1,label_offset=(0,48))])
    terminal=Shot('terminal','RT adds the terminal cytosines',[rna,tail_dna],
        [Event('tail','cdna','SuperScript II',persistent=True,enzyme_prebound=True)],
        'At the RNA end, reverse transcriptase adds the cytosines used in this template-switch illustration. The RNA remains present; it has not become the other DNA strand.',
        action=3,pairings=rt.pairings,source=source,scope='Model template-switch boundary · RNA retained')
    tso=Strand('tso',tso_seq,(x(N),Y-DY),(x(F),Y-DY),'adapter','',polymer='mixed',
               nucleotide_ids=list(derived['tso_ids']),sugars=[(0,T-3,'DNA'),(T-3,T-1,'RNA'),(T-1,T,'LNA')],
               regions=[Region(0,H,'TSO / ISPCR handle','adapter',-1,label_offset=(0,-80)),Region(H,T-3,'ACAT · 2014 form','primer',-1,label_offset=(0,-45)),Region(T-3,T,'rG rG +G','bio',-1,label_offset=(0,-115))])
    switched=replace(first_dna,sequence=full,end=(x(N),Y),initially_visible=U,nucleotide_ids=list(derived['full_ids']),regions=regions+[Region(U,N,'Copied TSO','adapter',1,label_offset=(0,48))])
    switch=Shot('switch','RT switches onto the template-switch oligo',[rna,shifted(tso,dx=1800),switched],
        [Event('anneal','tso',target=(-1800,0),fraction_end=.3),bind('tso_site','cdna',.3,.4),
         Event('extend','cdna','SuperScript II',template='tso',template_start=T-4,primer_length=U,primer_binding_length=3,exact=True,persistent=True,enzyme_prebound=True,fraction_start=.4,enzyme_label_offset=(0,125))],
        'The terminal guanosines of the template-switch oligo pair with the cytosines. The same reverse transcriptase copies its handle while RNA remains in the hybrid. This is the longer TSO form printed in the 2014 protocol.',
        action=7,pairings=[*rt.pairings,bp('cdna',F,'tso',T-1,3,'tso_site',source,False)],source=source,scope='2014 TSO form · two riboG and one LNA G · RNA retained')
    parent=replace(switched,initially_visible=None)
    minus=derived['first_reverse']
    reverse=Strand('reverse',minus,(x(N),Y+DY),(x(0),Y+DY),'new','',initially_visible=H,nucleotide_ids=list(derived['first_reverse_ids']),
                   regions=[Region(0,H,'ISPCR primer','primer',1,label_offset=(0,45))])
    firstcopy=Shot('firstcopy','The PCR primer makes the second DNA strand',[rna,tso,parent,shifted(reverse,dx=1800)],
        [Event('unpair','rna',pairing='hybrid',fraction_end=.15),Event('unpair','tso',pairing='tso_site',fraction_end=.15),
         Event('separate','rna',target=(-2000,-300),fraction_start=.15,fraction_end=.45),Event('separate','tso',target=(2000,-300),fraction_start=.15,fraction_end=.45),
         Event('anneal','reverse',target=(-1800,0),fraction_start=.15,fraction_end=.45),bind('first_primer','reverse',.45,.55),
         Event('extend','reverse','KAPA HiFi',template='cdna',template_start=N-H-1,primer_length=H,exact=True,fraction_start=.55)],
        'Heat separates RNA and the template-switch oligo from cDNA. An ISPCR primer enters from the edge, binds the copied handle and primes the first complementary DNA strand. Polymerase extends it; no DNA strand appears without synthesis.',
        action=8,pairings=[*rt.pairings,bp('cdna',F,'tso',T-1,3,'tso_site',source),bp('reverse',0,'cdna',N-1,N,'first_primer',source,False)],source=source,scope='First complementary DNA · exact primer and copy mapping')
    # PCR keeps old DNA parents while daughter strands start only as primers.
    oldminus=replace(reverse,initially_visible=None,regions=[])
    oldplus=replace(parent,regions=[Region(0,H,'ISPCR site','primer',-1,label_offset=(0,-42)),Region(N-H,N,'ISPCR site','adapter',-1,label_offset=(0,-42))])
    daughter_minus=shifted(replace(reverse,id='daughter_minus',nucleotide_ids=list(derived['preamp_r2_ids'])),dx=1800,dy=-100)
    daughter_plus=Strand('daughter_plus',full,(x(0)-1800,Y+148),(x(N)-1800,Y+148),'new','',initially_visible=H,nucleotide_ids=list(derived['preamp_f2_ids']),
                         regions=[Region(0,H,'ISPCR primer','primer',1,label_offset=(0,50))])
    pcr=Shot('pcr','Both primer sites drive amplification',[oldplus,oldminus,daughter_minus,daughter_plus],
        [Event('unpair','cdna',pairing='parents',fraction_end=.15),Event('separate','cdna',target=(0,-100),fraction_start=.15,fraction_end=.35),Event('separate','reverse',target=(0,100),fraction_start=.15,fraction_end=.35),
         Event('anneal','daughter_minus',target=(-1800,0),fraction_start=.2,fraction_end=.5),Event('anneal','daughter_plus',target=(1800,0),fraction_start=.2,fraction_end=.5),
         bind('minus_site','daughter_minus',.5,.58),bind('plus_site','daughter_plus',.5,.58)],
        'The duplex melts. Primers enter from opposite edges and bind the two amplification sites. Each polymerase extends a primer while retaining its parent template. The paired strands remain close together.',
        action=8,pairings=[bp('cdna',0,'reverse',N-1,N,'parents',source),bp('daughter_minus',0,'cdna',N-1,N,'minus_site',source,False),bp('daughter_plus',0,'reverse',N-1,N,'plus_site',source,False)],source=source,scope='Later PCR cycle · retained templates and new daughters')
    tag,repair,tagged_top,tagged_bottom,insert=reaction(full,m.nx.ADAPTOR_S5,m.nx.ADAPTOR_S7,m.nx.ME,derived['cut_left'],derived['cut_right'],cell=CELL,x0=X0,y=480,derived=derived)
    shots=[anchor,rt,terminal,switch,firstcopy,pcr,tag,repair]
    index_shots,final_components=index_pcr(m,tagged_top,tagged_bottom,insert,source,derived)
    shots.extend(index_shots)
    for id,key in [('rt','rt'),('terminal','terminal'),('switch','switch'),('firstcopy','firstcopy')]:
        shot=next(s for s in shots if s.id==id)
        event=next(e for e in shot.events if e.kind in {'extend','tail'})
        event.molecular_operation=derived[key].to_dict()
    from .simulation_projection import configure,extensions
    extensions(pcr,derived['pcr_before'],derived['pcr_after'],derived['preamp_operations'],fraction_start=.58)
    configure(shots,derived)
    from .simulation_projection import binding_spans
    for shot in shots:
        if shot.id in {'index','index_copy','index_duplex'}:binding_spans(shot)
    # Each PCR group has a dedicated label lane outside its duplex and primer
    # annotations. Enzymes stay at their growing ends through the result hold.
    for shot in shots:
        if shot.id not in {'index','index_copy','index_duplex'}:continue
        tracks={s.id:s for s in shot.strands}
        for e in shot.events:
            if e.kind=='extend':
                e.persistent=True
                y=tracks[e.strand].start[1]+sum(x.target[1]for x in shot.events if x.strand==e.strand and x.kind=='separate')
                e.enzyme_label_offset=(0,-140 if y<500 else 145)
    steps={s.id:[{'anchor':'reverse_transcription','rt':'reverse_transcription','terminal':'reverse_transcription','switch':'template_switch','firstcopy':'preamp','pcr':'preamp','tagment':'tagment','fill':'gap_fill','index':'index_pcr','index_copy':'index_pcr','index_duplex':'index_pcr','library':'index_pcr'}[s.id]]for s in shots}
    return qualify(m,slug,'SMART-seq2',shots,steps,limits=[
        'Synthetic short transcript and cut coordinates illustrate a molecular example, not a measured insert or claimed yield.',
        'The ACAT TSO follows2014 Nature Protocols; added AT is not the defining2013 TSO sequence.',
        'Source-era May2012 vendor letter and15029964RevC support the adapter orders and8base widths. N501/N701 is an illustrative available kit pair, not an assigned study sample.',
        'Short nontransferred ME fate is illustrated as loss during repair; displacement versus digestion not assigned.'],
        extra={'source_dt_pattern':m.OLIGO_DT,'source_ispcr':h,'source_tso_order':m.TSO,'template':rseq,
               'index_example':{'i5_name':'N501','i5_oligo':derived['i5_oligo'],'i7_name':'N701','i7_oligo':derived['i7_oligo'],'source':'Illumina customer sequence letter May24,2012 p3; sourcekit FC1211011','notice':'Oligonucleotide sequences ©2007–2012 Illumina, Inc. All rights reserved.'},'full_length_cdna':full,'final_components':final_components,'selected_insert':insert,'presentation_version':'all-extendable-ends-v2','molecular_simulation':derived['trace'],'chemistry_authority':'Reactions generate nucleotide identities and backbone/duplex connectivity; no expected product supplied.'})


def index_pcr(m,top,bottom,insert,source,derived):
    """Both 5prime primer tails are retained; a later copy makes their complements."""
    nx=m.nx
    f=derived['index_f'];r=derived['index_r']
    # One right-to-left copy of the opposite end produces a complete template
    # for the subsequent cycle; endpoints are derived, not hand-typed sequences.
    core=top.sequence
    plus=derived['final_plus']
    minus=derived['final_minus']
    left_overhang=len(f)-len(nx.S5);right_overhang=len(r)-len(nx.S7)
    cell=8.3;leftx=110
    total=len(plus);endx=leftx+cell*total
    tx=lambda i:leftx+cell*i
    # Place parent templates exactly under the eventual binding columns. Primer
    # tails stay collinear; unpaired status, not bending, identifies the overhang.
    parent_top=replace(top,id='tagged_top',start=(tx(left_overhang),400),end=(tx(left_overhang+len(core)),400),regions=[],lineage=[])
    parent_bottom=replace(bottom,id='tagged_bottom',start=(tx(left_overhang+len(core)),424),end=(tx(left_overhang),424),regions=[],lineage=[])
    outplus=Strand('index_plus',derived['index_f1'],(leftx,564),(tx(len(derived['index_f1'])),564),'new','',initially_visible=len(f),schematic=False,nucleotide_ids=list(derived['index_f1_ids']),
                   regions=[Region(0,len(m.P5),'P5','primer',1,label_offset=(0,52)),Region(len(m.P5),left_overhang,'i5 · N501 example','bio',1,label_offset=(0,85)),Region(left_overhang,len(f),'s5 · binding','adapter',1,label_offset=(0,52))])
    outminus=Strand('index_minus',derived['index_r1'],(endx,424),(tx(left_overhang),424),'new','',initially_visible=len(r),schematic=False,nucleotide_ids=list(derived['index_r1_ids']),
                    regions=[Region(0,len(m.P7),'P7','primer',-1,label_offset=(0,-70)),Region(len(m.P7),right_overhang,'i7 · N701 example','bio',-1,label_offset=(0,-105)),Region(right_overhang,len(r),'s7 · binding','adapter',-1,label_offset=(0,-70))])
    # Both participants can extend: the incoming oligo and the recessed original
    # strand. The reaction engine discovers these ends without primer-role selection.
    plus_stop=len(derived['index_f1']);minus_stop=len(derived['index_r1'])
    first=Shot('index','Index primers bind through s5 and s7',[parent_top,parent_bottom,shifted(outplus,dx=-1800),shifted(outminus,dx=1800)],
        [Event('unpair','tagged_top',pairing='index_parents',fraction_end=.15),Event('separate','tagged_bottom',target=(0,116),fraction_start=.15,fraction_end=.35),Event('anneal','index_plus',target=(1800,0),fraction_start=.2,fraction_end=.45),Event('anneal','index_minus',target=(-1800,0),fraction_start=.2,fraction_end=.45),bind('index_f','index_plus',.45,.55),bind('index_r','index_minus',.45,.55)],
        'The index primers bind through s5 and s7. Each duplex has two extendable three-prime ends. One polymerase extends the added primer along the original strand. Another extends the original strand along the primer five-prime overhang, copying its index and flow-cell sequence. Four polymerases act across the two duplexes.',
        action=8,pairings=[bp('tagged_top',0,'tagged_bottom',len(core)-1,len(core),'index_parents',source),bp('index_plus',left_overhang,'tagged_bottom',len(core)-1,len(nx.S5),'index_f',source,False,False),bp('index_minus',right_overhang,'tagged_top',len(core)-1,len(nx.S7),'index_r',source,False,False)],source=source,scope='Source-era Nextera components · example N501/N701 pair')
    from .simulation_projection import extensions
    extensions(first,derived['index_before'],derived['index_after'],derived['index_operations'])
    # These two single strands now each contain one added primer tail. In the
    # next cycle the opposite primer supplies the missing tail and copies the
    # first tail already present in its template.
    partial_plus=replace(outplus,id='first_plus',sequence=derived['index_f1'],end=(tx(plus_stop),564),start=(leftx,564),initially_visible=None,regions=[])
    partial_minus=replace(outminus,id='first_minus',sequence=derived['index_r1'],start=(endx,424),end=(tx(left_overhang),424),initially_visible=None,regions=[])
    next_plus=replace(outplus,sequence=derived['index_f2'],nucleotide_ids=list(derived['index_f2_ids']),start=(leftx,448),end=(endx,448))
    next_minus=replace(outminus,sequence=derived['index_r2'],nucleotide_ids=list(derived['index_r2_ids']),start=(endx,588),end=(leftx,588))
    second=Shot('index_copy','The next cycle copies both added ends',[replace(parent_top,id='older_top'),shifted(replace(parent_bottom,id='older_bottom'),dy=116),partial_plus,partial_minus,shifted(next_plus,dx=-1800),shifted(next_minus,dx=1800)],
        [Event('unpair','first_plus',pairing='old_f',fraction_end=.15),Event('unpair','first_minus',pairing='old_r',fraction_end=.15),Event('separate','older_top',target=(-2200,0),fraction_start=.15,fraction_end=.35),Event('separate','older_bottom',target=(2200,0),fraction_start=.15,fraction_end=.35),Event('anneal','index_plus',target=(1800,0),fraction_start=.2,fraction_end=.45),Event('anneal','index_minus',target=(-1800,0),fraction_start=.2,fraction_end=.45),bind('next_f','index_plus',.45,.55),bind('next_r','index_minus',.45,.55)],
        'After melting, the opposite index primers bind the selected first-cycle products. Both ends of each duplex extend again. The new primers copy the template ends, while the original strands copy the remaining primer overhangs. These four extensions produce complete strands with both indexed ends.',
        action=8,pairings=[bp('first_plus',left_overhang,'older_bottom',len(core)-1,len(core),'old_f',source,exact=False),bp('first_minus',right_overhang,'older_top',len(core)-1,len(core),'old_r',source,exact=False),bp('index_plus',left_overhang,'first_minus',len(partial_minus.sequence)-1,total-left_overhang,'next_f',source,False,False),bp('index_minus',right_overhang,'first_plus',len(partial_plus.sequence)-1,total-right_overhang,'next_r',source,False,False)],source=source,scope='Subsequent index-PCR cycle · every copied tail has a template')
    extensions(second,derived['index_copy_before'],derived['index_copy_after'],derived['index_copy_operations'])
    third_plus=replace(next_plus,id='third_parent_plus',initially_visible=None,start=(leftx,448),end=(endx,448),regions=[])
    third_minus=replace(next_minus,id='third_parent_minus',initially_visible=None,start=(endx,588),end=(leftx,588),regions=[])
    dplus=replace(outplus,sequence=derived['index_f3'],nucleotide_ids=list(derived['index_f3_ids']),start=(leftx,612),end=(endx,612))
    dminus=replace(outminus,sequence=derived['index_r3'],nucleotide_ids=list(derived['index_r3_ids']),start=(endx,472),end=(leftx,472))
    duplex=Shot('index_duplex','The existing tails now pair and get copied',[partial_plus,partial_minus,third_plus,third_minus,shifted(dplus,dx=-1800),shifted(dminus,dx=1800)],
        [Event('unpair','third_parent_plus',pairing='old_f2',fraction_end=.15),Event('unpair','third_parent_minus',pairing='old_r2',fraction_end=.15),Event('separate','first_plus',target=(-2200,0),fraction_start=.15,fraction_end=.35),Event('separate','first_minus',target=(2200,0),fraction_start=.15,fraction_end=.35),Event('anneal','index_plus',target=(1800,0),fraction_start=.2,fraction_end=.45),Event('anneal','index_minus',target=(-1800,0),fraction_start=.2,fraction_end=.45),bind('third_f','index_plus',.45,.55),bind('third_r','index_minus',.45,.55)],
        'The simulated substrates now contain both indexed ends. Primers can pair through their full matching sequence, including the former overhangs. Extension produces fully paired complete duplexes. Further source PCR cycles reuse these products through the same reaction rule.',
        action=8,pairings=[bp('third_parent_plus',left_overhang,'first_minus',len(partial_minus.sequence)-1,len(partial_minus.sequence),'old_f2',source,exact=False),bp('third_parent_minus',right_overhang,'first_plus',len(partial_plus.sequence)-1,len(partial_plus.sequence),'old_r2',source,exact=False),bp('index_plus',0,'third_parent_minus',total-1,total,'third_f',source,False,False),bp('index_minus',0,'third_parent_plus',total-1,total,'third_r',source,False,False)],source=source,scope='Substrate-derived full primer binding · complete duplex formation')
    extensions(duplex,derived['index_duplex_before'],derived['index_cycle3_after'],derived['index_cycle3_operations'])
    # Use ALL canonical declared final components, resolving only the explicitly
    # synthetic insert. The sequence/feature inventory is independently compared.
    components=[]
    names={'illumina.p5':'Illumina P5','index.i5':'i5','nextera.s5':'s5','nextera.me':'ME','insert':'cDNA','nextera.s7':'s7','index.i7':'i7','illumina.p7':'Illumina P7'}
    state=derived['final_state']
    for id in derived['final_plus_ids']:
        nt=state.atoms[id]
        if len(nt.features)!=1 or nt.features[0]not in names:raise ValueError('unannotated generated final-library nucleotide')
        feature=nt.features[0]
        if not components or components[-1]['feature']!=feature:
            components.append({'name':names[feature],'feature':feature,'sequence':'','placeholder':bool(nt.symbol)})
        components[-1]['sequence']+=nt.base
    # The independently declared final-library model is an assertion only. It
    # cannot supply or repair the sequence generated by the reaction engine.
    canonical=''.join(insert if seg.name=='cDNA'else derived['i5_oligo']if seg.name=='i5'else rc(derived['i7_oligo'])if seg.name=='i7'else seg.top for seg in m.final_library().segments)
    if canonical!=plus:raise ValueError('simulated product disagrees with source-declared canonical library')
    expected=''.join(c['sequence']for c in components)
    if expected!=plus:raise ValueError('animated final molecule differs from protocol final_library component assembly')
    # Wrap at a component boundary so every nucleotide remains legible. Chunks
    # are explicitly continued views of one molecule, not additional DNA species.
    cut=sum(len(c['sequence'])for c in components[:5])
    strands=[];pairs=[];events=[];cursor=0
    for row,(lo,hi) in enumerate([(0,cut),(cut,total)]):
        y=370+row*240;step=11;origin=130
        seq=plus[lo:hi];regions=[]
        for j,c in enumerate(components):
            start=sum(len(x['sequence'])for x in components[:j]);end=start+len(c['sequence'])
            if lo<=start<hi:
                label=c['name']
                if label in {'i5','i7'}:label+=' · '+('N501'if label=='i5'else'N701')+' example'
                color='bio'if c['name']in{'i5','i7'}else 'dna'if c['name']=='cDNA'else 'adapter'if c['name']in{'ME','s5','s7'}else 'primer'
                regions.append(Region(start-lo,end-lo,label,color,-1,label_offset=(0,-65 if j%2==0 else -105)))
        upper=Strand(f'library_{row}',seq,(origin,y),(origin+len(seq)*step,y),'dna','',schematic=False,regions=regions,continued_5prime=lo>0,continued_3prime=hi<total,nucleotide_ids=list(derived['final_plus_ids'][lo:hi]))
        lower=Strand(f'library_bottom_{row}',rc(seq),(origin+len(seq)*step,y+24),(origin,y+24),'dna','',schematic=False,continued_5prime=hi<total,continued_3prime=lo>0,nucleotide_ids=list(derived['final_minus_ids'][total-hi:total-lo]),
                     bands=[(len(seq)-r.end,len(seq)-r.start,r.color)for r in regions])
        strands.extend([upper,lower]);pairs.append(bp(upper.id,0,lower.id,len(seq)-1,len(seq),f'library_pair_{row}',source,exact=False))
        for reg in regions:events.append(Event('highlight',upper.id,span=(reg.start,reg.end),fraction_start=cursor/len(components),fraction_end=(cursor+1)/len(components)));cursor+=1
    final=Shot('library','Complete library: every defined component and base',strands,events,
        'This is one molecule shown on two continued rows. It contains P5, the N501 example i5 index, s5, mosaic end, cDNA insert, reverse mosaic end, reverse s7, the N701 example i7 index and reverse P7. Both strands show every modeled nucleotide. The eight-base indexes are source-kit examples, not recovered study-well assignments.',
        action=10,hold=4,pairings=pairs,source=source,scope='One molecule, continued below · synthetic insert · source-era N501/N701 example indexes')
    return [first,second,duplex,final],components
