"""Selected Gibson, ATAC and RNase-H mechanisms from explicit molecular inputs.

Source-specific choices live here; no desired polymerase product is supplied.
"""
from dataclasses import replace
from .model import Movie,Shot,Event,Region,Pairing,Junction,complement
from .protocols import load,qualify
from .simulator import Simulator,Polymerase
from .generated_views import strand,sequence
from .simulation_adapters import generated,evidence


def gibson_simulation():
    from .other_simulations import gibson
    m=load('gibson-assembly__10.1038+nmeth.1318','gibson_assembly.py')
    return m,gibson(m)


def gibson():
    m,r=gibson_simulation();s=r['simulator'];source=m.SOURCE+'; '+m.EXAMPLE
    tracks=[];events=[]
    for row,key in enumerate(['a','b']):
        for side,direction,yy in [('top',1,380+row*210),('bottom',-1,404+row*210)]:
            ids=r['original'][key][0 if side=='top'else 1]
            t=strand(key+'_'+side,r['before'],ids,400 if direction==1 else 840,yy,cell=20,direction=direction)
            if side=='top':t.regions=[Region(6,16,'Part '+key.upper(),'dna',-1,label_offset=(0,-65)),Region(16,22,'Overlap','bio',-1,label_offset=(40,-35))]
            op=r['resect_ops'][row*2+(side=='bottom')]
            events.append(Event('resect',t.id,'5′ exonuclease',amount=len(op.details['removed']),enzyme_label_offset=(0,-65 if side=='top'else 70),molecular_operation=op.to_dict()))
            tracks.append(t)
    resect=Shot('resect','Expose complementary 3′ overlaps',tracks,events,
        'Five-prime exonuclease removes the selected terminal bases. The three-prime overlaps and duplex bodies remain. Removed nucleotides remain recorded as degradation products. The sequences and stopping positions are explicitly synthetic.',source=source,scope='INFERRED toy sequences and resection boundaries',action=7)
    generated(resect,r['before'])
    if hasattr(resect,'product_pool'):resect.product_pool=True
    binding=Pairing('a_top',16,'b_bottom',21,6,id='overlap',initially_bound=False,source=source)
    anneal=Shot('anneal','Bring the exposed overlaps together',[replace(t)for t in tracks],
        [Event('anneal','b_top',target=(320,-210),fraction_end=.7),Event('anneal','b_bottom',target=(320,-210),fraction_end=.7),Event('pair','a_top',pairing='overlap',fraction_start=.7,molecular_operation=r['anneal'].to_dict())],
        'Complementary three-prime overlaps anneal while each fragment retains its duplex body. In this exact-end example polymerase reaches the already paired downstream body immediately: no nucleotide gap exists, so no bases are added.',source=source,scope='One selected joint · exact-overlap example has no gap to fill',action=7,transition=1.8,pairings=[binding])
    generated(anneal,r['anneal_before'])
    if hasattr(anneal,'product_pool'):anneal.product_pool=True
    upper=strand('joined_top',r['seal_before'],r['joined_top'],520,380,cell=20)
    lower=strand('joined_bottom',r['seal_before'],r['joined_bottom'],1040,404,cell=20,direction=-1)
    upper.regions=[Region(0,10,'Part A','dna',-1,label_offset=(0,-70)),Region(10,16,'One shared overlap','bio',-1,label_offset=(10,-35)),Region(16,26,'Part B','dna',-1,label_offset=(0,-70))]
    seals=[]
    for name,ids,op in [('joined_top',r['joined_top'],r['seals'][0]),('joined_bottom',r['joined_bottom'],r['seals'][1])]:
        a,b=next(iter(set(op.after.next.items())-set(op.before.next.items())));bond=list(ids).index(b)
        (upper if name=='joined_top'else lower).sealed_marks=[bond]
        seals.append(Event('seal',name,'DNA ligase',bond=bond,junction=Junction(bond,source=source),fraction_start=.1 if name=='joined_top'else.5,fraction_end=.45 if name=='joined_top'else.85,enzyme_label_offset=(0,-100 if name=='joined_top'else 110),molecular_operation=op.to_dict()))
    seal=Shot('seal','Seal the remaining nicks',[upper,lower],seals,
        'DNA ligase seals both actual three-prime hydroxyl to five-prime phosphate junctions. Pink marks the joined bases and bond. The outer cohesive ends are retained: this is one selected joint, not a completed circular vector.',source=source,scope='Pink marks sealed bonds · outer ends retained',action=7,transition=1.8)
    generated(seal,r['seal_before'])
    if hasattr(seal,'product_pool'):seal.product_pool=True
    final=Shot('product','One overlap joins the two fragments',[replace(upper),replace(lower)],[],
        'Both backbones are continuous across the selected joint. Each overlap is present once, and the original fragment identities are conserved. A second vector junction lies outside this selected demonstration.',source=source,scope='One selected joint · no invented gap synthesis',action=2,hold=3,transition=1.8)
    generated(final,s.snapshot())
    if hasattr(final,'product_pool'):final.product_pool=True
    final.camera=(800,390,1.35)
    movie=qualify(m,'gibson-assembly__10.1038+nmeth.1318','Gibson assembly',[resect,anneal,seal,final],{'resect':['a_process','b_process'],'anneal':['anneal','repair'],'seal':['seal'],'product':['seal']},limits=[m.EXAMPLE,'One joint of the source-owned exact-overlap toy. Gap filling adds zero bases because the other fragment is already paired; no expected product is supplied.'])
    return evidence(movie,s)


def rnaseh_simulation():
    source='Lee et al.2022 doi:10.1093/nar/gkab1064; illustrative blunt RNA-DNA hybrid and cut sites'
    s=Simulator();rna=s.add('rna','AUGCAGUACCGUAGCUAGCUACGAUCGGAUCCAUGCAGUACCGUAGCU',sugar='RNA',source=source)
    dna=s.add('dna',complement(sequence(s.snapshot(),rna)[::-1]),source=source)
    s.duplex('rna','dna',source=source);before=s.snapshot()
    cuts=[s.rnase_h1([rna[i]],source=source)for i in [15,31]]
    cleaved=s.snapshot();release=s.unpair(rna[:16],source=source+' explicitly illustrative dissociation of first fragment; no RNase-driven melting claim')
    return dict(simulator=s,source=source,rna=rna,dna=dna,before=before,cuts=cuts,cleaved=cleaved,release=release,after=s.snapshot())


def rnaseh():
    r=rnaseh_simulation();s=r['simulator'];n=len(r['rna']);x=800-n*12/2
    a=strand('rna',r['before'],r['rna'],x,450,cell=12,regions=[Region(0,n,'RNA','rna',-1,label_offset=(0,-70))])
    b=strand('dna',r['before'],r['dna'],x+n*12,474,cell=12,direction=-1,regions=[Region(0,n,'DNA · retained','dna',1,label_offset=(0,75))])
    cut=Shot('cleave','RNase H1 cuts RNA in the hybrid',[a,b],
        [Event('cleave','rna','RNase H1',bond=16,fraction_start=.05,fraction_end=.4,enzyme_label_offset=(0,-80),molecular_operation=r['cuts'][0].to_dict()),Event('cleave','rna','RNase H1',bond=32,fraction_start=.5,fraction_end=.85,enzyme_label_offset=(0,-80),molecular_operation=r['cuts'][1].to_dict())],
        'RNase H1 recognizes consecutive RNA nucleotides paired to DNA and cleaves internal RNA backbone bonds. The DNA strand remains unchanged. These two cut positions are illustrative, not predictions of enzyme specificity.',source=r['source'],scope='Supplied blunt hybrid · illustrative H1 cleavage positions',action=8)
    generated(cut,r['before'])
    frag=replace(a,regions=[Region(0,16,'RNA fragment · 3′ OH','rna',-1,label_offset=(-20,-55)),Region(16,32,'5′ P · RNA · 3′ OH','rna',-1,label_offset=(10,-85)),Region(32,n,'5′ P · RNA','rna',-1,label_offset=(30,-55))])
    release=Shot('fragments','Cleavage makes fragments; DNA stays intact',[frag,replace(b)],
        [Event('unpair','rna',pairing='first_fragment',fraction_start=.1,fraction_end=.4,molecular_operation=r['release'].to_dict()),Event('displace','rna',span=(0,16),target=(-100,-155),fraction_start=.4,fraction_end=.9)],
        'Cleavage leaves RNA fragments with new three-prime hydroxyl and five-prime phosphate ends. The selected first fragment then loses its hybrid contacts and drifts away. This dissociation is an explicit illustrative choice, not automatic destruction of RNA by a cleavage event. The other fragments remain paired.',source=r['source'],scope='Illustrative fragment dissociation · no RNA converted to DNA or erased',action=7,hold=4,transition=1.8,pairings=[Pairing('rna',0,'dna',n-1,16,id='first_fragment',source=r['source'])])
    generated(release,r['cleaved'])
    return Movie('RNase H1',[cut,release],provenance={'representation':'illustrative_example','source':r['source'],'molecular_simulation':s.trace(),'limits':['Synthetic hybrid, selected endocuts and explicit dissociation; not a native specificity, rate or processivity model.','H1 four-RNA recognition rule only; H2 and commercial CEL-Seq enzyme identity are not inferred.']}).validate()


def atac_simulation(*,target=None,cuts=(28,88)):
    m=load('atac-seq__10.1038+nmeth.2688','atac_seq.py')
    nx=m.nx
    source='Buenrostro et al.2013 doi:10.1038/nmeth.2688 Online Methods; canonical Nextera ends; synthetic accessible-DNA sequence and selected cut positions'
    s=Simulator();s.supply_nt('DNA','ACGT',source=source)
    target=target or 'ACGTTGCAGATCTAGCGTACGATCGGACTGCA'*4
    top=s.add('target_top',target,source=source);bottom=s.add('target_bottom',complement(target[::-1]),source=source)
    s.duplex('target_top','target_bottom',source=source)
    r=dict(simulator=s,source=source,target=target,target_top_ids=top,target_bottom_ids=bottom,reagents={})
    for name,adapter in [('left_top_adapter',nx.ADAPTOR_S5),('left_bottom_adapter',nx.ADAPTOR_S5),('right_top_adapter',nx.ADAPTOR_S7),('right_bottom_adapter',nx.ADAPTOR_S7)]:
        r['reagents'][name]=s.add(name,adapter,source=source+' supplied canonical transferred end')
        short=name+'_nontransfer';r['reagents'][short]=s.add(short,complement(nx.ME[::-1]),five='unassigned',three='unassigned',source=source+' retained short nontransfer ME complement')
        s.anneal(short,name,binding=len(nx.ME),source=source)
    r['tag_before']=s.snapshot();a,b=cuts;r['cut_left']=a;r['cut_right']=b
    for name,position in [('left',a),('right',b)]:
        r['transfer_'+name]=s.tn5(name+'_top_adapter',name+'_bottom_adapter',top[position],bottom[len(bottom)-position-9],stagger=9,complex_id='tn5_'+name,source=source)
    for name in ['left','right']:r['release_'+name]=s.release_protein('tn5_'+name,source=source+' MinElute workup; generic protein-removal boundary')
    for name,position in [('left',a),('right',b)]:r[name+'_bridge_melt']=s.unpair(top[position:position+9],source=source+' illustrative separation of neighboring fragments; nucleotides retained')
    r['tag_top_ids']=s.chain(top[a]);r['tag_bottom_ids']=s.chain(bottom[len(bottom)-b-9]);r['tag_after']=s.snapshot();r['repair_before']=s.snapshot()
    pol=Polymerase('DNA polymerase',('DNA',),'DNA')
    for name,anchor,short in [('top',top[a],'right_bottom_adapter_nontransfer'),('bottom',bottom[len(bottom)-b-9],'left_top_adapter_nontransfer')]:
        r[name+'_gap']=s.extend(anchor,pol,source=source+' 72C extension boundary')
        r[name+'_me_release']=s.unpair(s.chain(short),source=source+' explicitly conditional ME-prime contact loss; displacement/digestion path unassigned')
        r[name+'_adapter_copy']=s.extend(anchor,pol,source=source+' retained 3prime end continues on available template')
    r['repair_after']=s.snapshot();r['repaired_top_ids']=s.chain(top[a]);r['repaired_bottom_ids']=s.chain(bottom[len(bottom)-b-9])
    for name in ['top','bottom']:r['repaired_'+name]=sequence(s.snapshot(),r['repaired_'+name+'_ids'])
    r['insert']=target[a:b+9]
    return r



def retain_offscreen(state,visible,previous):
    """Keep actual nonselected chains at their preceding offscreen coordinates."""
    from .evaluate import evaluate
    old,_=evaluate(previous,previous.timing()[2])
    locations={id:(st,j)for st in old.values()for j,id in enumerate(st.strand.nucleotide_ids)}
    covered={id for t in visible for id in t.nucleotide_ids}
    missing=set(state.atoms)-state.released-covered
    output=[]
    for head in sorted(id for id in missing if state.previous.get(id)not in missing):
        ids=[];id=head
        while id in missing:
            ids.append(id);id=state.next.get(id)
        track=strand('retained_'+str(len(output)),state,ids,0,0)
        route=[]
        for k in range(2*len(ids)+1):
            q=k/2;j=min(int(q),len(ids)-1);prior,at=locations[ids[j]]
            route.append((q,*prior.segment_point(at+q-j,at)))
        track.path=route;track.start=tuple(route[0][1:]);track.end=tuple(route[-1][1:]);output.append(track)
    return output

def atac():
    m=load('atac-seq__10.1038+nmeth.2688','atac_seq.py');nx=m.nx
    from .tagmentation import reaction
    r=atac_simulation();s=r['simulator'];m=load('atac-seq__10.1038+nmeth.2688','atac_seq.py')
    tag,fill,top,bottom,_=reaction(r['target'],nx.ADAPTOR_S5,nx.ADAPTOR_S7,nx.ME,r['cut_left'],r['cut_right'],cell=8,x0=200,derived=r)
    tag.scope='Synthetic accessible DNA · canonical Nextera ends · 9 bp stagger'
    for t in fill.strands:
        for region in t.regions:region.label=region.label.replace('cDNA insert','Genomic insert')
    tag.source=fill.source=r['source'];tag.narration='Two loaded Tn5 homodimers transfer adapter three-prime ends to target five-prime ends, nine bases apart. After the source purification boundary, neighboring DNA fragments slide aside while their nucleotide identities remain. We follow one productive heterologous fragment; neither target three-prime end has yet been extended.'
    # Straight transfer strands and their actual antiparallel ME′ partners
    # share one docking geometry. Flattening each strand independently would
    # preserve the old angled short-strand endpoint and stretch the duplex.
    adapter_tracks={track.id:track for track in tag.strands}
    for name,track in adapter_tracks.items():
        if '_adapter' not in name or name.endswith('_nontransfer'):continue
        track.start=(track.start[0],track.end[1])
        short=adapter_tracks[name+'_nontransfer']
        yy=track.end[1]+(24 if '_top_'in name else -24)
        short.start=(short.start[0],yy);short.end=(short.end[0],yy)
    generated(tag,r['tag_before']);tracks={t.id:t for t in tag.strands}
    for e in tag.events:
        if e.kind=='cleave':
            ids=tracks[e.strand].nucleotide_ids;bond=(ids[e.bond-1],ids[e.bond])
            e.molecular_operation=next(op.to_dict()for op in [r['transfer_left'],r['transfer_right']]if bond in set(op.before.next.items())-set(op.after.next.items()))
        elif e.kind=='strand_transfer':e.molecular_operation=r['transfer_left'if e.strand.startswith('left')else'transfer_right'].to_dict()
        elif e.kind=='release':e.molecular_operation=r['release_left'if e.strand.startswith('left')else'release_right'].to_dict()
        elif e.kind=='unpair':e.molecular_operation=r['left_bridge_melt'if e.pairing=='left_bridge'else'right_bridge_melt'].to_dict()
    fill.strands.extend(retain_offscreen(r['repair_before'],fill.strands,tag))
    generated(fill,r['repair_before']);fill.transition=1.8
    fill.scope='9 nt fill, then conditional ME′ release · chemistry of release unassigned'
    fill.narration='Each original recessed target three-prime end copies nine target bases and stops at the bound short mosaic-end complement. We then explicitly release those short oligo contacts, without assigning displacement versus digestion. The same three-prime ends continue copying the opposite attached adapter. No ligase or completed complementary strand is invented.'
    for e in fill.events:
        if e.kind=='unpair':e.molecular_operation=r['bottom_me_release'if e.pairing=='short_left'else'top_me_release'].to_dict()
    top.regions=[Region(0,len(nx.ADAPTOR_S5),'Transferred left end','adapter',-1,label_offset=(0,-65)),Region(len(nx.ADAPTOR_S5),len(nx.ADAPTOR_S5)+len(r['insert']),'Genomic insert','dna',-1,label_offset=(0,-65)),Region(len(nx.ADAPTOR_S5)+len(r['insert']),len(top.sequence),'Copied right end','new',-1,label_offset=(0,-65))]
    bottom.regions=[Region(0,len(nx.ADAPTOR_S7),'Transferred right end','adapter',1,label_offset=(0,75)),Region(len(bottom.sequence)-len(nx.ADAPTOR_S5),len(bottom.sequence),'Copied left end','new',1,label_offset=(0,75))]
    final=Shot('product','Both original 3′ ends have copied the opposite tag',[top,bottom],[],
        'The retained insert and its transferred tags are now complemented. This is an amplifiable pre-PCR product, not the final indexed sequencing library. Original target nucleotides, original adapters and newly copied bases retain their identities.',source=r['source'],scope='Pre-PCR duplex · subsequent indexed amplification outside this focused clip',action=2,hold=3,transition=1.8)
    final.strands.extend(retain_offscreen(r['repair_after'],final.strands,fill))
    generated(final,r['repair_after'])
    return evidence(qualify(m,'atac-seq__10.1038+nmeth.2688','ATAC-seq',[tag,fill,final],{'tagment':['tagment','purify_tagmented'],'fill':['extend'],'product':['extend']},limits=['Synthetic accessible target and selected cut positions; no nucleosome/accessibility/yield simulation.','Loaded adapter strands are canonical Nextera components, not a claim of complete disclosure of historical commercial transposome stocks.','ME-prime contact loss is explicitly conditional; displacement versus digestion and exact polymerase identity in the commercial master mix remain unassigned.','Only selected productive pre-PCR chemistry is animated; no flow-cell or sample-index sequence appears without later PCR.']),s)
