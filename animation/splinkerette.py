"""Native splinkerette animation projected from the protocol's actual reactions."""
from dataclasses import replace
from .protocols import load,qualify
from .model import Shot,Event,Pairing,Region,Junction
from .generated_views import strand
from .simulation_adapters import generated

SLUG='splinkerette-pcr__10.1371+journal.pone.0010168'
CELL=8.0
LEFT=170.0
TOP_Y=450.0
BOTTOM_Y=474.0


def movie():
    m=load(SLUG,'splinkerette.py');r=m.simulate();s=r['simulator'];final=s.snapshot();source=m.SOURCE
    N=len(r['linked']);G=len(r['genome_top']);x=lambda col:LEFT+CELL*col
    # A base's horizontal coordinate is inherited from its actual template.
    # Only source molecules receive initial coordinates.
    columns={i:j for j,i in enumerate(r['parent_top'])}
    columns.update({i:len(r['parent_top'])-1-j for j,i in enumerate(r['parent_bottom'])})
    columns.update({i:G+j for j,i in enumerate(r['adapter_top'])})
    columns.update({i:N-1-j for j,i in enumerate(r['adapter_bottom'])})
    for op in s.history:
        for atom in op.added:
            nt=op.after.atoms[atom]
            if nt.template:columns[atom]=columns[nt.template]
        # Primers obtain position from their *binding* partners, not names.
        if op.kind=='anneal':
            for a,b in op.to_dict()['pairing_added']:
                if a in columns:columns[b]=columns[a]
                elif b in columns:columns[a]=columns[b]

    def track(name,ids,snapshot,y,direction=1,*,dx=0,initial=None,label='',regions=()):
        ids=list(ids);col=columns[ids[0]]+(1 if direction<0 else 0)
        t=strand(name,snapshot,ids,x(col)+dx,y,cell=CELL,direction=direction,initial=initial,label=label,regions=regions)
        # Retained/copied template positions keep their visual identity across
        # generations; region labels remain local presentation choices.
        for j,i in enumerate(ids):
            tag='primer' if columns[i]<65 else 'bio' if columns[i]<G else 'adapter'
            if t.bands and t.bands[-1][2]==tag:t.bands[-1]=(t.bands[-1][0],j+1,tag)
            else:t.bands.append((j,j+1,tag))
        return t

    def gt(snapshot,y=TOP_Y):
        return track('genome_top',r['genome_top'],snapshot,y,regions=[
            Region(0,65,'Known insertion end · T1 / T2 / T3','primer',-1,label_offset=(-20,-95)),
            Region(65,G-1,'Unknown genomic flank','bio',-1,label_offset=(-75,-50))])
    def gb(snapshot,y=BOTTOM_Y):return track('genome_bottom',r['genome_bottom'],snapshot,y,-1)
    def linked(snapshot,y=BOTTOM_Y):
        t=track('linked',r['linked'],snapshot,y,-1)
        t.sealed_marks=[len(m.BOT)]
        t.regions=[Region(0,28,'S1 sequence · no complementary site yet','adapter',1,label_offset=(0,105)),
                   Region(28,53,'S2','adapter',1,label_offset=(0,55))]
        return t

    foldroute=[(0,x(G),TOP_Y),(27,x(G+27),TOP_Y),(30,x(G+30),365),
               (37,x(G+37),365),(39,x(G+39),353),(41,x(G+37),341),(48,x(G+30),341)]
    def adapter(snapshot,folded=True,dx=0):
        t=track('adapter_top',r['adapter_top'],snapshot,TOP_Y,dx=dx)
        t.regions=[Region(0,4,'5′ OH · GATC','adapter',-1,label_offset=(-50,-150)),
                   Region(30,48,'TOP hairpin','bio',-1,label_offset=(50,-70))]
        if folded:t.path=[(i,a+dx,b)for i,a,b in foldroute]
        return t

    def attach(shot,op,kind,a,b,*,start=0,end=.15):
        """Derive the visual validation anchor from a real pairing delta."""
        byid={t.id:t for t in shot.strands};aa=byid[a];bb=byid[b]
        pairs=op.to_dict()['pairing_removed'if kind=='unpair'else'pairing_added']
        hit=next(((aa.nucleotide_ids.index(i),bb.nucleotide_ids.index(j))
                   for i,j in pairs if i in aa.nucleotide_ids and j in bb.nucleotide_ids),None)
        if hit is None:
            hit=next(( (aa.nucleotide_ids.index(j),bb.nucleotide_ids.index(i))
                        for i,j in pairs if j in aa.nucleotide_ids and i in bb.nucleotide_ids),None)
        if hit is None:raise ValueError('reaction has no visible anchor')
        key=shot.id+'_'+str(len(shot.pairings))
        shot.pairings.append(Pairing(a,hit[0],b,hit[1],1,id=key,initially_bound=kind=='unpair',source=source))
        shot.events.append(Event(kind,a,pairing=key,fraction_start=start,fraction_end=end,molecular_operation=op.to_dict()))

    def extend(shot,op,*,start=.55,end=1):
        places={i:(t,j)for t in shot.strands for j,i in enumerate(t.nucleotide_ids)}
        first=op.added[0];t,j=places[first];template,k=places[op.after.atoms[first].template]
        p=len(t.sequence)-len(op.added)
        t.initially_visible=p
        shot.events.append(Event('extend',t.id,'Phusion',template=template.id,template_start=k,
            primer_length=p,primer_binding_length=1,exact=True,fraction_start=start,fraction_end=end,
            enzyme_class='polymerase',enzyme_label='Phusion',molecular_operation=op.to_dict()))

    def introduce(shot,stock,strand):
        shot.events.append(Event('introduce',strand,fraction_end=.01,
            molecular_operation=r['reagent_operations'][stock].to_dict()))

    def shot(id,title,tracks,snapshot,narration,*,action=8,events=(),scope=None):
        q=Shot(id,title,tracks,list(events),narration,action=action,hold=2.5,transition=1.8,
               source=source,scope=scope or'Source oligos · synthetic insertion locus · selected chemistry')
        generated(q,snapshot);return q

    parent=[track('genome_top',r['parent_top'],r['genome_initial'],TOP_Y),
            track('genome_bottom',r['parent_bottom'],r['genome_initial'],BOTTOM_Y,-1)]
    parent[0].regions=[*gt(r['genome_initial']).regions,
                       Region(G-1,G+5,'BstYI · A/GATCC','new',1,label_offset=(80,85))]
    digest=shot('digest','A restriction cut supplies a genomic 5′ phosphate',parent,r['genome_initial'],
        'BstYI cuts the synthetic locus to expose a four-base GATC overhang. The known transposon end and unknown flank remain on the same selected fragment. A cut inside the transposon is not required.',action=6)
    digest.transition=0
    for j,t in enumerate(parent):
        bond=next(t.nucleotide_ids.index(b) for a,b in r['digest'].to_dict()['backbone_removed'] if a in t.nucleotide_ids and b in t.nucleotide_ids)
        digest.events.append(Event('cleave',t.id,'BstYI' if j==0 else '',bond=bond,fraction_start=.2,fraction_end=.7,
                                   enzyme_class='nuclease',enzyme_label='BstYI' if j==0 else '',molecular_operation=r['digest'].to_dict()))
    other_t=track('other_top',r['off_top'],r['digested'],TOP_Y)
    other_b=track('other_bottom',r['off_bottom'],r['digested'],BOTTOM_Y,-1)
    separate=shot('separate','Expose the compatible GATC end',[gt(r['digested']),gb(r['digested']),other_t,other_b],r['digested'],
        'The cut fragments separate. The other fragment moves out of view but remains in the reaction inventory. This movement does not degrade DNA.',action=5)
    attach(separate,r['bridge'],'unpair','other_top','genome_bottom',end=.3)
    separate.events.extend([Event('transfer',name,target=(1800,0),fraction_start=.3)for name in ['other_top','other_bottom']])

    at=adapter(r['adapter_unfolded'],False,dx=1800)
    ab=track('adapter_bottom',r['adapter_bottom'],r['adapter_unfolded'],BOTTOM_Y,-1,dx=1800)
    ab.regions=[Region(0,28,'BOT · S1 arm','adapter',1,label_offset=(0,95)),Region(28,53,'S2','adapter',1,label_offset=(0,50))]
    folding=shot('adapter','An unmodified adapter folds into an asymmetric splinkerette',
        [gt(r['adapter_unfolded']),gb(r['adapter_unfolded']),at,ab],r['fragment_ready'],
        'Both adapter oligos enter with five-prime hydroxyl ends. TOP forms its terminal hairpin and pairs with the recessed end of BOT. The long BOT arm remains single stranded. Folding is a selected complementary geometry, not an energy prediction.',action=9)
    introduce(folding,'adapter_top','adapter_top');introduce(folding,'adapter_bottom','adapter_bottom')
    folding.events.extend([Event('anneal',name,target=(-1800,0),fraction_end=.3)for name in ['adapter_top','adapter_bottom']])
    # Routes are local to the moving strand; translation is evaluated separately.
    folding.events.append(Event('fold','adapter_top',route=[(i,a+1800,b)for i,a,b in foldroute],fraction_start=.3,fraction_end=.65,molecular_operation=r['fold'].to_dict()))
    attach(folding,r['core'],'pair','adapter_top','adapter_bottom',start=.65,end=.95)

    seal=shot('ligate','Ligase seals one strand; the opposite nick stays open',
        [gt(r['adapter_ready']),adapter(r['adapter_ready']),linked(r['adapter_ready'])],r['adapter_ready'],
        'The GATC overhang anneals to the restriction end. T4 DNA ligase uses ATP to join the recessed BOT three-prime hydroxyl to the genomic five-prime phosphate. The opposite junction stays nicked because TOP has a five-prime hydroxyl.',action=8)
    attach(seal,r['dock'],'pair','adapter_top','linked',start=.05,end=.35)
    seal.events.append(Event('seal','linked','T4 ligase',bond=len(m.BOT),junction=Junction(len(m.BOT),source=source),
        fraction_start=.4,fraction_end=.85,enzyme_class='ligase',enzyme_label='T4 ligase',molecular_operation=r['seal'].to_dict()))

    t1=track('first_top',r['first_t1_primer'],r['first_before'],TOP_Y,dx=-1800,
             regions=[Region(0,len(m.T1),'T1 · insertion-specific primer','primer',-1,label_offset=(0,-55))])
    opening=shot('select','The insertion-specific primer starts the selective copy',
        [gt(r['first_before']),adapter(r['first_before']),linked(r['first_before']),t1],r['ligated'],
        'Heat separates the original strands and opens the hairpin. The ligated BOT and genomic bottom remain one covalent template. T1 enters from the left and binds its known insertion site. S1 still has no complementary site on this template.',action=7)
    introduce(opening,'first_t1','first_top');introduce(opening,'first_s1','linked')
    attach(opening,r['melt1'],'unpair','genome_top','linked',end=.15)
    opening.events.extend([Event('transfer','genome_top',target=(-1800,-300),fraction_start=.15,fraction_end=.45),
        Event('transfer','adapter_top',target=(1800,-300),fraction_start=.15,fraction_end=.45),
        Event('anneal','first_top',target=(1800,0),fraction_start=.25,fraction_end=.65)])
    attach(opening,r['t1_bind'],'pair','first_top','linked',start=.65,end=.9)
    # The unpaired hairpin opens continuously before its free strand leaves.
    opening.events.append(Event('fold','adapter_top',route=[(0,x(G),TOP_Y),(48,x(G+48),TOP_Y)],fraction_end=.15))

    t1copy=track('first_top',r['first_top'],r['first_copied'],TOP_Y,initial=len(m.T1),
        regions=[Region(0,len(m.T1),'T1','primer',-1,label_offset=(0,-55)),
                 Region(N-len(m.BOT),N,'Copied BOT · creates adapter-primer sites','new',-1,label_offset=(0,-90))])
    first=shot('copy','T1 copies through the junction into the BOT arm',
        [linked(r['first_bound']),t1copy],r['first_bound'],
        'Phusion extends the real T1 three-prime end through the unknown flank and the sealed junction, then copies BOT to its physical five-prime end. The new DNA now contains the complement needed for adapter priming. TOP is not incorporated.',action=9)
    extend(first,r['first_extend'],start=.05,end=.9)

    old=linked(r['s1_before']);old.regions=[]
    top=track('first_top',r['first_top'],r['s1_before'],TOP_Y)
    bottom=track('round1_bottom',r['round1_bottom'],r['round1'],BOTTOM_Y,-1,dx=1800,initial=len(m.S1),
                 regions=[Region(0,len(m.S1),'S1 · now has a binding site','primer',1,label_offset=(0,65))])
    reverse=shot('reverse','Only the newly copied strand lets S1 prime back',
        [top,old,bottom],r['s1_before'],
        'The next thermal cycle separates the first copy from its ligated template. S1 enters from the right, binds the newly created complement and copies back to the T1 end. The original ligated template remains in the inventory.',action=10)
    # s1_before is already melted; initialize from the preceding paired state
    # while keeping the stock S1 primer from the same molecular inventory.
    generated(reverse,r['first_copied'])
    attach(reverse,r['melt2'],'unpair','first_top','linked',end=.15)
    reverse.events.extend([Event('transfer','linked',target=(0,950),fraction_start=.15,fraction_end=.4),
        Event('anneal','round1_bottom',target=(-1800,0),fraction_start=.15,fraction_end=.4)])
    attach(reverse,r['s1_bind'],'pair','round1_bottom','first_top',start=.4,end=.5)
    extend(reverse,r['s1_extend'],start=.5)

    def cycle(number,topkey,bottomkey,topname,bottomname,top_y,bottom_y,previous_top_y,previous_bottom_y):
        before=r[f'nested{number}_before'];after=r[f'nested{number}_after']
        fkey=f'nested_t{number}_ids';qkey=f'nested_s{number}_ids'
        # T2 copies the lower antisense template. S2 copies the upper sense
        # template. Layout is independent of which generation owns that strand.
        sense=track(topname,r[topkey],before,previous_top_y)
        antisense=track(bottomname,r[bottomkey],before,previous_bottom_y,-1)
        newtop=track('nested_t'+str(number),r[fkey],after,bottom_y-24,dx=-1800,initial=len(m.T2),
                     regions=[Region(0,len(m.T2),'T2','primer',-1,label_offset=(0,-55))])
        newbottom=track('nested_s'+str(number),r[qkey],after,top_y+24,-1,dx=1800,initial=len(m.S2),
                        regions=[Region(0,len(m.S2),'S2','adapter',1,label_offset=(0,55))])
        q=shot('nested'+str(number),'Nested primers copy both physical templates' if number==1 else 'The next cycle establishes both nested boundaries',
            [sense,antisense,newtop,newbottom],r['round1'] if number==1 else r['nested1_after'],
            ('T2 and S2 bind inside the first-round product. Both three-prime ends extend. Each first nested copy continues to its actual template end, so one outer-primer boundary remains on each molecule.' if number==1 else
             'The first nested products become templates. Fresh T2 and S2 primers bind and both extend again. Copying these physical templates now creates the complementary products bounded by T2 and S2.'),action=11)
        introduce(q,'nested_t'+str(number),newtop.id);introduce(q,'nested_s'+str(number),newbottom.id)
        if number==1:
            attach(q,r['nested1_melt'],'unpair',topname,bottomname,end=.15)
            q.events.extend([Event('separate',topname,target=(0,top_y-previous_top_y),fraction_start=.15,fraction_end=.4),
                             Event('separate',bottomname,target=(0,bottom_y-previous_bottom_y),fraction_start=.15,fraction_end=.4)])
        else:
            # Keep the two old full-length parents visible until melting lets
            # them leave. Their strands are retained, never destroyed.
            oldtop=track('first_top',r['first_top'],before,320)
            oldbottom=track('round1_bottom',r['round1_bottom'],before,650,-1)
            q.strands.extend([oldtop,oldbottom])
            attach(q,r['nested2_melt'],'unpair',topname,'round1_bottom',end=.15)
            q.events.extend([Event('transfer','first_top',target=(0,-950),fraction_start=.15,fraction_end=.4),
                             Event('transfer','round1_bottom',target=(0,950),fraction_start=.15,fraction_end=.4)])
        q.events.extend([Event('anneal',newtop.id,target=(1800,0),fraction_start=.15,fraction_end=.45),
                         Event('anneal',newbottom.id,target=(-1800,0),fraction_start=.15,fraction_end=.45)])
        attach(q,r[f'nested_t{number}_bind'],'pair',newtop.id,bottomname,start=.45,end=.55)
        attach(q,r[f'nested_s{number}_bind'],'pair',newbottom.id,topname,start=.45,end=.55)
        for op in r[f'nested{number}_operations']:extend(q,op,start=.55)
        return q

    nested1=cycle(1,'first_top','round1_bottom','first_top','round1_bottom',320,650,TOP_Y,BOTTOM_Y)
    nested2=cycle(2,'nested_t1_ids','nested_s1_ids','nested_t1','nested_s1',626,344,626,344)
    # Select complementary endpoints by reassociation; not by fabricating or
    # truncating first-cycle molecules. Their earlier partners leave view.
    finaltop=track('final_top',r['final_top'],r['final_before'],320)
    finalbottom=track('final_bottom',r['final_bottom'],r['final_before'],650,-1)
    oldtop=track('nested_t1',r['nested_t1_ids'],r['final_before'],626)
    oldbottom=track('nested_s1',r['nested_s1_ids'],r['final_before'],344,-1)
    reassociate=shot('reassociate','The bounded strands reassociate into the final amplicon',
        [finaltop,finalbottom,oldtop,oldbottom],r['final_before'],
        'A representative endpoint duplex forms when complementary bounded strands reassociate. Earlier templates move out of view but remain in the model. The final amplicon contains copied adapter sequence, not the original hairpin oligo.',action=7)
    attach(reassociate,r['final_melt'],'unpair','final_top','nested_s1',end=.15)
    reassociate.events.extend([Event('transfer','nested_t1',target=(0,950),fraction_start=.15,fraction_end=.55),
        Event('transfer','nested_s1',target=(0,-950),fraction_start=.15,fraction_end=.55),
        Event('transfer','final_top',target=(0,TOP_Y-320),fraction_start=.15,fraction_end=.65),
        Event('transfer','final_bottom',target=(0,BOTTOM_Y-650),fraction_start=.15,fraction_end=.65)])
    attach(reassociate,r['final_pair'],'pair','final_top','final_bottom',start=.65,end=1)
    a=track('final_top',r['final_top'],r['final'],TOP_Y)
    b=track('final_bottom',r['final_bottom'],r['final'],BOTTOM_Y,-1)
    lo=len(m.T2);t3=lo+4;flank=t3+len(m.T3)+4;junction=G-columns[r['final_top'][0]]
    a.regions=[Region(0,lo,'T2 boundary','primer',-1,label_offset=(0,-65)),
               Region(lo,t3,'Spacer','bio',1,label_offset=(0,210)),
               Region(t3,t3+len(m.T3),'T3 · Sanger reads toward flank →','primer',-1,label_offset=(40,-125)),
               Region(t3+len(m.T3),flank,'Spacer','bio',1,label_offset=(0,130)),
               Region(flank,junction,'Synthetic genomic flank','bio',1,label_offset=(0,75)),
               Region(junction,junction+4,'GATC junction','new',1,label_offset=(10,210)),
               Region(junction+4,len(a.sequence)-len(m.S2),'Copied core','adapter',-1,label_offset=(-25,-65)),
               Region(len(a.sequence)-len(m.S2),len(a.sequence),'S2 boundary','adapter',-1,label_offset=(70,-125))]
    product=shot('product','A flanking amplicon ready for nested-primer Sanger sequencing',
        [a,b],r['final'],
        'T3 reads outward from the known transposon into the recovered genomic flank. The final molecule retains the restriction junction and the copied S2 adapter boundary. The source cleans the PCR product before Sanger sequencing. No Illumina library or sample index is implied.',action=4)
    product.hold=6
    shots=[digest,separate,folding,seal,opening,first,reverse,nested1,nested2,reassociate,product]
    steps={'digest':['digest'],'separate':['digest'],'adapter':['fold','anneal_adapter'],'ligate':['ligate'],
           'select':['selective_copy'],'copy':['selective_copy'],'reverse':['adapter_copy'],
           'nested1':['nested_pcr'],'nested2':['nested_pcr'],'reassociate':['nested_pcr'],'product':['nested_pcr']}
    return qualify(m,SLUG,'Splinkerette PCR',shots,steps,
        limits=[m.EXAMPLE,'Hairpin is a selected complementary geometry; folding thermodynamics, background suppression and whole-vessel PCR populations are not simulated.',
                'One insertion-facing adapter junction is expanded. Other fragments and parental PCR strands remain in the molecular inventory when outside the viewport.',
                'Formal endpoint is an amplicon for Sanger analysis; source cleanup and sequencing are described, not animated.'],
        extra={'representation':'illustrative_example','molecular_simulation':s.trace(),
               'chemistry_authority':'The protocol page and animation consume one shared nucleotide/backbone reaction trajectory.'})
