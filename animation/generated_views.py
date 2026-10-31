"""Small molecular views whose sequences/connectivity come only from reactions."""
from dataclasses import replace
from .model import Strand,Shot,Event,Region,Pairing,Junction
from .simulation_adapters import generated,evidence


def sequence(state,ids):return ''.join(state.atoms[i].base for i in ids)


def strand(name,state,ids,x,y,*,cell=10,direction=1,label='',initial=None,regions=()):
    ids=list(ids);types=[state.atoms[i].sugar for i in ids];sugars=[]
    for j,t in enumerate(types):
        if sugars and sugars[-1][2]==t:sugars[-1]=(sugars[-1][0],j+1,t)
        else:sugars.append((j,j+1,t))
    polymer=types[0]if len(set(types))==1 else'mixed'
    nicks=[j for j in range(1,len(ids))if ids[j]in state.atoms and state.previous.get(ids[j])!=ids[j-1]]
    # Only initial physical connectivity determines nick/circle status. For a
    # growing strand the newly added nucleotides are supplied by the after-state;
    # generated molecular frames override their future backbone display.
    return Strand(name,sequence(state,ids),(x,y),(x+direction*cell*len(ids),y),'rna'if polymer=='RNA'else'dna',label,
                  polymer=polymer,five_prime=state.five.get(ids[0],'none'),three_prime=state.three.get(ids[-1],'none'),
                  initially_visible=initial,nicks=nicks,regions=list(regions),nucleotide_ids=ids,sugars=sugars if polymer=='mixed'else[])


def extension(id,title,op,output_ids,template_ids,*,source,enzyme_label=None,narration='',scope='',cell=10,include_template_partner=False):
    before,after=op.before,op.after;output_ids=list(output_ids);template_ids=list(template_ids)
    initial=len(output_ids)-len(op.added);T=len(template_ids)
    first=after.atoms[op.added[0]].template if op.added else ''
    start=template_ids.index(first)if first else 0
    # A real paired primer base anchors the new strand. RNA transcription has
    # no primer and begins at its declared physical initiation nucleotide.
    output_left=120+(T-start-1-initial)*cell
    template_left=120
    shift=max(0,100-min(template_left,output_left));template_left+=shift;output_left+=shift
    template=strand('template',before,template_ids,template_left+T*cell,430,cell=cell,direction=-1)
    out=strand('product',after,output_ids,output_left,454,cell=cell,initial=initial)
    tracks=[template,out]
    if include_template_partner:
        partner=[before.paired.get(i)for i in reversed(template_ids)]
        if all(partner):tracks.insert(0,strand('coding',before,partner,template_left,406,cell=cell))
    kind='transcribe'if op.kind=='transcribe'else'rca'if op.kind=='rca'else'extend'
    e=Event(kind,'product',enzyme_label or op.details['enzyme']['name'],template='template',template_start=start,primer_length=initial,primer_binding_length=1 if initial else 0,
            exact=not(set(sequence(after,output_ids)+sequence(before,template_ids))-set('ACGTU')),molecular_operation=op.to_dict())
    shot=Shot(id,title,tracks,[e],narration or title,source=source,scope=scope or'Generated from actual substrate connectivity',action=7)
    generated(shot,before);return shot


def rnaseh():
    from .simulator import Simulator
    from .model import Movie
    source='Lee et al.2022 doi:10.1093/nar/gkab1064, H1-class RNA:DNA hybrid substrate; synthetic sequence and illustrative cut positions'
    s=Simulator();rna=s.add('rna','AUGCAGUACCGUAGCUAGCUACGAUCGGAUCC'*3,sugar='RNA',source=source)
    # These are the two explicitly supplied strands of the starting hybrid,
    # not a claimed reverse-transcription product appearing without synthesis.
    from .model import complement
    dna=s.add('dna',complement(''.join(s.atoms[i].base for i in rna)[::-1]),source=source)
    s.duplex('rna','dna',source=source);before=s.snapshot()
    cut1=s.rnase_h1([rna[31]],source=source);cut2=s.rnase_h1([rna[63]],source=source)
    a=strand('rna',before,rna,224,440,cell=12,regions=[Region(0,len(rna),'RNA','rna',-1,label_offset=(0,-60))])
    b=strand('dna',before,dna,224+len(dna)*12,464,cell=12,direction=-1,regions=[Region(0,len(dna),'DNA · retained','dna',1,label_offset=(0,60))])
    shot=Shot('cleave','RNase H cuts the RNA within an RNA–DNA hybrid',[a,b],
        [Event('cleave','rna','RNase H1',bond=32,fraction_start=.1,fraction_end=.4,molecular_operation=cut1.to_dict()),Event('cleave','rna','RNase H1',bond=64,fraction_start=.5,fraction_end=.8,molecular_operation=cut2.to_dict())],
        'RNase H1 cleaves internal RNA backbone bonds in an RNA DNA hybrid. DNA is retained. The cut RNA fragments remain RNA and stay paired here; cleavage does not automatically dissolve their hybrid contacts. The resulting RNA ends carry three-prime hydroxyl and five-prime phosphate groups.',
        source=source,scope='Illustrative H1 cleavage sites',action=9,hold=4)
    generated(shot,before)
    return Movie('RNase H', [shot],provenance={'representation':'illustrative_example','molecular_simulation':s.trace(),'source':source,'limits':['Selected H1 cleavage sites, not a sequence-specificity or processivity prediction. H2 is not simulated.','RNA–DNA hybrid supplied as input; no native CEL-Seq DNA-polymerase identity is inferred.']}).validate()


def celseq(m,qualify,slug):
    from .other_simulations import celseq as simulate
    r=simulate(m);s=r['simulator'];source='CEL-Seq source T7-bearing primer and illustrative insert; '+getattr(m,'AUTHORITY','doi:10.1016/j.celrep.2012.08.003')
    first=extension('ivt','T7 starts within the promoter’s terminal GGG',r['op1'],r['rna1'],r['bottom'],source=source,include_template_partner=True,
                    narration='T7 RNA polymerase initiates at the promoter terminal GGG, then follows the actual DNA template toward its five-prime end. It makes RNA from supplied ribonucleotides without a primer. Both DNA strands remain present.',cell=10)
    # Promoter annotation is attached to its physical nucleotide interval.
    coding=next(t for t in first.strands if t.id=='coding')
    offset=list(r['top']).index(r['promoter'][0]);coding.regions=[Region(offset,offset+len(r['promoter']),'T7 promoter','bio',-1,label_offset=(0,-55))]
    second=extension('repeat','The retained DNA can initiate another transcript',r['op2'],r['rna2'],r['bottom'],source=source,include_template_partner=True,
                     narration='The same unchanged promoter-bearing DNA can initiate another RNA molecule. The earlier transcript remains a separate RNA chain. This reaction reuses DNA rather than converting it into RNA.',cell=10)
    old=strand('earlier_rna',r['second_before'],r['rna1'],first.strands[-1].start[0],650,cell=10,label='Earlier RNA')
    old.label_position=(800,730);second.strands.append(old)
    movie=qualify(m,slug,'CEL-Seq',[first,second],{'ivt':['second_strand','ivt'],'repeat':['ivt']},limits=['Synthetic short insert and barcode variable; source RT primer supplies T7 initiation sequence.','This focused transcription view starts from a prepared duplex. Native second-strand kit enzyme identity remains bounded in the source review.'])
    return evidence(movie,s)


def gibson(m,qualify,slug):
    from .other_simulations import gibson as simulate
    r=simulate(m);s=r['simulator'];source=m.SOURCE+'; '+m.EXAMPLE
    tracks=[];events=[]
    for row,key in enumerate(['a','b']):
        for suffix,ids,direction,y in [('top',r['original'][key][0],1,340+row*230),('bottom',r['original'][key][1],-1,364+row*230)]:
            x=270 if direction==1 else 270+len(ids)*30
            tracks.append(strand(key+'_'+suffix,r['before'],ids,x,y,cell=30,direction=direction))
    for ids,op in zip([r['original']['a'][0],r['original']['a'][1],r['original']['b'][0],r['original']['b'][1]],r['resect_ops']):
        name=next(t.id for t in tracks if t.nucleotide_ids==list(ids))
        events.append(Event('resect',name,'5′ exonuclease',amount=len(op.details['removed']),molecular_operation=op.to_dict()))
    resect=Shot('resect','Resect actual 5′ ends to expose 3′ overlaps',tracks,events,'Five-prime exonuclease removes nucleotides from the two supplied duplexes. The exact stopping budgets are invented for this toy example. Removed nucleotide identities remain in the product pool.',source=source,scope='Invented DNA parts · explicit illustrative resection budgets',action=7)
    generated(resect,r['before'])
    # Display the same processed fragments, not newly fabricated single strands.
    tracks=[]
    for key,dy,dx in [('a',0,0),('b',210,0)]:
        for suffix,ids,direction in [('top',r['processed'][key][0],1),('bottom',r['processed'][key][1],-1)]:
            x=(380+len(m.LEFT)*20)if direction==1 else 380+(len(m.LEFT)+len(m.A))*20
            tracks.append(strand(key+'_'+suffix,r['anneal_before'],ids,x,420+dy+(24 if direction<0 else 0),cell=20,direction=direction))
    # Move both strands of part B as one duplex. Pairing is added only from the
    # generated annealing operation, independent of its approach coordinates.
    approach=Shot('anneal','Anneal complementary 3′ overlaps',[*tracks],
        [Event('anneal','b_top',target=(len(m.A)*20,-210),fraction_end=.7),Event('anneal','b_bottom',target=(len(m.A)*20,-210),fraction_end=.7),Event('pair','a_top',pairing='overlap',fraction_start=.7,molecular_operation=r['anneal'].to_dict())],
        'The exposed complementary overlaps anneal while each fragment retains its original duplex body. This exact-overlap toy has no nucleotide gap: a polymerase encounter adds zero bases because the downstream body is already paired.',source=source,scope='Same retained substrate nucleotides · zero-gap result generated by polymerase',action=7,pairings=[Pairing('a_top',len(m.A),'b_bottom',len(r['processed']['b'][1])-1,len(m.RIGHT),id='overlap')])
    generated(approach,r['anneal_before'])
    upper=strand('joined_top',s.snapshot(),r['joined_top'],400,440,cell=20)
    lower=strand('joined_bottom',s.snapshot(),r['joined_bottom'],400+(len(m.A)+len(m.RIGHT)+len(m.B))*20,464,cell=20,direction=-1)
    upper.nicks=[len(r['processed']['a'][0])];lower.nicks=[len(r['processed']['b'][1])]
    seals=[]
    for name,ids,op in [('joined_top',r['joined_top'],r['seals'][0]),('joined_bottom',r['joined_bottom'],r['seals'][1])]:
        a,b=next(iter(set(op.after.next.items())-set(op.before.next.items())));bond=list(ids).index(b)
        seals.append(Event('seal',name,'DNA ligase',bond=bond,junction=Junction(bond,source=source),fraction_start=.1 if name=='joined_top'else.55,fraction_end=.5 if name=='joined_top'else 1,molecular_operation=op.to_dict()))
    seal=Shot('seal','Seal both nicks at the selected joint',[upper,lower],seals,'Ligase seals the actual free three-prime hydroxyl and five-prime phosphate ends on each strand. No nucleotide is inserted. This clip follows one joint; the remaining vector junction is outside this view.',source=source,scope='One selected joint · remaining cohesive ends retained · no invented gap fill',action=7)
    generated(seal,r['seal_before'])
    movie=qualify(m,slug,'Gibson assembly',[resect,approach,seal],{'resect':['a_process','b_process'],'anneal':['anneal'],'seal':['repair','seal']},limits=[m.EXAMPLE,'One selected junction is followed, not a completed circular-vector simulation.'])
    return evidence(movie,s)
