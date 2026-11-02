"""Selected, source-qualified transcription, capture/RT and bisulfite mechanisms.

Display inputs are explicit examples. Simulator operations determine births,
backbone connectivity and pairing; no expected product is supplied to an enzyme.
"""
from dataclasses import replace
from .model import Event,Shot,Region,Pairing,Decoration,complement
from .simulator import Simulator,Polymerase
from .generated_views import strand,sequence
from .simulation_adapters import generated,evidence
from .protocols import load,qualify

CEL='cel-seq__10.1016+j.celrep.2012.08.003'
DROP='drop-seq__10.1016+j.cell.2015.05.002'
RRBS='scrrbs__10.1101+gr.161679.113'
rc=lambda seq:complement(seq[::-1])


def celseq_simulation(body='ACGTGACCTGATCGTAGC'):
    m=load(CEL,'native.py');s=Simulator();source='Hashimshony2012 published mmc4 CEL-Seq Protocol2012-07-16, primer#1 and IVT sections'
    parts=m.rt_parts();bases=[];features=[];at=0
    for part in parts:
        value='CATCACGC'if part.name=='cell barcode'else'A'if part.name=='V anchor'else part.top
        bases.append(value)
        if part.name in {'cell barcode','T7 promoter'}:features.append((at,at+len(value),part.name))
        at+=len(value)
    coding=''.join(bases)+body
    top=s.add('coding',coding,features=features,source=source+'; source primer#1, V=A and illustrative cDNA insert')
    bottom=s.add('template',rc(coding),source=source+' prepared second-strand input; native enzyme not reconstructed')
    s.duplex('coding','template',source=source)
    # Propagate identifiers to the supplied complementary strand explicitly.
    for a,b in zip(top,reversed(bottom)):
        s.atoms[b]=replace(s.atoms[b],features=s.atoms[a].features)
    s.supply_nt('RNA','ACGU',source=source)
    promoter_span=next((a,b)for a,b,label in features if label=='T7 promoter');promoter=top[promoter_span[0]:promoter_span[1]]
    start_coding=promoter_span[1]-len(m.C.T7_TRANSCRIBED_G)
    before=s.snapshot();pol=Polymerase('T7 RNA polymerase',('DNA',),'RNA',de_novo=True)
    rna1,op1=s.transcribe('template',len(bottom)-1-start_coding,pol,promoter=promoter,source=source)
    between=s.snapshot();rna2,op2=s.transcribe('template',len(bottom)-1-start_coding,pol,promoter=promoter,source=source)
    return dict(simulator=s,module=m,source=source,top=top,bottom=bottom,before=before,between=between,after=s.snapshot(),rna1=rna1,rna2=rna2,op1=op1,op2=op2,promoter_span=promoter_span,start=start_coding,features=features)


def celseq():
    r=celseq_simulation();s=r['simulator'];src=r['source'];cell=10;left=170;y=430;n=len(r['top'])
    def dna(state):
        top=strand('coding',state,r['top'],left,y,cell=cell,regions=[Region(*r['promoter_span'],'T7 promoter','bio',-1,label_offset=(0,-65))])
        bottom=strand('template',state,r['bottom'],left+n*cell,y+24,cell=cell,direction=-1)
        return [top,bottom]
    anatomy=Shot('promoter','A paired T7 promoter recruits RNA polymerase',dna(r['before']),
        [Event('highlight','coding',span=r['promoter_span'])],
        'This prepared cDNA duplex carries the source T7 promoter, sequencing handle and cell barcode. We follow transcription of a short illustrative insert. The native second-strand preparation precedes this view.',source=src,scope='Prepared duplex · source primer #1 · illustrative insert',action=4,hold=2)
    anatomy.strands[0].regions.extend([Region(a,b,'Cell barcode','primer',-1,label_offset=(0,-65))for a,b,label in r['features']if label=='cell barcode'])
    generated(anatomy,r['before'])
    x=left+r['start']*cell
    def transcript(name,ids,initial,y=520):
        out=strand(name,r['after'],ids,x,y,cell=cell,initial=initial)
        out.color='rna';out.label='RNA transcript';out.label_position=(1050,y+75)
        out.bands=[(j,j+1,'primer')for j,id in enumerate(ids)if 'cell barcode'in r['after'].atoms[id].features]
        return out
    first=Shot('ivt','T7 makes RNA while retaining both DNA strands',dna(r['before'])+[transcript('rna',r['rna1'],0)],
        [Event('transcribe','rna','T7 RNA polymerase',template='template',template_start=len(r['bottom'])-1-r['start'],primer_length=0,actor_offset=(0,-65),enzyme_label_offset=(-65,-70),fraction_end=.75,molecular_operation=r['op1'].to_dict()),Event('separate','rna',target=(0,110),fraction_start=.8)],
        'T7 initiates within the promoter terminal GGG and makes RNA without a primer. The RNA inherits the barcode sequence from its DNA template. The unpaired transcript moves away, while both DNA strands remain available.',source=src,scope='De novo RNA synthesis · DNA retained',action=10,hold=2,transition=1.8)
    generated(first,r['before'])
    second=Shot('repeat','The same DNA produces another RNA molecule',dna(r['between'])+[transcript('earlier_rna',r['rna1'],len(r['rna1']),630),transcript('rna',r['rna2'],0)],
        [Event('transcribe','rna','T7 RNA polymerase',template='template',template_start=len(r['bottom'])-1-r['start'],primer_length=0,actor_offset=(0,-65),enzyme_label_offset=(-65,-70),molecular_operation=r['op2'].to_dict())],
        'The retained promoter supports another initiation. The earlier RNA is still present as a separate chain. Reusing DNA is linear amplification, rather than converting the DNA substrate into RNA.',source=src,scope='Two illustrative transcripts · no yield prediction',action=9,hold=3,transition=1.8)
    generated(second,r['between'])
    return evidence(qualify(r['module'],CEL,'CEL-Seq',[anatomy,first,second],{'promoter':['second_strand'],'ivt':['ivt'],'repeat':['ivt']},limits=['Source primer#1 barcode; V=A and short cDNA body are selected examples, not an observed sample molecule.','Starts from the prepared duplex; native second-strand-kit enzyme identity and subsequent fragmentation/library preparation are outside this focused clip.','T7 initiation includes the promoter terminal GGG. RNA departure is illustrative geometry, not an elongation-complex structure.']),s)


def dropseq_simulation(body='AUGCAGUACCGUAGCUAGCUACGC',barcode='ACGTACGATGCA',umi='GACTACGT'):
    m=load(DROP,'drop_seq.py');s=Simulator();src='Macosko2015 Cell published TableS6 SeqB and Supplemental Experimental Procedures reverse transcription'
    values=[];features=[];regions=[];at=0
    for part in m.bead_primer():
        value=barcode if part.name=='cell barcode'else umi if part.name=='UMI'else part.top
        values.append(value)
        if part.name in {'cell barcode','UMI'}:features.append((at,at+len(value),part.name))
        regions.append((at,at+len(value),part.name));at+=len(value)
    primer_sequence=''.join(values)
    primer=s.add('primer',primer_sequence,five='bead',features=features,source=src+' synthetic selected cell-barcode/UMI members')
    rna=s.add('rna',body+'A'*30,sugar='RNA',source=src+' illustrative RNA body and30-base tail window')
    s.supply_nt('DNA','ACGT',source=src)
    before=s.snapshot();bind=s.anneal('primer','rna',binding=30,source=src)
    bound=s.snapshot();pol=Polymerase('Maxima H− reverse transcriptase',('RNA',),'DNA')
    if s.extension_sites(pol,['primer','rna'])!=[s.tail('primer')]:raise ValueError('capture must supply exactly one eligible DNA3prime')
    extend=s.extend('primer',pol,source=src);after=s.snapshot()
    return dict(simulator=s,module=m,source=src,primer=primer,rna=rna,before=before,bound=bound,bind=bind,extend=extend,after=after,product=s.chain('primer'),regions=regions,body=body)


def dropseq():
    r=dropseq_simulation();src=r['source'];cell=11;left=190;y=500;pl=len(r['primer']);q=pl-30
    def tracks(state,product=False):
        p=strand('cdna',r['after']if product else state,r['product']if product else r['primer'],left,y,cell=cell,initial=pl if product else None)
        p.color='primer';p.regions=[]
        for a,b,label in r['regions']:
            names={'SMART handle':'PCR handle','cell barcode':'Cell barcode','UMI':'UMI','dT30':'oligo(dT)'}
            if label not in names:continue
            off={'SMART handle':(0,95),'cell barcode':(-20,-130),'UMI':(10,95),'dT30':(40,95)}[label]
            p.regions.append(Region(a,b,names[label],'bio'if label in {'cell barcode','UMI'}else'primer',1,label_offset=off))
        a=strand('rna',state,r['rna'],left+(pl+len(r['body']))*cell,y-24,cell=cell,direction=-1,regions=[Region(len(r['body']),len(r['rna']),'poly(A) tail','rna',-1,label_offset=(0,-70))])
        a.continued_3prime=True
        return [p,a]
    initial=tracks(r['before']);initial[1].start=(initial[1].start[0]+1800,initial[1].start[1]);initial[1].end=(initial[1].end[0]+1800,initial[1].end[1])
    bead=Decoration('bead',left-35,y,32,'Bead','primer')
    capture=Shot('capture','RNA binds the bead primer’s oligo(dT)',initial,
        [Event('anneal','rna',target=(-1800,0),fraction_end=.65),Event('pair','cdna',pairing='polya',fraction_start=.65,molecular_operation=r['bind'].to_dict())],
        'After droplet lysis, a polyadenylated RNA binds the oligo dT end of a bead primer. The full source SeqB primer is shown. Its cell barcode and UMI already exist before reverse transcription.',source=src,scope='One bead primer · synthetic barcode, UMI and RNA body',action=8,hold=2,decorations=[bead],pairings=[Pairing('cdna',q,'rna',len(r['rna'])-1,30,id='polya',initially_bound=False,exact=True,source=src)])
    generated(capture,r['before'])
    rt=Shot('rt','Extend the bead primer; inherit its identifiers',tracks(r['bound'],True),
        [Event('extend','cdna','Maxima H− RT',template='rna',template_start=len(r['body'])-1,primer_length=pl,primer_binding_length=30,exact=True,enzyme_label_offset=(-40,-100),molecular_operation=r['extend'].to_dict())],
        'Maxima H minus reverse transcriptase extends the existing bead-bound three-prime end along RNA. The barcode and UMI remain at the DNA five-prime end with their original nucleotide identities. Only the sequence copied from RNA is newly synthesized.',source=src,scope='Selected first-strand step · original RNA and bead retained',action=9,hold=3,transition=1.8,decorations=[bead])
    generated(rt,r['bound'])
    final=Shot('product','The cDNA retains the bead barcode and UMI',tracks(r['after'],True),
        [Event('highlight','cdna',span=(q-20,q))],
        'The same primer has become a longer cDNA. The cell barcode identifies this bead and the UMI identifies this captured molecule. The RNA is retained here. Template switching and subsequent amplification occur later and are outside this selected step.',source=src,scope='Retained bead–primer linkage · selected RNA–DNA hybrid',action=4,hold=3,transition=1.8,decorations=[bead])
    final.strands[0].initially_visible=len(r['product']);generated(final,r['after'])
    return evidence(qualify(r['module'],DROP,'Drop-seq',[capture,rt,final],{'capture':['partition_capture'],'rt':['reverse_transcribe'],'product':['reverse_transcribe']},limits=['Exact source SeqB primer scaffold with explicitly synthetic barcode/UMI members; no claim these values occurred in the study.','Illustrative RNA body and poly(A) window; droplet statistics, transcript occupancy and yields are not predicted.','Focused capture and first-strand copying only; template switching, exonuclease cleanup and downstream library construction are not silently omitted reactions within these shots.']),r['simulator'])


def scrrbs_simulation(seq='ACGCCGATCCGACGTCACGTACGA',protected=4):
    m=load(RRBS,'scrrbs.py');src='Guo2013 Genome Research doi:10.1101/gr.161679.113, Methods: Construction of single-cell RRBS sequencing libraries'
    s=Simulator();original=s.add('original',seq,modifications={protected:'5mC'},source=src+' synthetic genomic window and selected methylation occupancy')
    initial=s.snapshot();sites=[id for id in original if s.atoms[id].base=='C'and s.atoms[id].modification!='5mC']
    conversion=s.convert(sites,from_base='C',to_base='U',protected=('5mC',),source=src+' bisulfite plus desulfonation endpoint')
    converted=s.snapshot();current=sequence(converted,original)
    primer1=s.add('primer1',rc(current[-5:]),source=src+' synthetic genomic priming window; not source adapter primer')
    intake1=s.history[-1];primed1=s.snapshot()
    bind1=s.anneal('primer1','original',binding=5,template_site=len(original)-1,source=src+' illustrative selected binding site')
    bound1=s.snapshot();s.supply_nt('DNA','ACGT',source=src)
    pol=Polymerase('PfuTurbo Cx polymerase',('DNA',),'DNA')
    copy1=s.extend('primer1',pol,source=src+' selected first complementary copy; uracil-tolerant source polymerase')
    dna1=s.chain('primer1');after1=s.snapshot()
    melt=s.unpair(original,source=src+' illustrative next-cycle denaturation')
    melted=s.snapshot();primer2=s.add('primer2',rc(sequence(melted,dna1)[-5:]),source=src+' synthetic opposing priming window; not source adapter primer')
    intake2=s.history[-1];primed2=s.snapshot();bind2=s.anneal('primer2','primer1',binding=5,template_site=len(dna1)-1,source=src+' illustrative next-cycle binding site')
    bound2=s.snapshot();copy2=s.extend('primer2',pol,source=src+' subsequent complementary copying, not a claim about round2 timing')
    dna2=s.chain('primer2')
    return dict(simulator=s,module=m,source=src,original=original,initial=initial,sites=sites,conversion=conversion,converted=converted,primer1=primer1,intake1=intake1,primed1=primed1,bind1=bind1,bound1=bound1,copy1=copy1,dna1=dna1,after1=after1,melt=melt,melted=melted,primer2=primer2,intake2=intake2,primed2=primed2,bind2=bind2,bound2=bound2,copy2=copy2,dna2=dna2,after=s.snapshot(),protected=protected)


def scrrbs():
    r=scrrbs_simulation();src=r['source'];n=len(r['original']);cell=30;left=430;right=left+n*cell;y=450
    def original(state,y=y):
        a=strand('original',state,r['original'],left,y,cell=cell,label='Original DNA')
        a.label_position=(left+n*cell/2,y-70);a.marks[str(r['protected'])]='5mC'
        a.regions=[Region(r['protected'],r['protected']+1,'5mC retained','bio',-1,label_offset=(0,-130))]
        return a
    conversion=Shot('conversion','Bisulfite converts C to U; selected 5mC remains',
        [original(r['initial'])],[Event('convert','original',positions=[list(r['original']).index(i)for i in r['sites']],replacement='U',molecular_operation=r['conversion'].to_dict())],
        'This synthetic genomic window starts after source adapter ligation and denaturation. Bisulfite treatment and desulfonation convert susceptible cytosines into uracils in the original DNA. The selected five-methylcytosine remains unchanged. This is chemical modification, not copying.',source=src,scope='Prepared single-strand genomic window · illustrative methylation',action=8,hold=3)
    generated(conversion,r['initial'])
    p1=strand('copy1',r['after1'],r['dna1'],right+1800,y+24,cell=cell,direction=-1,initial=0)
    p1.color='new';p1.label='New complementary DNA';p1.label_position=(800+1800,y+100)
    introduce1=Shot('primer','Add and anneal a primer after bisulfite recovery',[original(r['converted']),p1],
        [Event('introduce','copy1',fraction_start=.05,molecular_operation=r['intake1'].to_dict()),Event('anneal','copy1',target=(-1800,0),fraction_start=.05,fraction_end=.7),Event('pair','copy1',pairing='copy1',fraction_start=.7,molecular_operation=r['bind1'].to_dict())],
        'After conversion and recovery, a stock DNA primer enters intact and anneals to the converted template. The synthetic five-base genomic priming window makes copying visible; it is not presented as the published adapter PCR primer.',source=src,scope='Illustrative genomic priming window · intact stock primer',action=6,hold=2,transition=1.8,pairings=[Pairing('copy1',0,'original',n-1,5,id='copy1',initially_bound=False,exact=True,source=src)])
    generated(introduce1,r['converted'])
    p1copy=strand('copy1',r['after1'],r['dna1'],right,y+24,cell=cell,direction=-1,initial=5,label='First complementary DNA')
    p1copy.color='new';p1copy.label_position=(800,y+110)
    copy1=Shot('copy','PfuTurbo Cx places A opposite the original U',[original(r['bound1']),p1copy],
        [Event('extend','copy1','PfuTurbo Cx',template='original',template_start=n-6,primer_length=5,primer_binding_length=5,exact=True,enzyme_label_offset=(0,180),molecular_operation=r['copy1'].to_dict())],
        'The source uracil-tolerant polymerase extends the existing DNA primer. It reads the physical U-bearing original and incorporates A opposite U. The original uracil is retained; it has not become thymine. Five-methylcytosine templates G.',source=src,scope='Chemical product U · complementary copy A',action=8,hold=3,transition=1.8)
    generated(copy1,r['bound1'])
    # Separate complete strands only after the physical duplex contacts melt.
    p1full=replace(p1copy,initially_visible=n)
    p2=strand('copy2',r['after'],r['dna2'],left-1800,y+24+100+24,cell=cell,initial=0)
    p2.color='primer';p2.label='Next complementary DNA';p2.label_position=(800-1800,y+24+100+120)
    nextcycle=Shot('next_cycle','Separate the strands and anneal the next primer',[original(r['after1']),p1full,p2],
        [Event('unpair','original',pairing='hybrid',fraction_end=.25,molecular_operation=r['melt'].to_dict()),Event('separate','copy1',target=(0,100),fraction_start=.25,fraction_end=.5),Event('introduce','copy2',fraction_start=.5,molecular_operation=r['intake2'].to_dict()),Event('anneal','copy2',target=(1800,0),fraction_start=.5,fraction_end=.8),Event('pair','copy2',pairing='copy2',fraction_start=.8,molecular_operation=r['bind2'].to_dict())],
        'A later thermal cycle separates the same original and daughter strands. A second intact stock primer enters and anneals to the new DNA. The uracil-bearing original remains separately visible, and no molecule is replaced by a readout label.',source=src,scope='Illustrative next cycle · original DNA retained',action=7,hold=2,transition=1.8,pairings=[Pairing('original',0,'copy1',n-1,n,id='hybrid',initially_bound=True,exact=True,source=src),Pairing('copy2',0,'copy1',n-1,5,id='copy2',initially_bound=False,exact=True,source=src)])
    generated(nextcycle,r['after1'])
    orig=original(r['bound2']);orig.label_position=(800,340)
    daughter=replace(p1full,start=(right,y+124),end=(left,y+124),label_position=(800,y+70))
    p2copy=strand('copy2',r['after'],r['dna2'],left,y+148,cell=cell,initial=5,label='T-bearing DNA copy');p2copy.color='primer';p2copy.label_position=(800,y+245)
    copy2=Shot('report','Copying A produces T in the next DNA strand',[orig,daughter,p2copy],
        [Event('extend','copy2','PfuTurbo Cx',template='copy1',template_start=n-6,primer_length=5,primer_binding_length=5,exact=True,enzyme_label_offset=(0,-210),molecular_operation=r['copy2'].to_dict())],
        'When the complementary A-bearing DNA is copied, polymerase incorporates T opposite A. This explains how bisulfite sequencing reports C to T while the original chemical change was C to U. The new T, the earlier A and the original U belong to three distinct retained molecules.',source=src,scope='Original U retained · first copy A · next copy T',action=8,hold=4,transition=1.8)
    generated(copy2,r['bound2'])
    return evidence(qualify(r['module'],RRBS,'scRRBS',[conversion,introduce1,copy1,nextcycle,copy2],{'conversion':['bisulfite_conversion'],'primer':['column_recovery','library_pcr'],'copy':['library_pcr'],'next_cycle':['library_pcr'],'report':['library_pcr']},limits=['Selected synthetic genomic window and methylation occupancy, not a sample-derived molecule or methylation estimate.','Starts after MspI digestion, end preparation, source premethylated-adapter ligation and denaturation; these preparative reactions are not reconstructed here.','Five-base synthetic genomic primers illustrate complementary copying; they are not source PCR oligos and do not establish native adapter geometry.','PfuTurbo Cx is used for both illustrated first-round PCR cycles. This is not a claim that the later Phusion round begins at the second strand shown.']),r['simulator'])
