"""Source-qualified adapters for selected protocol mechanisms, not whole-population simulations.

Exact values come from existing protocol constructors. Geometric fixtures, hidden
population members, sugar uncertainty and selected-success scope remain explicit.
"""
from __future__ import annotations
from dataclasses import replace
import hashlib
import importlib.util
import math
from pathlib import Path
import re
import sys
from .model import Region, Decoration, Event, Junction, Movie, Pairing, Shot, Strand, complement

ROOT=Path(__file__).resolve().parents[1]


def load(slug,filename):
    path=ROOT/slug/'tools'/filename
    original_path=sys.path[:]
    sys.path[:0]=[str(path.parent),str(ROOT/'lib'),str(ROOT)]
    name='animation_protocol_'+hashlib.sha256(str(path).encode()).hexdigest()[:12]
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    try:spec.loader.exec_module(module)
    finally:sys.path[:]=original_path
    return module


def arc(n,cx=820,cy=480,r=135,start=math.pi/2):
    return [(n*i/48,cx+r*math.cos(start+math.tau*i/48),cy+r*math.sin(start+math.tau*i/48))for i in range(49)]


def straight(n,x0,y0,x1,y1):return [(0,x0,y0),(n,x1,y1)]


def qualify(module,slug,title,shots,steps,*,limits=(),extra=None):
    old=sys.path[:];sys.path[:0]=[str(ROOT/slug/'tools'),str(ROOT/'lib'),str(ROOT)]
    try:g=module.protocol_graph()
    finally:sys.path[:]=old
    for shot in shots:
        if shot.id not in steps or not all(k in g.steps for k in steps[shot.id]):
            raise ValueError(f'{slug}/{shot.id}: animation has no corresponding declared protocol step')
        if not shot.source:shot.source=re.sub('<[^>]*>','',g.authority)+'; selected steps: '+', '.join(steps[shot.id])
    identities={}
    for shot in shots:
        identities.update(shot.inventory)
        for s in shot.strands:
            if s.lineage:
                for _,_,origin,_ in s.lineage:
                    if origin not in shot.inventory and origin in identities:shot.inventory[origin]=identities[origin]
                continue
            if any(i!=0 for i in s.nicks):continue
            initial=len(s.sequence)if s.initially_visible is None else s.initially_visible
            if initial:
                origin=shot.id+'_'+s.id+'_input'
                shot.inventory[origin]=s.sequence[:initial]
                s.lineage=[(0,initial,origin,0)]
        identities.update(shot.inventory)
    path=Path(module.__file__)
    provenance={'representation':'source-qualified selected mechanism',
                'protocol_module':str(path.relative_to(ROOT)),'module_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                'authority':g.authority,'steps_used':steps,
                'limits':['Selected successful event, not yield or kinetics. Motion and layout are illustrative.',
                          'Prepared-source inputs and workups outside these focused scenes are not reconstructed.',*limits]}
    if extra:provenance.update(extra)
    return Movie(title,shots,slug,provenance).validate()


def closure_event(strand,source,enzyme='Ligase',start=.5):
    return Event('seal',strand,enzyme,bond=0,junction=Junction(0,source=source,mode='splinted'),fraction_start=start)


def crispr_mip():
    slug='crispr-mip__10.1101+2024.03.28.587082';m=load(slug,'crisprmip.py')
    probe=m.probe_construct().top();fill=m.capture_fill();full=probe+fill
    e,l=len(m.EXT_ARM),len(m.LIG_ARM);n=len(probe);gap=len(fill)
    target=complement((m.EXT_ARM+fill+m.LIG_ARM)[::-1])
    # Actual nucleotide intervals determine arm/gap positions. Backbone is a routed loop.
    target_start=350;scale=800/(e+gap+l)
    ext_end=target_start+e*scale;lig_start=ext_end+gap*scale
    route=[(0,lig_start,495),(l,1150,495),(l+20,1240,380),
           (n-e-20,270,380),(n-e,350,495),(n,ext_end,495),
           (len(full),lig_start,495)]
    probe_route=route[:-1]
    probe_strand=Strand('probe',probe,(350,370),(1150,370),'primer','One MIP probe',five_prime='P',
                        path=straight(n,350,360,1150,360),label_position=(800,300))
    probe_strand.regions=[Region(0,l,'Ligation arm · 5′ P','primer',1),Region(n-e,n,'Extension arm · 3′ OH','bio',1)]
    template=Strand('target',target,(1150,545),(350,545),'dna','Guide-cassette target',label_position=(800,630))
    pairing=[Pairing('probe',0,'target',l-1,l,source='Model capture ligation arm'),Pairing('probe',n-e,'target',len(target)-1,e,source='Model capture extension arm')]
    fold=Shot('arms','Wrap around the target',[probe_strand,template],
              [Event('fold','probe',route=probe_route)],
              'One continuous probe binds the two flanking arms. Its backbone loops away from the target, leaving the guide cassette between the extendable end and the phosphorylated end.',
              scope='Illustrative cassette · published probe arms',action=6,pairings=pairing)
    extended=Strand('probe',full,(0,0),(0,0),'primer','Probe + copied gap',five_prime='P',initially_visible=n,
                    path=route,nicks=[0],bands=[(0,n,'primer')],label_position=(800,300))
    extended.regions=[Region(0,l,'5′ P · ligation arm','primer',1),Region(n-e,n,'3′ extension arm','bio',1),Region(n,len(full),'Copied target gap','new',-1,(800,440))]
    gapshot=Shot('gap','Fill the gap',[extended,template],
                [Event('extend','probe','Phusion',template='target',template_start=gap+l-1,primer_length=n,primer_binding_length=e,exact=True)],
                'Phusion extends the probe three-prime end across the guide cassette. The original target remains intact. The extended probe is still linear until its final nick is sealed.',
                scope='Selected capture · gap copied from model target',action=6,pairings=pairing)
    filled=replace(extended,initially_visible=None)
    seal=Shot('closure','Seal the circle',[filled,template],
              [closure_event('probe','Source probe5P/3OH; modeled_capture and filled_probe_construct.','Ampligase',0)],
              'Ampligase joins the filled probe three-prime hydroxyl to its original five-prime phosphate. Closing this last nick creates the selected covalent capture circle.',
              scope='Selected successful closure · target retained',action=4,pairings=pairing)
    # Check copying against the actual supplied template independently of the long
    # primer backbone, which is intentionally not complementary to the gap template.
    assert fill==complement(target[l:l+gap][::-1])
    c=m.modeled_capture()
    assert full[n-e:]+full[:n-e]==c.circle.seq
    return qualify(m,slug,'CRISPR-MIP',[fold,gapshot,seal],{'arms':['hybridize'],'gap':['extend_gap'],'closure':['seal_circle']},
                   limits=['Guide cassette and leading G convention are illustrative, not every screened clone.',
                           'Source probe126nt versus Methods134bp discrepancy retained in protocol notes; TableS2 constructor governs probe.'],
                   extra={'probe_sequence_5to3':probe,'copied_gap_5to3':fill,'circle_rotation_checked':True})


def boloramis():
    from .boloramis import movie
    return movie()


def atac():
    slug='atac-seq__10.1038+nmeth.2688';m=load(slug,'atac_seq.py')
    # Existing constructor knows end components; native target bases remain unknown.
    import nextera
    tag=tagmentation_shot()
    template=Strand('template','N'*40,(1270,440),(330,440),'dna','Tagged template',schematic=True)
    copy=Strand('fill','N'*40,(330,495),(1270,495),'adapter','Extended DNA',initially_visible=30,schematic=True)
    extend=Shot('fill','Extend the tagged ends',[template,copy],
                [Event('extend','fill','DNA polymerase',template='template',template_start=9,primer_length=30)],
                'After transposition and the selected workup, polymerase extends the tagged ends. The source then uses limited indexed PCR to make the productive library.',
                scope='Source extension boundary · whole early products unresolved',action=4)
    return qualify(m,slug,'ATAC-seq',[tag,extend],{'transfer':['tagment'],'fill':['purify_tagmented','extend','indexed_pcr']},
                   limits=['Staggered cuts are qualitative; no unvalidated genomic cut-coordinate claim.',
                           'Canonical mosaic-end sequence imported; exact native transposome loading and strand-transfer intermediates are not simulated.'])


def smallseq():
    slug='small-seq__10.1038+nbt.3701';m=load(slug,'smallseq.py')
    insert='N'*22;ra3=m.RA3+'C';first=insert+ra3;n=len(first)
    track=Strand('rna_adapter',first,(330,495),(1270,495),'rna','RNA + DNA adapter',polymer='mixed',three_prime='ddC',
                 nicks=[len(insert)],schematic=True,bands=[(len(insert),n,'adapter')],sugars=[(0,len(insert),'RNA'),(len(insert),n,'DNA')])
    track.pieces=[(len(insert),n,150,0)]
    track.regions=[Region(0,len(insert),'Small RNA · 3′ OH','rna',-1,(530,380)),Region(len(insert),n,'5′ App adapter · 3′ ddC','adapter',1,(1090,630))]
    lig=Shot('lig3','Add the 3′ adapter',[track],
             [Event('displace','rna_adapter',span=(len(insert),n),target=(-150,0),fraction_end=.55),Event('seal','rna_adapter','T4 RNA ligase 2 KQ',bond=len(insert),fraction_start=.6,junction=Junction(len(insert),donor_5prime='App',mode='single-strand',source=m.SOURCE+' RA3 preadenylated; no ATP.'))],
             'The truncated ligase joins an accessible RNA three-prime hydroxyl to the preadenylated DNA adapter. The adapter terminal ddC blocks extension. No ATP is added in this ligation.',
             scope='Selected RNA3′OH–5′App adapter join · unknown insert',action=4)
    ra5=m.ra5_oligo().replace('T','U');whole=ra5+insert+ra3
    track2=Strand('both',whole,(330,495),(1270,495),'rna','Doubly ligated small RNA',polymer='mixed',five_prime='NH2',three_prime='ddC',nicks=[len(ra5)],schematic=True,label_position=(800,625),
                  bands=[(0,len(ra5),'primer'),(len(ra5)+len(insert),len(whole),'adapter')],
                  sugars=[(0,len(ra5)+len(insert),'RNA'),(len(ra5)+len(insert),len(whole),'DNA')])
    track.lineage=[(0,len(insert),'rna_input',0),(len(insert),n,'ra3_input',0)]
    lig.inventory={'rna_input':insert,'ra3_input':ra3}
    track2.lineage=[(0,len(ra5),'ra5_input',0),(len(ra5),len(whole),'rna_ra3_retained',0)]
    track2.pieces=[(0,len(ra5),-120,0)]
    track2.regions=[Region(0,len(ra5)-m.RA5_UMI_LEN,'5′ adapter','primer',-1,(430,360)),Region(len(ra5)-m.RA5_UMI_LEN,len(ra5),'UMI','bio',1,(650,680)),Region(len(ra5),len(ra5)+len(insert),'Small RNA','rna',-1,(950,360)),Region(len(ra5)+len(insert),len(whole),'3′ adapter · ddC','adapter',1,(1160,680))]
    lig5=Shot('lig5','Add the 5′ adapter',[track2],
              [Event('displace','both',span=(0,len(ra5)),target=(120,0),fraction_end=.55),Event('seal','both','T4 RNA ligase 1',bond=len(ra5),fraction_start=.6,junction=Junction(len(ra5),mode='single-strand',source=m.SOURCE+' RA5RNA3OH–RNA5P, ATP.'))],
              'After free three-prime adapter processing, the second ligase joins the five-prime adapter to the RNA. Its eight-base molecular identifier is carried by the adapter, not added by polymerase.',
              scope='Selected RA5–RNA join · free-adapter workup separate',action=4)
    lig5.inventory={'ra5_input':ra5,'rna_ra3_retained':first}
    mask=m.MASK_58S
    conditional=m.MASK_58S_TARGET.replace('T','U')
    masking=Shot('mask','Occlude the rRNA end',[
        Strand('conditional_rrna',conditional,(1270,440),(330,440),'rna','Conditional rRNA target','RNA',label_position=(800,345)),
        Strand('mask',mask,(330,650),(1270,650),'primer','3′-biotinylated mask',three_prime='biotin',label_position=(800,755))],
        [Event('anneal','mask',target=(0,-160))],
        'The masking oligo occludes the ribosomal RNA end rather than pulling the RNA out of the sample. The complementary target shown here is conditional; its reference coordinates remain unverified in the model.',
        scope='INFERRED complementary target · source mask ordering',action=4,
        pairings=[Pairing('mask',0,'conditional_rrna',len(mask)-1,len(mask),source='Source mask and conditional reverse-complement query')])
    free=Strand('free_adapter',ra3,(330,440),(1270,440),'adapter','Free RA3',five_prime='P',three_prime='ddC')
    rtp=Strand('rtp',m.RTP,(1270,490),(330+940/len(ra3),490),'primer','Protected RTP',polymer='mixed',five_prime='biotin',
               sugars=[(0,len(m.RTP)-1,'DNA'),(len(m.RTP)-1,len(m.RTP),'RNA')],label_position=(800,680))
    remove=Shot('free_adapter','Digest the unligated adapter',[free,rtp],
        [Event('resect','free_adapter','Lambda exonuclease',amount=len(ra3)-1,note='RA3 already deadenylated to5′P; RTP5′ biotin protects its end')],
        'Free RA3 is deadenylated and paired to RTP. Lambda exonuclease then degrades the phosphorylated adapter strand while RTP is protected by its five-prime biotin. This acts on free adapter, not the successfully ligated RNA.',
        scope='Free-adapter branch · degradation products unexpanded',action=5,
        pairings=[Pairing('free_adapter',0,'rtp',len(m.RTP)-1,len(m.RTP),source='Existing free_ra3_rtp_duplex constructor')])
    return qualify(m,slug,'Small-seq',[masking,lig,remove,lig5],{'mask':['mask_rrna'],'lig3':['ligate_3p'],'free_adapter':['deadenylate','digest_free_adapter'],'lig5':['ligate_5p']},
                   limits=['Unknown small-RNA sequence/length is represented by an illustrative22nt class, not a sample value.',
                           'Masking and adapter digestion are separate source steps, not poly(A) tailing.',
                           'The extra RA3 terminal ddC is physical; not an omitted nonnucleotide label.'],extra={'ra5_ordering':m.RA5_MODS,'ra3_ordering':m.RA3_MODS,'umi_length':m.RA5_UMI_LEN})


def gibson():
    slug='gibson-assembly__10.1038+nmeth.1318';m=load(slug,'gibson_assembly.py')
    parent,topspan,bottomspan=m.precursor(m.example_fragments()[0]);seq=parent.top();n=len(seq)
    top=Strand('top',seq,(350,435),(1250,435),'dna','Toy duplex',schematic=False)
    bottom=Strand('bottom',complement(seq[::-1]),(1250,485),(350,485),'primer','',schematic=False)
    resect=Shot('resect','Expose complementary overlaps',[top,bottom],
                [Event('resect','top','5′ exonuclease',amount=topspan[0]),Event('resect','bottom','5′ exonuclease',amount=bottomspan[0])],
                'Five-prime exonuclease recession exposes complementary three-prime overlaps. Removed bases become degradation products. These exact stopping points are invented for the repository toy example.',
                scope='INFERRED toy bases and resection boundaries',action=5)
    before=m.A+m.RIGHT;after=complement((m.RIGHT+m.B)[::-1])
    overlap=Shot('anneal','Anneal the exposed overlaps',[
        Strand('upstream',before,(350,440),(830,440),'dna','Part A',label_position=(510,335)),
        Strand('downstream',after,(1130,650),(650,650),'primer','Part B',label_position=(1100,725))],
        [Event('anneal','downstream',target=(0,-160))],
        'The exposed three-prime ends are complementary. They anneal without incorporating the overlap twice. The source mixture contains polymerase, but this exact-end toy already has no gap to fill.',
        scope='INFERRED toy overlap · exact supplied complementarity',action=4,
        pairings=[Pairing('upstream',len(m.A),'downstream',len(after)-1,len(m.RIGHT),source=m.EXAMPLE)])
    overlap.strands[0].regions=[Region(len(m.A),len(before),'Complementary overlap','bio',-1,(850,340))]
    overlap.strands[1].regions=[Region(len(after)-len(m.RIGHT),len(after),'Complementary overlap','bio',1,(850,720))]
    overlap.events.extend([Event('highlight','upstream',span=(len(m.A),len(before))),Event('highlight','downstream',span=(len(after)-len(m.RIGHT),len(after)))])
    # The existing exact-end toy has no gaps to fill. Do not manufacture synthesis.
    geom=m.protocol_graph().materials['joint'].chemical_state
    seqout=geom.top() if hasattr(geom,'top')else ''.join(x.top for x in geom)
    seam=len(m.LEFT)+len(m.A)
    joined=Strand('joined',seqout,(350,495),(1250,495),'dna','Paired toy assembly',nicks=[seam],label_position=(800,610))
    seal=Shot('seal','Seal the remaining nick',[joined],
              [Event('seal','joined','DNA ligase',bond=seam,junction=Junction(seam,source=m.SOURCE+' illustrative source-compatible phosphorylated nick; no gap in exact-end toy.'))],
              'Compatible overlaps anneal. Polymerase repair is part of the source method, but this exact-end toy has no gap to copy. Ligase seals the illustrated remaining nick without inventing added nucleotides.',
              scope='Illustrative compatible nick · zero-gap toy',action=4)
    return qualify(m,slug,'Gibson assembly',[resect,overlap,seal],{'resect':['a_process'],'anneal':['anneal'],'seal':['anneal','repair','seal']},
                   limits=[m.EXAMPLE,'One nick shown at magnification; not a complete constructed-vector movie.',
                           'Other protocol examples exercise real gap filling; none is invented in this zero-gap toy.'])


def hic():
    slug='in-situ-hi-c__10.1016+j.cell.2014.11.021';m=load(slug,'insitu_hic.py')
    seq='N'*14+m.MBOI.site+'N'*14;n=len(seq)
    a=Strand('a',seq,(330,435),(1270,435),'dna','Chromatin DNA',schematic=True)
    b=Strand('b',complement(seq[::-1]),(1270,485),(330,485),'primer','',schematic=True)
    cut=Shot('cut','Create staggered ends',[a,b],
             [Event('cleave','a','MboI',bond=14+m.MBOI.top_cut),Event('cleave','b',bond=n-14-m.MBOI.bottom_cut)],
             'MboI cleaves its GATC site on both strands, leaving cohesive ends. The original genomic flanks are represented as unknown regions.',
             scope='Source cut coordinates · illustrative genomic flanks',action=4)
    ca=14+m.MBOI.top_cut;cb=n-14-m.MBOI.bottom_cut
    cut.events.extend([Event('displace','a',span=(0,ca),target=(-100,0),fraction_start=.6),Event('displace','a',span=(ca,n),target=(100,0),fraction_start=.6),Event('displace','b',span=(0,cb),target=(100,0),fraction_start=.6),Event('displace','b',span=(cb,n),target=(-100,0),fraction_start=.6)])
    cut.events[0].fraction_end=.55;cut.events[1].fraction_end=.55
    a.regions=[Region(14,18,'MboI recognition site · GATC','bio',-1,(800,330))]
    tail=m.FILLED.left_new_top
    seq='N'*10+tail;n=len(seq)
    filling=Strand('filled',seq,(400,495),(1200,495),'dna','Filled end',initially_visible=10,schematic=True,marks={str(10+i):'biotin'for i in m.FILLED.biotin_left_top})
    tmpl=Strand('template',complement(seq[::-1]),(1200,445),(400,445),'primer','Template',schematic=True)
    fill=Shot('fill','Fill and label the ends',[tmpl,filling],
              [Event('extend','filled','Klenow',template='template',template_start=len(tail)-1,primer_length=10,note='Biotin-dATP incorporated during fill-in')],
              'Klenow fills the recessed ends while incorporating biotin-bearing nucleotide. Biotin marks the new material near the future contact junction; it is not a change to the base letter reported by sequencing.',
              scope='Source-derived fill sequence · selected labeled end',action=4)
    j=m.FILLED.junction;joined='N'*12+j+'N'*12;nick=12+m.MBOI.bottom_cut
    top=Strand('contact',joined,(360,495),(1240,495),'dna','Locus A — junction — locus B',schematic=True,nicks=[nick],bands=[(nick,len(joined),'primer')],marks={str(12+m.FILLED.junction_biotin_top[0]):'biotin'},label_position=(800,625))
    top.regions=[Region(0,12,'Locus A','dna',-1,(500,360)),Region(12,12+len(j),'Contact junction','bio',1,(800,700)),Region(12+len(j),len(joined),'Locus B','primer',-1,(1100,360))]
    top.label=''
    lig=Shot('contact','Join proximal DNA ends',[top],
             [Event('seal','contact','T4 DNA ligase',bond=nick,junction=Junction(nick,source=m.CITATION+' source filled-end proximity ligation.'))],
             'Ligation joins two filled ends that were close in the fixed nucleus. This selected junction records proximity. It does not mean the genomic loci were adjacent in the original sequence.',
             scope='Selected proximity junction · no contact probability claim',action=4)
    cap=Shot('capture','Recover biotin-bearing fragments',[
        replace(top,nicks=[],start=(360,375),end=(1240,375),regions=[r for r in top.regions if r.label!='Contact junction'])],
        [Event('capture','contact',target=(0,160))],
        'After crosslink reversal and fragmentation, streptavidin enriches biotin-bearing fragments. Binding and the subsequent workup are distinct from covalent ligation; recovery is not assumed complete.',
        scope='Source physical enrichment · schematic bead/fragment',action=4,decorations=[Decoration('bead',800,595,62,'Streptavidin','bio')])
    return qualify(m,slug,'in situ Hi-C',[cut,fill,lig,cap],{'cut':['digest'],'fill':['fill'],'contact':['contact_ligation'],'capture':['reverse_crosslinks','fragment','capture']},
                   limits=['One strand/end shown at magnification; complete population, opposite-strand seals and yields are not expanded.',m.LIBRARY_CAVEAT])


def dropseq():
    slug='drop-seq__10.1016+j.cell.2015.05.002';m=load(slug,'drop_seq.py')
    bead=m.bead_primer();data=[{'name':s.name,'sequence':s.top,'role':getattr(s.feature,'role',None)}for s in bead]
    dt=next(s.top for s in bead if s.name=='dT30')
    # Zoomed region is actual dT tract; barcode values remain unspecified upstream.
    template=Strand('rna','N'*20+'A'*len(dt),(1050,385),(755,385),'rna','RNA','RNA',schematic=True,label_position=(1050,305))
    primer=Strand('primer',dt,(755,565),(932,565),'primer','Bead-bound oligo(dT)',five_prime='bead',label_position=(930,685))
    capture=Shot('capture','Capture RNA on one bead',[template,primer],
                 [Event('anneal','rna',target=(0,130))],
                 'A droplet brings a cell and a barcoded bead together. Released polyadenylated RNA binds oligo dT primers tethered to the bead. The cell barcode and UMI are already present on those primers.',
                 scope='Schematic droplet · zoomed dT binding region',action=5,
                 decorations=[Decoration('droplet',850,490,245,'','new'),Decoration('bead',700,565,55,'Bead','primer'),Decoration('cell',1000,435,48,'','rna')])
    parts=list(bead);prefix=''.join('N'*len(s.top) if s.placeholder else s.top for s in parts[:-1])
    primer_seq=prefix+dt;pl=len(primer_seq);q=len(prefix)
    xpos=lambda i:300+8*i
    regions=[];cursor=0
    for seg in parts:
        color='bio' if seg.name in {'cell barcode','UMI'} else 'primer'
        name={'SMART handle':'PCR handle','dT30':'oligo(dT) · TTTTT…'}.get(seg.name,seg.name)
        if seg.name not in {'bead linker','bead constant'}:
            placements={'SMART handle':(440,390),'cell barcode':(665,700),'UMI':(790,390),'dT30':(960,700)}
            regions.append(Region(cursor,cursor+len(seg.top),name,color,1,placements[seg.name]))
        cursor+=len(seg.top)
    fullprimer=Strand('cdna',primer_seq+'N'*20,(300,565),(xpos(pl+20),565),'dna','',five_prime='bead',initially_visible=pl,schematic=True,regions=regions)
    fullrna=Strand('rna','N'*20+'A'*30,(xpos(pl+20),515),(xpos(q),515),'rna','',polymer='RNA',schematic=True,
                   regions=[Region(20,50,'poly(A) tail','rna',-1,(990,325))])
    pairing=[Pairing('cdna',q,'rna',49,30,exact=False,source='Source bead oligo dT30; illustrative RNA body.')]
    anatomy=Shot('primer','Zoom into the captured bead primer',[fullrna,fullprimer],
                 [Event('highlight','cdna',span=(0,pl))],
                 'The bead primer already contains a PCR handle, a cell barcode, a UMI and oligo dT. Only oligo dT binds the RNA poly A tail. The other regions remain as the five-prime part of the primer.',
                 scope='Full primer anatomy · unspecified barcode / UMI values',action=5,pairings=pairing,
                 decorations=[Decoration('bead',250,565,50,'Bead','primer')])
    rtshot=Shot('rt','Copy RNA; retain the barcode and UMI',[fullrna,fullprimer],
        [Event('extend','cdna','Reverse transcriptase',template='rna',template_start=19,primer_length=pl,primer_binding_length=30),Event('highlight','cdna',span=(q-20,q))],
        'Reverse transcriptase extends the existing bead primer. The cell barcode and UMI remain attached to the growing cDNA; they are inherited from the primer, rather than created by polymerase.',
        scope='Selected tethered copy · barcode / UMI values unspecified',action=5,pairings=pairing,decorations=[Decoration('bead',250,565,50,'Bead','primer')])
    return qualify(m,slug,'Drop-seq',[capture,anatomy,rtshot],{'capture':['partition_capture'],'primer':['partition_capture'],'rt':['reverse_transcribe']},
                   limits=['Droplet/cell layout is illustrative; occupancy statistics and yields are not simulated.',
                           'Partition view zooms into oligo(dT); subsequent full-primer view retains the upstream barcode/UMI handles.',
                           'RNA body and the depicted target count are schematic, not actual sample sequence.'],extra={'source_bead_primer_segments':data})


def scrrbs():
    slug='scrrbs__10.1101+gr.161679.113';m=load(slug,'scrrbs.py')
    # Explicit synthetic sequence to explain the chemical/readout distinction;
    # source mouse/lambda methylation occupancy is not manufactured.
    seq='ACGCCGATCCGACGTCACGTACGA';protected=3
    track=Strand('original',seq,(350,470),(1250,470),'dna','Original DNA',marks={str(protected):'5mC',str(1):'C → U'},label_position=(800,630))
    sites=[i for i,b in enumerate(seq)if b=='C'and i!=protected]
    fillseq='N'*10+'CGA'
    endprep=Shot('ends','Fill the recessed end, then add A',[
        Strand('end',fillseq,(350,495),(1250,495),'dna','Prepared end',initially_visible=10,schematic=True,label_position=(800,650)),
        Strand('overhang_template',complement(fillseq[:-1][::-1]),(1181,445),(350,445),'primer','Template',schematic=True,label_position=(800,340))],
        [Event('extend','end','Klenow exo−',template='overhang_template',template_start=1,primer_length=10,copy_end=12,fraction_end=.65),
         Event('tail','end','Klenow exo−',primer_length=12,fraction_start=.65,note='Source extra3′A · no template for terminal addition')],
        'Klenow fills the recessed MspI end and adds the source extra three-prime A. The final A is a terminal addition, distinct from copying the overhang. Premethylated adapter ligation follows in the source route.',
        scope='Source CG fill/A-tail · illustrative genomic flank',action=5)
    chem=Shot('conversion','Convert selected unmethylated C sites',[track],
              [Event('convert','original',positions=sites,replacement='U',note='Bisulfite + desulfonation; chemical treatment, not an enzyme')],
              'Bisulfite treatment followed by desulfonation converts susceptible cytosine to uracil in the original DNA. The marked methylated cytosine is retained in this illustrative example. This is chemical conversion, not DNA copying.',
              scope='Illustrative sequence and occupancy · chemical DNA U',action=6)
    converted=''.join('U'if i in sites else b for i,b in enumerate(seq));daughter=complement(converted[::-1])
    # Copy the U-bearing original, not a preemptively C→T-rewritten template.
    template=Strand('converted',converted,(1250,440),(350,440),'dna','U-bearing original',marks={str(protected):'5mC'},label_position=(800,340))
    copy=Strand('copy',daughter,(350,495),(1250,495),'primer','New complementary DNA',initially_visible=5,label_position=(800,640))
    pcr=Shot('copy','Copy the converted template',[template,copy],
             [Event('extend','copy','Uracil-tolerant polymerase',template='converted',template_start=len(seq)-6,primer_length=5,exact=True)],
             'A uracil-tolerant polymerase copies the converted original. It places A opposite U. A later complementary copy can report T at that position; the chemical product and the sequencing readout are different claims.',
             scope='Illustrative copied strand · U→A complementary relation',action=5)
    return qualify(m,slug,'scRRBS',[endprep,chem,pcr],{'ends':['end_prepare','adapter_ligation'],'conversion':['bisulfite_conversion','column_recovery'],'copy':['library_pcr']},
                   limits=['Illustrative bases and selected methylation state; no quantitative conversion or native adapter occupancy predicted.',m.CAVEAT],extra={'illustrative_sequence':seq,'selected_protected_c':protected})


def celseq():
    slug='cel-seq__10.1016+j.celrep.2012.08.003';m=load(slug,'native.py')
    promoter=next(s.top for s in m.rt_parts()if s.name=='T7 promoter')
    body=40;p=len(promoter);n=p+body
    seq=promoter+'N'*body
    non=Strand('nontemplate',seq,(350,475),(1250,475),'dna','',schematic=True,
               path=[(0,350,475),(p,550,475),(n,1250,475)],regions=[Region(0,p,'T7 promoter','bio',-1,(440,310)),Region(p,n,'cDNA insert','dna',-1,(950,330))])
    template=Strand('template',complement(seq[::-1]),(1250,430),(350,430),'dna','',schematic=True,
                    path=[(0,1250,430),(body,550,430),(n,350,430)],bands=[(body,n,'bio')])
    rna=Strand('rna','N'*body,(550,540),(1250,540),'rna','New antisense RNA','RNA',initially_visible=0,schematic=True,label_position=(930,660))
    first=Shot('ivt','Transcribe the DNA template',[template,non,rna],
               [Event('transcribe','rna','T7 RNA polymerase',template='template',template_start=body-1,primer_length=0,note='rATP · rCTP · rGTP · rUTP')],
               'A paired T7 promoter recruits RNA polymerase to the double-stranded cDNA. Transcription makes antisense RNA with ribonucleotides. The DNA template is retained rather than converted into RNA.',
               scope='Source T7 amplification · transcript geometry schematic',action=5,pairings=[Pairing('nontemplate',0,'template',n-1,p,exact=False,source=m.AUTHORITY if hasattr(m,'AUTHORITY') else promoter)])
    first.events.append(Event('highlight','nontemplate',span=(0,p)))
    old=replace(rna,id='older_rna',initially_visible=None,start=(550,625),end=(1250,625),label='Earlier transcript',label_position=(800,750))
    second=Shot('repeat','Reuse the same DNA template',[template,non,old,replace(rna,label='Another RNA copy')],
                [Event('transcribe','rna','T7 RNA polymerase',template='template',template_start=body-1,primer_length=0)],
                'The same promoter-bearing DNA can initiate another transcript while earlier RNA copies remain. This is linear amplification from a reusable template, rather than PCR doubling through repeated strand separation.',
                scope='Illustrative successive initiations · no rate/yield claim',action=5)
    return qualify(m,slug,'CEL-Seq',[first,second],{'ivt':['second_strand','ivt'],'repeat':['ivt']},
                   limits=[m.CAVEAT,'Full source transcript lengths and promoter initiation complex are unassembled; region-level motion only.'],extra={'source_t7_promoter':promoter,'source_rt_primer_segments':[{'name':x.name,'sequence':x.top}for x in m.rt_parts()]})


def smartseq():
    from .smartseq import movie
    return movie()


def tagmentation_shot(id='transfer'):
    """Two loaded homodimers bracket one selected fragment; stagger graphic only."""
    from .model import Region
    top=Strand('top','N'*140,(220,465),(1340,465),'dna','',schematic=True, regions=[Region(30,105,'Selected cDNA / DNA fragment','dna',1,(800,670))])
    bottom=Strand('bottom','N'*140,(1340,525),(220,525),'dna','',schematic=True)
    left=Strand('adapter_left','N'*15,(340,300),(460,300),'adapter','',schematic=True,
                regions=[Region(0,15,'Loaded adapter DNA','adapter',-1,(390,240))])
    right=Strand('adapter_right','N'*15,(1252,380),(1132,380),'adapter','',schematic=True,
                 regions=[Region(0,15,'Loaded adapter DNA','adapter',-1,(1200,240))])
    outside_left=Strand('outside_left','N'*15,(652,360),(532,360),'adapter','',schematic=True)
    outside_right=Strand('outside_right','N'*15,(940,320),(1060,320),'adapter','',schematic=True)
    events=[Event('transfer','outside_left',target=(0,165),fraction_end=.55),
            Event('transfer','outside_right',target=(0,145),fraction_end=.55),
            Event('separate','outside_left',target=(-200,0),fraction_start=.8),
            Event('separate','outside_right',target=(200,0),fraction_start=.8),
            Event('transfer','adapter_left','Tn5 homodimer',target=(0,165),fraction_end=.55,persistent=True,complex='homodimer',actor_offset=(36,30)),
            Event('transfer','adapter_right','Tn5 homodimer',target=(0,145),fraction_end=.55,persistent=True,complex='homodimer',actor_offset=(-36,-30))]
    events.extend([Event('release','adapter_left','Tn5 homodimer',complex='homodimer',actor_offset=(36,30),fraction_start=.65,fraction_end=.8),Event('release','adapter_right','Tn5 homodimer',complex='homodimer',actor_offset=(-36,-30),fraction_start=.65,fraction_end=.8)])
    for strand,bonds,intervals in [('top',[30,105],[(0,30,-200),(105,140,200)]),('bottom',[26,101],[(0,26,200),(101,140,-200)])]:
        events.extend(Event('cleave',strand,bond=b,fraction_start=.4,fraction_end=.65)for b in bonds)
        events.extend(Event('displace',strand,span=(a,b),target=(dx,0),fraction_start=.8)for a,b,dx in intervals)
    return Shot(id,'Tn5 dimers cut and transfer adapter DNA',[top,bottom,left,right,outside_left,outside_right],events,
        'Each Tn5 complex has two identical subunits and carries adapter DNA. Two insertion events bracket a selected fragment. Staggered cleavage and strand transfer attach adapters; after complex release, the cut flanking DNA pieces slide apart to reveal the fragment.',
        action=8,hold=3,scope='Two loaded homodimers · selected fragment, qualitative stagger')
