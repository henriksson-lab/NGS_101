"""Selected source-qualified capture and ligation mechanisms, driven by substrates."""
from dataclasses import replace
import math
from .model import Event,Shot,Movie,Region,Pairing,Junction
from .protocols import load,qualify
from .simulator import Simulator,Ligase
from .simulation_adapters import generated,evidence
from .generated_views import strand
from .other_simulations import mip as simulate_mip,hic as simulate_hic

MIP='crispr-mip__10.1101+2024.03.28.587082'
SMALL='small-seq__10.1038+nbt.3701'
HIC='in-situ-hi-c__10.1016+j.cell.2014.11.021'


def seal_event(name,ids,op,enzyme,source,start=.1,end=.85):
    added=set(op.after.next.items())-set(op.before.next.items());a,b=next(iter(added));bond=list(ids).index(b)
    return Event('seal',name,enzyme,bond=bond,junction=Junction(bond,source=source,mode=op.details['mode'],donor_5prime=op.before.five[b]),fraction_start=start,fraction_end=end,molecular_operation=op.to_dict())


def finish(m,slug,title,shots,steps,sim,limits):
    result=qualify(m,slug,title,shots,steps,limits=limits)
    for shot in result.shots[1:]:shot.transition=shot.establish
    return evidence(result,sim)


def crispr_mip():
    m=load(MIP,'crisprmip.py');r=simulate_mip(m);s=r['simulator'];source='CRISPR-MIP source v1 Methods lines458–493 and Table S2; selected illustrative cassette'
    n=len(r['probe']);N=len(r['filled']);l=len(m.LIG_ARM);e=len(m.EXT_ARM);gap=N-n;cell=7.4;x=220;y=515
    ext=x+e*cell;lig=ext+gap*cell;right=x+len(r['target'])*cell
    route=[(0,lig,y),(l,right,y),(l+18,right,325),(n-e-18,x,325),(n-e,x,y),(n,ext,y),(N,lig,y)]
    def tracks(state,growing=False,folded=True):
        ids=r['filled']if growing else r['probe']
        probe=strand('probe',r['filled_state'] if growing else state,ids,x,325,cell=cell,initial=n if growing else None,
            regions=[Region(0,l,'5′ ligation arm','primer',1,label_offset=(55,80)),Region(n-e,n,'3′ extension arm','bio',1,label_offset=(-50,80))])
        probe.color='primer';probe.sealed_marks=[0]
        if folded:probe.path=route if growing else route[:-1]
        target=strand('target',state,r['target'],right,y+24,cell=cell,direction=-1,
            regions=[Region(l,l+gap,'Guide cassette to capture','dna',1,label_offset=(0,95))])
        free=strand('free',state,r['unused'],230,720,cell=8,label='Unreacted probe');free.color='primer';free.label_position=(750,780)
        return [probe,target,free]
    start=tracks(r['initial'],folded=False)
    start[1].start=(start[1].start[0]+1800,start[1].start[1]);start[1].end=(start[1].end[0]+1800,start[1].end[1])
    arms=Shot('arms','One probe brackets the guide cassette',start,
        [Event('anneal','target',target=(-1800,0),fraction_end=.4),Event('fold','probe',route=route[:-1],fraction_start=.2,fraction_end=.65),
         Event('pair','probe',pairing='lig_arm',fraction_start=.65,fraction_end=.8,molecular_operation=r['bind5'].to_dict()),Event('pair','probe',pairing='ext_arm',fraction_start=.8,fraction_end=.95,molecular_operation=r['bind3'].to_dict())],
        'The phosphorylated probe binds both flanks of the selected guide cassette. Its free three-prime end faces the gap. The central backbone loops away from the target. An excess probe represents the unreacted linear pool.',source=source,scope='Source probe · illustrative guide cassette · representative excess probe',action=9,
        pairings=[Pairing('probe',0,'target',l-1,l,id='lig_arm',initially_bound=False),Pairing('probe',n-e,'target',len(r['target'])-1,e,id='ext_arm',initially_bound=False)])
    generated(arms,r['initial'])
    tr=tracks(r['before'],True);tr[0].regions.append(Region(n,N,'Copied cassette','new',-1,label_offset=(0,-50)))
    gapshot=Shot('gap','Phusion copies the gap up to the other arm',tr,
        [Event('extend','probe','Phusion polymerase',template='target',template_start=l+gap-1,primer_length=n,primer_binding_length=e,exact=True,molecular_operation=r['extend'].to_dict())],
        'Phusion extends the actual probe three-prime hydroxyl along the target. Synthesis stops at the occupied ligation arm. Copying fills the gap but does not close the backbone.',source=source,scope='Extension and ligation share the source incubation · operations shown separately',action=10)
    generated(gapshot,r['before'])
    tr=tracks(r['filled_state'],True);tr[0].initially_visible=None;tr[0].nicks=[0];tr[0].regions.append(Region(n,N,'Copied cassette','new',-1,label_offset=(0,-50)))
    seal=Shot('seal','Ampligase closes the capture circle',tr,[seal_event('probe',r['filled'],r['seal'],'Ampligase ligase',source)],
        'Ampligase seals the adjacent three-prime hydroxyl and five-prime phosphate. The target is a separate molecule, retained beside the circle. Pink marks show the new covalent junction.',source=source,scope='One actual closure bond · target and unused probe retained',action=5)
    generated(seal,r['filled_state'])
    tr=tracks(r['exo_before'],True);tr[0].initially_visible=None;tr[0].circular=True;tr[0].five_prime=tr[0].three_prime='none'
    exo=Shot('exonucleases','Accessible linear ends are digested',tr,
        [Event('resect3','free','Exonuclease I',amount=len(r['unused_digest'].details['removed']),fraction_end=.8,molecular_operation=r['unused_digest'].to_dict()),
         Event('resect3','target','Exonuclease III',amount=len(r['target_digest'].details['removed']),fraction_start=.15,molecular_operation=r['target_digest'].to_dict())],
        'Exonuclease I walks from the accessible single-stranded probe three-prime end; its source-qualified model retains a five-prime dinucleotide. Exonuclease III follows an accessible duplex three-prime end of this selected target fragment. The closed circle has no free entry end for these activities. Released nucleotides remain products.',source=source,scope='Selected accessible substrates · released nucleotide pool retained · no clearance yield',action=12,hold=4)
    exo.product_pool=True
    generated(exo,r['exo_before'])
    return finish(m,MIP,'CRISPR-MIP',[arms,gapshot,seal,exo],{'arms':['hybridize'],'gap':['extend_gap'],'seal':['seal_circle'],'exonucleases':['digest_linear']},s,
        ['Illustrative lentiCRISPRv2 guide-cassette member; source126nt versus Methods134bp discrepancy retained in source protocol.','Both enzymes act in one source incubation. Separate frames distinguish extension from sealing, not separate experimental additions.','Exonuclease activities and selected accessible ends are explicit; this is not complete genomic clearance or universal circle immunity.'])


def smallseq_simulation(m):
    from .simulator import Exonuclease
    source='Small-seq2018 detailed protocol, adapter ligation and cleanup; source orders'
    s=Simulator();rna=s.add('rna','N'*22,sugar='RNA',five='P',source=source+' illustrative22nt insert')
    ra3=s.add('ra3',m.RA3+'C',five='App',three='ddC',source=source)
    free=s.add('free',m.RA3+'C',five='App',three='ddC',source=source+' representative unligated excess adapter')
    rtp=s.add('rtp',m.RTP,sugars=[(0,len(m.RTP)-1,'DNA'),(len(m.RTP)-1,len(m.RTP),'RNA')],five='biotin',source=source)
    ra5=s.add('ra5',m.ra5_oligo().replace('T','U'),sugar='RNA',five='NH2',features=[(len(m.ra5_oligo())-10,len(m.ra5_oligo())-2,'umi.H8')],source=source)
    initial=s.snapshot();lig3=s.ligate(rna[-1],ra3[0],Ligase('T4 RNA ligase2 truncated KQ',('RNA',),('DNA',),('App',),(),()),mode='single-strand',cofactors=(),source=source)
    ligated=s.snapshot();dead=s.deadenylate('free',source=source)
    pair=s.anneal('rtp','free',binding=len(rtp),template_site=len(m.RTP)-1,source=source)
    digest_before=s.snapshot();digest=s.resect5('free',required_five=('P',),require_duplex=True,source=source+' lambda exonuclease; explicit nucleotide-level substrate boundary')
    lig5_before=s.snapshot();lig5=s.ligate(ra5[-1],rna[0],Ligase('T4 RNA ligase1',('RNA',),('RNA',),('P',),('ATP',),()),mode='single-strand',cofactors=('ATP',),source=source)
    return locals()


def smallseq():
    m=load(SMALL,'smallseq.py');r=smallseq_simulation(m);s=r['s'];source=r['source'];a=len(r['ra5']);b=len(r['rna']);c=len(r['ra3']);cell=13;x=300;y=420
    allids=list(r['ra5'])+list(r['rna'])+list(r['ra3'])
    def tracks(state):
        # The final continuous order is also a projection of the earlier separate
        # molecules: actual nicks/end groups come from the supplied state.
        main=strand('main',state,allids,x,y,cell=cell,label='Selected small RNA and adapters')
        main.initially_visible=None;main.color='rna';main.label='';main.label_position=(800,340)
        main.regions=[Region(0,a-10,'5′ adapter','primer',-1,label_offset=(-10,-80)),Region(a-10,a-2,'UMI · H8','bio',1,label_offset=(0,85)),Region(a,a+b,'Small RNA','rna',-1,label_offset=(0,-60)),Region(a+b,a+b+c,'3′ DNA adapter · ddC','adapter',1,label_offset=(0,80))]
        free=strand('free',state,r['free'],500,625,cell=cell,regions=[Region(0,c,'Free 3′ adapter','adapter',-1,label_offset=(230,-95))])
        rtp=strand('rtp',state,r['rtp'],500+len(r['rtp'])*cell,649,cell=cell,direction=-1,regions=[Region(0,len(r['rtp']),'5′ biotin RTP','primer',1,label_offset=(0,65))])
        return [main,free,rtp]
    tr=tracks(r['initial']);tr[2].start=(tr[2].start[0]+1800,tr[2].start[1]);tr[2].end=(tr[2].end[0]+1800,tr[2].end[1]);tr[0].pieces=[(0,a,-1800,0),(a+b,a+b+c,1800,0)];tr[1].start=(500,625);tr[1].end=(500+c*cell,625)
    lig3=Shot('lig3','Join the preadenylated 3′ adapter',tr,
        [Event('displace','main',span=(a+b,a+b+c),target=(-1800,0),fraction_end=.55),seal_event('main',allids,r['lig3'],'T4 RNA ligase 2 KQ',source,.6,.9)],
        'An eligible small RNA three-prime hydroxyl joins the preadenylated DNA adapter. The adapter ends in a physical dideoxycytidine. No ATP is supplied. An unreacted adapter remains in the lower branch.',source=source,scope='Source adapter orders · unknown22nt insert · no ATP in first ligation',action=7)
    next(e for e in lig3.events if e.kind=='seal').enzyme_label_offset=(0,-140)
    generated(lig3,r['initial'])
    tr=tracks(r['ligated']);tr[0].pieces=[(0,a,-1800,0)];tr[2].start=(tr[2].start[0]+1800,tr[2].start[1]);tr[2].end=(tr[2].end[0]+1800,tr[2].end[1])
    cleanup=Shot('cleanup','Prepare only the free adapter for digestion',tr,
        [Event('highlight','free','5′ deadenylase',span=(0,c),fraction_end=.4,enzyme_label_offset=(-100,-120),molecular_operation=r['dead'].to_dict()),Event('anneal','rtp',target=(-1800,0),fraction_start=.4,fraction_end=.7),Event('pair','rtp',pairing='free_adapter',fraction_start=.7,molecular_operation=r['pair'].to_dict())],
        'Deadenylase removes the free adapter five-prime adenyl group, leaving phosphate. The biotin-blocked RT primer anneals to this adapter. The successfully ligated small RNA remains intact above.',source=source,scope='Free-adapter branch · App → P · biotin protects RTP',action=7,
        pairings=[Pairing('rtp',0,'free',len(r['rtp'])-1,len(r['rtp']),id='free_adapter',initially_bound=False)])
    generated(cleanup,r['ligated'])
    tr=tracks(r['digest_before']);tr[0].pieces=[(0,a,-1800,0)]
    digest=Shot('digest','Lambda exonuclease removes phosphorylated adapter',tr,
        [Event('resect','free','Lambda exonuclease',amount=len(r['digest'].details['removed']),enzyme_label_offset=(-80,-115),molecular_operation=r['digest'].to_dict())],
        'Lambda exonuclease walks from the free adapter phosphorylated five-prime end. It stops at the unpaired terminal dideoxycytidine in this explicit substrate model. The protected RT primer and the ligated RNA are retained. Removed nucleotides remain in the released pool.',source=source,scope='Explicit substrate boundary · retained terminal ddC · released products retained',action=7)
    digest.product_pool=True
    generated(digest,r['digest_before'])
    tr=tracks(r['lig5_before']);tr[0].pieces=[(0,a,-1800,0)]
    lig5=Shot('lig5','Add the RNA 5′ adapter and its UMI',tr,
        [Event('displace','main',span=(0,a),target=(1800,0),fraction_end=.55),seal_event('main',allids,r['lig5'],'T4 RNA ligase 1',source,.6,.9)],
        'The RNA adapter arrives with its pre-existing eight-base H-class UMI. T4 RNA ligase one uses ATP to join its three-prime hydroxyl to the small RNA five-prime phosphate. The resulting molecule contains RNA and a terminal DNA adapter, not a newly synthesized DNA strand.',source=source,scope='ATP-dependent 5′ ligation · H means A/C/U · source amino block',action=7,hold=4)
    next(e for e in lig5.events if e.kind=='seal').enzyme_label_offset=(0,-140)
    generated(lig5,r['lig5_before'])
    return finish(m,SMALL,'Small-seq',[lig3,cleanup,digest,lig5],{'lig3':['ligate_3p'],'cleanup':['deadenylate'],'digest':['digest_free_adapter'],'lig5':['ligate_5p']},s,
        ['Selected adapter-ligation and free-adapter branch only; rRNA masking, RT and library PCR outside this movie.','Unknown source insert represented by22N; source H8 UMI remains a class, not a fabricated barcode value.','Lambda stopping at the terminal unpaired ddC is a strict declared substrate-boundary illustration, not a measured efficiency or residual-fragment claim.'])


def hic():
    m=load(HIC,'insitu_hic.py');r=simulate_hic(m);s=r['simulator'];source='Rao et al.2014 published in-situ Hi-C; Hi-C Protocol steps12–17; synthetic two-locus substrates'
    cell=15;y=470;orig=r['originals'];a,ab=orig['a'];b,bb=orig['b'];fillstate=r['seal_before'];At=r['left_after'];Bt=r['joined'][len(At):];Bb=list(r['other'])[:18];Ab=list(r['other'])[18:]
    # Source cuts leave four physical duplex fragments. The two partners not
    # selected for this junction remain represented in the same reaction vessel.
    leftovers=[('a_off_top',a[14:],450,280,1),('a_off_bottom',ab[:14],720,304,-1),('b_off_top',b[:14],920,680,1),('b_off_bottom',bb[14:],1190,704,-1)]
    def tracks(state,filled=False,joined=False):
        if joined:
            top=strand('top',state,r['joined'],180,y,cell=cell,label='Selected A–B contact');bottom=strand('bottom',state,r['other'],180+36*cell,y+24,cell=cell,direction=-1)
            top.label_position=(450,365);top.sealed_marks=[18];bottom.sealed_marks=[18]
            top.bands=[(0,18,'dna'),(18,36,'primer')]
            tracks=[top,bottom]
        else:
            top=strand('a_top',fillstate,At,180,y,cell=cell,initial=14 if not filled else None,
                regions=[Region(0,14,'Locus A','dna',-1,label_offset=(0,-75)),Region(14,18,'Fill-in','bio',1,label_offset=(-35,75))])
            bottom=strand('a_bottom',state,Ab,450,y+24,cell=cell,direction=-1)
            bt=strand('b_top',state,Bt,1050,y,cell=cell,regions=[Region(4,18,'Locus B','primer',-1,label_offset=(0,-75))])
            bot=strand('b_bottom',fillstate,Bb,1320,y+24,cell=cell,direction=-1,initial=14 if not filled else None)
            tracks=[top,bottom,bt,bot]
        tracks +=[strand(name,state,ids,x,yy,cell=cell,direction=direction)for name,ids,x,yy,direction in leftovers]
        for t in tracks:
            t.marks={str(i):'biotin' for i,id in enumerate(t.nucleotide_ids)if id in fillstate.atoms and fillstate.atoms[id].modification=='biotin'}
            lower=t.id in {'bottom','b_bottom'}
            t.mark_label_offsets={key:((65,85)if lower else(-65,-80)if joined else(90,-85))for key in t.marks}
        return tracks
    # Restriction view uses intact input strands; nucleotide positions agree with
    # the following split-fragment view after explicit displacement.
    tr=[]
    for key,top,bottom,x in [('a',a,ab,180),('b',b,bb,1050-14*cell)]:
        up=strand(key+'_top',r['before'],top,x,y,cell=cell,regions=[Region(14,18,'MboI · GATC','bio',-1,label_offset=(0,-75))])
        down=strand(key+'_bottom',r['before'],bottom,x+32*cell,y+24,cell=cell,direction=-1)
        tr.extend([up,down])
    events=[];pairs=[]
    for key,op,unpair in zip(['a','b'],r['cuts'],r['exposed']):
        events.append(Event('cleave',key+'_top','MboI restriction enzyme',bond=14,fraction_start=.05,fraction_end=.25,molecular_operation=op.to_dict()))
        pairs.append(Pairing(key+'_top',14,key+'_bottom',17,4,id=key+'_cohesive'))
        events.append(Event('unpair',key+'_top',pairing=key+'_cohesive',fraction_start=.3,fraction_end=.45,molecular_operation=unpair.to_dict()))
    events.extend([Event('displace','a_top',span=(14,32),target=(60,-190),fraction_start=.5),Event('displace','a_bottom',span=(0,14),target=(60,-190),fraction_start=.5),Event('displace','b_top',span=(0,14),target=(80,210),fraction_start=.5),Event('displace','b_bottom',span=(14,32),target=(80,210),fraction_start=.5)])
    cut=Shot('cut','MboI cuts both strands at GATC',tr,events,
        'MboI cuts two supplied duplexes and leaves four-base five-prime overhangs. The newly cut fragments separate after their cohesive contacts dissociate. Two ends are selected to illustrate a proximity junction; the other two fragments remain in the scene.',source=source,scope='Two synthetic loci · source GATC geometry · selected contact',action=8,pairings=pairs)
    generated(cut,r['before'])
    tr=tracks(r['fill_before'])
    fill=Shot('fill','Klenow fills both selected ends with biotin-dATP',tr,
        [Event('extend','a_top','Klenow polymerase',template='a_bottom',template_start=3,primer_length=14,primer_binding_length=14,exact=True,molecular_operation=r['fill_a'].to_dict()),Event('extend','b_bottom','Klenow polymerase',template='b_top',template_start=3,primer_length=14,primer_binding_length=14,exact=True,molecular_operation=r['fill_b'].to_dict())],
        'Klenow extends each selected recessed three-prime end across its own four-base overhang. Biotin-dATP is incorporated at adenine positions during this fill-in. The filled ends are blunt; their labelled bases are not pre-existing methylation or a change in sequence identity.',source=source,scope='Both selected ends filled · biotin-dATP supplied · unused partners retained',action=7)
    for e in fill.events:
        e.actor_offset=(0,20);e.enzyme_label_offset=(-200,-120)
    generated(fill,r['fill_before'])
    tr=tracks(fillstate,filled=True)
    approach=Shot('approach','Bring the selected blunt ends together',tr,
        [Event('anneal','b_top',target=(-600,0)),Event('anneal','b_bottom',target=(-600,0))],
        'The selected filled fragment ends approach each other. This illustrates proximity in the fixed nucleus, not base-pair-directed pairing of blunt ends. The other cut partners remain separate.',source=source,scope='Illustrative proximity · no contact-frequency prediction',action=5)
    generated(approach,fillstate)
    tr=tracks(fillstate,filled=True,joined=True)
    seal=Shot('contact','T4 DNA ligase seals both strands of the contact',tr,
        [seal_event('top',r['joined'],r['seals'][0],'T4 DNA ligase',source,.1,.5),seal_event('bottom',r['other'],r['seals'][1],'T4 DNA ligase',source,.55,.95)],
        'T4 DNA ligase joins the selected blunt ends, then seals the opposite-strand nick. Both original loci persist in a single duplex with the GATCGATC ligation junction. The biotin-bearing fill-in bases remain beside the junction. Pink marks identify both newly sealed bonds.',source=source,scope='Source GATCGATC junction · two covalent joins · biotin retained',action=8,hold=4)
    for e in seal.events:
        e.actor_offset=(0,20);e.enzyme_label_offset=(230,-120)
    seal.camera=(560,470,1.2)
    generated(seal,fillstate)
    return finish(m,HIC,'in situ Hi-C',[cut,fill,approach,seal],{'cut':['digest'],'fill':['fill'],'approach':['contact_ligation'],'contact':['contact_ligation']},s,
        ['Synthetic32bp loci and one selected junction; no chromatin-contact frequency or complete-library simulation.','Crosslink formation/reversal, fragmentation and streptavidin capture are outside this focused restriction/fill/ligation movie.','Two non-selected cut duplex fragments remain in view. Biotin is a physical modification on newly incorporated A, not a distinct reported base.'])
