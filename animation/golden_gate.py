"""Selected source mechanism on protocol-owned invented Golden Gate parts."""
import math
from dataclasses import replace
from .protocols import load,qualify
from .simulator import Simulator,Ligase
from .generated_views import strand
from .simulation_adapters import generated
from .model import Shot,Event,Region,Pairing,Junction

SLUG='golden-gate__10.1371+journal.pone.0003647'
SOURCE='Engler et al.2008 doi:10.1371/journal.pone.0003647, Results basic strategy and Figure1; protocol-owned invented linear parts, not study plasmids'


def simulate(*,phosphate=True,cofactors=('ATP',),fusion_mutation=False):
    m=load(SLUG,'golden_gate.py');s=Simulator();parents={};retained={};cutops=[]
    for key,f in zip(['a','b'],m.example_fragments()):
        p,t,b=m.precursor(f);seq=p.top()
        if fusion_mutation and key=='b':seq=seq[:7]+'A'+seq[8:]
        top=s.add(key+'_top',seq,source=SOURCE);bottom=s.add(key+'_bottom',m.revcomp(seq),source=SOURCE)
        s.duplex(key+'_top',key+'_bottom',source=SOURCE);parents[key]=(top,bottom);retained[key]=(top[t[0]:t[1]],bottom[b[0]:b[1]])
    initial=s.snapshot();offset=len(m.BsaI.site)
    for key in ['a','b']:
        for ref in [key+'_top',key+'_bottom']:
            cutops.append(s.restrict(ref,0,m.BsaI.site,offset+m.BsaI.top_cut,offset+m.BsaI.bottom_cut,source=SOURCE))
    cut=s.snapshot();release=[]
    for key in ['a','b']:
        kept=set(sum(retained[key],()))
        cross=[id for ids in retained[key]for id in ids if s.paired.get(id)not in kept]
        release.append(s.unpair(cross,source=SOURCE+' illustrative end fraying/exchange before cohesive annealing; not hydrolytic removal'))
    ready=s.snapshot();a,ab=retained['a'];b,bb=retained['b']
    pair1=s.bind_spans(a[:4],tuple(reversed(bb[:4])),source=SOURCE)
    pair2=s.bind_spans(b[:4],tuple(reversed(ab[:4])),source=SOURCE)
    paired=s.snapshot();lig=Ligase('T4 DNA ligase',('DNA',),('DNA',),('P',),('ATP',),('DNA',))
    if not phosphate:s.five[b[0]]='OH'
    seals=[]
    for acceptor,donor in [(a[-1],b[0]),(bb[-1],ab[0]),(b[-1],a[0]),(ab[-1],bb[0])]:
        seals.append(s.ligate(acceptor,donor,lig,mode='nick',cofactors=cofactors,source=SOURCE))
    return dict(module=m,simulator=s,parents=parents,retained=retained,initial=initial,cut=cut,cutops=cutops,release=release,ready=ready,pairs=[pair1,pair2],paired=paired,seals=seals,final=s.snapshot(),top=a+b,bottom=bb+ab)


def parent_tracks(r,state):
    tracks=[]
    for key,x in [('a',260),('b',920)]:
        top,bottom=r['parents'][key];n=len(top)
        for ids,suffix,xx,y,d in [(top,'top',x,420,1),(bottom,'bottom',x+n*12,444,-1)]:
            t=strand(key+'_'+suffix,state,ids,xx,y,cell=12,direction=d)
            t.regions=[Region(0,6,'BsaI site','bio',-1,label_offset=(0,-65)),Region(n-6,n,'BsaI site','bio',-1,label_offset=(0,-65))]if suffix=='top'else[]
            if suffix=='top':t.label='Invented part '+key.upper();t.label_position=(x+n*6,610)
            tracks.append(t)
    return tracks


def assembled_tracks(r,state):
    n=len(r['top']);tracks=[]
    for name,ids,radius,phase,direction in [('top',r['top'],135,0,1),('bottom',r['bottom'],159,4,-1)]:
        t=strand(name,state,ids,0,0)
        t.path=[(i/2,800+radius*math.cos(math.tau*(phase+direction*i/2)/n),470+radius*math.sin(math.tau*(phase+direction*i/2)/n))for i in range(2*n+1)]
        t.start=tuple(t.path[0][1:]);t.end=tuple(t.path[-1][1:]);t.nicks=[0,14];t.five_prime='P';t.three_prime='OH';t.sealed_marks=[0,14]
        t.bands=[(0,14,'adapter'),(14,n,'new')];tracks.append(t)
    kept=set(r['top']+r['bottom'])
    for t in parent_tracks(r,state):
        original=t.nucleotide_ids
        for j,(lo,hi) in enumerate([(0,7),(21,len(original))]):
            ids=original[lo:hi]
            # Exact remnants still contain paired and single-stranded ends.
            x=100 if t.id.startswith('a')else 1320
            flank=0 if (('top'in t.id and lo==0)or('bottom'in t.id and lo==21))else 1
            y=300+280*flank+(24 if 'bottom'in t.id else 0)
            origin=0 if flank==0 else 21
            position=lo if 'top'in t.id else len(original)-lo
            xx=x+(position-origin)*10
            t2=strand(t.id+'_off_'+str(j),state,ids,xx,y,cell=10,direction=1 if 'top'in t.id else-1)
            tracks.append(t2)
    return tracks


def movie():
    r=simulate();m=r['module'];tracks=parent_tracks(r,r['initial']);cuts=[]
    for j,op in enumerate(r['cutops']):
        actor=True
        for t in tracks:
            for left,right in op.to_dict()['backbone_removed']:
                if left in t.nucleotide_ids and right in t.nucleotide_ids:
                    cuts.append(Event('cleave',t.id,'BsaI'if actor else'',bond=t.nucleotide_ids.index(right),fraction_start=.1+j*.18,fraction_end=.24+j*.18,enzyme_class='nuclease',enzyme_label='BsaI',molecular_operation=op.to_dict()))
                    actor=False
    cut=Shot('cut','BsaI cuts outside its recognition sites',tracks,cuts,
        'BsaI recognizes its site but cuts outside it, making four-base five-prime overhangs. Both strands are cut at each flank. These short invented DNA parts illustrate the chemistry, rather than the source plasmids.',source=SOURCE,scope='Invented linear parts · source Type IIS mechanism',action=9)
    generated(cut,r['initial'])
    # One continuous view exposes the cohesive ends and retains every cut piece.
    separate=Shot('exchange','Free the compatible four-base ends',parent_tracks(r,r['cut']),[],
        'The terminal base pairs exchange so that matching cohesive ends can find one another. All cut flanks remain present. No DNA bases are degraded in this view.',source=SOURCE,scope='Illustrative end exchange · flanks retained',action=6,transition=1.8)
    for j,op in enumerate(r['release']):
        key='a'if j==0 else'b'
        for suffix in ['top','bottom']:
            t=next(t for t in separate.strands if t.id==key+'_'+suffix)
            # The central selected strand stays fixed; both flanks move outwards.
            left,right=((-90,-75),(90,75))if suffix=='top'else((90,75),(-90,-75))
            separate.events.extend([Event('displace',t.id,span=(0,7),target=left,fraction_start=.5),Event('displace',t.id,span=(21,32),target=right,fraction_start=.5)])
        # Pairing IDs are required for the explicit unpairing event; generated
        # frame contacts, not this singleton declaration, determine rendering.
        pair=Pairing(key+'_top',7,key+'_bottom',24,1,id=key+'_exchange',exact=True)
        separate.pairings.append(pair)
        separate.events.append(Event('unpair',key+'_top',pairing=pair.id,fraction_start=.1,fraction_end=.5,molecular_operation=op.to_dict()))
    generated(separate,r['cut'])
    assembled=assembled_tracks(r,r['ready'])
    anneal=Shot('anneal','Matching ends specify the circular order',assembled,[],
        'The two distinct complementary overhangs select the illustrated orientation. They pair without synthesizing DNA. Both circular-looking backbones still have two nicks, so neither strand is a covalently closed circle yet.',source=SOURCE,scope='Two compatible seams · four nicks remain',action=6,transition=1.8)
    for j,op in enumerate(r['pairs']):
        pair=Pairing('top',0 if j==0 else 14,'bottom',3 if j==0 else 17,4,id='fusion'+str(j))
        anneal.pairings.append(pair);anneal.events.append(Event('pair','top',pairing=pair.id,fraction_start=.15+j*.35,fraction_end=.45+j*.35,molecular_operation=op.to_dict()))
    generated(anneal,r['ready'])
    sealed=assembled_tracks(r,r['paired'])
    seal=Shot('seal','T4 DNA ligase seals all four nicks',sealed,[],
        'T4 DNA ligase uses ATP to join each three-prime hydroxyl to the adjacent five-prime phosphate. Both strands must be sealed at both seams. Pink marks the new bonds; neither recognition site is retained in the selected circular product.',source=SOURCE,scope='Both strands sealed · recognition sites discarded with flanks',action=9,transition=1.8)
    for j,op in enumerate(r['seals']):
        name='top'if j%2==0 else'bottom';bond=14 if j<2 else 0
        seal.events.append(Event('seal',name,'T4 DNA ligase',bond=bond,junction=Junction(bond,source=SOURCE,mode='nick'),fraction_start=.05+j*.22,fraction_end=.24+j*.22,enzyme_label='T4 DNA ligase',molecular_operation=op.to_dict()))
    generated(seal,r['paired'])
    final=Shot('product','The selected circle has lost the BsaI sites',assembled_tracks(r,r['final']),[],
        'The product retains both invented bodies and one copy of each fusion sequence. Both strands are covalently closed. The source mechanism enriches products that lose the restriction sites; this selected example does not predict competing products or yield.',source=SOURCE,scope='Invented construct · flanking pieces retained · no yield prediction',action=4,hold=4,transition=1.8)
    for t in final.strands[:2]:t.circular=True;t.five_prime=t.three_prime='none';t.nicks=[]
    generated(final,r['final'])
    return qualify(m,SLUG,'Golden Gate',[cut,separate,anneal,seal,final],{'cut':['a_process','b_process'],'exchange':['a_process','b_process'],'anneal':['anneal'],'seal':['seal'],'product':['seal']},limits=['Invented short linear parts from the protocol definition; not source plasmid sequences or a working vector.','End-fraying and selected assembly are illustrative; alternative ligations, kinetics and temperature cycles are not predicted.'],extra={'molecular_simulation':r['simulator'].trace()})
