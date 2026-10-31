"""Presentation choreography for original SMART-seq; functional unknowns stay schematic."""
from dataclasses import replace
from .model import Strand, Region, Shot, Event, Pairing, complement


def movie():
    from .protocols import load, qualify, tagmentation_shot
    slug='smart-seq__10.1038+nbt.2282';m=load(slug,'chem.py')
    handle=m.DT_SEQUENCE[:-32];h=len(handle)
    # One illustrative VN member (GA), transcript body and short terminal boundary.
    # N graphic units do not assign native TSO sequence, tail length or footprint.
    rseq='N'*38+'UC'+'A'*30
    primer=handle+'T'*30+'GA';p=len(primer)
    first=handle+complement(rseq[::-1]);f=len(first)
    tail=first+'NNN';t=len(tail);full=tail+'N'*h;n=len(full)
    x=lambda i:220+8*i
    rna=Strand('rna',rseq,(x(f),405),(x(h),405),'rna','',polymer='RNA',schematic=True,
               regions=[Region(0,38,'RNA transcript','rna',-1,(840,310)),Region(40,70,'poly(A) tail · AAAAA…','rna',-1,(540,345))])
    regs=[Region(0,h,'IS-PCR handle','primer',1,(315,570)),
          Region(h,h+30,'oligo(dT) · TTTTT…','bio',1,(550,570))]
    cdna=Strand('cdna',first,(220,465),(x(f),465),'dna','',schematic=True,initially_visible=p,regions=regs)
    pairs=[Pairing('cdna',h,'rna',69,32,exact=False,source=m.SOURCE+' CDS oligo(dT)30VN; illustrative GA anchor.')]
    anchor=Shot('anchor','RT primer binds the poly(A) tail',[rna,replace(cdna,start=(220,650),end=(x(f),650))],
        [Event('anneal','cdna',target=(0,-185),fraction_end=.75),Event('highlight','rna',span=(40,70)),Event('highlight','cdna',span=(0,p))],
        'The oligo dT reverse transcription primer binds the RNA poly A tail. Its five-prime amplification handle remains unpaired. That overhang provides an IS PCR binding site in the copied products.',
        action=5,pairings=pairs,scope='Original SMART-seq · representative transcript / functional kit boundaries')
    rt=Shot('rt','Reverse transcriptase copies the RNA',[rna,cdna],
        [Event('extend','cdna','Reverse transcriptase',persistent=True,template='rna',template_start=37,primer_length=p,primer_binding_length=32)],
        'Reverse transcriptase extends the bound primer toward the RNA five-prime end. The poly A tail and the amplification handle remain visible. RNA stays paired with its complementary DNA.',action=6,pairings=[Pairing('cdna',h,'rna',69,70,exact=False,source=m.SOURCE)],scope='RNA retained throughout first-strand synthesis')
    added=replace(cdna,sequence=tail,end=(x(t),465),initially_visible=f)
    terminal=Shot('terminal','Make the template-switching boundary',[rna,added],
        [Event('tail','cdna','Reverse transcriptase',persistent=True),Event('highlight','cdna',span=(f,t))],
        'At the RNA end, the retained reverse transcriptase makes a short untemplated boundary for template switching. The small drawn tail is schematic; its native sequence and length are not assigned.',action=3,scope='Terminal additions schematic · RNA remains paired')
    tso=Strand('tso','N'*(h+3),(x(n),405),(x(f),405),'adapter','',schematic=True,
        regions=[Region(0,h,'Template-switch oligo','adapter',-1,(1160,305))])
    switched=replace(cdna,sequence=full,end=(x(n),465),initially_visible=t,
        regions=regs+[Region(t,n,'Copied PCR handle','adapter',1,(1160,570))])
    switch=Shot('switch','Copy the template-switch oligo',[rna,replace(tso,start=(x(n),280),end=(x(f),280)),switched],
        [Event('anneal','tso',target=(0,125),fraction_end=.3),Event('extend','cdna','Reverse transcriptase',persistent=True,template='tso',template_start=h-1,primer_length=t,primer_binding_length=3,fraction_start=.35),Event('highlight','tso',span=(0,h))],
        'The same enzyme switches from RNA to the bound template-switch oligo. It copies the second amplification handle onto cDNA while the original RNA remains in place. The handles now bracket the transcript.',action=6,pairings=[*rt.pairings,Pairing('cdna',f,'tso',h+2,3,exact=False,source=m.SOURCE)],scope='Native TSO bases schematic · RNA–cDNA hybrid retained')
    finalcdna=replace(switched,initially_visible=None)
    minus=complement(full[::-1])
    reverse=Strand('reverse',minus,(x(n),650),(220,650),'new','',initially_visible=h,schematic=True,
        regions=[Region(0,h,'IS-PCR primer','primer',1,(1160,725))])
    firstcopy=Shot('firstcopy','Bind the PCR primer; make the second DNA strand',[rna,tso,finalcdna,reverse],
        [Event('separate','rna',target=(0,-110),fraction_end=.3),Event('separate','tso',target=(0,-110),fraction_end=.3),
         Event('anneal','reverse',target=(0,-125),fraction_start=.2,fraction_end=.48),
         Event('highlight','cdna',span=(t,n)),
         Event('extend','reverse','DNA polymerase',template='cdna',template_start=n-h-1,primer_length=h,fraction_start=.5)],
        'Heating separates the original templates from cDNA. An IS PCR primer binds the template-switch-derived handle. DNA polymerase extends it to create the complementary DNA strand; that strand does not appear spontaneously.',
        action=7,pairings=[Pairing('reverse',0,'cdna',n-1,h,exact=False,source=m.SOURCE)],scope='Functional first-copy illustration · exact native primer footprint schematic')
    # PCR parents remain present while two daughter strands grow from short primers.
    top=replace(finalcdna,regions=[Region(0,h,'PCR binding site','primer',-1,(320,310)),Region(t,n,'PCR binding site','adapter',-1,(1160,310))])
    bottom=replace(reverse,initially_visible=None,start=(x(n),525),end=(220,525),regions=[])
    newminus=replace(reverse,id='daughter_minus',start=(x(n),300),end=(220,300),regions=[Region(0,h,'IS-PCR primer','adapter',-1,(1160,365))])
    newplus=replace(finalcdna,id='daughter_plus',start=(220,760),end=(x(n),760),initially_visible=h,regions=[Region(0,h,'IS-PCR primer','primer',1,(320,740))])
    pcr=Shot('pcr','Both handles support PCR amplification',[top,bottom,newminus,newplus],
        [Event('separate','cdna',target=(0,-70),fraction_end=.2),Event('separate','reverse',target=(0,100),fraction_end=.2),
         Event('anneal','daughter_minus',target=(0,155),fraction_start=.2,fraction_end=.42),Event('anneal','daughter_plus',target=(0,-75),fraction_start=.2,fraction_end=.42),
         Event('highlight','cdna',span=(t,n)),Event('highlight','reverse',span=(n-h,n)),
         Event('extend','daughter_minus','DNA polymerase',template='cdna',template_start=n-h-1,primer_length=h,fraction_start=.45),
         Event('extend','daughter_plus','DNA polymerase',template='reverse',template_start=n-h-1,primer_length=h,fraction_start=.45)],
        'The next cycle separates the two DNA templates. Primers bind their highlighted amplification handles. Two new strands grow from the primers while the original DNA strands remain as templates.',action=7,
        pairings=[Pairing('daughter_minus',0,'cdna',n-1,h,exact=False,source=m.SOURCE),Pairing('daughter_plus',0,'reverse',n-1,h,exact=False,source=m.SOURCE)],scope='Later PCR cycle · template and daughter strands distinguished')
    tag=tagmentation_shot('tagment')
    # Functional tagged ends: no modern i5/i7/P5/P7 ordering is asserted for original kit.
    # Separate partial-copy step makes gap completion visible before index PCR.
    tagged='N'*90
    template=Strand('tagged_template',tagged,(1220,440),(380,440),'dna','',schematic=True,
        regions=[Region(0,15,'Transferred adapter','adapter',-1,(1150,335)),Region(75,90,'Transferred adapter','adapter',-1,(450,335))])
    fill=Strand('tagged_copy',tagged,(380,495),(1220,495),'dna','',initially_visible=75,schematic=True,
        regions=[Region(0,15,'Adapter end','adapter',1,(450,600))])
    gapfill=Shot('fill','Complete the tagged ends',[template,fill],
        [Event('extend','tagged_copy','DNA polymerase',template='tagged_template',template_start=14,primer_length=75,primer_binding_length=15)],
        'The tagged ends are completed before amplification. This view follows a selected productive fragment. The adapter boundaries are functional regions rather than an assignment of original Epicentre oligo sequences.',action=4,scope='Selected productive fragment · early end geometry schematic')
    # A primer includes an unpaired 5prime overhang plus a paired 3prime landing.
    product='N'*110
    copy=Strand('indexed_copy',product,(220,665),(1220,620),'new','',initially_visible=35,schematic=True,
        path=[(0,220,685),(12,320,685),(20,380,620),(35,520,620),(110,1220,620)],
        regions=[Region(0,12,'5′ overhang','bio',1,(245,720)),Region(12,20,'Sample index','ligase',1,(470,720)),Region(20,35,'3′ binding region','primer',1,(700,720))])
    opposite=Strand('opposite_copy',product,(1380,305),(380,370),'new','',initially_visible=35,schematic=True,
        path=[(0,1380,305),(12,1280,305),(20,1220,370),(35,1080,370),(110,380,370)],
        regions=[Region(0,20,'Opposite 5′ overhang','adapter',-1,(1320,270)),Region(20,35,'3′ binding region','primer',-1,(1030,270))])
    forward_template=replace(fill,initially_visible=None,regions=[],bands=[(0,15,'adapter'),(75,90,'adapter')])
    indexed=Shot('index','Index PCR: binding sites and 5′ overhangs',[replace(template,regions=[],bands=[(0,15,'adapter'),(75,90,'adapter')]),copy,forward_template,opposite],
        [Event('separate','tagged_template',target=(0,-50),fraction_end=.2),Event('separate','tagged_copy',target=(0,105),fraction_end=.2),
         Event('anneal','indexed_copy',target=(0,-170),fraction_start=.2,fraction_end=.4),Event('anneal','opposite_copy',target=(0,170),fraction_start=.2,fraction_end=.4),
         Event('highlight','tagged_template',span=(75,90)),Event('highlight','tagged_copy',span=(75,90)),
         Event('highlight','indexed_copy',span=(0,35)),Event('highlight','opposite_copy',span=(0,35)),
         Event('extend','indexed_copy','DNA polymerase',template='tagged_template',template_start=74,primer_length=35,primer_binding_length=15,fraction_start=.45),
         Event('extend','opposite_copy','DNA polymerase',template='tagged_copy',template_start=74,primer_length=35,primer_binding_length=15,fraction_start=.45)],
        'Each PCR primer binds through its three-prime region. Its five-prime overhang remains unpaired. Polymerase copies the insert, retaining those primer tails in the daughters. The sample index is carried by a primer; subsequent cycles copy the added ends into duplex library DNA.',
        action=7,scope='Functional primer tails · original kit sequence / index placement schematic',
        pairings=[Pairing('indexed_copy',20,'tagged_template',89,15,exact=False,source=m.SOURCE),Pairing('opposite_copy',20,'tagged_copy',89,15,exact=False,source=m.SOURCE)])
    # Final selected library: functional architecture with unknown ends plainly scoped.
    final=Strand('library','N'*120,(220,465),(1380,465),'dna','',schematic=True,
        regions=[Region(0,12,'Flow-cell attachment','primer',-1,(260,325)),Region(12,22,'Sample index','ligase',1,(430,610)),Region(22,35,'Read primer site','adapter',-1,(570,325)),Region(35,85,'cDNA insert','dna',1,(850,610)),Region(85,105,'Opposite adapter','adapter',-1,(1170,325)),Region(105,120,'Flow-cell attachment','primer',1,(1330,610))])
    finalbottom=Strand('library_complement','N'*120,(1380,515),(220,515),'dna','',schematic=True,
        bands=[(0,15,'primer'),(15,35,'adapter'),(85,98,'adapter'),(98,108,'ligase'),(108,120,'primer')])
    anatomy=Shot('library','The library that goes to the sequencer',[final,finalbottom],
        [Event('highlight','library',span=(a,b),fraction_start=i/6,fraction_end=(i+1)/6) for i,(a,b) in enumerate([(0,12),(12,22),(22,35),(35,85),(85,105),(105,120)])],
        'The selected sequencing library combines a cDNA insert with adapter ends, a read-primer binding site, flow-cell attachment regions and a sample index. These are functional labels; the original kit sequence, widths and precise end ordering remain unassigned.',action=8,hold=4,scope='Functional library anatomy · end ordering / widths schematic, not a modern dual-index recipe')
    shots=[anchor,rt,terminal,switch,firstcopy,pcr,tag,gapfill,indexed,anatomy]
    steps={s.id:[{'anchor':'reverse_transcription','rt':'reverse_transcription','terminal':'reverse_transcription','switch':'template_switch','firstcopy':'preamp','pcr':'preamp','tagment':'tagment','fill':'index_pcr','index':'index_pcr','library':'index_pcr'}[s.id]] for s in shots}
    return qualify(m,slug,'SMART-seq',shots,steps,limits=[m.CAVEAT,
        'Illustrative transcript, GA anchor member and graphic terminal tail; not native sequence/length assignments.',
        'First complementary strand and index-PCR anatomy are functional illustrations; native primer footprints and original Epicentre end geometry require source curation.',
        'Final adapter subregion ordering is a functional teaching diagram, not a validated sequence construct.'],
        extra={'source_cds_primer_5to3':m.DT_SEQUENCE,'original_protocol':'Ramskold2012; not SMART-seq2','omitted_presentation_steps':['strip_cleanup']})
