"""Source-specific selected reaction recipes for the complementary examples.

Only molecular inputs/capabilities are declared; products and intermediate bond
changes come from Simulator. Source-defined representative inputs stay local.
"""
from .simulator import Simulator,Polymerase,Ligase,Exonuclease
from .model import complement
rc=lambda s:complement(s[::-1])


def mip(m):
    source='CRISPR-MIP source probe and illustrative lentiCRISPRv2 cassette'
    s=Simulator();s.supply_nt('DNA','ACGT',source=source)
    probe=s.add('probe',m.probe_construct().top(),five='P',source=source)
    spacer=m.EXAMPLE_SPACER
    cassette=m.EXT_ARM+(''if spacer.startswith('G')else m.U6_PLUS1)+spacer+m.SCAFFOLD_V1+m.CAPTURE_3PRIME_FLANK+m.LIG_ARM
    target=s.add('target',rc(cassette),source=source)
    unused=s.add('unused_probe',m.probe_construct().top(),five='P',source=source+' one illustrative unligated excess-probe representative')
    initial=s.snapshot();l=len(m.LIG_ARM);e=len(m.EXT_ARM)
    bind5=s.bind_spans(probe[:l],tuple(reversed(target[:l])),source=source)
    bind3=s.anneal('probe','target',binding=e,template_site=len(target)-1,source=source)
    before=s.snapshot();extend=s.extend('probe',Polymerase('Phusion',('DNA',),'DNA'),source=source)
    if extend.details['stop']!='occupied_template':raise ValueError('MIP extension did not reach the occupied ligation arm')
    filled=s.chain('probe');filled_state=s.snapshot()
    ligase=Ligase('Ampligase',('DNA',),('DNA',),('P',),('NAD',),('DNA',))
    seal=s.ligate(s.tail('probe'),probe[0],ligase,mode='splinted',cofactors=('NAD',),source=source+' source buffer/NAD dependency')
    exo_before=s.snapshot()
    exoi=Exonuclease('Exonuclease I EN0581','3to5',('DNA',),('OH',),'single',2)
    exoiii=Exonuclease('Exonuclease III M0206','3to5',('DNA',),('OH',),'duplex',0,maximum_3prime_overhang=3)
    circle_exoi=s.digest('probe',exoi,source=source+' declared3prime-exonuclease activity only')
    circle_exoiii=s.digest('probe',exoiii,source=source+' declared3prime-exonuclease activity only; noAPactivity enabled')
    unused_digest=s.digest('unused_probe',exoi,source=source+' sourceEN0581 retains5prime dinucleotide')
    target_digest=s.digest('target',exoiii,source=source+' selected accessible linear target fragment; no quantitative genomic clearance prediction')
    return dict(simulator=s,unused=unused,exo_before=exo_before,circle_exoi=circle_exoi,circle_exoiii=circle_exoiii,unused_digest=unused_digest,target_digest=target_digest,probe=probe,target=target,initial=initial,bind5=bind5,bind3=bind3,before=before,extend=extend,filled=filled,filled_state=filled_state,seal=seal,closed=s.snapshot())



def smallseq(m):
    source='Small-seq source RA3/RTP/RA5 orders; unknown eligible small-RNA molecule'
    s=Simulator();rna=s.add('rna','N'*22,sugar='RNA',five='P',source=source+' illustrative insert length')
    ra3=s.add('ra3',m.RA3+'C',five='App',three='ddC',source=source)
    lig3_before=s.snapshot();rnl2=Ligase('T4 RNA ligase2 KQ',('RNA',),('DNA',),('App',),(),())
    lig3=s.ligate(rna[-1],ra3[0],rnl2,mode='single-strand',cofactors=(),source=source)
    first=s.chain('rna')
    free=s.add('free_adapter',m.RA3+'C',five='App',three='ddC',source=source)
    rtp=s.add('rtp',m.RTP,sugars=[(0,len(m.RTP)-1,'DNA'),(len(m.RTP)-1,len(m.RTP),'RNA')],five='biotin',source=source)
    s.anneal('rtp','free_adapter',binding=len(rtp),template_site=len(m.RTP)-1,source=source)
    deadenylate=s.deadenylate('free_adapter',source=source+' deadenylase branch')
    digest_before=s.snapshot()
    digest=s.resect5('free_adapter',required_five=('P',),require_duplex=True,source=source+' lambda exonuclease; walk until unpaired/non-DNA substrate')
    ra5=s.add('ra5',m.ra5_oligo().replace('T','U'),sugar='RNA',five='NH2',source=source)
    lig5_before=s.snapshot();rnl1=Ligase('T4 RNA ligase1',('RNA',),('RNA',),('P',),('ATP',),())
    lig5=s.ligate(ra5[-1],rna[0],rnl1,mode='single-strand',cofactors=('ATP',),source=source)
    return dict(simulator=s,rna=rna,ra3=ra3,first=first,free=free,rtp=rtp,ra5=ra5,lig3_before=lig3_before,lig3=lig3,deadenylate=deadenylate,digest_before=digest_before,digest=digest,lig5_before=lig5_before,lig5=lig5,product=s.chain('ra5'))


def celseq(m,body='ACGT'*8):
    source='CEL-Seq source RT primer; illustrative cDNA insert; initiation includes promoter-terminalGGG'
    s=Simulator();s.supply_nt('RNA','ACGU',source=source)
    parts=m.rt_parts();primer=''.join('N'*len(p.top)if p.placeholder else p.top for p in parts)
    # A concrete V=A member only in this synthetic molecule; barcode remains a
    # symbolic source variable whose identity is propagated into every transcript.
    if parts[-1].name=='V anchor':primer=primer[:-1]+'A'
    coding=primer+body
    top=s.add('coding',coding,source=source)
    bottom=s.add('template',rc(coding),source=source)
    s.duplex('coding','template',allow_ambiguous=True,source=source+' paired symbolic barcode copies')
    offset=sum(len(p.top)for p in parts[:next(i for i,p in enumerate(parts)if p.name=='T7 promoter')])
    prom_len=len(next(p.top for p in parts if p.name=='T7 promoter'))
    promoter=top[offset:offset+prom_len]
    # RNA starts at the terminal GGG within the printed promoter, not after it.
    initiation=offset+prom_len-len(m.C.T7_TRANSCRIBED_G)
    enzyme=Polymerase('T7 RNA polymerase',('DNA',),'RNA',de_novo=True)
    before=s.snapshot();rna1,op1=s.transcribe('template',len(bottom)-1-initiation,enzyme,promoter=promoter,source=source)
    second_before=s.snapshot();rna2,op2=s.transcribe('template',len(bottom)-1-initiation,enzyme,promoter=promoter,source=source)
    return dict(simulator=s,top=top,bottom=bottom,promoter=promoter,initiation=initiation,before=before,second_before=second_before,rna1=rna1,rna2=rna2,op1=op1,op2=op2)


def conversion(sequence,protected):
    source='scRRBS illustrative input sequence and selected5mC occupancy; not sample methylation'
    s=Simulator();s.supply_nt('DNA','ACGT',source=source)
    original=s.add('original',sequence,modifications={protected:'5mC'},source=source)
    before=s.snapshot();sites=[i for i in original if s.atoms[i].base=='C'and s.atoms[i].modification!='5mC']
    op=s.convert(sites,from_base='C',to_base='U',protected=('5mC',),source=source)
    converted=s.snapshot()
    # PCR/extension copies physical U into A; original sugar and5mC stay put.
    current=''.join(s.atoms[i].base for i in original)
    primer=s.add('primer',rc(current[-5:]),source=source+' illustrative primer complementary to converted terminus')
    s.anneal('primer','original',binding=len(primer),source=source)
    copy_before=s.snapshot()
    copy=s.extend('primer',Polymerase('uracil-tolerant polymerase',('DNA',),'DNA'),source=source)
    return dict(simulator=s,original=original,before=before,conversion=op,converted=converted,sites=sites,primer=primer,copy_before=copy_before,copy=copy,product=s.chain('primer'))



def gibson(m):
    source=m.EXAMPLE+'; source5prime isothermal mechanism; selected first junction'
    s=Simulator();s.supply_nt('DNA','ACGT',source=source);original={};processed={};ops=[]
    for key,fragment in zip(['a','b'],m.example_fragments()):
        parent,span_top,span_bottom=m.precursor(fragment);sequence=parent.top()
        a=s.add(key+'_top',sequence,source=source);b=s.add(key+'_bottom',rc(sequence),source=source)
        s.duplex(key+'_top',key+'_bottom',source=source)
        original[key]=(a,b)
    before=s.snapshot()
    for key,fragment in zip(['a','b'],m.example_fragments()):
        _,ts,bs=m.precursor(fragment)
        ops.append(s.resect5(key+'_top',ts[0],required_five=('OH','P'),source=source+' illustrative resection budget'))
        ops.append(s.resect5(key+'_bottom',bs[0],required_five=('OH','P'),source=source+' illustrative resection budget'))
        processed[key]=(s.chain(key+'_top'),s.chain(key+'_bottom'))
    anneal_before=s.snapshot()
    anneal=s.anneal('a_top','b_bottom',binding=len(m.RIGHT),template_site=len(s.chain('b_bottom'))-1,source=source)
    # The zero-gap toy genuinely stalls immediately at the other fragment's
    # already-paired body. The simulator adds zero nucleotides, without being
    # told to make a zero-gap product.
    gap=s.extend('a_top',Polymerase('DNA polymerase',('DNA',),'DNA'),source=source)
    if gap.added:raise ValueError('source exact-overlap toy unexpectedly required synthesis')
    seal_before=s.snapshot();lig=Ligase('DNA ligase',('DNA',),('DNA',),('P',),('ATP',),('DNA',))
    seal1=s.ligate(s.tail('a_top'),s.chain('b_top')[0],lig,mode='nick',cofactors=('ATP',),source=source)
    seal2=s.ligate(s.tail('b_bottom'),s.chain('a_bottom')[0],lig,mode='nick',cofactors=('ATP',),source=source)
    return dict(simulator=s,original=original,processed=processed,before=before,resect_ops=ops,anneal_before=anneal_before,anneal=anneal,gap=gap,seal_before=seal_before,seals=[seal1,seal2],joined_top=s.chain('a_top'),joined_bottom=s.chain('b_bottom'))


def hic(m):
    source='in-situ Hi-C source MboI/fill/ligation; two synthetic loci, no observed contact probability'
    s=Simulator();s.supply_nt('DNA','ACGT',modifications={'A':'biotin'},source=source+' biotin-dATP supply')
    originals={}
    for key,flank in [('a','ACGTACGTACGTAC'),('b','TGCATGCATGCATG')]:
        seq=flank+m.MBOI.site+flank
        top=s.add(key+'_top',seq,source=source);bottom=s.add(key+'_bottom',rc(seq),source=source)
        s.duplex(key+'_top',key+'_bottom',source=source);originals[key]=(top,bottom)
    before=s.snapshot();cuts=[]
    for key in ['a','b']:cuts.append(s.restrict(key+'_top',14,m.MBOI.site,m.MBOI.top_cut,m.MBOI.bottom_cut,source=source))
    a_top,a_bottom=originals['a'];b_top,b_bottom=originals['b']
    exposed=[]
    for top,_ in originals.values():exposed.append(s.unpair(top[14:18],source=source+' selected accessible cohesive ends; not backbone-induced DNA removal'))
    fill_before=s.snapshot()
    left_anchor=a_top[0];right_anchor=b_bottom[0]
    pol=Polymerase('Klenow',('DNA',),'DNA')
    fill_a=s.extend(left_anchor,pol,source=source);fill_b=s.extend(right_anchor,pol,source=source)
    a_selected=s.chain(left_anchor);b_selected=s.chain(b_top[14])
    seal_before=s.snapshot();lig=Ligase('T4 DNA ligase',('DNA',),('DNA',),('P',),('ATP',),('DNA',))
    seal_a=s.ligate(s.tail(left_anchor),b_selected[0],lig,mode='blunt',cofactors=('ATP',),source=source)
    seal_b=s.ligate(s.tail(right_anchor),s.chain(a_bottom[len(a_bottom)-18])[0],lig,mode='nick',cofactors=('ATP',),source=source)
    joined=s.chain(left_anchor);other=s.chain(right_anchor)
    capture_before=s.snapshot();capture=s.capture('streptavidin',left_anchor,required_modification='biotin',source=source)
    return dict(simulator=s,originals=originals,before=before,cuts=cuts,exposed=exposed,fill_before=fill_before,fill_a=fill_a,fill_b=fill_b,left_before=a_top[:14],left_after=a_selected,seal_before=seal_before,seals=[seal_a,seal_b],joined=joined,other=other,capture_before=capture_before,capture=capture)
