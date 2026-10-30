from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
class B:
    def __init__(self,title,authority,start,kind=K.BIOLOGICAL):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',start,kind,primary=True)
    def stock(self,key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind,inputs=(),geometry=None,product='Source mixture including incomplete/unreacted/off-target classes; occupancy/recovery unassigned',reason=None):
        y=self.g.product(key+'_product',product,K.MOLECULAR);kw={}
        if kind==S.MAGIC:kw['reason']=reason or 'Positive source operation; complete occupied molecular inputs/products and historical formulation unassembled, not a missing-information conclusion.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label,extra=()):
        y=self.g.product(key+'_kept','Actual source retained physical fraction; yield/carryover unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual separated source workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label))+tuple(extra),(y,w));self.g.discard(w,'Actual workup fraction excluded; no perfect specificity or molecular yield claim.');self.x=y
    def split(self,key,label,selected,remainder,disposition=Disposition.UNRESOLVED):
        y=self.g.product(key+'_selected',selected,K.MOLECULAR);r=self.g.product(key+'_remainder',remainder,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r));self.g.mark_terminal(Terminal(r.id,disposition,remainder));self.x=y;return y

def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('source products unassembled','[unassembled]',placeholder=True)]);return s
AUTHORITY='AssociatedAMNIHMS2201910/PMC13509204main238–286, oneoriginal; TableS7SD8gated. No olderSpatialATACRNA/final-body inheritance.'
def mix(b,key,label,parts):
    stock=b.g.input(key+'_stock',label,K.MOLECULAR);return components(b,key,stock,parts)
def components(b,key,stock,parts):
    out=[];ins=[]
    for suffix,kind,role,text in parts:
        x=b.g.product(key+'_'+suffix,text,kind,state=unknown_oligo_scene(text)if kind==K.OLIGO else None,chemical_state={'spatial_identifier':feature('spatial_hic_rna.axis_'+key[-1],'spatial_barcode','unknown',group='spatial_hic_rna.coordinate',part=key[-1]),'whole_coordinate':'Unassembled, same per-channel array'}if key.startswith('channel_')and suffix=='barcode'else{});out.append(x);ins.append((x,role))
    b.g.react(key+'_components','CopresentcomponentsfromONEactualphysicalpremix, notfreshadditions',S.GENERIC,((stock,R.TEMPLATE),),tuple(out),changes_molecule=False);return tuple(ins)
def pcrmix(b,key,label):
    return mix(b,key,label,[('enzyme',K.CATALYST,R.CATALYST,'Same-source-premixDNApolymerase'),('nt',K.MOLECULAR,R.STRUCTURAL,'Same-source-premixdNTPs, amountsunassembled'),('buffer',K.CONDITION,R.CONDITION,'Same-source-premixbuffer, fullformulationunassembled'),('forward',K.OLIGO,R.STRUCTURAL,'Sourcefirstexternalprimerofsamepreparedmix; bodyunassembled'),('reverse',K.OLIGO,R.STRUCTURAL,'Sourcesecondexternalprimerofsamepreparedmix; bodyunassembled')])
def unknown_oligo_scene(label):
    s=Scene();s.strand('source',[Segment(label,'[whole oligo unassembled]',placeholder=True)]);return s
def opaque():
    s=Scene();s.strand('source',[Segment('P5splintoligos/targetgeometryunassembled','[source P5 splint preparation]',placeholder=True)]);return s

def protocol_graph():
    b=B('Spatial Hi-C-RNA wholeROIandphysicalDNA/RNApartition',AUTHORITY,'Actual10umfrozentissue; adjacentH&Esectiondifferent/outsideproduction');g=b.g
    b.op('fa','Thaw37C1min/1percentFA+inhibitor5min',S.CROSS_LINKING,(b.cond('fa_condition','ActualsourceFA/inhibitor'),))
    b.op('qfa','1.25Mglycine5min',S.QUENCHING,(b.cond('glycine1','Actualglycine'),))
    b.op('dsg','3mMDSG/inhibitor40min',S.CROSS_LINKING,(b.cond('dsg_condition','ActualsourceDSG/inhibitor'),))
    b.op('qdsg','1.25Mglycine5min',S.QUENCHING,(b.cond('glycine2','Actualglycine'),));b.clean('fix_wash','Waterwash/10xbrightfield')
    b.op('permeabilize','NP40lysis15min/DPBSstop/.5percentSDS65C10min/1percentTriton5min',S.MAGIC,(b.cond('lysis','SourceIGEPAL/Tris/NaCl/proteaseinhibitor/RNaseinhibitor'),));b.clean('perm_wash','.5xDPBS/inhibitor5min')
    rt=mix(b,'initial_rt','OneactualRTpremix12.5ul5xbuffer/6.2ulMaxima/20.3ulDPBS/9.6ulwater/.8ulinhibitor/10ul100uMprimer; NTinventoryunassembled',[('enzyme',K.CATALYST,R.CATALYST,'ActualMaximaHMinus'),('primer',K.OLIGO,R.STRUCTURAL,'ActualRTprimer, body/biotinpositionunassembled'),('inhibitor',K.CATALYST,R.CATALYST,'ActualRNaseinhibitor'),('buffer',K.CONDITION,R.CONDITION,'Sourcebuffer/DPBS/water; noextraNTdoseinvented')])
    b.op('rt','RT30min/42C90min BEFOREcontactchemistry',S.MAGIC,rt,reason='SourceRTpremix positive; wholeprimer/NT/occupiedproducts/biotinpositionunassembled.');b.clean('rt_wash','NEBr3.1/inhibitor5min')
    b.op('digest','200ulDpnII/DdeI/HinfIeach1Uul37Cnight',S.DIGESTION,(b.stock('dpnii','DpnII'),b.stock('ddei','DdeI'),b.stock('hinfi','HinfI'),b.cond('digest_buffer','NEBr3.1/Triton/inhibitor')));b.clean('digest_wash','NEBr2.1/inhibitor5min')
    b.op('fill','DNApolI.2Uul/100uMdNTPRT2h',S.MAGIC,(b.stock('pol_i','ActualDNApolI'),b.stock('fill_nt','ActualdNTPmix, NOsourcebiotindATP',K.MOLECULAR,R.STRUCTURAL),b.cond('fill_buffer','NEBr2.1')),reason='Sourcebluntfill, exactends/productsunassembled; no canonicalHiCbiotinfill.');b.clean('fill_wash','T4buffer/inhibitor5min')
    b.op('contact','16UulT4ligaseRTnight BEFORETn5',S.MAGIC,(b.stock('contact_ligase','ActualT4ligase'),b.cond('contact_buffer','T4buffer/Triton/inhibitor')))
    b.op('post_sds','.1percentSDS55C1h/Triton5min',S.MAGIC,(b.cond('post_sds_condition','ActualsourceSDS/Tritonprogramme'),))
    wash=mix(b,'pre_tn5_washmix','ONEactualwashbuffer1percentBSA/Tris/NaCl/MgCl2/Tween/inhibitor',[('bsa',K.MOLECULAR,R.STRUCTURAL,'Copresent1percentBSA ofsamewashbuffer'),('buffer',K.CONDITION,R.CONDITION,'Copresentsalts/detergent/inhibitorofsamewashbuffer')])
    keep=g.product('pre_tn5_wash_kept','Actualwashedtissue, BSA/carryoverunassigned',K.MOLECULAR);waste=g.product('pre_tn5_wash_other','Actualremovedwashfractions',K.WASTE)
    g.react('pre_tn5_wash','Source2xwash5min withBSA-containingbuffer',S.CLEANUP,((b.x,R.TEMPLATE),)+wash,(keep,waste));g.discard(waste,'Actualwashfraction, noperfectBSAremovalclaimed.');b.x=keep

    loaded=mix(b,'loaded_tn5','ONEactual5ulB-onlyloadedtransposomestock; componentsremainloadedcomplex, notfreshfreeadditions',[('enzyme',K.CATALYST,R.CATALYST,'CopresentTn5catalystinloadedcomplex'),('rev',K.OLIGO,R.STRUCTURAL,'SourceMErevcomponentofSAMEloadedcomplex, bodyunassembled'),('b',K.OLIGO,R.STRUCTURAL,'SourceME-BcomponentofSAMEloadedcomplex, bodyunassembled')])
    b.op('tn5','50ulmix/5ulMErev-ME-BONLYtransposome55C1h',S.MAGIC,(*loaded,b.cond('tn5_buffer','Sourcewater/DPBS/2xbuffer/digitonin/Tween/inhibitor')),reason='B-onlyHiC source; ATACA+Bbranchnotmerged; exactends/gaps/productsunassembled.')
    b.op('tn5_stop','40mMEDTA/inhibitor5min',S.QUENCHING,(b.cond('edta','ActualEDTA'),));b.clean('tn5_wash','2xDPBS/inhibitor5min')
    for axis in ('a','b'):
        prep=g.input('barcode_'+axis+'_stock','ActualarrayofSEPARATE5.5ulpreparedchannelmixtures, NOTpooledchannels',K.MOLECULAR,chemical_state={'spatial_identifier':feature('spatial_hic_rna.axis_'+axis,'spatial_barcode','unknown',group='spatial_hic_rna.coordinate',part=axis),'source_recipe':'1.5ulannealedbarcode+4ulmix; source100channelrecipe288.2ulwater/23.2ulNEBr3.1/54ulT4buffer/1.8ulinhibitor/10.8ul5percentTriton/22ulligase','source_anneal':'1:1barcode/linker/2xannealbuffer, printed100umunitretained','whole_oligos':'TableS7unassembled'})
        load=g.product('barcode_'+axis+'_loaded','Actualarrayof4.5ulperchannel-loadingaliquots',K.MOLECULAR);rest=g.product('barcode_'+axis+'_rest','Actualpreparedchannelstockremainders',K.MOLECULAR)
        g.react('select_barcode_'+axis,'Select4.5ulperinletfrom5.5ulpreparedstock',S.SPLIT,((prep,R.TEMPLATE),),(load,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Actualuntransferredchannelstocks; retentionunassigned.'))
        channel=components(b,'channel_'+axis,load,[('barcode',K.OLIGO,R.STRUCTURAL,'ActualspatialbarcodecomponentofSAME4.5ulperchannelarray; bodyunassembled'),('linker',K.OLIGO,R.STRUCTURAL,'ActualligationlinkercomponentofSAME4.5ulperchannelarray; bodyunassembled'),('enzyme',K.CATALYST,R.CATALYST,'ActualT4ligasecomponentofSAMEchannelmixtures'),('buffer',K.CONDITION,R.CONDITION,'Actualcopresentchannelbuffer/Triton/inhibitor/water')])
        b.op('barcode_'+axis,'Orthogonalchannel'+axis+'37C30min/images',S.MAGIC,channel,reason='Actualchannelarrayapplication; no pooled80/200solutions orphysicalpixelsplit. Wholeligationgeometryunassembled.');b.clean('wash_barcode_'+axis,'Removechip/clamp/waterwash, retainROI')
    b.op('dapi','300nMDAPI10min/2xwaterwash/images',S.GENERIC,(b.stock('dapi_stock','ActualDAPIstain',K.MOLECULAR,R.STRUCTURAL),))
    b.op('roi_lyse','100ul.4mgmlPK/EDTA/Tris/NaCl/SDS58C2h/collect60Cnight',S.MAGIC,(b.stock('pk','ActualsourcePK'),b.cond('roi_lysis','SourcewholeROIlysisbuffer')));b.clean('lysate_purify','ZymoDNAcleanconcentrator5/elute100ulwater')
    cdna=g.product('bound_cdna','ActualC1-boundcDNA, exactbiotinlocation/yieldunassembled',K.MOLECULAR);gdna=g.product('gdna_supernatant','ActualgDNAsupernatant forHiC, NOTdiscardedwaste',K.MOLECULAR)
    g.react('partition','40ulC1in100ul2xB&W+100ulDNA RT1h/magnet',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),b.stock('c1','Actual40ulC1',K.SELECTION,R.SELECTION),b.cond('c1_buffer','Source2xB&W/inhibitor')),(cdna,gdna))
    b.x=gdna;b.clean('gdna_purify','ActualsupernatantZymopurify/elute20ul');b.op('gdna_denature','96C5min/ice2min',S.THERMAL_TREATMENT,(b.cond('dna_denature_condition','Actualheatprogramme'),))
    adapter=g.input('p5_splint_stock','ActualSL_P5_rc/S5MEA_N10P5splintstock12uMthen.75uMworkingstock',K.MOLECULAR,chemical_state={'source_recipe':'30.6ulwater/6ul100uMSL_P5_rc/8.4ul100uMS5MEA_N10/5ul10xRNA-ligasebuffer;95C1min/cool.1Cs to10C/5min','whole_bodies_ends':'Unassembled; N10notUMIproof'})
    selected=g.product('p5_splint_10','Actual10ul.75uMadapterinput',K.MOLECULAR);remaining=g.product('p5_splint_rest','Actualuntransferredadapterworkingstock',K.MOLECULAR);g.react('select_adapter','Source10ulselectedworkingstock',S.SPLIT,((adapter,R.TEMPLATE),),(selected,remaining));g.mark_terminal(Terminal(remaining.id,Disposition.UNRESOLVED,'Untransferredstock.'))
    lig=mix(b,'p5_master','ONEsource80ulP5ligationmasterdescriptor40ulPEG/12.5ulSCR/10ulDTT/10mMATP/1.25ulPNK/6.25ulT4ligase; ATPvolume/closureunassembled',[('pnk',K.CATALYST,R.CATALYST,'Actualsame-masterPNK'),('ligase',K.CATALYST,R.CATALYST,'Actualsame-masterT4ligase'),('atp',K.MOLECULAR,R.STRUCTURAL,'Actualsame-masterATP'),('buffer',K.CONDITION,R.CONDITION,'Actualsame-masterSCR/PEG/DTT/buffer')])
    ad=components(b,'p5_adapter',selected,[('sl',K.OLIGO,R.STRUCTURAL,'ActualSL_P5_rccomponentofSAMEselected10uladapterstock; bodyunassembled'),('mea',K.OLIGO,R.STRUCTURAL,'ActualS5MEA_N10componentofSAMEselected10uladapterstock; bodyunassembled'),('buffer',K.CONDITION,R.CONDITION,'CopresentbufferofSAMEselected10uladapterstock')])
    b.op('p5_ligate','P5splintligation37C45min/65C20min',S.MAGIC,ad+lig,geometry=AdapterLigationGeometry('unavailable','unavailable',(opaque(),),'SourceP5splintreactionpositive, wholeadapter/DNAends/junctionunassembled.'),reason='ActualPNK/T4ligasestock; no borrowedwholeadapter/N10/primerconstructor.')
    b.op('gdna_pcr','200ulmix150ul2xHighFidelitymaster/5uleachN70X-N50X/40ulwater;95C5min/10–18cycles95C30s56C30s72C60s/72C5min',S.MAGIC,pcrmix(b,'gdna_pcr_mix','ONEactual200ulpreparedPCRmixturewith150ul2xmaster'),reason='SourcePCR same-masterenzyme/NT/buffer/primerinventory, wholelanding/productsunassembled.');b.clean('hic_library_clean','.7xSPRI/elute20ul');hic=b.x
    b.x=cdna;b.clean('cdna_wash','2x400ulB&W/Tween/inhibitor;Tris/Tween/inhibitor;water')
    tso=mix(b,'tso','ONEsourceTSOpremixdescribed200ul, listed22+44+44+88+5.5+11+5.5=220ul; exactdoseunassembled',[('enzyme',K.CATALYST,R.CATALYST,'Actual11ulMaximaHMinus'),('nt',K.MOLECULAR,R.STRUCTURAL,'Actual22ul10mMdNTP'),('oligo',K.OLIGO,R.STRUCTURAL,'Actual5.5ul100uMTSO, bodyunassembled'),('inhibitor',K.CATALYST,R.CATALYST,'Actual5.5ulinhibitor'),('buffer',K.CONDITION,R.CONDITION,'Actual44ul5xbuffer/44ulFicoll/88ulwater')])
    b.op('tso','On-beadRT30min/42C90min',S.MAGIC,tso,reason='Source200ulversus220ulretained; finalformulation/wholeTSOlandingunassembled.');b.clean('tso_wash','400ulTris/Tweenthenwater')
    b.op('prepare_cdna_pcr','220ulbeadPCRstock110ul2xKAPAHiFiHotStart/8.8ul10uMprimer1and2/92.4ulwater',S.GENERIC,pcrmix(b,'cdna_pcr_mix','ONEactualsource220ulKAPApreparedPCRstock;8.8uleachprimerinterpretationunassembled'))
    b.split('cdna_50ul','Source50ulbeadPCRmixturepertube, selectedone','Actual50ulPCRinput','ActualotherbeadPCRallocations/stock, number/poolingscheduleunassembled')
    b.op('cdna_pcr5','95C3min/5cycles98C20s65C45s72C3min',S.MAGIC,reason='Actualselectedtube withcarriedKAPA/NT/primers, nofreshaddition; wholeproductsunassembled.')
    free=g.product('free_cdna_pcr','ActualPCRproductseparatedfrombeads, carriedreactioncomponents',K.MOLECULAR);beads=g.product('post_pcr_beads','Actualseparatedbeads/residualclasses',K.MOLECULAR);g.react('release_pcr','SourcePCRproductseparatedfrombeads BEFOREqPCR',S.SPLIT,((b.x,R.TEMPLATE),),(free,beads));g.mark_terminal(Terminal(beads.id,Disposition.UNRESOLVED,'BeadsnotimplicitlyreturnedtoqPCR.'));b.x=free
    b.op('sybr','AddSYBRGreen1xfinal',S.GENERIC,(b.stock('sybr_stock','ActualSYBRGreen',K.MOLECULAR,R.STRUCTURAL),))
    b.op('cdna_qpcr','95C3min/up-to15cycles98C20s65C20s72C3min/stopplateau/72C5min',S.MAGIC,reason='Actualpost-beadPCRproduct/carriedpremixplusSYBR, nofreshKAPA/NT.');b.clean('cdna_amplicon_clean','.6xSPRI/elute20ul')
    b.split('xt_2ng','Actual2ngcDNAstockaliquot/diluteto5ul','Actual2ngXTinput','ActualuntransferredcDNAstock/otherPCRallocations')
    b.op('xt_tagment','10ulTagmentDNAbuffer/5ulAmpliconTagmentmix55C5min',S.MAGIC,(b.stock('xt_mix','ActualhistoricalXT5ulAmpliconTagmentmix; oligos/endinventoryunassembled',K.MOLECULAR,R.STRUCTURAL),b.cond('xt_buffer','Actual10ulTagmentDNAbuffer')),reason='ActualsourceXT2ngboundary, no canonicalvendoradapter/gapconstructor.')
    b.op('xt_stop','5ulNTbufferRT5min',S.QUENCHING,(b.cond('xt_nt','ActualNTbuffer'),))
    b.op('xt_pcr','25ulmix15ulmaster/1uleachN50X-N70X/8ulwater;12cycles95C10s55C30s72C30s/72C5min',S.MAGIC,pcrmix(b,'xt_pcr_mix','ONEsource25ulXTpreparedPCRmixturewith15ulmaster'),reason='Source12cyclePCR, fullhistoricalmaster/primers/productsunassembled.');b.clean('rna_library_clean','.7xSPRI finalRNA library');rna=b.x
    for name,start in (('hic',hic),('rna',rna)):
        b.x=start;b.split(name+'_qc','ActualD5000analyticalallocation, amountunassembled','ActualproductionstockafterQCallocation','ActualQCaliquotnotreturned',Disposition.QC)
        b.split(name+'_loading','ActualNGSloadingallocation, amount/poolingunassembled','ActualNGScontribution','Actualuntransferredlibrary stock')
        lib=g.product(name+'_library','Actual'+name+'libraryclassfromphysicalpartition; noDNA/RNAmergeclaimed',K.LIBRARY);g.react(name+'_load','SourceNovaSeqXPlusPE150',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(read_lengths={'R1':150,'R2':150},limitation='Wholeprimers/barcodes/biotinpositions/indexmask/replicatepoolingunassembled. gDNAsupernatantandcDNAbeadsphysicalseparate; spatialassignmentcomputed.'))
    g.add_source_item(SourceItem('mixture_limits','RNA-TSOdescribed200ulbutcomponents220ul; barcodeprepared5.5ulversusloaded4.5ul retainedasstock/aliquot; RTandP5recipeNT/volumeclosureunassembled.','No silentunit/volumeclosure or publicationerrorclaim; TableS7SD8gated, N10suffixnotUMIproof.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();c('RTbeforecontactdigest',g.steps['digest'].inputs[0].material,'rt_wash_kept');c('ContactbeforeTn5',g.steps['post_sds'].inputs[0].material,'contact_product');c('ActualC1twofractions',set(g.steps['partition'].outputs),{'bound_cdna','gdna_supernatant'});c('HiCusesactualsupernatant',g.steps['gdna_purify'].inputs[0].material,'gdna_supernatant');c('RNAusesactualbeads',g.steps['cdna_wash'].inputs[0].material,'bound_cdna');c('PCRusesselected50ul',g.steps['cdna_pcr5'].inputs[0].material,'cdna_50ul_selected');c('SYBRafterbeadseparation',g.steps['sybr'].inputs[0].material,'free_cdna_pcr');c('XTusesselected2ng',g.steps['xt_tagment'].inputs[0].material,'xt_2ng_selected');c('Barcodeloadselectedstock',g.steps['channel_a_components'].inputs[0].material,'barcode_a_loaded');c('TwophysicalNGSclasses',set(g.sequencing),{'hic_library','rna_library'});c('OneactualKAPAstock',g.steps['cdna_pcr_mix_components'].inputs[0].material,'cdna_pcr_mix_stock');c('NofreshKAPAintoqPCR',len(g.steps['cdna_qpcr'].inputs),1)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'j.cell.2026.07.039_PMC13509204.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1016/j.cell.2026.07.039')):
        t=p.read_text();c('OriginalgDNAfraction','supernatant (containing gDNA)'in t,True);c('OriginalBonlytransposome','For Hi-C, Tn5 was assembled with Tn5MErev/ME-B only'in t,True);c('Originalbeadproductseparation','product was separated from beads'in t,True)
    c('P5oligosfromactualselected10ul',g.steps['p5_adapter_components'].inputs[0].material,'p5_splint_10');c('BSAfromonewashbufferstock',g.steps['pre_tn5_washmix_components'].inputs[0].material,'pre_tn5_washmix_stock');c('LoadedTn5componentsfromonephysicalstock',g.steps['loaded_tn5_components'].inputs[0].material,'loaded_tn5_stock');c('NamedRTTSOandP5inventorytypedoligo',all(g.materials[k].kind==K.OLIGO for k in ('initial_rt_primer','tso_oligo','p5_adapter_sl','p5_adapter_mea')),True);c('Partialnotpublicready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
