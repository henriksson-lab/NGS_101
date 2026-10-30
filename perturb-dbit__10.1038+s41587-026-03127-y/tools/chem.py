from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition

def unknown_scene(label='Unassembled whole molecular product'):
    sc=Scene();sc.strand('unassembled',[Segment('unassembled','[unassembled]',placeholder=True)]);return sc

def add(g,key,label,kind=K.MOLECULAR,role=R.STRUCTURAL,**kw):
    return (g.input(key,label,kind,**kw),role)
def condition(g,key,label):return add(g,key,label,K.CONDITION,R.CONDITION)
def opaque(g,key,label,x,ins=(),kind=S.MAGIC):
    y=g.product(key+'_product',label+'; whole molecular identities/junctions unassembled',K.MOLECULAR)
    kw={'reason':'Source-positive operation; whole kit composition, substrate context or product geometry unassembled, not publication absence.'} if kind==S.MAGIC else {}
    g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(ins),(y,),**kw);return y
def cleanup(g,key,label,x):
    y=g.product(key+'_kept','Source retained preparation; recovery/carryover unassigned',K.MOLECULAR);w=g.product(key+'_removed','Source workup fraction',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),condition(g,key+'_condition',label)),(y,w));g.discard(w,'Physical workup, not perfect molecular selectivity.');return y
def aliquot(g,key,label,x,remainder_label='Untransferred source stock; later fate unassigned'):
    y=g.product(key+'_chosen',label,K.MOLECULAR);r=g.product(key+'_remainder',remainder_label,K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,remainder_label));return y
def qc(g,key,label,x):
    a=g.product(key+'_analytical',label,K.MOLECULAR);k=g.product(key+'_kept','Retained production stock',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(a,k));g.mark_terminal(Terminal(a.id,Disposition.QC,'Analytical aliquot not returned; amount unassigned.'));return k
from chemdraw import feature
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from illumina import P5,P7
from nextera import ADAPTOR_S5,ADAPTOR_S7
TITLE='Perturb-DBiT: published selected polyadenylation-capture branch'
AUTHORITY='PublishedNaturemain1091–1107, publishedSI1computationalscope, actualMOESM4TableS3wholeoligos/representativeA1B1 andTableS6headers originalZIP/XML. Three originals; DC/vendorwholechemistrynotinherited.'
RT_HANDLE='CATCGGCGTACGACT'
RT_ORDER='/5Phos/'+RT_HANDLE+'N'*10+'/iBiodT/'+'T'*15+'VN'
TSO_PREFIX='AAGCAGTGGTATCAACGCAGAGTGAAT'
TSO_ORDER=TSO_PREFIX+'rGrG+G'
PCR1='CAAGCGTTGGCTTCTCGCATCT';PCR2='AAGCAGTGGTATCAACGCAGAGT'
A_HEAD='AGGCCAGAGCATTCG';A_TAIL='GTGGCCGATGTTTCG';B_TAIL='ATCCACGTGCTTGAG';BC1='AACGTGAT'
A1=A_HEAD+BC1+A_TAIL;B1=PCR1+BC1+B_TAIL
LINK1='AGTCGTACGCCGATGCGAAACATCGGCCAC';LINK2='CGAATGCTCTGGCCTCTCAAGCACGTGGAT'
N501=P5+'TAGATCGC'+ADAPTOR_S5
N701=P7+'TCGCCTTA'+ADAPTOR_S7+PCR1
SOURCE_ORDERINGS={'RT':RT_ORDER,'TSO':TSO_ORDER,'PCR1':PCR1,'PCR2':PCR2,'Link1':LINK1,'Link2':LINK2,'A1':'/5Phos/'+A1,'B1':B1,'N501':N501,'N701':N701,'N501direct':N501+PCR2,'N701CRISPR':P7+'TCGCCTTA'+ADAPTOR_S7+'GTCAGATGTGTATAAGAGACAG'}
def rt_parts(prefix='rt'):
    return[Segment(prefix+'_handle',RT_HANDLE),Segment(prefix+'_umi','N'*10,placeholder=True,feature=feature('perturb.umi','umi','random')),Segment(prefix+'_bioT','T',note='Source /iBiodT/ is one EXTRA thymidine nucleotide bearing internal biotin'),Segment(prefix+'_polyT','T'*15),Segment(prefix+'_VN','VN',placeholder=True,note='Anchoredprimingdegeneracy,notidentifier')]
def rt_scene():
    s=Scene();s.strand('rt',rt_parts(),mod5='P');return s
def oligo(seq,key='source',mod5=''):
    s=Scene();s.strand(key,[Segment(key,seq)],mod5=mod5);return s
def tso_scene():
    s=Scene();s.strand('tso',[Segment('DNAhandle',TSO_PREFIX),Segment('RNAterminal','GGG',note='Source rG,rG,+G; finalG LNA-modified, allthreephysicalresidues')]);return s
def barcode_parts(which):
    head,tail=(A_HEAD,A_TAIL)if which=='A'else(PCR1,B_TAIL)
    return[Segment(which+'_head',head),Segment(which+'_bc',BC1,feature=feature('perturb.'+which.lower(),'spatial_barcode','fixed',group='perturb.pixel'),note='Representative sourcebarcode1; whole50/100channelallocationunassembled'),Segment(which+'_tail',tail)]
def end_scene(which, include_acceptor=True):
    s=Scene();s.strand('barcode',barcode_parts(which),mod5='P'if which=='A'else'')
    link=LINK1 if which=='A'else LINK2
    s.anneal('link',[Segment('acceptor_comp',link[:15]),Segment('adapter_comp',link[15:])],to='barcode',pair=('adapter_comp',which+'_tail'),unpaired=('acceptor_comp',))
    if not include_acceptor:return s
    if which=='A':parts=rt_parts('acceptor');pair='acceptor_handle'
    else:parts=[Segment('acceptor_head',A_HEAD),Segment('previous_A','[A barcode and RT/cDNA context]',placeholder=True)];pair='acceptor_head'
    s.anneal('acceptor',parts,to='link',pair=(pair,'acceptor_comp'),unpaired=tuple(p.name for p in parts if p.name!=pair),mod5='P')
    s.ligation_end('barcode',which+'_tail',kind='splinted',polarity='3-prime',owner='adapter',text='Sourcebarcode3OH localinterface; remaining5Pacceptor/cDNA/polymercontextunassembled')
    s.ligation_end('acceptor',pair,kind='splinted',polarity='5-prime',owner='substrate',text='SourceRT5P orpreviousbarcodeA5P; nofullsealedproductclaim');return s

def pcr_from_mixture(g,key,mix,label,fwd= PCR1,rev= PCR2):
    template=g.product(key+'_template','Same-tubeDNAtemplatecomponent; wholeampliconunassembled',K.MOLECULAR);cat=g.product(key+'_enzyme','SamepreparedPCRenzymecomponent',K.CATALYST);nt=g.product(key+'_nt','SamepreparedPCRnucleotidecomponent',K.MOLECULAR)
    f=g.product(key+'_f','Samepreparedsourceforwardprimercomponent',K.OLIGO,state=oligo(fwd,key+'_f')if fwd else None);r=g.product(key+'_r','Samepreparedsourcereverseprimercomponent',K.OLIGO,state=oligo(rev,key+'_r')if rev else None);buf=g.product(key+'_buffer','SamepreparedPCRbuffer',K.CONDITION)
    g.react(key+'_components','Samephysicalmixturecomponents, notnewstocks/orphysicalfractionation',S.GENERIC,((mix,R.TEMPLATE),),(template,cat,nt,f,r,buf),changes_molecule=False)
    y=g.product(key+'_product','SourceamplifiedDNApopulation; wholeproducts/unreactedclassesunassembled',K.MOLECULAR)
    g.react(key,label,S.PCR,((template,R.TEMPLATE),(cat,R.CATALYST),(nt,R.STRUCTURAL),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(buf,R.CONDITION)),(y,),geometry=PCRGeometry.unavailable(unknown_scene(),reason='WholeRNA/cDNA/TSOproducts/sourceprimerlandingsunassembled, notsourceabsence.'));return y

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('prepared_tissue','SourcePACtissuesectionROI; upstreamfixation/deparaffinizationunassembled, sourcegRNAexpressionconstructsoutsidecaptureboundary',K.BIOLOGICAL,primary=True)
    x=opaque(g,'permeabilize','Source1percentTritonDPBS20minRT→DPBS-RIwash/airdry/reservoir',x,(condition(g,'perm_buffer','SourceTriton/DPBS-RI'),))
    x=cleanup(g,'pap_equilibrate','Source100ulpolyAwashbuffer5minRTthenREMOVEbeforeenzymaticmix',x)
    x=opaque(g,'polyadenylate','Source60ulPAPreaction37C25min',x,(add(g,'pap','Source6ul5U/ulEcoliPAP',K.CATALYST,R.CATALYST),add(g,'pap_atp','Source6ul10mMATP structuraladenylateprecursor'),condition(g,'pap_conditions','Source6ul10xbuffer/water/inhibitors')),S.TAILING)
    x=cleanup(g,'post_pap','SourceDPBSshakewash beforeRT',x)
    x=opaque(g,'rt','Source60.2ulPACRT30minRT+90min42C',x,(add(g,'rt_primer','Source10ul25uMtranscriptomeRTprimer; noDCSCaffoldprimer',K.OLIGO,R.STRUCTURAL,state=rt_scene(),chemical_state={'source_ordering':RT_ORDER,'internal_biotin_thymidine':'One physicalTbetweenUMI10andT15VN, notnonntspacer'}),add(g,'maxima_rt','Source6ul200U/ulMaximaHMinusRT',K.CATALYST,R.CATALYST),add(g,'rt_nt','Source3ul10mMdNTP,perbase/totalconventionunassigned'),condition(g,'rt_conditions','Source12ul5xRTbuffer/DPBS-RI/water/RNaseinhibitors')),S.REVERSE_TRANSCRIPTION)
    x=cleanup(g,'post_rt','Source5min50mlDPBSshakewash',x)
    for which,raw in [('A',A1),('B',B1)]:
        q=which.lower()
        stock=g.input('channel_'+q+'_stocks','Sourceprepared50or100barcode/linkerplate; equimolar100uMstocks/annealing95C5min→20C-.1C/s→12C3min; exactperchannelstockdistributionunassembled',K.ADAPTER,state=end_scene(which,False),chemical_state={'representative_source':'/5Phos/'+raw if which=='A'else raw,'linker':LINK1 if which=='A'else LINK2,'selected_channel_scope':'Representativefamilyinterface,notwholeplateallocation'})
        chosen=g.product(q+'_chosen_barcodes','Source1ulbarcodealiquots selectedperchannel; full50/100channelallocationsunassembled',K.ADAPTER,state=end_scene(which,False));remaining=g.product(q+'_barcode_remainder','Untransferredpreparedbarcodeplate stocks',K.ADAPTER)
        g.react(q+'_barcode_aliquot','Sourcechannelbarcodestockselection, amountdistributedperchannelknown/totalstockunassigned',S.SPLIT,((stock,R.TEMPLATE),),(chosen,remaining));g.mark_terminal(Terminal(remaining.id,Disposition.UNRESOLVED,'Untransferredbarcodestocknotassumedallconsumed.'))
        y=g.product('barcode_'+q+'_product','Source'+which+'barcodedcDNAmixture; completechannelproducts/junctions/recoveryunassembled',K.MOLECULAR)
        g.react('barcode_'+q,'Sourceorthogonal'+which+'channels:each5ul=4ulligationmix+1ulselectedbarcode,37C35min; totalpreparedmixes/individualchannelallocationunassembled',S.LIGATION,((x,R.TEMPLATE),(chosen,R.STRUCTURAL),add(g,q+'_ligase','ActualT4DNAligasecomponentofsourcechannelmix; originalmix15ulstock, preparedmixcountunassigned',K.CATALYST,R.CATALYST),condition(g,q+'_ligation_buffer','SourceNEBuffer3.1/T4buffer/water/Triton/inhibitors; perchannel4ulselectedmixedpreparation')), (y,),geometry=AdapterLigationGeometry.known('splinted',end_scene(which),statement='Source-ownedlocalbarcode3OH/acceptor5Psplintedinterface. FullRNA/cDNA/channelproducts andsealedjunctionsunassembled; notwholeligationproof.'));x=y
        x=cleanup(g,'post_'+q+'_wash','SourceDPBSshakewash/device removal after'+which, x)
    x=cleanup(g,'channel_wash','Source5x500ulwaterflowwash afterB/finalscan',x)
    x=opaque(g,'lyse','Source110ulPKlysis55C2h; collect110ulcDNA-richlysate',x,(add(g,'pk','Source10ul20ug/ulProteinaseK',K.CATALYST,R.CATALYST),condition(g,'lysis_buffer','Source50ulDPBS+50ul2x20mMTris/400mMNaCl/100mMEDTA/4.4percentSDS')),S.DEPROTEINIZATION)
    selected=g.product('lysate_100','Source100ulselectedlysateforcapture',K.MOLECULAR);rest=g.product('lysate_10','Source10uluntransferredlysatestock',K.MOLECULAR);g.react('lysate_aliquot','Actual100of110ullysate',S.SPLIT,((x,R.TEMPLATE),),(selected,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Untransferredlysate.'))
    bound=g.product('bound_cdna','SourceC1-retainedbiotincDNAenrichedmixture; occupancyunassigned',K.MOLECULAR);flow=g.product('capture_flowthrough','Sourcenonretainedcapturefraction',K.MOLECULAR)
    g.react('capture','Source40ulC1beads/100ul2xBWbufferRT60min',S.AFFINITY_CAPTURE,((selected,R.TEMPLATE),add(g,'c1','SourceDynabeadsC1',K.SELECTION,R.SELECTION),condition(g,'bw','SourceBWbuffer')),(bound,flow));g.mark_terminal(Terminal(flow.id,Disposition.UNRESOLVED,'Nonretainedfractionnotimplicitlysequenced.'));x=cleanup(g,'capture_wash','Source2xBW/.05Tween+2xTris/.1Tweenwashes',bound)
    x=opaque(g,'template_switch','Source200ulTSOmix30minRT+90min42C; alreadyprimedcDNA/RNAcontextretained',x,(add(g,'tso','Source10ul100uMTSOacceptor',K.OLIGO,R.STRUCTURAL,state=tso_scene(),chemical_state={'source_ordering':TSO_ORDER,'terminal_sugars':'rG,rG,+GphysicalGresidues'}),add(g,'maxima_ts','Source10ul200U/ulfreshMaximaHMinusRT',K.CATALYST,R.CATALYST),add(g,'ts_nt','Source20ul10mMdNTP'),condition(g,'ts_conditions','Source40ul5xRTbuffer/40ul20percentFicoll/5ulinhibitor/75ulwater')),S.TEMPLATE_SWITCHING)
    x=cleanup(g,'post_ts','SourceTris/.1Tweenwashthenwaterwash',x)
    mixed=opaque(g,'pcr_mix','Sourcebeadsresuspended200ulPACPCRmix:100ul2xKAPA/84ulwater/8ul10uMPCR1/8ul10uMPCR2',x,(add(g,'pcr_stock','Actualone200ulpreparedPACPCRstock, enzyme/nucleotide/primercomponentsco-present'),),S.GENERIC)
    tubes=[g.product('pcr_tube_'+str(i),'SourceoneoffourdistributedPCRtubemixtures',K.MOLECULAR)for i in range(4)];g.react('four_tubes','Source200ulmixturephysicallydistributedfourstrip tubes',S.SPLIT,((mixed,R.TEMPLATE),),tuple(tubes))
    products=[pcr_from_mixture(g,'pcr'+str(i),t,'Source13cycles98C20s/63C45s/72C3min')for i,t in enumerate(tubes)]
    wta=g.product('wta_fraction','Source0.7xpelletWTAenrichedfraction',K.MOLECULAR);small=g.product('sgrna_fraction','Sourcesaved0.7xsupernatantsgRNAenrichedfraction',K.MOLECULAR)
    g.react('pcr_fractionate','SourcefullPCRproduct0.7xSPRI; exactper-tubeversuscombinedworkupscheduleunassembled',S.MAGIC,tuple((p,R.TEMPLATE)for p in products)+(condition(g,'fractionation','Source0.7xSPRI; bothfractionsretained'),),(wta,small),reason='Sourcephysicalfractionspositive; timingofcombiningfourPCRtubesunassembled,nofakepool.')
    wta=cleanup(g,'wta_elute','Sourcepelletpurification/elute20ulwater',wta);small=cleanup(g,'sgrna_select','Sourcesupernatantsecond1.2xSPRI/elute20ulwater',small)
    for name,part in [('wta',wta),('sgrna',small)]:
        part=qc(g,name+'_initial_qc','SourceTapeStationD5000analyticalaliquot',part);part=aliquot(g,name+'_depletion_input','Source20ngselectedcDNAforribodepletion',part)
        for n in (1,2):part=opaque(g,name+'_ribo'+str(n),'SourceSEQuoiaribodepletionround'+str(n),part,(add(g,name+'_ribo_stock'+str(n),'Actualvendorribodepletionpreparation; reagentidentities/physicalclassesunassembled'),))
        part=aliquot(g,name+'_final_pcr_input','Sourceapproximately42ulprocessedsolutionfor100ulfinalPCR',part)
        mix=opaque(g,name+'_final_mix','Source100ulfinalPCRmixture:50ul2xKAPA/~42ulsolution/4ul10uMP5/4ul10uMP7',part,(add(g,name+'_final_pcr_stock','ActualpreparedP5/P7PCRstock; TableS3plain/extendedN501familyroleassignment/sampleindicesunassembled'),),S.GENERIC)
        # P5/P7 whole-family assignment stays unknown, rather than guessing plainN501 versus N501direct.
        part=pcr_from_mixture(g,name+'_final_pcr',mix,'Sourcesevencyclesusingaforementionedprogramme, PCRattachmentnotanextraligasestep',None,None)
        part=cleanup(g,name+'_final_clean','Source'+('0.7xSPRI'if name=='wta'else'0.6xthen1.2xSPRI; intermediatefractionselectionunassembled'),part);part=qc(g,name+'_final_qc','SourceTapeStationlibraryQCaliquot',part)
        lib=g.product(name+'_library','Source'+name+'PACsequencinglibrary',K.LIBRARY);g.react(name+'_load','Sourcelibraryloadingallocationunassembled',S.MAGIC,((part,R.TEMPLATE),),(lib,),reason='SourceNGSpositive; actualloadingaliquot/fullmoleculesunassembled,notabsence.');g.sequence(lib,SequencingPlan(limitation='PublishedPACNovaSeq6000PE150. SourceRNA/cDNARead1/UMI+spatialRead2overview; completeproducedmolecules/primersites/runlandingsunassembled. DistinctDCbranch/experimentalindicesnotinherited.'))
    g.add_source_item(SourceItem('unassembled','ExactwholecDNA/TSOsealedproducts/allchannelallocations/vendorSEQuoiadepletion/actualP5familymapping andrunlandingunassembled.','LocalA/Bsplintedinterfacesknown;notwholeligationorcompletechemistryproof.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    import re
    c=Check();g=protocol_graph();c('PAPbeforeRT',list(g.steps).index('polyadenylate')<list(g.steps).index('rt'),True);c('RTbiotinTphysicalnucleotide',sum(len(p.top)for p in rt_parts()),43);c('TSOstepnotfreshprimingRT',g.steps['template_switch'].kind,S.TEMPLATE_SWITCHING);c('Actual110to100lysateremainder','lysate_10'in g.terminals,True);c('FouractualPCRtubes',len(g.steps['four_tubes'].outputs),4);c('No inventedPCRpooltiming',g.steps['pcr_fractionate'].kind,S.MAGIC);c('Bothsourcefractionlibraries',set(g.sequencing),{'wta_library','sgrna_library'});c('Firstsplintlocalpairingcomputes',len(end_scene('A').strands),3);c('Secondsplintlocalpairingcomputes',len(end_scene('B').strands),3)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'s41587-026-03127-y_41587_2026_3127_MOESM4_ESM.xlsx.txt'
    if have(c,Source(p,'https://doi.org/10.1038/s41587-026-03127-y')):
        t=re.sub(r'\s+','',p.read_text())
        for k,v in SOURCE_ORDERINGS.items():c('Originalwholeordering '+k,re.sub(r'\s+','',v)in t,True)
    c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
