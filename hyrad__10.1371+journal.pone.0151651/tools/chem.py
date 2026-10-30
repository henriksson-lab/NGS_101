"""Private partial hyRAD donor-probe factory and butterfly target-library capture."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'lib'))
from chemdraw import Scene,Segment,feature,revcomp
from illumina import P5,P7
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem,SequencingPlan
AUTHORITY='PublishedPMC4801390 main149–250/Table3; actualsameoriginalPDFp7/p8. Closed1original; no updatedGitHub/vendor/predecessor recipe inherited.'
ORDERS={'RAD-P1.1': 'ACACTCTTTCCCTACACGACGCTCTTCCGATCTxxxxxxCCTGCA', 'RAD-P1.2': 'GGxxxxxxAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT', 'RAD-P2.1': 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT', 'RAD-P2.2': 'TAAGATCGGAAGAGCGAGAACAA', 'P1.1': 'ACACTCTTTCCCTACACGACGCTCTTCCGATCTxxxxxx', 'P1.2': 'xxxxxxAGATCGGAAGAGC', 'P2': 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCCCCC', 'ILLPCR1': 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTT', 'ILLPCR2_01': 'CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGC', 'ILLPCR2_02': 'CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCAGACGTGTGC', 'ILLPCR2_03': 'CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCAGACGTGTGC', 'ILLPCR2_04': 'CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCAGACGTGTGC', 'ILLPCR2_05': 'CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCAGACGTGTGC', 'ILLPCR2_06': 'CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCAGACGTGTGC', 'ILLPCR2_07': 'CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCAGACGTGTGC', 'ILLPCR2_08': 'CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCAGACGTGTGC', 'ILLPCR2_09': 'CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCAGACGTGTGC', 'ILLPCR2_10': 'CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCAGACGTGTGC', 'ILLPCR2_11': 'CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCAGACGTGTGC', 'ILLPCR2_12': 'CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCAGACGTGTGC', 'BO1.P5.F': 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT', 'BO2.P5.R': 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT', 'BO3.P7.F': 'CAAGCAGAAGACGGCATACGAGATIIIIIIGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT', 'BO4.P7.R': 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIATCTCGTATGCCGTCTTCTGCTTG'}
SOURCE_OLIGO_EXCEPTIONS={k:'PublishedTable3actualPDFp7 prints six literalx positions for an unassigned inline samplebarcode; generic source guard normalizes modeled window toN6 but source parser doesnot retainx6. Whole printed sourceorder fixture plus actual Table3 originalimage verify allflanks/windowwidth. Narrow sourceformattingexception only; noselectedbarcode or molecularproduct asserted.'for k in ('rad_p1_a','rad_p1_b','target_p1_a','target_p1_b')}
INDEXES=tuple(ORDERS['ILLPCR2_'+str(i).zfill(2)][len(P7):len(P7)+6]for i in range(1,13))
def segments(name):
    text=ORDERS[name]
    if 'xxxxxx'in text:
        a,b=text.split('xxxxxx');return([Segment('source5arm',a)]if a else[])+[Segment('samplebarcode','xxxxxx',placeholder=True,feature=feature('hyrad.inline_sample','sample_index','unknown',note='SourceX6inlinebarcode, assignedsampleunknown; notUMI'))]+([Segment('source3arm',b)]if b else[])
    if name.startswith('ILLPCR2_'):
        return[Segment('canonicalP7',P7),Segment('PCRsampleindex',text[len(P7):len(P7)+6],feature=feature('hyrad.pcr_sample','sample_index','whitelist',whitelist='PublishedTable3ILLPCR2_01–12')),Segment('sourcePCR3body',text[len(P7)+6:])]
    if 'IIIIII'in text:
        a,b=text.split('IIIIII');return[Segment('source5body',a),Segment('sixphysicalinosines','IIIIII',placeholder=True,length_bp=6,note='SourceI=inosine nucleotide, notUMI/randomN'),Segment('source3body',b)]
    if text.startswith(P5):return[Segment('canonicalP5',P5),Segment('wholeprintedsource3body',text[len(P5):])]
    return[Segment('wholeprintedsourcebody',text)]
def oligo_scene(name):
    sc=Scene();sc.strand('oligo',segments(name));return sc
def seq(sc):return ''.join(s.top for st in sc.strands.values()for s in st.segs)
def adapter_scene(which='rad_p1'):
    sc=Scene()
    if which=='rad_p1':
        top=ORDERS['RAD-P1.1'];bot=ORDERS['RAD-P1.2'];n=len(top.split('xxxxxx')[0])
        sc.strand('top',[Segment('pairedprefix',top[:n]),Segment('barcode','xxxxxx',placeholder=True,feature=feature('hyrad.inline_sample','sample_index','unknown')),Segment('pairedCC','CC'),Segment('3TGCA','TGCA')])
        sc.anneal('bottom',[Segment('pairedGG','GG'),Segment('barcode','xxxxxx',placeholder=True,feature=feature('hyrad.inline_sample','sample_index','unknown')),Segment('prefixcomplement',bot[8:])],to='top',pair=('prefixcomplement','pairedprefix'))
        sc.ligation_end('top','3TGCA',kind='cohesive',polarity='3-prime',text='SourceSbfI-compatible terminal3TGCA; exactoccupiedDNA/adaptor5Pstatesunassembled')
    elif which=='rad_p2':
        top=ORDERS['RAD-P2.1'];bot=ORDERS['RAD-P2.2'];core=revcomp(bot[2:15])
        if not top.endswith(core):raise ValueError('SourceRADP2pairedcoredoesnotmatch')
        sc.strand('top',[Segment('Y5arm',top[:-13]),Segment('pairedcore',core)])
        sc.anneal('bottom',[Segment('5TA','TA'),Segment('corecomplement',bot[2:15]),Segment('Y3arm',bot[15:])],to='top',pair=('corecomplement','pairedcore'),unpaired=('Y3arm',))
        sc.ligation_end('bottom','5TA',kind='cohesive',polarity='5-prime',text='SourceMseI-compatible5TA; wholeligationproductunassembled')
    elif which=='shotgun_p1':
        top=ORDERS['P1.1'];bot=ORDERS['P1.2'];core=revcomp(bot[6:])
        sc.strand('top',[Segment('Y5arm',top[:-19]),Segment('pairedcore',core),Segment('barcode','xxxxxx',placeholder=True,feature=feature('hyrad.inline_sample','sample_index','unknown'))])
        sc.anneal('bottom',[Segment('barcode','xxxxxx',placeholder=True,feature=feature('hyrad.inline_sample','sample_index','unknown')),Segment('corecomplement',bot[6:])],to='top',pair=('corecomplement','pairedcore'))
        sc.ligation_end('top','barcode',kind='blunt',text='Sourceadapterpairedbluntclass; sourceDNA5P/bluntpreparationknown, wholeendmods/productunassembled')
    else:raise ValueError('Selectsourceadapter')
    return sc
def unknown_scene():
    sc=Scene();sc.strand('unassembled',[Segment('wholeoccupiedmolecules','[unassembled]',placeholder=True)]);return sc
def index_scene(index=None):
    if index is not None and index not in INDEXES:raise ValueError('Indexmustbesource-listed')
    sc=Scene();body=ORDERS['ILLPCR2_01'];sc.strand('oligo',[Segment('canonicalP7',P7),Segment('sampleindex',index or'NNNNNN',placeholder=index is None,feature=feature('hyrad.pcr_sample','sample_index','whitelist',whitelist='PublishedTable3ILLPCR2_01–12')),Segment('sourcePCR3body',body[len(P7)+6:])]);return sc

def protocol_graph():
    g=ProtocolGraph('hyRAD sourcebutterfly intacttarget/probefactory preparation',authority=AUTHORITY,representation='published_instance')
    x=g.input('fresh_donor','ActualfreshdonorgDNAforRADprobes; exactgenome/individualunassigned',K.BIOLOGICAL,primary=True)
    def inp(k,label,kind=K.CONDITION,role=R.CONDITION,**kw):return g.input(k,label,kind,**kw),role
    def op(k,label,ports=(),kind=S.MAGIC,product='Actualsourceheterogeneousmixture; wholeproducts/endclasses/occupancy/recoveryunassembled',**kw):
        nonlocal x
        y=g.product(k+'_product',product,K.MOLECULAR)
        if k in ('rad_pcr','rad_other_pcr','target_pcr','capture_pcr'):
            kind=S.PCR;kw['geometry']=PCRGeometry.unavailable(unknown_scene(),reason='Wholeprintedprimersretained; occupiedILLPCR1possiblegappedlandingandcompleteampliconsunassembled.')
        if kind==S.MAGIC:kw.setdefault('reason','Positive sourceoperation; fulloccupiedmolecularproducts/endchemistry/yieldsunassembled, noabsenceclaim.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        g.react(k,label,kind,((x,R.TEMPLATE),)+tuple(ports),(y,),**kw);x=y
    def clean(k,label):
        nonlocal x
        y=g.product(k+'_kept','Actualretainedphysicalfraction; recovery/carryoverunassigned',K.MOLECULAR);w=g.product(k+'_waste','Actualremovedworkupfraction',K.WASTE)
        g.react(k,label,S.CLEANUP,((x,R.TEMPLATE),inp(k+'_conditions',label)),(y,w));g.discard(w,'Actualworkupfractions; no perfectpurityclaim');x=y
    def alloc(k,label,rest='Untransferredstock; quantity/fateunassigned',qc=False):
        nonlocal x
        y=g.product(k+'_used',label,K.MOLECULAR);r=g.product(k+'_rest',rest,K.MOLECULAR)
        g.react(k,label,S.SPLIT,((x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(y.id if qc else r.id,Disposition.QC if qc else Disposition.UNRESOLVED,label if qc else rest));x=r if qc else y
    def supplies(k,label,spec,prepared=False):
        st=g.input(k,label,K.MOLECULAR);out=[];ports=[]
        for item in spec:
            suffix,kind,role,text=item[:4];m=g.product(k+'_'+suffix,text,kind,state=item[4]if len(item)>4 else None);out.append(m);ports.append((m,role))
        g.react(k+'_components','CopresentcomponentsofONEactualpreparedstock'if prepared else'Actualrecipeinventorygrouping; packaging/stockcountunassigned',S.GENERIC,((st,R.TEMPLATE),),tuple(out),changes_molecule=False);return tuple(ports)
    def anneal(k,one,two,which):
        a=g.input(k+'_a','SourceTable3 '+one,K.OLIGO,state=oligo_scene(one));b=g.input(k+'_b','SourceTable3 '+two,K.OLIGO,state=oligo_scene(two));p=g.product(k+'_stock','Actualannealedsourceadapterstock; occupancynotperfect',K.ADAPTER,state=adapter_scene(which))
        g.react(k+'_anneal','Source95Cheat/gradualcoolannealing',S.HYBRIDIZATION,((a,R.TEMPLATE),(b,R.STRUCTURAL),inp(k+'_conditions','Sourceannealingconditions; amount/untransferredstockunassigned')),(p,))
        used=g.product(k+'_used','Actualadapterallocation; samestockcomponentclasses, amountperselectedreaction',K.ADAPTER);rest=g.product(k+'_rest','Unusedannealedadapterstock; amountunassigned',K.ADAPTER)
        g.react(k+'_split','Allocateactualannealedstock, notrepeatfreshstrands',S.SPLIT,((p,R.TEMPLATE),),(used,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Unusedannealedadapterstock'));return used
    def pcr_supplies(k,label,prepared=False):
        return supplies(k,label,[('f',K.OLIGO,R.STRUCTURAL,'ActualwholeprintedILLPCR1, no canonicalrepair',oligo_scene('ILLPCR1')),('r',K.OLIGO,R.STRUCTURAL,'ActualsourceILLPCR2family; studyindexunassigned',index_scene()),('enzyme',K.CATALYST,R.CATALYST,'ActualsourceQ5HotStartenzymecomponent'),('nts',K.MOLECULAR,R.STRUCTURAL,'ActualsourcedNTPcomponent'),('buffer',K.CONDITION,R.CONDITION,'ActualsourceQ5buffer/remainingconditions')],prepared)
    alloc('donor_input','Actual6uldonorgDNAfor9uldigest')
    op('donor_digest','Source9ul/1UMseI/2USbfI-HF/CutSmart37C3h',(inp('msei','Actual1UMseI',K.CATALYST,R.CATALYST),inp('sbfi','Actual2USbfI-HF',K.CATALYST,R.CATALYST),inp('donor_buffer','Actual1xCutSmart')))
    clean('donor_clean','Source2xAMPure/elute10ulTris')
    p1=anneal('rad_p1','RAD-P1.1','RAD-P1.2','rad_p1');p2=anneal('rad_p2','RAD-P2.1','RAD-P2.2','rad_p2')
    op('rad_ligate','Source20ul/10DNA/.5uMeachRADadapter/400UT4ligase16C3h',((p1,R.STRUCTURAL),(p2,R.STRUCTURAL),inp('rad_ligase','Actual400UT4DNAligase',K.CATALYST,R.CATALYST),inp('rad_lig_buffer','Actual1xT4ligasebuffer')),geometry=AdapterLigationGeometry.known('cohesive',adapter_scene('rad_p1'),adapter_scene('rad_p2'),statement='Table3sourceSbfI3TGCA/MseI5TAadapterendclasses; fullrestrictionfragmentmolecules/barcoderegisters/5Pstates/sealedproductsunassembled.'))
    clean('rad_clean','Source1xAMPure/Triselution')
    clean('rad_size','SourcePippinPrep270bptightselection; noperfectsizepurity')
    first_template=g.product('rad_first_template','Actualtemplatealiquotforone10ulRADPCRreaction; volumeunassigned',K.MOLECULAR)
    other_templates=g.product('rad_other_templates','ActualtemplatealiquotsforOTHERseparateRADPCRreactions; number/volumesempiricalandunassigned',K.MOLECULAR)
    unused_template=g.product('rad_unused_template','AnyunallocatedRADtemplate stock; existence/amountunassigned',K.MOLECULAR)
    g.react('rad_template_split','AllocateactualdistincttemplatecontributionsBEFOREseparatePCRreactions; nofixedreplicatecount',S.SPLIT,((x,R.TEMPLATE),),(first_template,other_templates,unused_template));g.mark_terminal(Terminal(unused_template.id,Disposition.UNRESOLVED,'AnyunallocatedRADtemplateisnotautomaticallyamplified'));x=first_template
    op('rad_pcr','SourceONE10ultotalPCR/.2mMdNTP/.6uMeachprimer/.2UQ5/30cycles98C20s60C30s72C40s',pcr_supplies('rad_pcr_supplies','ActualsuppliesforONEfinal10ulPCRrecipe; packagingunassigned'))
    alloc('rad_qubit','ActualQubitQCaliquotfromfirstPCR; exactamountunassigned',qc=True)
    alloc('rad_fragment_qc','ActualFragmentAnalyzerQCaliquotfromfirstPCR; exactamountunassigned',qc=True)
    first_product=x;x=other_templates
    op('rad_other_pcr','SourceOTHERactualseparate10ulPCRreactions: each.2mMdNTP/.6uMprimers/.2UQ5/same30cycleprogramme; reactioncountempirical',pcr_supplies('rad_other_pcr_supplies','ActualsuppliesforOTHERseparatePCRreactionarray; exactnumber/totalstockquantitiesunassigned'),product='ActualOTHERindependentlyamplifiedRADPCRproducts fromtheirallocatedtemplates; count/yields/wholeproductsunassembled')
    alloc('rad_other_qubit','ActualQubitQCaliquotsfromotherPCRproducts; exactamounts/countunassigned',qc=True)
    alloc('rad_other_fragment_qc','ActualFragmentAnalyzerQCaliquotsfromotherPCRproducts; exactamounts/countunassigned',qc=True)
    other_product=x
    pool=g.product('rad_pooled','ActualpooledDISTINCTRADPCRreactionproducts',K.MOLECULAR);g.react('rad_pool','PhysicallypooldistinctPCRproductsAFTERtheiractualQCaliquots',S.POOL,((first_product,R.TEMPLATE),(other_product,R.TEMPLATE)),(pool,));x=pool
    clean('rad_pool_clean','Source1xAMPureafterreplicatepool')
    ref=g.product('rad_reference_aliquot','ActualRADlibraryaliquotforMiSeqSE300reference; amountunassigned',K.MOLECULAR);probe=g.product('rad_probe_rest','ActualremainingRADlibrarystockforprobefactory',K.MOLECULAR)
    g.react('rad_branches','Actualaliquotsequenced/restconvertedtoprobes',S.SPLIT,((x,R.TEMPLATE),),(ref,probe));ref_lib=g.product('rad_reference_library','ActualsourceRADreferenceNGSlibrary; wholemoleculesunassembled',K.LIBRARY);g.react('rad_reference_load','SourcebutterflyRADreferenceMiSeqSE300',S.GENERIC,((ref,R.TEMPLATE),),(ref_lib,),changes_molecule=False);g.sequence(ref_lib,SequencingPlan(limitation='OptionalsourceRADreferenceMiSeqSE300branch only; wholeprimers/readlandingunassembled.'));x=probe
    alloc('probe_redigest_input','Actual30ulofremainingRADstockper50ulredigest; replicatecountunassigned')
    op('adapter_strip','Source50ul/5UMseI/10USbfI37C3h65C20min',(inp('strip_msei','Actualfresh5UMseI',K.CATALYST,R.CATALYST),inp('strip_sbfi','Actualfresh10USbfI',K.CATALYST,R.CATALYST),inp('strip_conditions','ActualCutSmart/sourceconditions')))
    clean('strip_clean','Source1.5xAMPure')
    bn=supplies('bionick_supplies','ActualBioNickDNALabelingSysteminventory; packaging/formulationunassembled',[('enzyme',K.CATALYST,R.CATALYST,'PositiveDNA-nicklabelingenzymeclass; exactenzymes/dosesunassembled'),('label',K.MOLECULAR,R.STRUCTURAL,'Actualbiotinlabelingnucleotideclass; exactnucleotide/linkage/amountunassembled'),('other',K.MOLECULAR,R.STRUCTURAL,'RemainingactualBioNicksuppliesunassembled')])
    op('bionick','SourceBioNickbiotinnick-labelDNAbaits, notRNAbaits',bn)
    clean('probe_clean','Source1.5xAMPurebiotinDNAprobes')
    alloc('probe_capture_allocation','Actual500–1000ngbiotinDNAprobesforcapture');bait=x
    target=g.input('target_gdna','SeparatesourcemuseumbutterflygDNA; exactDNAenddamages/identityunassigned',K.BIOLOGICAL,primary=True)
    intact=g.product('target_intact','Actualsourceintact/nonsonicatedaliquotselected',K.MOLECULAR);sonic=g.product('target_sonication_input','ActualotheraliquotforparallelCovaris300bpbranch',K.MOLECULAR)
    g.react('target_split','SourcesplitsindividualDNAintointactandsonicatedparallelpreparations',S.SPLIT,((target,R.TEMPLATE),),(intact,sonic));sonicated=g.product('other_sonicated','ActualotherCovaris~300bpaliquot; parallelfinishingoutside selectedtrace',K.MOLECULAR);g.react('other_sonication','SourceCovaris300bpbranch',S.MAGIC,((sonic,R.TEMPLATE),inp('covaris','Actualmechanicalsonicationconditions')),(sonicated,),reason='Positiveparallelbranch, exactendstates/otherlibraryconstructionoutside selectedintacttrace.');g.mark_terminal(Terminal(sonicated.id,Disposition.UNRESOLVED,'OtheractualsonicatedDNAstocknotmergedintoselectedintactlibrary'));x=intact
    alloc('target_input','Actual8ulDNAforPNK')
    op('target_pnk','Source95C10minice/10ulPNKreaction10U37C30min65C20min',(inp('pnk','Actual10UT4PNK',K.CATALYST,R.CATALYST),inp('pnk_conditions','Actual1xPNKbuffer; sourceexplicitly5primephosphorylatesDNA')),product='Source5P-treatedDNAclassplusunreactedclasses; exactoccupiedendsunassembled')
    clean('pnk_clean','Source2xAMPure/elute10ulTris')
    op('tail','Source95C10minice/20ulTdT10U/4mMGTP/.25mMCoCl2/NEB4/37C30min70C10min',(inp('tdt','Actual10UTdT',K.CATALYST,R.CATALYST),inp('gtp','Actual4mMGTPASPRINTED; tailribose/deoxyriboseformunassembled, nodGTPsubstitution',K.MOLECULAR,R.STRUCTURAL),inp('tail_conditions','ActualNEB4/cobaltconditions')),product='Sourceguanine-tailedDNAclassinCARRIED20ulTdTmix; exacttailresidues/sugar/lengthunassembled')
    second=supplies('second_mix','ONEactual10ulsecondstrandmixaddedtoENTIRE20ulTdTmix',[('primer',K.OLIGO,R.STRUCTURAL,'Actual1ulP2atprinted15mM',oligo_scene('P2')),('polymerase',K.CATALYST,R.CATALYST,'Actual2ul5UulKlenowexo-minus'),('nts',K.MOLECULAR,R.STRUCTURAL,'Actual.6ul25mMeachdNTP'),('buffer',K.CONDITION,R.CONDITION,'Actual1ul10xNEB4/5.4ulwater')],prepared=True)
    op('second_strand','Sourcecarry20ulTdTmix+10ulsecondmix/23C3h75C20min',second)
    blunt=supplies('blunt_mix','Actual5ulbluntingmixaddedtocarriedsecondstrandreaction',[('polymerase',K.CATALYST,R.CATALYST,'Actual.2ul3UulT4DNApolymerase'),('bsa',K.MOLECULAR,R.STRUCTURAL,'Actual.35ul10mgmlBSAcarrier'),('conditions',K.CONDITION,R.CONDITION,'Actual.5ul10xNEB4/3.95ulwater')],prepared=True)
    op('blunt','SourceT4polymeraseblunting12C15min',blunt,product='Source5P/blunt-treatedDNAclassplusotherclasses; whole5Pdistribution/productstrandsunassembled')
    clean('blunt_clean','Source2xAMPure/elute10ulTris')
    tp1=anneal('target_p1','P1.1','P1.2','shotgun_p1')
    op('target_ligate','Source20ul/10DNA/1ul25uMP1/400UT4ligase16Covernight',((tp1,R.STRUCTURAL),inp('target_ligase','Actualfresh400UT4DNAligase',K.CATALYST,R.CATALYST),inp('target_lig_buffer','Actual1xT4ligasebuffer')),geometry=AdapterLigationGeometry.known('blunt',adapter_scene('shotgun_p1'),statement='Source5P-treatedbluntDNArecipientandlocallypairedP1endclass; adapterchemicalendmods/fullsealedproductunassembled.'))
    clean('target_lig_clean','Source1xAMPure/elute20ulTris')
    op('bst_fill','Source40ul/20DNA/12UBst/.25mMdNTPeach/ThermoPol37C20min',(inp('bst','Actual12UBstpolymerase',K.CATALYST,R.CATALYST),inp('bst_nts','Actual.25mMdNTPeach',K.MOLECULAR,R.STRUCTURAL),inp('bst_conditions','Actual1xThermoPolbuffer')))
    alloc('target_pcr_input','Actual10ulof40ulBstfillstock','Actualnominal30ulBstfillstocknottransferred; fateunassigned')
    op('target_pcr','Source10ultemplate+15ulpreparedmix:5buffer/.2dNTP/2.5primerpair/.5UQ5;25cycles98C10s60C20s72C25s',pcr_supplies('target_pcr_mix','ONEactual15ulpreparedQ5PCRmix; exactlistedcomponentvolumes/balanceunassigned',True))
    alloc('target_gel_qc','ActualgelQCaliquotofPCR; amountunassigned',qc=True)
    clean('target_pcr_clean','Source.7xAMPure')
    alloc('target_pool_allocation','Actualequimolarsamplelibrarycontribution; amountunassigned')
    otherlibs=g.input('other_sample_libraries','Otheractualindependentlypreparedshotgunlibrarycontributions; quantities/identitiesunassigned',K.MOLECULAR)
    lp=g.product('shotgun_pool','Actualequimolarshotgunlibrarypool',K.MOLECULAR);g.react('shotgun_pooling','Physicallypoolactualsamplelibrarycontributions',S.POOL,((x,R.TEMPLATE),(otherlibs,R.TEMPLATE)),(lp,));x=lp
    alloc('capture_library_allocation','Actual500–1000ngshotgunlibrarypoolforcapture')
    blocks=[]
    for n in ('BO1.P5.F','BO2.P5.R','BO3.P7.F','BO4.P7.R'):
        blocks.append(inp(n.lower().replace('.','_'),'Actualsource2uM '+n,K.OLIGO,R.STRUCTURAL,state=oligo_scene(n),chemical_state={'source_ordering':ORDERS[n],'physical_inosine_positions':tuple(i for i,b in enumerate(ORDERS[n])if b=='I'),'I_meaning':'Sourcephysicalinosinenucleotide; Sceneplaceholderflagpermitsunsupportedletterdisplay, notuncertainsequence/randomUMI'}))
    op('hybridize','Source40ul6xSSC/50mMEDTA/1percentSDS/2xDenhardt/2uMeachBO/500–1000ngbaits+library/95C10min65C48h',((bait,R.STRUCTURAL),)+tuple(blocks)+(inp('hybridization_conditions','Actualbutterflycaptureconditions; grasshopperCot1notselected'),),kind=S.HYBRIDIZATION)
    bound=g.product('capture_bound','ActualM280-boundbiotinDNAprobe/librarycomplexes; incomplete/offtargetclassesunassigned',K.MOLECULAR);free=g.product('capture_supernatant','Actualremovedhybridizationsupernatant',K.WASTE)
    g.react('capture','Source10ulM280washed3x200TEN/resuspend200TEN/add40hybmix/RT30min',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),inp('m280','ActualpreparedM280beads',K.SELECTION,R.SELECTION)),(bound,free));g.discard(free,'Actualsourceunboundfractionremoved');x=bound
    clean('capture_wash','Source4x200ul65C15min:1xSSC.1SDS twice/.5xonce/.1xonce')
    eluate=g.product('target_eluate','Actual30ul80CwatereluatecontainsTARGETLIBRARYwashedOFFprobes',K.MOLECULAR);resid=g.product('probe_bead_residue','Actualprobe/beadresidue; laterfateunassigned',K.MOLECULAR)
    g.react('elute','Source30ul80Cwater10min targetlibraryrelease',S.MAGIC,((x,R.TEMPLATE),inp('elution_conditions','Actualsourcehotwaterelution')),(eluate,resid),reason='Sourcetargetlibraryelutesoffprobes; notclaimedbiotin-streptavidinbondcleavage, wholehybrid/offproductfatesunassembled.');g.mark_terminal(Terminal(resid.id,Disposition.UNRESOLVED,'Residualprobe/beadstocknotcapturePCRinput'));x=eluate
    alloc('capture_pcr_input','Actual15ulof30ulcaptureeluate','Actualremainingnominal15ulcaptureeluate; fateunassigned')
    op('capture_pcr','Source50ulPCR/15template/1UQ5/.2mMdNTP/.5uMeachprimer/25cycles98C10s60C20s72C25s',pcr_supplies('capture_pcr_supplies','Actualfinal50ulQ5PCRrecipeinventory; packagingunassigned'))
    clean('capture_pcr_clean','Source1xAMPure')
    alloc('final_pool_allocation','Actualequimolarcapture-librarycontribution; amountunassigned')
    others=g.input('other_capture_libraries','Otheractualcapture-enrichedlibrarycontributions; identities/amountsunassigned',K.MOLECULAR);fp=g.product('final_pool','Actualequimolarcapture-librarysequencingpool',K.MOLECULAR);g.react('final_pooling','Sourceequimolarlibrarypoolforsequencing',S.POOL,((x,R.TEMPLATE),(others,R.TEMPLATE)),(fp,));x=fp
    alloc('loading','ActualMiSeqloadingallocation; volumeunassigned')
    lib=g.product('captured_library','ActualbutterflyhyRADcapturelibraryNGSmaterial; wholegeometryunassembled',K.LIBRARY);g.react('load','SourcebutterflyMiSeqPE150; grasshopperPE100separate',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='MiSeqPE150butterflycapture. Exactwholeamplicons/ILLPCR1possiblegapped3primelanding/readsites/barcodecyclesunassembled; no primerbodyrepair orUMIassignment.'))
    g.add_source_item(SourceItem('source_literals','Main234/PDFp8prints4mMGTP;main236/PDFp8printsP2primer15mM;Table3wholeILLPCR1asprinted.','SourcevaluesretainedwithoutdGTP/uM/canonicalprimersubstitution; tailactualsugar/concentrationinterpretation/wholeprimerlandingunassembled.'))
    g.add_source_item(SourceItem('scope','Published2016selectedbutterflyintacttarget+RADprobefactory;twoactualsequencingbranches.','No2016equivalencetoupdatedGitHubchiasto/hyRAD; BioNick/vendor/predecessororiginalsunread. X6/indexessamplemeaning, BOsixIphysicalinosinesnotUMI.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();p=sources_root()/HERE.parent.name/'journal.pone.0151651_PMC4801390.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1371/journal.pone.0151651')):
        text=p.read_text();t=''.join(text.split())
        for name,order in ORDERS.items():c('WholeTable3ordering '+name,order in t,True)
        c('SourceGTPnotnormalized','4 mM GTP'in text,True)
        c('SourceP2literal15mM','P2 oligonucleotide (15 mM)'in text,True)
        c('SourceIphysicalinosine','I = inosine'in text,True)
    for which,top,bottom in [('rad_p1','RAD-P1.1','RAD-P1.2'),('rad_p2','RAD-P2.1','RAD-P2.2'),('shotgun_p1','P1.1','P1.2')]:
        sc=adapter_scene(which);c('Adapterwholeorders '+which,tuple(''.join(s.top for s in st.segs)for st in sc.strands.values()),(ORDERS[top],ORDERS[bottom]))
    try:index_scene('AAAAAA');bad=False
    except ValueError:bad=True
    c('Rejectunlistedindex',bad,True)
    c('Unknownstudyindexnotdegeneratesynthesis',index_scene().strands['oligo'].segs[1].placeholder,True)
    c('SecondstrandcarriesTdTmixture',g.steps['second_strand'].inputs[0].material,'tail_product')
    c('Secondstrandprimerfromsame10ulmix',g.producer_of('second_mix_primer').id,'second_mix_components')
    c('BluntBSAfromsameactualmix',g.producer_of('blunt_mix_bsa').id,'blunt_mix_components')
    c('TargetPCRuses10of40allocation',g.steps['target_pcr'].inputs[0].material,'target_pcr_input_used')
    c('Bstremainingstockretained',g.terminals['target_pcr_input_rest'].disposition,Disposition.UNRESOLVED)
    c('Probeconversionusesrestnotreference',g.steps['probe_redigest_input'].inputs[0].material,'rad_probe_rest')
    c('OtherPCRhasdistinctallocatedtemplate',g.steps['rad_other_pcr'].inputs[0].material,'rad_other_templates')
    c('FirstPCRhasdistinctallocatedtemplate',g.steps['rad_pcr'].inputs[0].material,'rad_first_template')
    c('PooltakesdistinctPCRproductsafterQC',set(p.material for p in g.steps['rad_pool'].inputs),{'rad_fragment_qc_rest','rad_other_fragment_qc_rest'})
    c('Actualreplicatesphysicallypooled',g.steps['rad_pool'].kind,S.POOL)
    c('CapturePCRuses15of30allocation',g.steps['capture_pcr'].inputs[0].material,'capture_pcr_input_used')
    c('Targetelutesoffbeadsnotprobesubstitution',g.steps['capture_pcr_input'].inputs[0].material,'target_eluate')
    c('ResidualprobesnotPCRinput',g.terminals['probe_bead_residue'].disposition,Disposition.UNRESOLVED)
    c('Parallelsonicatedbranchseparate',g.terminals['other_sonicated'].disposition,Disposition.UNRESOLVED)
    c('TwoactualNGSbranches',set(g.sequencing),{'rad_reference_library','captured_library'})
    c('Blockerinosinehasnomoleculeidentifier',next(s for s in oligo_scene('BO3.P7.F').strands['oligo'].segs if s.name=='sixphysicalinosines').feature,None)
    c('ActualI6physicaldescriptor',len(g.materials['bo3_p7_f'].chemical_state['physical_inosine_positions']),6)
    c('Wholecapturelibraryunassembled',g.materials['captured_library'].state,None)
    c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
