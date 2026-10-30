"""SPACE-seq associated-author-manuscript selected Tigre route; private partial model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature
from illumina import TRUSEQ_READ1
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
TITLE='SPACE-seq: associated-manuscript selected Tigre spatial library'
AUTHORITY='PMC13401910.1/NIHMS2190436 associated-author-manuscript Methods246–254/analysis264, DOI10.1016/j.stem.2026.05.017. Final published body/SI not verified; no complete-version equivalence. One selected original.'
PRINTED_BEAD='TTT_PC_GCCGGTAATACGACTCACTATAGGGCTACACGACGCTCTTCCGATCTJJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJNNNNNNNNVT30'
PREFIX='GCCGGTAATACGACTCACTATAGGG'+TRUSEQ_READ1[11:]
SPACER='TCTTCAGCGTTCCCGAGA'
NOTICE='WITHIN ASSOCIATED AUTHOR MANUSCRIPT: Methods246 printed bead has8N but prose says7nt UMI; Methods246 prints8+6J, while analysis264 prints8+7J. Local inventory preserves Methods246 literal8N/8+6J. Final-body/SI identity not verified; this is not a claim of a final-publication error.'

def inventory():
    sc=Scene();sc.strand('bead',[Segment('source_start','TTT'),Segment('photocleavable_linker','[PC]',placeholder=True),Segment('source_handle',PREFIX),Segment('spatial_barcode_first','J'*8,placeholder=True,feature=feature('space_bead','spatial_barcode','unknown')),Segment('source_split_barcode_spacer',SPACER),Segment('spatial_barcode_second','J'*6,placeholder=True,feature=feature('space_bead','spatial_barcode','unknown')),Segment('printed_UMI','N'*8,placeholder=True,feature=feature('space_umi','umi','random'),note='Literal8N; prose7nt conflictingwithinassociatedMS'),Segment('pre_polyT_anchor','V',placeholder=True,note='PrintedVprecedesT30;notaddedtoUMI'),Segment('polyT','T'*30)],mod5='bead attachment; detailed linker chemistry unassembled')
    return sc

def unknown():
    sc=Scene();sc.strand('unknown',[Segment('unknown','[unassembled]',placeholder=True)]);return sc

def capture_end_inventory():
    sc=Scene();sc.strand('capture_oligo',[Segment('capture_sequence','[CS1/CS2/CS2L unassembled]',placeholder=True)],mod5='P')
    sc.strand('splint',[Segment('splint_sequence','[Splint Oligo unassembled]',placeholder=True)],label='Independent source splint inventory; no pairing assigned')
    return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance');x=g.input('section','Source10umcryosectionfromengineeredDARLINtissue; genotyping/editingupstreamnotassembled',K.BIOLOGICAL,primary=True)
    def terminal(m,why,disp=Disposition.STORED):g.mark_terminal(Terminal(m.id,disp,why))
    def clean(key,label):
        nonlocal x
        b=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','Sourcephysicallyretainedmaterial;recovery/carryoverunassigned',K.MOLECULAR);w=g.product(key+'_waste','Actualsourcewash/cleanupunretainedfraction',K.WASTE);g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(b,R.CONDITION)),(y,w));g.discard(w,'Physicalfractionnotperfectmolecularselection');x=y
    def mix(key,label,primers=(),enz=True):
        stock=g.input(key+'_bulk',label,K.MOLECULAR);b=g.product(key+'_buffer','Buffercomponentofsamebulkaddition',K.CONDITION);parts=[b];ps=[];cat=None
        for i,p in enumerate(primers):
            m=g.product(key+'_primer_'+str(i),p+';samebulkcomponent;wholeorderingunassembled',K.OLIGO);parts.append(m);ps.append(m)
        if enz:cat=g.product(key+'_enzyme','Enzymecomponentofsamebulkaddition;sourceidentityinbulklabel',K.CATALYST);parts.append(cat)
        g.react(key+'_components','Representonebulkadditioncomponents;notphysicalseparation',S.GENERIC,((stock,R.TEMPLATE),),parts,changes_molecule=False)
        return b,ps,cat
    def pcr(key,label,primernames):
        nonlocal x
        b,ps,e=mix(key,'SourceKAPAHiFiHotStartReadyMix079589350001/actualPCRprimerpair',primernames)
        y=g.product(key+'_product','Sourceampliconreaction;wholecovalentproduct/primerlandingunassembled',K.MOLECULAR);g.react(key,label,S.PCR,((x,R.TEMPLATE),(b,R.CONDITION),(e,R.CATALYST),*((p,R.STRUCTURAL)for p in ps)),(y,),geometry=PCRGeometry.unavailable(unknown(),reason='Sourcepairnamesknown;TableS1gated,wholeordering/landingunassembled.'));x=y
    puck=g.input('puck','SourcepredecodedSlide-seqbarcodedbeadpuck;fullpuckassembly/versionunassembled',K.MOLECULAR,state=inventory(),chemical_state={'literal_ordering':PRINTED_BEAD,'PC_after_first_bases':3,'scope':NOTICE})
    cs=g.input('capture_oligos','Source5primephosphorylatedCS1/CS2/CS2Lmixture;wholeorderingunassembled',K.OLIGO,chemical_state={'five_prime':'P'})
    spl=g.input('splint','SourceSplintOligo;wholeorderingunassembled',K.OLIGO);lig=g.input('capture_ligase','SourceT4DNAligase37C1hM0202',K.CATALYST);targetpuck=g.product('targeted_puck','Sourcecapture-modifiedpuckwithretainedpolyTandtarget-specificcapturepopulation;fractionunassigned',K.MOLECULAR)
    g.react('modify_puck','AssociatedMS248sourcecaptureligations',S.LIGATION,((puck,R.TEMPLATE),(cs,R.STRUCTURAL),(spl,R.STRUCTURAL),(lig,R.CATALYST)),(targetpuck,),geometry=AdapterLigationGeometry('unavailable','unavailable',(capture_end_inventory(),),'Source5primephosphorylatedcaptureoligosandsplintnamed;TableS1gated,fullcapture/bead3primejunctionunassembled.No sourceabsenceconclusion.'))
    b=g.input('hyb_buffer','Source6xSSC/RNaseinh2Uul/minimum15minRT',K.CONDITION);cap=g.product('captured','Sourcepuck-associatedRNA/tissuepreparation;RNAsequences/occupancyunassembled',K.MOLECULAR);g.react('capture_rna','SourceRNAhybridizationoncapturepuck',S.HYBRIDIZATION,((x,R.TEMPLATE),(targetpuck,R.STRUCTURAL),(b,R.CONDITION)),(cap,));x=cap
    clean('rt_dip','Sourcepuckdippedin1xMaximaRTbuffer')
    template=g.product('rna_component','SamepuckRNAtemplatecomponent;notphysicalfraction',K.MOLECULAR);primer=g.product('retained_bead_primers','Co-retainedbeadpolyT/customcaptureprimercomponent;notnewstockaddition',K.OLIGO);g.react('bead_components','Representretainedpucktemplate/primercomponents',S.GENERIC,((x,R.TEMPLATE),),(template,primer),changes_molecule=False);x=template
    b,ps,e=mix('rt','SourceMaximaHMinus10Uul/RTbuffer/dNTP1mM/NxGenRNaseinh1Uul/TSO2.5uM',('TemplateSwitchOligo',))
    y=g.product('rt_product','SourceRT/TSreactiononpuck;RNA/cDNAanduncopiedcompanions;wholejunctionunassembled',K.MOLECULAR);g.react('rt','Source30minRTthen52C90mincombinedRT/templateswitching',S.MAGIC,((x,R.TEMPLATE),(primer,R.STRUCTURAL),(ps[0],R.STRUCTURAL),(e,R.CATALYST),(b,R.CONDITION)),(y,),reason='SourcecombinedRT/TSwithbeadRTprimersandTSO;wholeRNA/cDNAandTSOjunctionunassembled,TSOnotmislabelledinitialRTprimer.');x=y
    pk=g.input('clearing_mix','Sourceequalvolume100mMTris/200mMNaCl/2percentSDS/5mMEDTA/16UmlPK;stocknotfinaldoubledconcentrations',K.CATALYST);y=g.product('cleared','Sourceproteolyzedtissue/puckRTmixture',K.MOLECULAR);g.react('clear_tissue','Source37C30minPKclearing',S.DEPROTEINIZATION,((x,R.TEMPLATE),(pk,R.CATALYST)),(y,));x=y
    b=g.input('bead_release_condition','Sourcevigorouspipetting/washes/centrifugation',K.CONDITION);beads=g.product('recovered_beads','Sourceindividualbeadsrecoveredfrompuck;boundcDNAsretained',K.MOLECULAR);res=g.product('residual_puck','Residualpuck/tissuematerial;composition/laterfateunassigned',K.BIOLOGICAL);wash=g.product('release_washes','Actualbead-recoverywash/supernatantfractions',K.WASTE);g.react('recover_beads','Sourcephysicalbeadreleaseandcollection',S.EXTRACTION,((x,R.TEMPLATE),(b,R.CONDITION)),(beads,res,wash));terminal(res,'Residualmaterialnotdeclaredempty',Disposition.UNRESOLVED);g.discard(wash,'Actualcollectionwashes');x=beads
    b,ps,e=mix('second','SourceKlenow125Uml/dNTP1mM/MaximaRTbuffer/dN-SMRT10uM',('dN-SMRT',));y=g.product('second_product','Sourcebead-associatedsecondstrand/cDNAmixture;wholeproductsunassembled',K.MOLECULAR);g.react('second','SourceKlenowsecond-strandsynthesis;time/temperaturedelegated23notassigned',S.EXTENSION,((x,R.TEMPLATE),(ps[0],R.STRUCTURAL),(e,R.CATALYST),(b,R.CONDITION)),(y,));x=y;clean('recover_second','Sourcecentrifugationrecoversbeads')
    b,ps,e=mix('amp','SourceTerraDirectPCRbuffer/Terra25Uml/forward2uM/reverse2uM;selectedSMART-PCR-Ronly',('TruSeqPCRhandleforward','SMART-PCR-R'))
    master=g.product('amp_master','Sourcebeadswithcompleteamplificationmixturebefore50uldivision',K.MOLECULAR);g.react('amp_mix','SourcecDNAamplificationmixtureaddition',S.GENERIC,((x,R.TEMPLATE),(e,R.CATALYST),(ps[0],R.STRUCTURAL),(ps[1],R.STRUCTURAL),(b,R.CONDITION)),(master,),changes_molecule=False)
    selected=g.product('amp_selected50','Sourceone50ulproductionreaction;totalmastervolumenotassigned',K.MOLECULAR);others=g.product('amp_parallel','Other50ulparallelproductionreactions;numberandpoolingnotexpanded',K.MOLECULAR);g.react('amp_allocate','Source50ulreactiondivision',S.SPLIT,((master,R.TEMPLATE),),(selected,others));terminal(others,'Otheractualparallelproductionreactionsnotreconstructed',Disposition.UNRESOLVED);x=selected
    y=g.product('amp_product','SourceamplifiedcDNA;wholeamplicon/retainedreagentsunassembled',K.MOLECULAR);g.react('amplify','SourceTerraamplification;cycleprogrammeunassignedinthisselectedscope',S.MAGIC,((x,R.TEMPLATE),),(y,),reason='Source50ulaliquotalreadycontainsactualenzymeandbothprimers;program/wholeproductunassembled,noextraenzymeaddition.');x=y;clean('amp_clean','Source1.2xSPRIselect')
    target=g.product('target10ng','Source10ngcDNAallocatedtoselectedTigreDARLINroute',K.MOLECULAR);rna=g.product('transcriptome_allocation','SourceallocationforwholetranscriptomeSlide-seqV2library;amountunassigned',K.MOLECULAR);rem=g.product('other_cdna','Otherlocus/libraryallocationsorstock;amount/fatesunexpanded',K.MOLECULAR);g.react('split_libraries','SourcepostcleanupcDNAsplit',S.SPLIT,((x,R.TEMPLATE),),(target,rna,rem));terminal(rna,'Transcriptomebranchdelegated23;actualfullrouteunassembled',Disposition.UNRESOLVED);terminal(rem,'Otheractualallocationsnotfabricated',Disposition.UNRESOLVED);x=target
    pcr('darlin1','Source10cyclefirstspatialDARLINPCR',('P5-CARLIN-F1','TA-RA-Spatial-PCR1-F'));clean('darlin1_clean','Source1.2xSPRIselect')
    pcr('darlin2','SourceTigresecondDARLINPCR;cyclecountunassigned',('P5-PR1','Tigre_Bulk_PCR2_R'));clean('darlin2_clean','Source1.2xSPRIselect')
    half=g.product('final_half','SourcehalfpurifiedsecondPCRproduct',K.MOLECULAR);other=g.product('other_half','Sourceotherhalfpurifiedproduct',K.MOLECULAR);g.react('half_for_index','Sourcehalfproductallocation',S.SPLIT,((x,R.TEMPLATE),),(half,other));terminal(other,'Untransferredhalf;lateruseunassigned');x=half
    pcr('index','Source10cyclefinalindexPCR',('P5','Spatial_Lib_Prep_i7_index;8ntsamplecode'));clean('size_select','Source.5x/1.2xdouble-sidedSPRIselect;retainedfraction/sideclassesnotsequence-assigned')
    qc=g.product('qc','SourceQubit/BioAnalyzeranalyticalallocations;volumesunassigned',K.QC);prod=g.product('production','Sourceproductionlibraryafteranalyticalallocations',K.MOLECULAR);g.react('qc_split','Actualsourceanalyticalallocations',S.SPLIT,((x,R.TEMPLATE),),(prod,qc));terminal(qc,'Actualanalyticalmaterialnotreturned',Disposition.QC);x=prod
    phix=g.input('phix','SourcePhiXControlv3DNAsequencing-libraryspike10percent;normalizationvolumeunassigned',K.MOLECULAR);lib=g.product('library','SourceDARLINlibraryplusPhiX;wholeDNA/readlandinggeometryunassembled',K.LIBRARY);g.react('sequence_mix','SourcephysicalPhiXspike',S.GENERIC,((x,R.TEMPLATE),(phix,R.SPIKE_IN)),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='AssociatedMS254selectedspatialDARLIN: R1=50/I1i7=8/R2=350;MiSeqv2500orNextSeq1k/2kP1600;10percentPhiX. Transcriptome50/8/50andbulk250/8/2505percentPhiXareseparatebranches. '+NOTICE))
    g.add_source_item(SourceItem('selected',AUTHORITY,NOTICE));return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph();sc=inventory()
    c('Wholeliteralprintedbead','TTT_PC_'+PREFIX+'J'*8+SPACER+'J'*6+'N'*8+'VT30',PRINTED_BEAD)
    c('PrintedUMIkeeps8notprose7',len(sc.strands['bead'].segs[6].top),8)
    c('Printedsecondspatialblockkeeps6notanalysis7',len(sc.strands['bead'].segs[5].top),6)
    c('Samefeatureidentityacrosssplitspatialcode',sc.strands['bead'].segs[3].feature,sc.strands['bead'].segs[5].feature)
    c('CustomcapturebeforeRNAhybridization',g.steps['capture_rna'].inputs[1].material,'targeted_puck')
    c('RTretainedprimercomponentnotnewbeadstock',[(i.material,i.role)for i in g.steps['rt'].inputs if i.role==R.STRUCTURAL][0],('retained_bead_primers',R.STRUCTURAL))
    c('Proteinaseclearingtypedproteolysis',g.steps['clear_tissue'].kind,S.DEPROTEINIZATION)
    c('ParallelPCRfatesexplicit',g.terminals['amp_parallel'].disposition,Disposition.UNRESOLVED)
    c('Transcriptomeallocationpreserved',g.terminals['transcriptome_allocation'].disposition,Disposition.UNRESOLVED)
    c('Finalhalfstockpreserved',g.terminals['other_half'].disposition,Disposition.STORED)
    c('PhiXphysicalDNArole',[(i.material,i.role)for i in g.steps['sequence_mix'].inputs if i.role==R.SPIKE_IN],[('phix',R.SPIKE_IN)])
    c('Knowncapture5Pdrawn',capture_end_inventory().strands['capture_oligo'].mod5,'P')
    from checks import Source,have
    from storage import sources_root
    source=Source(sources_root()/'space-seq__10.1016+j.stem.2026.05.017'/'j.stem.2026.05.017_PMC13401910.xml.txt','https://pmc-oa-opendata.s3.amazonaws.com/PMC13401910.1/PMC13401910.1.xml via tools/get_sources.py; Methods246')
    if have(c,source):c('Wholeprintedbeadagainstselectedoriginal',PRINTED_BEAD in source.path.read_text(),True)
    c('Privatepartialgraph',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
