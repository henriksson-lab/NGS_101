"""Rho-seq associated-manuscript chemical R+/R− pair; private partial graph."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature,revcomp
from illumina import TRUSEQ_READ1,TRUSEQ_READ2
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
TITLE='Rho-seq: selected S.pombe R+/R− endogenous-RNA preparation'
SOURCE='AssociatedauthorMS PMC8792297.1/NIHMS1756155 main260/268/288/keyadapters884–886. FinalpublishedbodyandgatedSIunverified; STAR2022companionordinarydeferred,notanotherprotocol. Oneoriginal;notfullsource/absenceaudit.'
RNA_ADAPTER_BODY=revcomp(TRUSEQ_READ1)[:19]
DNA_ADAPTER_BODY=revcomp(TRUSEQ_READ2)[:22]

SOURCE_ORDERINGS={arm+'_'+kind+'_adapter':{'source':'AssociatedMS keyresources '+('884'if kind=='rna'else'886'),'written':'5′-P-'+('AGATCGGAAGAGCGTCGTG'if kind=='rna'else'AGATCGGAAGAGCACACGTCTG')+'-ddC'}for arm in ('plus','minus')for kind in ('rna','dna')}
SOURCE_OLIGO_EXCEPTIONS={key:'AssociatedMS keyresource wholeordering ends in literal -ddC. Model retains the additional terminal modified C nucleotide; current source-text matching omits that terminal base. Exact whole body, phosphate and ddC are locally checked against selected original, not repaired or declared absent.'for key in SOURCE_ORDERINGS}

def adapter(name):
    if name not in {'rna','dna'}:raise ValueError('ExplicitRNA-targetorDNA-targetadapter')
    body=RNA_ADAPTER_BODY if name=='rna'else DNA_ADAPTER_BODY
    sc=Scene();sc.strand('adapter',[Segment('source_body',body),Segment('terminal_ddC','C',note='Source-ddCisadditionalterminalmodifiedcytidine')],mod5='P',mod3='ddC; no3OH');return sc

def written_ordering(name):
    st=adapter(name).strands['adapter']
    terminal='ddC'if st.segs[-1].top=='C'and st.mod3.startswith('ddC')else'INVALID_TERMINAL'
    return '5′-'+st.mod5+'-'+''.join(x.top for x in st.segs[:-1])+'-'+terminal

def ends(name):
    sc=adapter(name);sc.strand('substrate',[Segment('source_substrate','[RNA fragment unassembled]'if name=='rna'else'[cDNA unassembled]',placeholder=True)],mod3='OH')
    sc.ligation_end('adapter','source_body',through='terminal_ddC',kind='single-strand',polarity='5-prime',owner='adapter',text='Source5Pdonor;3ddCblocksoppositeend')
    sc.ligation_end('substrate','source_substrate',kind='single-strand',polarity='3-prime',owner='substrate',text='SourceRNAafterFastAP orcDNA3OH;wholeproductunassembled');return sc

def indexed_primer(arm):
    sc=Scene();sc.strand('barcoded_primer',[Segment('source_prefix','[PCR context unassembled]',placeholder=True),Segment('source_index','[sample index unassigned]',placeholder=True,feature=feature('rho_'+arm,'sample_index','unknown')),Segment('source_binding','[binding arm unassembled]',placeholder=True)]);return sc

def unknown():
    sc=Scene();sc.strand('unknown',[Segment('unassembled','[unassembled]',placeholder=True)]);return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance');stock=g.input('rna','SourcetotalRNAstockfromS.pombestudysample;WT/Δ4genotypeunassignedintheselectedchemistrypair;3ngsyntheticspikeadditiontimingnotreconstructed',K.MOLECULAR,primary=True)
    pair={a:g.product(a+'_input','Source30ugtotalRNAfor'+a+'chemicalcondition',K.MOLECULAR)for a in ('plus','minus')};rem=g.product('rna_remainder','UntransferredstartingRNAstockifany;amountunassigned',K.MOLECULAR);g.react('allocate_pair','SourceR+/R−chemistryallocations',S.SPLIT,((stock,R.TEMPLATE),),(*pair.values(),rem));g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,'Actualstockallocationbeyondselectedconditionpairunassembled'))
    for arm,x in pair.items():
        def key(k):return arm+'_'+k
        def op(k,label,kind,stocklabel,product,cat=False,extra=(),reason=None,structural=False):
            nonlocal x
            c=g.input(key(k+'_stock'),stocklabel,K.CATALYST if cat else K.MOLECULAR if structural else K.CONDITION);y=g.product(key(k+'_product'),product,K.MOLECULAR);kw={}
            if kind==S.MAGIC:kw['reason']=reason or'Sourceoperationpositive;fullstock/productgeometryunassembled.'
            if kind==S.GENERIC:kw['changes_molecule']=False
            if structural:
                conditions=g.input(key(k+'_conditions'),'Source37C90min dark; acidic formate mixture retained from preceding step',K.CONDITION)
                extra=(*extra,(conditions,R.CONDITION))
            g.react(key(k),label,kind,((x,R.TEMPLATE),(c,R.CATALYST if cat else R.STRUCTURAL if structural else R.CONDITION),*extra),(y,),**kw);x=y
        def clean(k,label):
            nonlocal x
            c=g.input(key(k+'_buffer'),label,K.CONDITION);y=g.product(key(k+'_kept'),'Sourcephysicallyretainedpreparation;yield/carryoverunassigned',K.MOLECULAR);w=g.product(key(k+'_waste'),'Actualsourcephysicalwash/organic/cleanupfraction',K.WASTE);g.react(key(k),label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(y,w));g.discard(w,'Actualphysicalfractionnotperfectmolecularselection');x=y
        def aliquot(k,label):
            nonlocal x
            take=g.product(key(k+'_selected'),label,K.MOLECULAR);left=g.product(key(k+'_remaining'),'Untransferredstockifany;amount/fateunassigned',K.MOLECULAR);g.react(key(k),'Sourceproductionstockallocation',S.SPLIT,((x,R.TEMPLATE),),(take,left));g.mark_terminal(Terminal(left.id,Disposition.UNRESOLVED,'Sourceallocationdoesnotestablishwholestockconsumptionorlaterfate'));x=take
        op('reduce','Source25C60mindark750rpm',S.CHEMICAL_CONVERSION,'Source400ul40mMTris7.4/'+('4mgNaBH4freshin1mMKOH'if arm=='plus'else'KOHalonewithoutNaBH4'),'Source'+arm+'chemicallyprocessedRNA;modificationselectivitynotperfect')
        op('quench','SourceglacialAcOHfinalapproximately.3M',S.QUENCHING,'Sourceaceticacid','Sourceacid-quenchedRNA/reagentmixture');clean('precip1','Sourceethanolprecipitationasmain250')
        op('resuspend','Source5ulwater+85ul100mMformatepH4/65C7min',S.THERMAL_TREATMENT,'Sourcewater/formatebuffer','SourceresuspendedRNA')
        op('rhodamine','SourceRho1102.2mM37C90mindarkBOTHconditions',S.CHEMICAL_CONVERSION,'SourceRhodamine110chloride;10xstockinmethanol','SourceRNAwithdifferentialRhoadductsandunmodifiedclasses;notD-specificwithoutgenotypecontrols',structural=True)
        op('neutralize','SourcepHraised7.5with.8MTris8.5',S.QUENCHING,'SourceTris','SourceneutralizedRhoRNAmixture');clean('phenol','Sourcephenolextraction65C60min/collectupperphase');clean('precip2','Sourceethanolprecipitationasmain250')
        aliquot('for_ribodepletion','Source2–5ugchemicallyprocessedRNAforlibraryribodepletion')
        op('ribodeplete','SourceS.pombemodifiedThermoA15020/K155005manufacturerprotocol',S.MAGIC,'Sourcehistoricalribodepletionkitconditions;exactguide/versionunread','SourceRNAafterribodepletion;exactphysicalsubsteps/wholeprobepoolunassembled')
        qc=g.product(key('rna_qc'),'SourceQubitRNAHS/BioAnalyzeranalyticalallocations;volumesunassigned',K.QC);prod=g.product(key('postrnaqc'),'SourceproductionribodepletedRNAstock',K.MOLECULAR);g.react(key('rna_qc_split'),'ActualRNAanalyticalallocation',S.SPLIT,((x,R.TEMPLATE),),(prod,qc));g.mark_terminal(Terminal(qc.id,Disposition.QC,'Analyticalmaterialnotreturned'));x=prod;aliquot('for_fragment','Source90–125ngribodepletedRNAforzincfragmentation')
        op('fragment','Sourcezincfragmentation70C15min/ThermoAM8740stoppermanufacturer',S.FRAGMENTATION,'Sourcezincfragmentationreagents/stopboundaryunassembled','SourcefragmentedRNA;end/lengthdistributionunassembled');clean('silane1','SourceMyOneSILANE/RLTbetweenlibraryoperations')
        dn=g.input(key('dnase'),'SourceTURBODNase2U',K.CATALYST);ap=g.input(key('fastap'),'SourceFastAP3U',K.CATALYST);b=g.input(key('end_buffer'),'Source37C30min/TrisMgCaK/TritonDTT/RNaseinh40U',K.CONDITION);y=g.product(key('end_processed'),'SourceDNA-degraded/dephosphorylatedRNAclasses;completeendoccupancy/unremovedcleavageproductsunassembled',K.MOLECULAR);g.react(key('end_process'),'SourcesimultaneousTURBODNase/FastAP',S.MAGIC,((x,R.TEMPLATE),(dn,R.CATALYST),(ap,R.CATALYST),(b,R.CONDITION)),(y,),reason='SourcecompoundDNAdigestionandRNAenddephosphorylation;reactionorderwithinmixtureandwholeendproductsunassembled.');x=y;clean('silane2','SourceMyOneSILANE/RLTbetweenlibraryoperations')
        a=g.input(key('rna_adapter'),'Source20pmolRNA-targetadapter;targetlabeldoesnotassignallsugarchemistry',K.ADAPTER,state=adapter('rna'),chemical_state={'five_prime':'P','three_prime':'ddC','polymer':'RNA-targetadapter;exactbodiesaswrittenwithT,sugarpartitionunassigned'})
        b=g.input(key('rna_adapter_heat'),'Source70C2min',K.CONDITION);templ=g.product(key('rna_template'),'RNAcomponentofsameheatedmixture',K.MOLECULAR);aret=g.product(key('retained_rna_adapter'),'Co-retained20pmoladaptercomponent;notnewaddition',K.ADAPTER,state=adapter('rna'));g.react(key('rna_adapter_prepare'),'Sourceadapteradditionanddenaturation',S.THERMAL_TREATMENT,((x,R.TEMPLATE),(a,R.STRUCTURAL),(b,R.CONDITION)),(templ,aret))
        lig=g.input(key('rna_ligase'),'SourceT4RNAligase39U/RNaseinh12U/Ligase1buffer/9percentDMSO/1mMATP/20percentPEG8000',K.CATALYST);y=g.product(key('rna_adapted'),'Source3adapter-ligated/unligatedRNAclasses;wholejunctionunassembled',K.MOLECULAR);g.react(key('rna_ligate'),'Source25C90min3endligation',S.LIGATION,((templ,R.TEMPLATE),(aret,R.STRUCTURAL),(lig,R.CATALYST)),(y,),geometry=AdapterLigationGeometry.known('single-strand',ends('rna'),statement='AssociatedMS288FastAPdephosphorylatesRNAends;key884adapter5P/body/terminalddCknown. FullRNA/productandsugarpartitionunassembled.'));x=y;clean('silane3','SourceMyOneSILANE/RLTbetweenlibraryoperations')
        p=g.input(key('rt_primer'),'Source10pmolDNAprimer;wholeorderingunassembled,notanalytical32Pprimer',K.OLIGO);b=g.input(key('rt_heat'),'Source70C2min',K.CONDITION);templ=g.product(key('rt_template'),'SourceRNAtemplatecomponentafterprimeraddition/heat',K.MOLECULAR);pret=g.product(key('retained_rt_primer'),'Same-tubeDNAprimercomponent;notnewfreeprimeraddition',K.OLIGO);g.react(key('rt_prepare'),'SourceRTprimeradditionanddenaturation',S.THERMAL_TREATMENT,((x,R.TEMPLATE),(p,R.STRUCTURAL),(b,R.CONDITION)),(templ,pret))
        rt=g.input(key('ssiii'),'SourceSSIIIasmain272RTconditions;not32Pprimer-labelingreaction',K.CATALYST);nts=g.input(key('rt_dntps'),'Sourceas-aboveRT150uMdNTP;per-base-versus-totalconventionunassigned',K.MOLECULAR);y=g.product(key('cdna'),'SourceRNA/cDNA/freeprimerRTmixture;exactRTstop/endstatesunassembled',K.MOLECULAR);g.react(key('rt'),'SourceSSIIIreverse-transcriptionstage',S.REVERSE_TRANSCRIPTION,((templ,R.TEMPLATE),(pret,R.STRUCTURAL),(nts,R.STRUCTURAL),(rt,R.CATALYST)),(y,));x=y
        op('exosap','SourceExoSAPIT2ul37C4minthen80C1min',S.DIGESTION,'SourceExoSAPITstock','Sourceprimer-digestedcDNA/RNAmixture;completeprimerclearancenotasserted',cat=True)
        op('hydrolyze_rna','SourceNaOH100mM70C10min',S.CHEMICAL_CONVERSION,'SourceNaOH','SourceRNAhydrolyzed/cDNAmixture');op('neutralize2','Sourceaceticacidneutralization',S.QUENCHING,'Sourceaceticacid','SourceneutralizedcDNAmixture');clean('silane6','SourceMyOneSILANE/RLTafterRT/primer-removal/RNAhydrolysisgroup;exactintermediatesubstepgroupingunassembled')
        a=g.input(key('dna_adapter'),'Source40pmolDNAadapter',K.ADAPTER,state=adapter('dna'),chemical_state={'polymer':'DNA','five_prime':'P','three_prime':'ddC'});b=g.input(key('dna_heat'),'Source75C2min',K.CONDITION);templ=g.product(key('dna_template'),'cDNAcomponentofheatedmixture',K.MOLECULAR);aret=g.product(key('retained_dna_adapter'),'Co-retained40pmolDNAadaptercomponent;notnewaddition',K.ADAPTER,state=adapter('dna'));g.react(key('dna_prepare'),'SourceDNAadapteradditionanddenaturation',S.THERMAL_TREATMENT,((x,R.TEMPLATE),(a,R.STRUCTURAL),(b,R.CONDITION)),(templ,aret))
        lig=g.input(key('dna_ligase'),'SourceT4RNAligase48U/Ligase1buffer/4percentDMSO/1mMATP/23.75percentPEG8000',K.CATALYST);y=g.product(key('dna_adapted'),'Source3adapter-ligatedandunligatedcDNAclasses;wholejunctionunassembled',K.MOLECULAR);g.react(key('dna_ligate'),'Source25CovernightDNA3endligation',S.LIGATION,((templ,R.TEMPLATE),(aret,R.STRUCTURAL),(lig,R.CATALYST)),(y,),geometry=AdapterLigationGeometry.known('single-strand',ends('dna'),statement='AssociatedMS288DNA3endligationandkey8865P/terminalddCadapter;fullcDNAsealedproductunassembled.'));x=y;clean('silane7','SourceMyOneSILANE/RLTbetweenlibraryoperations')
        mix=g.input(key('pcr_mix'),'Source50ulPCRmixwith25ulNEBNextHiFi2xM0541/1uMbarcodedanduniversalprimers',K.MOLECULAR);e=g.product(key('pcr_e'),'CatalystcomponentofsamePCRbulkaddition',K.CATALYST);f=g.product(key('pcr_f'),'BarcodedPCRprimercomponent;exactordering/indexunassembled',K.OLIGO,state=indexed_primer(arm));r=g.product(key('pcr_r'),'UniversalPCRprimercomponent;wholeorderingunassembled',K.OLIGO);b=g.product(key('pcr_b'),'Buffercomponentofsamebulkaddition',K.CONDITION);nt=g.product(key('pcr_nt'),'Molecular dNTP substrate component of same PCR bulk addition; concentration convention unassigned',K.MOLECULAR);g.react(key('pcr_components'),'RepresentoneactualPCRbulkaddition',S.GENERIC,((mix,R.TEMPLATE),),(e,f,r,b,nt),changes_molecule=False)
        y=g.product(key('pcr_product'),'Sourcebarcodedamplicons/primerandothercompanions;wholeproductunassembled',K.MOLECULAR);g.react(key('pcr'),'Source98C30s/5cycles98C10s67C30s72C30s/13cycles98C10s72C30s/final72C60s',S.PCR,((x,R.TEMPLATE),(e,R.CATALYST),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(nt,R.STRUCTURAL),(b,R.CONDITION)),(y,),geometry=PCRGeometry.unavailable(unknown(),reason='HistoricalPCRpairnamesknown;fullprimerorderings/landingandindexedampliconunassembled.LinkedSIgated,itsoligodisclosurescopeunverified.'));x=y;clean('size1','Source1.2xAMPure');clean('size2','Source1.1xAMPure')
        q=g.product(key('dna_qc'),'SourceQubit/DNAmicrochipanalyticalallocations;volumesunassigned',K.QC);prod=g.product(key('production'),'SourceproductionlibraryafterDNAQC',K.MOLECULAR);g.react(key('dna_qc_split'),'ActualDNAanalyticalallocation',S.SPLIT,((x,R.TEMPLATE),),(prod,q));g.mark_terminal(Terminal(q.id,Disposition.QC,'Analyticalmaterialnotreturned'));lib=g.product(key('library'),'Source'+arm+'Rholibrary;wholeends/readlandingunassembled',K.LIBRARY);g.react(key('identity'),'Physicallibraryidentity',S.GENERIC,((prod,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='AssociatedMS288HiSeq2500orNovaSeqS1/paired2x50nt;indexmasksunassigned. EndogenousRNApreparationcomponentonly:3ngsyntheticspikesequence/additiontimingunassembled. WT/Δ4controlcomparisonsrequiredforD-specificcalling;RNA/Riboprofilingcompanionlibrariesnotmerged.'))
    g.add_source_item(SourceItem('selected',SOURCE,'RNA/DNAadapterendsknown;wholeproducts,studyspikeentrypoint,exactSILANEsubstepgroupingsandfinalversionunassembled.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('WholeRNAadapterordering',RNA_ADAPTER_BODY,'AGATCGGAAGAGCGTCGTG')
    c('WholeDNAadapterordering',DNA_ADAPTER_BODY,'AGATCGGAAGAGCACACGTCTG')
    c('WholeRNAmodelmodordering',written_ordering('rna'),SOURCE_ORDERINGS['plus_rna_adapter']['written'])
    c('WholeDNAmodelmodordering',written_ordering('dna'),SOURCE_ORDERINGS['plus_dna_adapter']['written'])
    c('ddCisadditionalnucleotide',len(''.join(s.top for s in adapter('dna').strands['adapter'].segs)),len(DNA_ADAPTER_BODY)+1)
    c('BotharmsreceiveRhodamine',all(a+'_rhodamine'in g.steps for a in ('plus','minus')),True)
    c('MockhasnoNaBH4',g.materials['minus_reduce_stock'].label.count('withoutNaBH4'),1)
    c('Labelstockallocationexplicit',g.terminals['plus_for_ribodepletion_remaining'].disposition,Disposition.UNRESOLVED)
    c('Fragmentstockallocationexplicit',g.terminals['plus_for_fragment_remaining'].disposition,Disposition.UNRESOLVED)
    c('RTconsumesretainedprimer',[(i.material,i.role)for i in g.steps['plus_rt'].inputs if i.material=='plus_retained_rt_primer'],[('plus_retained_rt_primer',R.STRUCTURAL)])
    c('RNAadapterconsumesretainednotnewstock',[(i.material,i.role)for i in g.steps['plus_rna_ligate'].inputs if i.role==R.STRUCTURAL],[('plus_retained_rna_adapter',R.STRUCTURAL)])
    c('DNAadapterconsumesretainednotnewstock',[(i.material,i.role)for i in g.steps['plus_dna_ligate'].inputs if i.role==R.STRUCTURAL],[('plus_retained_dna_adapter',R.STRUCTURAL)])
    c('RTnucleotidesarestructuralsubstrate',[(i.material,i.role)for i in g.steps['plus_rt'].inputs if i.material=='plus_rt_dntps'],[('plus_rt_dntps',R.STRUCTURAL)])
    c('Knownlocalendgeometry',g.steps['plus_dna_ligate'].geometry.availability,'known')
    src=Source(sources_root()/'rho-seq__10.1016+j.molcel.2021.11.003'/'j.molcel.2021.11.003_PMC8792297.xml.txt','https://pmc-oa-opendata.s3.amazonaws.com/PMC8792297.1/PMC8792297.1.xml;keyadapters884–886')
    if have(c,src):
        t=src.path.read_text();c('WholeRNAadaptermodorderingagainstoriginal',written_ordering('rna')in t,True);c('WholeDNAadaptermodorderingagainstoriginal',written_ordering('dna')in t,True)
    c('Privatepartialgraph',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
