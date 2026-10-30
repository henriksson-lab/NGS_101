"""Private TaDRIM-seq source-owned loaded stock, bridge and selected cellular workflow."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp
from nextera import ME_RC
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
TITLE='TaDRIM-seq: selected mammalian targeted RNA/DNA proximity preparation'
AUTHORITY='PublishedPMC11811046main69–71/169–193 and SI1SupplementaryNote1 originalPDF13–16 selectedsourceorderings/stocktables/downstream200nginput. Twooriginals; no plant/comparator/modernkit inheritance.'
A=ME_RC+'T';B=revcomp(A)
F='CCTCGATAC'+'T'+'TGAGCTGAC';RBR='GTCAGC'+'T'+'CAAGTATCGAG'
SOURCE_ORDERINGS={'A':'5-phos-'+A,'B':'5-phos-'+B,'bridge_f':'[5APP]CCTCGATAC/iBIOdT/TGAGCTGAC','bridge_r':'[5Phos]GTCAGC/iBIOdT/CAAGTATCGAG'}
def oligo(seq,mod5='P'):
    s=Scene();s.strand('oligo',[Segment('source_whole',seq)],mod5=mod5);return s
def tte_adapter():
    s=Scene();s.strand('A',[Segment('A',A)],mod5='P');s.anneal('B',[Segment('B',B)],to='A',pair=('B','A'),mod5='P');return s
def bridge():
    s=Scene();s.strand('forward',[Segment('forward_overhang',F[:1]),Segment('forward_stem',F[1:],note='Internalbiotin-dT at whole-strand0-basedposition9; retained as actualT')],mod5='App')
    s.anneal('reverse',[Segment('reverse_stem',RBR,note='Internalbiotin-dT at whole-strand0-basedposition6; retained as actualT')],to='forward',pair=('reverse_stem','forward_stem'),mod5='P');return s
def rna_geometry(*,inferred=True):
    s=bridge();s.strand('rna',[Segment('rna','[RNAfragment]',placeholder=True,inferred=inferred,note='Intended3OHacceptorclassafterPNK; actualenddistributionunassembled')],mod3='OH')
    s.ligation_end('forward','forward_overhang',kind='single-strand',polarity='5-prime',owner='adapter',text='Sourcebridge5App/unpairedfirstC, notborrowedChRDPETterminalT')
    s.ligation_end('rna','rna',kind='single-strand',polarity='3-prime',owner='substrate',text='ConditionalRNA3OHclass;wholeligatedRNA/bridgeproductunassembled');return s

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance');x=g.input('cells','Sourceabout1millionK562cellsatselectedCTCFtargetingboundary; otherantibodies/celltypesnotmerged',K.BIOLOGICAL,primary=True)
    def inp(k,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(g.input(k,label,kind,**kw),role)
    def cond(k,label):return inp(k,label,K.CONDITION,R.CONDITION)
    def nt(k,label):return inp(k,label,K.MOLECULAR,R.STRUCTURAL)
    def op(k,label,kind,ins=(),geometry=None,product='SourcecopresentDNA/RNA/reagentmixture with residual/unprocessedclasses untilphysicalworkup'):
        nonlocal x
        y=g.product(k+'_product',product,K.MOLECULAR);kw={}
        if geometry is not None:kw['geometry']=geometry
        if kind==S.GENERIC:kw['changes_molecule']=False
        if kind==S.MAGIC:kw['reason']='Sourceoperationpositive; completehistoricalstockcomposition/priminggeometry/wholeproductsunassembled,notabsence.'
        g.react(k,label,kind,((x,R.TEMPLATE),)+tuple(ins),(y,),**kw);x=y
    def clean(k,label):
        nonlocal x
        y=g.product(k+'_kept','Sourcephysicallyretainedpreparation; recovery/carryoverunassigned',K.MOLECULAR);w=g.product(k+'_waste','Sourcephysicalworkupfractions',K.WASTE);g.react(k,label,S.CLEANUP,((x,R.TEMPLATE),cond(k+'_condition',label)),(y,w));g.discard(w,'Sourcephysicalwashing,notperfectmolecularspecificity.');x=y
    # Independent, actual source TTE loading; not fresh enzyme/adapter re-additions at tagmentation.
    a=g.input('adapter_a','Source100uM5P-'+A,K.OLIGO,state=oligo(A));b=g.input('adapter_b','Source100uM5P-'+B,K.OLIGO,state=oligo(B));ab=g.product('annealed_ab','SourceequimolarA/Bprimer mixture',K.ADAPTER,state=tte_adapter());g.react('anneal_ab','Source75C15min60C10min50C10min40C10min25C30min',S.HYBRIDIZATION,((a,R.TEMPLATE),(b,R.STRUCTURAL),cond('anneal_buffer','Sourceannealingbuffer/thermalprogramme')),(ab,))
    selected_ab=g.product('ab3_5','Source3.5ulannealedadapterselected',K.ADAPTER,state=tte_adapter());abrem=g.product('ab_remainder','Untransferredannealedadapterstock;volume/fateunassigned',K.MOLECULAR);g.react('ab_split','Source3.5ulstockaliquot',S.SPLIT,((ab,R.TEMPLATE),),(selected_ab,abrem));g.mark_terminal(Terminal(abrem.id,Disposition.UNRESOLVED,'Actualsourceadapterstockremainder.'))
    pg=g.input('unloaded_pg','Source20ulpG-Tn5at500ng/ul',K.MOLECULAR);loaded=g.product('tte_stock','Source37.5ulTTEstock,exactloadingstoichiometry/freecomponentsunassembled',K.MOLECULAR);g.react('load_tte','Source20ulpGTn5+3.5ulprimer mix+14ulcouplingbuffer30C1h',S.BINDING,((pg,R.TEMPLATE),(selected_ab,R.STRUCTURAL),cond('coupling_buffer','Sourcecouplingbuffer')),(loaded,))
    tte=g.product('tte5','Source5ulTTEaliquot',K.MOLECULAR);tr=g.product('tte_remainder','Sourceuntransferrednominal32.5ulTTEstock,-20Cstored;laterfateunassigned',K.MOLECULAR);g.react('tte_split','Source5ulselectedfrom37.5ulTTEstock',S.SPLIT,((loaded,R.TEMPLATE),),(tte,tr));g.mark_terminal(Terminal(tr.id,Disposition.UNRESOLVED,'Storedstocknotre-addedlater.'))
    # Bridge forward phosphate is adenylated before reverse-strand annealing.
    f=g.input('bridge_f_precursor','Source6ul100uMforward5Phos CCTCGATAC/iBIOdT/TGAGCTGAC',K.OLIGO,state=oligo(F),chemical_state={'internal_biotin_dT':9,'written':'[5Phos]CCTCGATAC/iBIOdT/TGAGCTGAC'})
    fa=g.product('adenylated_f','Sourceforward5App,internalbiotinTretained;conversionyieldunassigned',K.OLIGO,state=oligo(F,'App'),chemical_state={'internal_biotin_dT':9});g.react('adenylate','Source120ulreaction:12ulbuffer/12ul1mMATP/6ulforward/12ulMth/78ulwater65C1h85C5min',S.CHEMICAL_CONVERSION,((f,R.TEMPLATE),inp('mth','SourceMthRNAligase'),nt('adenylation_atp','SourceATPadenylatecontribution'),cond('adenylation_buffer','Source5DNAadenylationbuffer')),(fa,))
    rb=g.input('bridge_r','Source6ul100uMreverse5PhosGTCAGC/iBIOdT/CAAGTATCGAG',K.OLIGO,state=oligo(RBR),chemical_state={'internal_biotin_dT':6,'written':SOURCE_ORDERINGS['bridge_r']});bs=g.product('annealed_bridge','Sourceannealedbiotinbridge;fullproductsunassembled',K.ADAPTER,state=bridge(),chemical_state={'forward_biotin_dT':9,'reverse_biotin_dT':6,'forward_5prime':'App','reverse_5prime':'P'});g.react('anneal_bridge','Sourceaddreverse/staged95Cto20Canneal',S.HYBRIDIZATION,((fa,R.TEMPLATE),(rb,R.STRUCTURAL),cond('bridge_anneal','SourceSI14thermalprogramme')),(bs,))
    ready=g.product('bridge_ready','Sourcecolumnpurifiedbridgeeluted50ul;quantitativerecoveryunassigned',K.ADAPTER,state=bridge(),chemical_state={'forward_biotin_dT':9,'reverse_biotin_dT':6,'forward_5prime':'App','reverse_5prime':'P'});bw=g.product('bridge_column_waste','Sourcebridgecolumnworkupfractions',K.WASTE);g.react('bridge_clean','SourceD4060oligocolumn/elute50ul',S.CLEANUP,((bs,R.TEMPLATE),cond('bridge_column','Sourceoligocleanconcentrator')),(ready,bw));g.discard(bw,'Sourcecolumnworkupfractions.')
    op('fix','Source1percentformaldehyde10minRT',S.CROSS_LINKING,(cond('formaldehyde','Sourceformaldehyde'),));op('fix_quench','Source.2Mglycine5min',S.QUENCHING,(cond('glycine','Sourceglycine'),));clean('fix_wash','Source3xPBS')
    op('nuclei','Source1mlTris/NaCl/.2percentNP40/PI/RNaseinhibitorice10min',S.PERMEABILIZATION,(cond('lysis_buffer','Sourcesmallmammaliancelllysisbuffer'),));clean('nuclei_wash','Source2500g3min4C/2x300ulwashA')
    op('antibody','Source10ugCTCFantibody4C2h',S.BINDING,(inp('ctcf','SourceCTCFAbab188408',K.SELECTION,R.SELECTION),cond('antibody_buffer','Sourceantibodybuffer/digitonin/EDTA/BSA')));clean('antibody_wash','Source5xwashA')
    op('target_bind','Sourceadd5ulpreparedTTEin500ulincubationbuffer4C1h',S.BINDING,((tte,R.STRUCTURAL),));clean('target_wash','Source5x300ultransposaseincubationbuffer')
    t=g.product('tag_template','Same-tubenuclearDNA/RNAcomponent',K.MOLECULAR);cat=g.product('retained_pg','RetainedloadedpGTn5catalystcomponent,notnewstock',K.CATALYST);ad=g.product('retained_ab','RetainedloadedA/Badaptercomponent,notnewstock',K.ADAPTER,state=tte_adapter());g.react('tag_components','Same-tubecomponentdeclarations,notphysicalfractionation',S.GENERIC,((x,R.TEMPLATE),),(t,cat,ad),changes_molecule=False);x=t
    op('tag','Source500ulbuffer+5ul1MMgCl237C1h',S.TAGMENTATION,((cat,R.CATALYST),(ad,R.STRUCTURAL),cond('tag_buffer','SourceMg/transposasebuffer')))
    op('tag_stop','Source500ul40mMEDTA',S.QUENCHING,(cond('edta','SourceEDTA'),));clean('tag_wash','Source2xPBS')
    op('sds','Source200ul.5percentSDS62C10min',S.PERMEABILIZATION,(cond('sds_buffer','SourceSDS'),));op('sds_quench','Source580ulwater+100ul10percentTriton37C15min',S.QUENCHING,(cond('triton','SourceTritonX100'),));clean('rnase_prepare','Source2xwashB')
    op('rnase','Source1.6ultenfolddilutedRNaseIf/500ulCutSmart37C3minice5min',S.DIGESTION,(inp('rnase_if','SourceRNaseIf'),));clean('rnase_pellet','Sourcecentrifuge/retainnuclei')
    op('pnk','Source400ulPNKsolution/ATP/RNaseinhibitor37C40min',S.PHOSPHORYLATION,(inp('pnk_enzyme','SourceT4PNK'),nt('pnk_atp','Source1mMATPphosphatedonor')));clean('pnk_wash','Source2x500ulwashB')
    op('rna_ligate','Source400ulRNA-ligationmixture/50ulpreparedbridge/Rnl2truncated16Covernight',S.LIGATION,((ready,R.STRUCTURAL),inp('rnl2','SourceT4RNAligase2truncated'),cond('rna_lig_buffer','Source15percentPEG/RNaseinhibitor/Triton/ligationbuffer')),AdapterLigationGeometry.known('single-strand',rna_geometry(),statement='Source5AppfirstC/conditionalRNA3OHinterface;biotinatedduplexretained,wholeligatedsubstrateunassembled.'))
    clean('rna_lig_wash','Source5x500ulwashB')
    t=g.product('rt_template','Same-tubeligatedRNA/forwardbridge/nuclearcontext',K.MOLECULAR);p=g.product('retained_bridge_primer','Co-retainedreversebridgeDNAprimercomponent; occupancy/primingdistributionunassigned',K.OLIGO,state=oligo(RBR),chemical_state={'internal_biotin_dT':6,'written':SOURCE_ORDERINGS['bridge_r']});g.react('rt_components','Existingbridgeprimercomponent,notnewfreeRTstock',S.GENERIC,((x,R.TEMPLATE),),(t,p),changes_molecule=False);x=t
    op('rt','Source400ulMMLVRNaseHminus/.5mMdNTP42C1h',S.REVERSE_TRANSCRIPTION,((p,R.STRUCTURAL),inp('mmlv','SourceMMLVRNaseHminus'),nt('rt_nt','Source.5mMdNTPmixture,perbase/totalconventionunassigned')));clean('rt_wash','Source3x500ulwashB')
    op('proximal_ligate','SourceT4DNAligase/10mMATPbuffer16Covernight;SIalsooffers4hRTnotmixedhere',S.PROXIMITY_LIGATION,(inp('dna_ligase','SourceT4DNAligase'),cond('dna_lig_buffer','SourceATP/Triton/BSA/ligationbuffer')),AdapterLigationGeometry('unavailable','unavailable',(bridge(),),'SourceRNA-DNA/RNA-RNAproximaljoiningpositive;genomictranspositionends/wholechimerajunctionsunassembled,notabsence.'))
    op('reverse_crosslink','SourceproteinaseK/SDS65C3h800rpm',S.DEPROTEINIZATION,(inp('pk','SourceproteinaseK'),cond('pk_buffer','SourceTris/SDS/NaCl/EDTA')))
    clean('organic','SourceequivalentphenolchloroformIAA/MaxTractcollectupperaqueousfraction')
    op('carrier_add','SourceSIstep22add2ulGlycoBluecarrierbeforeprecipitation',S.GENERIC,(inp('glycoblue','Source2ulGlycoBluemolecularcarrier;concentrationandquantitativecarryoverunassigned',K.MOLECULAR,R.STRUCTURAL),),product='SourceaqueousDNA/RNApreparationplusGlycoBluecarrier;recovery/carryoverunassigned')
    clean('precipitate','Sourceglycoblue/acetate/isopropanol-80C30min/13000g15min2x75percentEtOH/elute50ul')
    op('gap_fill','SourceSI100ulgapfill:50ulDNA/10ulbuffer/8ul10mMdNTP/6ulT4pol/26ulwaterRT40min;main1mMlabelretainedasdiscrepancy',S.EXTENSION,(inp('t4pol','SourceT4DNApolymerase'),nt('gap_nt','SourceSI8ul10mMdNTP;main1mMlabeldiffers,conventionunassigned')));clean('gap_clean','SourceAMPure/elute50ul')
    op('second_strand','SourceE6111S8ulbuffer/4ulenzyme/50ulDNA/18ulwater16C1h',S.MAGIC,(inp('second_strand_stock','SourcehistoricalSecondStrandSynthesispreparation; enzyme/nucleotidestockcomponentsunassembled',K.MOLECULAR,R.STRUCTURAL),));clean('second_clean','SourceAMPure/elute50ulEB')
    used=g.product('library200ng','Source200ngDNA-cDNAcomplexesselected',K.MOLECULAR);rem=g.product('complex_remainder','UntransferredDNA-cDNAstock;mass/fateunassigned',K.MOLECULAR);g.react('library_aliquot','Source200ngaliquot',S.SPLIT,((x,R.TEMPLATE),),(used,rem));g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,'Untransferredstock.'));x=used
    op('trueprep','SourcehistoricalTruePrepTD501libraryconversionbeforebiotincapture;wholekitgeometryunassembled',S.MAGIC,(inp('trueprep_stock','SourceTruePrepDNA-librarypreparationstocks;wholeadapter/enzyme/nucleotidecomponentsunassembled',K.MOLECULAR,R.STRUCTURAL),))
    bound=g.product('biotin_bound','SourceM280retainedbiotin-bridgefragmentlibraryclasses',K.MOLECULAR);unbound=g.product('biotin_unbound','UnretainedM280capturefraction;laterfateunassigned',K.MOLECULAR);g.react('capture','SourcebiotinM280enrichmentBEFOREPCR',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),inp('m280','SourceM280streptavidinbeads',K.SELECTION,R.SELECTION)),(bound,unbound));g.mark_terminal(Terminal(unbound.id,Disposition.UNRESOLVED,'Sourceenrichmentnotcompletecapture.'));x=bound
    op('pcr','SourcePCRenrichment;completeindex/primer/cyclechoiceunassembled',S.MAGIC,(inp('pcr_preparation','SourcePCRpreparationstocksfromlibraryworkflow;wholecomponentsunassembled',K.MOLECULAR,R.STRUCTURAL),));clean('size_select','SourceAMPureXPsize-selection,exactratiosunassembled')
    lib=g.product('library','SourceTaDRIMDNAsequencinglibrary',K.LIBRARY);g.react('load','Sourcesequencerloading',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='SourceHiSeqXTen2x150.Wholeindexedconstruct/primerlanding/genomicjunctions/quantitativeyieldunassembled. ContactDNA-RNA/RNA-RNAreportersnotproofdirectbiologicalbinding.'))
    g.add_source_item(SourceItem('pending','SourceSIgapfill8ul10mMdNTP/100ulversusmain1mMlabelpreserved;completehistoricalkits/sealedchimeras/readlandingunassembled.','Twopublishedoriginals,noconfirmedabsence.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();c('FullpGTn5sourcepair',B,'AAGATGTGTATAAGAGACAG');c('BridgeFmodifiedTisbase',F,'CCTCGATACTTGAGCTGAC');c('BridgeRmodifiedTisbase',RBR,'GTCAGCTCAAGTATCGAG')
    c('NoChRDPETterminalTimport',F.endswith('GACT'),False)
    c('Loadedstockremainderretained',g.terminals['tte_remainder'].disposition,Disposition.UNRESOLVED)
    c('Tagmentusesproducedloadedcomponents',all(g.producer_of(p.material)is not None for p in g.steps['tag'].inputs if p.role in{R.CATALYST,R.STRUCTURAL}),True)
    c('RTprimercomesfromexistingbridge',g.producer_of('retained_bridge_primer').id,'rt_components')
    c('RNAacceptorclassdefaults inferred',rna_geometry().strands['rna'].segs[0].inferred,True)
    c('CapturebeforePCR',g.steps['pcr'].inputs[0].material,'biotin_bound')
    c('Actual200ngstockremainder',g.terminals['complex_remainder'].disposition,Disposition.UNRESOLVED)
    c('Molecularcarrierprecedesprecipitation',g.steps['precipitate'].inputs[0].material=='carrier_add_product'and any(p.material=='glycoblue'and p.role==R.STRUCTURAL for p in g.steps['carrier_add'].inputs),True)
    path=sources_root()/Path(__file__).resolve().parents[1].name/'s41467-024-53534-5_PMC11811046.xml.txt'
    if have(c,Source(path,'https://doi.org/10.1038/s41467-024-53534-5')):
        text=path.read_text();c('OriginalwholeA',A in text,True);c('OriginalwholeB',B in text,True);c('OriginalwholebridgeF',SOURCE_ORDERINGS['bridge_f']in text,True);c('OriginalwholebridgeR',SOURCE_ORDERINGS['bridge_r']in text,True)
    c('Privateincompletegraph',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
