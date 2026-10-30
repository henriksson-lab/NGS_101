"""Published2014 SSB-Seq nick-translation/DIG capture branch, private partial model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from illumina import P7,TRUSEQ_READ1
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
FORWARD=TRUSEQ_READ1[:22]
REVERSE=P7[:22]+'C'  # source terminalC; not a full canonicalP7 or a silent correction.
AUTHORITY='PublishedPMC4139894main54–72/94–120; originalFigure1aPDFp3. PubMed25056547identity only. No later2018chapter inheritance.'
def oligo(seq):
    s=Scene();s.strand('primer',[Segment('whole source primer',seq)]);return s

def a_tail_scene():
    # One source-known local A-tailed end, not an assembled complete fragment.
    s=Scene();s.strand('top',[Segment('fragment body','XXXX',placeholder=True,note='Schematic body only; length/sequence unassigned'),Segment('source protruding dA','A')],mod3='OH')
    s.anneal('bottom',[Segment('body complement','xxxx',placeholder=True,note='Complement of same unknown schematic body')],to='top',pair=('body complement','fragment body'))
    s.overhang('top','source protruding dA',polarity='3-prime',owner='substrate',text='Source3primeA at this local end; opposite end/full substrate unassembled');return s

def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('whole source products','[unassembled]',placeholder=True)]);return s

def protocol_graph():
    g=ProtocolGraph('SSB-Seq2014 selected DIG branch',authority=AUTHORITY,representation='published_instance')
    x=g.input('cells','SourceHCT1168x10^7cells; cultureDMEM/10percentFCS outside selected molecular boundary',K.BIOLOGICAL,primary=True)
    def inp(k,label,kind=K.CONDITION,role=R.CONDITION,**kw):return(g.input(k,label,kind,**kw),role)
    def op(k,label,kind,inputs=(),product='Source mixture including incomplete/off-target/unreacted classes; recovery/carryover unassigned',**kw):
        nonlocal x
        y=g.product(k+'_product',product,K.MOLECULAR)
        if kind==S.GENERIC:kw['changes_molecule']=False
        g.react(k,label,kind,((x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);x=y
    def clean(k,label):
        nonlocal x
        y=g.product(k+'_kept','Actual source retained workup fraction; no perfect purification assumed',K.MOLECULAR);w=g.product(k+'_removed','Actual separated workup fractions',K.WASTE)
        g.react(k,label,S.CLEANUP,((x,R.TEMPLATE),inp(k+'_conditions',label)),(y,w));g.discard(w,'Source removed fractions, molecular recovery/carryover unassigned');x=y
    def split(k,label,branches):
        nonlocal x
        yy=tuple(g.product(k+'_'+name,desc,K.MOLECULAR)for name,desc,disp in branches);g.react(k,label,S.SPLIT,((x,R.TEMPLATE),),yy)
        for y,(_,desc,disp)in zip(yy,branches):
            if disp is not None:g.mark_terminal(Terminal(y.id,disp,desc))
        x=yy[0];return yy
    op('etoposide','SourceMethods3.1:20uMetoposide37C20min; Figure5caption5min is a different source statement',S.CHEMICAL_CONVERSION,(inp('etoposide_condition','Sourceetoposidetreatment; actual accessibleTop2associatedbreakclasses unassembled'),),product='Source treatedcells; heterogeneous DNA/protein complexes, not a pure set of all genomic SSBs')
    clean('cell_wash','SourcetwiceicecoldPBS')
    op('lyse','Source10ml10mMTris8/100mMEDTA/.5percentSDS',S.MAGIC,(inp('lysis_conditions','Sourcelysisbuffermix'),),reason='SourceSDSlysispositive; fullcellularsoluble/insolublecomponentscopresentandunassembled, notaseparatedpureDNAproduct.')
    op('pk1','Source200ugmlProteinaseK52Covernight',S.DEPROTEINIZATION,(inp('pk1_enzyme','SourceRocheProteinaseK',K.CATALYST,R.CATALYST),inp('pk1_conditions','Source52Covernight')))
    clean('extract1','Source2xphenol/1xphenolchloroform; ethanolprecipitationwith2Mammoniumacetate')
    op('rnase','Source5ugpancreaticRNase37C1h',S.DIGESTION,(inp('rnase_enzyme','SourceRochepancreaticRNase',K.CATALYST,R.CATALYST),inp('rnase_conditions','Source37C1h')))
    op('pk2','Sourceadjust.5percentSDS/200ugmlProteinaseK55C1h',S.DEPROTEINIZATION,(inp('pk2_enzyme','SourceRocheProteinaseK secondactualaddition',K.CATALYST,R.CATALYST),inp('pk2_conditions','Source.5percentSDS55C1h')))
    clean('extract2','Source2xphenol/ethanol2Mammoniumacetate; resuspend1mlTE/RT12hgentlerotation')
    split('qc','Actualgel/NanoDropmeasurements separate fromproductionDNAstock',[('production','RetainedsourceHMWDNAstock',None),('gel','ActualgelsizeQCaliquot; amountunassigned',Disposition.QC),('nanodrop','ActualNanoDropND1000QCaliquot; amountunassigned',Disposition.QC)])
    split('label_input','Selectsource500ugDNAforDIGlabelreaction; controlreaction separate',[('selected','Source500ugDIG-branchDNA',None),('stock','UntransferredHMWDNAstock, mass/fateunassigned',Disposition.UNRESOLVED)])
    op('nick_translation','Source1000U E.coliDNApolI16C40s withDIG/chainterminatornucleotides',S.MAGIC,
       (inp('pol1','SourceNEBE.coliDNApolymeraseI',K.CATALYST,R.CATALYST),inp('canonical_nts','Source200uM dATP/dGTP/dCTP; printed per-base/total convention unassigned',K.MOLECULAR,R.STRUCTURAL),inp('dig_dutp','Source20uMDigoxigenin11dUTP, covalentnucleotide label',K.MOLECULAR,R.STRUCTURAL),inp('dd_nts','Source117uM ddATP/ddGTP/ddCTP, actualchainterminators; printed per-base/total convention unassigned',K.MOLECULAR,R.STRUCTURAL),inp('nick_conditions','Source16C40s, fullbuffer/unreactedfractionsunassembled')),
       reason='Source nicktranslation positive; occupied3OHnick geometry, DIGpatch lengths/DDtermination products and residue classes notassembled. Labeling does not report every SSB or prove quantitativeTop2specificity.',product='SourceDNApopulations with DIG-labelled patches/dd-terminated and unreacted molecules; exact product ends/polymers unassembled')
    op('stop_label','Source50uMEDTAstop; retain sourceunit literal',S.QUENCHING,(inp('stop_conditions','Source50uMEDTA'),))
    clean('label_clean','Sourcephenol then2xethanolprecipitation with2.5Mammoniumacetate')
    op('sonicate','SourceBioruptor30x30smedium/30siceintervals; mean200–400bp',S.FRAGMENTATION,(inp('sonication','SourceBioruptor physicalfragmentation'),))
    op('anti_dig','Source10ugantiDIGovernight4Cgentlerotation',S.BINDING,(inp('anti_dig_ab','SourceRocheantiDigoxigenin antibody, physicalboundcomponent',K.SELECTION,R.SELECTION),))
    op('protein_g','Source60ulProteinGSepharose4h4C',S.BINDING,(inp('protein_g_beads','SourceRocheProteinGSepharose physicalcapturematerial',K.SELECTION,R.SELECTION),))
    clean('capture_wash','Sourcebeadretention:1xPBS/3xNP40buffer/2xTE;10minagitate+1500g4minperwash; resuspend200ulTE')
    op('release','Source.5percentSDS/200ugmlProteinaseK65Covernight',S.DEPROTEINIZATION,(inp('release_pk','SourceProteinaseKactualaddition',K.CATALYST,R.CATALYST),inp('release_conditions','SourceSDS65Covernight')))
    clean('capture_dna','SourceQIAquickPCRpurification/quantify; exactmanufacturerinternalstocks unassembled')
    op('end_repair','SourceEpicentreENDRepair45minRT',S.END_REPAIR,(inp('end_enzyme','Source1ulT4DNApolymerase+T4PNKmix',K.CATALYST,R.CATALYST),inp('end_nts','Source.25mMeachdNTP+1mMATP',K.MOLECULAR,R.STRUCTURAL),inp('end_conditions','Source33mMTrisacetate7/66mMpotassiumacetate/10mMMgacetate/.5mMDTT')),
       product='Sourceend-repairedDNApopulation; intendedbluntclass, fullends/DDterminatorprocessing/yieldunassigned')
    clean('end_clean','Sourcepurification afterENDRepair; exactworkupunassembled')
    op('a_tail','Source15UKlenowexoMinus37C30min/.2mMdATP',S.TAILING,(inp('klenow','SourceKlenowexoMinus',K.CATALYST,R.CATALYST),inp('tail_nt','Source.2mMdATP',K.MOLECULAR,R.STRUCTURAL),inp('tail_conditions','Source37C30min')),
       product='SourceA-tailedfragmentclasses pluspartial/unreacted substrates; one local3A endrepresented separately, fullfragment unassembled')
    op('ligate','Source.1ulIlluminaadaptoroligomix/1000UT4ligaseRT30min',S.LIGATION,(inp('adapter','SourceactualIlluminaAdaptoroligomix; historicalwholeoligos/endsunassembled',K.ADAPTER,R.STRUCTURAL),inp('ligase','SourceT4DNAligase1000U',K.CATALYST,R.CATALYST),inp('ligation_conditions','SourceRT30min; fullbufferunassembled')),
       geometry=AdapterLigationGeometry('unavailable','unavailable',(a_tail_scene(),),'Sourceknown substrate3A retained; adapteroligobody/endpolarityandfulltargetgeometry unassembled, notconfirmedabsent.'))
    clean('ligation_clean','SourceQIAquickPCRpurification')
    clean('pre_pcr_gel','Source2percentEgel/200–400bpbandexcise/MinEluteelute12ul; retainedbandandactualrejectedgel separate')
    pcrstock=g.input('pcr_stock','Sourceprepared18cyclePCRreactionstock; exactstockquantitiesunassembled',K.MOLECULAR)
    f=g.product('pcr_forward','SourcewholeFwprimer',K.OLIGO,state=oligo(FORWARD));rev=g.product('pcr_reverse','SourcewholeRvprimer',K.OLIGO,state=oligo(REVERSE));e=g.product('pcr_enzyme','SourcePCRpolymerasecomponent; identityunassembled',K.CATALYST);nts=g.product('pcr_nts','SourcedNTPcomponentofsamePCRstock',K.MOLECULAR);buf=g.product('pcr_buffer','SourceconditioncomponentofsamePCRstock',K.CONDITION)
    g.react('pcr_components','Representcomponentsofsamepreparedstock, notfreshduplicateadditions',S.GENERIC,((pcrstock,R.TEMPLATE),),(f,rev,e,nts,buf),changes_molecule=False)
    op('pcr','Source18cycles98C30s/65C30s/72C30s',S.PCR,((f,R.STRUCTURAL),(rev,R.STRUCTURAL),(e,R.CATALYST),(nts,R.STRUCTURAL),(buf,R.CONDITION)),geometry=PCRGeometry.unavailable(unknown_scene(),reason='Whole historicaladapter/primerlanding and amplifiedDNAproducts unassembled; sourcewholeprimersretained.'))
    clean('post_pcr_gel','Source2.5percentagarose/220–500bpbandexcise/Qiagengelpurification')
    split('load_aliquot','Sourceclusterloadingallocation; amountunassigned',[('selected','ActualGenomeAnalyzerloadingaliquot',None),('remainder','Untransferredsourcepurifiedlibrarystock; fateunassigned',Disposition.UNRESOLVED)])
    lib=g.product('library','SelectedSSB-Seq2014productionlibrary',K.LIBRARY);g.react('load','SourceIlluminaGenomeAnalyzerclustergeneration/sequencing',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='Sourcealigned36bpreads; precisephysicalrunmask/primerlanding/indexstructure unassigned. NoUMI inferred fromcomputationalduplicate removal. DIGlabelpatchesnotbaseexactbreakpositions. SeparateDSBTdT/biotin/S1route anddTTPnegativecontrolnotmerged.'))
    g.add_source_item(SourceItem('source_variant','MainMethods3.1 etoposide20min versusFigure5caption5min; graphselectsMethods, intendedperfigureexperiment unresolved.','PublishedPMC4139894main94/92; no silentdose reconciliation.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();c('SourceDIGbeforephysicalshear',g.steps['sonicate'].inputs[0].material,'label_clean_kept');c('AntiDIGbindsbeforeProteinG',g.steps['protein_g'].inputs[0].material,'anti_dig_product');c('AffinitywashprecedescaptureDNArelease',g.steps['release'].inputs[0].material,'capture_wash_kept');c('PCRusesactualgelband',g.steps['pcr'].inputs[0].material,'pre_pcr_gel_kept');c('Source500ugselectedafterQC',g.steps['label_input'].inputs[0].material,'qc_production');c('WholeRvsource23nt',len(REVERSE),23)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'ijms150713111_PMC4139894.xml.txt'
    if have(c,Source(p,'https://doi.org/10.3390/ijms150713111')):
        t=''.join(p.read_text().split());c('PublishedwholeFwsource',FORWARD in t,True);c('PublishedwholeRvsource',REVERSE in t,True);c('SourceDIG/terminatorinventory','digoxigenin-11-dUTP' in t and 'ddATP' in t and '1000units' in t,True)
    c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
