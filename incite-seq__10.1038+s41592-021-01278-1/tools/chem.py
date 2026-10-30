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
TITLE='inCITE-seq: associated-manuscript nuclear protein/RNA preparation'
AUTHORITY='AssociatedNIHMS1745895/PMC8532076.1 selectedMethods173–209; publishedMOESM1 reporting-summaryidentityonly. No final-body equivalence or whole modern10X/CITE recipe inherited.'
def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('hela_cells','SourceHeLa sample forselectednuclearp65/HTO/RNAlibrary branch; priorculture/stimulationoutsideboundary',K.BIOLOGICAL,primary=True)
    x=opaque(g,'nucleus_extract','Source2mlEZlysis/RNaseinhibitor,24dounceA+24B/add3ml/ice5min/500g5min4C;5mlrepeatlysis/ice/pellet;1mlPBS3mMMg/RNaseinhibitor/35umfilter',x,(add(g,'extraction_preparation','SourceactualEZlysis/nucleuspreparationstocks; wholeformulation/physicalsubfractionsunassembled'),))
    x=opaque(g,'fix_permeabilize','Source3ml1.33percentFA-NTstock addedto1mlnuclei,10min4Crocking',x,(condition(g,'fa_nt','Source1.33percentformaldehyde/.2percentNP40/.1percentTween/3ulglacialaceticacid stockdescriptor; no silently substituted finalconcentrations'),),S.CROSS_LINKING)
    x=opaque(g,'quench','Source3ul1Mglycine then20umfilter',x,(condition(g,'glycine','Sourceactual3ul1Mglycine, notborrowed125mM'),),S.QUENCHING)
    x=cleanup(g,'fix_pellet','Source850g5min4C pellet/resuspend500ulblockingbuffer',x)
    x=opaque(g,'block','Source500uloptimizedblockingbuffer15min4Crocking',x,(add(g,'block_stock','Sourceactualprepared500ulblockingstock,co-present1percentBSA/FcX1:100/dextransulfate.05percent/Tween.2percent/RNaseinhibitor/PBS; noextraBSAaddition',chemical_state={'co_present_carriers':'Source1percentBSAandFcX inoneactualblockingbufferaddition'}),))
    x=cleanup(g,'block_pellet','Sourcepelletbeforefreshstainingstock',x)
    x=opaque(g,'stain','Source200ulblockingbuffer withp65TotalSeq-A1:400 andnucleushash1:200,1h4Crocking',x,(add(g,'p65_antibody_tag','Sourcep65DNA-conjugatedTotalSeq-Aformatantibody; completeprotein/DNA/linker/featurebarcodeorderingsunassembled'),add(g,'nucleus_hash','Sourcesample-specificBioLegend682213/682215 nucleus-hashingDNA-antibodypreparation; tagvalue/fullorderingunassembled'),add(g,'stain_buffer_stock','Sourceactualfresh200uloptimizedblockingstock includes1percentBSA/FcX/dextransulfate/Tween/RNaseinhibitor; notanotherdoseintotheearlier500ulstock',chemical_state={'co_present_carriers':'Source1percentBSA/FcXcopresentinfreshstainingbuffer'})),S.BINDING)
    x=cleanup(g,'stain_wash','Source2x500ul.2percentPBST+RNaseinhibitor/5min/re-pellet',x)
    x=opaque(g,'resuspend','Source300ulPBS for10X thenPBS3mMMg/10umfilter/count; flowcytometrysecondary/DAPIbranchNOTincluded',x,(condition(g,'resuspend_buffer','SourcePBS/Mg/RNasehandling'),),S.GENERIC)
    own=aliquot(g,'nucleus_load','Source10000selectedstainedhashednuclei foroneHeLav3channel',x)
    other=g.input('other_hashed_nuclei','Actual10000other-treatmentstainedhashedHeLanuclei fromseparatepreparation, notcopyofselectedsample; upstreamstockallocationunassembled',K.BIOLOGICAL,primary=True)
    pool=g.product('nucleus_pool','Sourceonechannel10000NT+10000TNFhashednuclei',K.MOLECULAR);g.react('pool_nuclei','Actualtwo-sampleHeLa nucleuspool',S.POOL,((own,R.TEMPLATE),(other,R.TEMPLATE)),(pool,));x=pool
    x=opaque(g,'gem','SourceChromium3primev3GEM;53C45minreversecrosslink/RNARTandDNA-tagcopyingthen85C5min;GEMrecovery',x,(add(g,'gem_preparation','Sourceactualhistorical3primev3beads/RT/structuralnucleotidesandothermolecularstocks; fullformulation/oligosunassembled'),))
    x=opaque(g,'cdna_amp','Sourcev3cDNAamplificationwithHTOadditive.1uM/ADTadditive.2uM',x,(add(g,'amp_preparation','Sourcev3amplificationpreparedstock+HTO/ADTadditiveprimers; fullenzyme/NT/wholeprimersunassembled'),))
    wta=g.product('wta_beads','Source.6xSPRIbead-boundWTAenrichedfraction; purityunassigned',K.MOLECULAR);tags=g.product('tag_supernatant','SourcesupernatantHTO/ADT-derivedDNAcopyenrichedfraction; purityunassigned',K.MOLECULAR)
    g.react('size_split','Source.6xSPRI5minRT physicalbeads/supernatant separation',S.SPLIT,((x,R.TEMPLATE),condition(g,'size_conditions','SourceSPRIselect')), (wta,tags))
    wta=cleanup(g,'wta_wash','SourceWTAbeads2x80percentethanol/elute40ulEB',wta)
    a=g.product('wta_ten_ul','Source10ulof40ulWTAforGEXconstruction',K.MOLECULAR);r=g.product('wta_thirty_ul','Source30ulWTAstockremaining',K.MOLECULAR);g.react('wta_aliquot','ActualWTA10of40ulaliquot',S.SPLIT,((wta,R.TEMPLATE),),(a,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,'UntransferredWTAstock.'))
    a=opaque(g,'gex_library','Sourcehistoricalv3enzymaticfragmentation/adapterligation/sampleindex/elute30ulEB',a,(add(g,'gex_preparation','Sourceactualhistoricalv3GEXmolecularconstructionstocks; fullends/primersunassembled'),))
    tags=cleanup(g,'tag_cleanup','Source1.4xSPRI5minRT/5minmagnet/2x80percentethanol/elute24ulwater',tags)
    h=g.product('hto_six_ul','Source6ulof24ultagstockforHTOindexPCR',K.MOLECULAR);d=g.product('adt_six_ul','Source6ulof24ultagstockforADTindexPCR',K.MOLECULAR);r=g.product('tag_twelve_ul','Source12ultagstockremainder',K.MOLECULAR)
    g.react('tag_aliquots','Actual6ulHTO+6ulADT aliquots fromsame24ulstock',S.SPLIT,((tags,R.TEMPLATE),),(h,d,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,'Actualunusedtagstock.'))
    libs=[('gex',a)]
    for name,part in [('hto',h),('adt',d)]:
        part=opaque(g,name+'_index_pcr','Source6ultemplate,unique'+name+'indexprimermix/NEBNext2xMastermixM0541L;21cycles98C2s/72C15s',part,(add(g,name+'_pcr_stock','SourceactualNEBNextMastermix+indexprimerpreparedreaction; whole molecularcomponents/indicesunassembled'),))
        part=cleanup(g,name+'_cleanup','Source2xSPRI10min/5minmagnet/2xethanol/elute20ulEB',part);libs.append((name,part))
    chosen=[]
    for name,part in libs:
        part=qc(g,name+'_qc','SourceQubit/TapeStationanalyticalaliquots,amountunassigned',part);chosen.append(aliquot(g,name+'_pool_input','Selected'+name+'libraryaliquotforactualNGSpool;volumeunassigned',part))
    pool=g.product('library_pool','ActualGEX+HTO+ADTselectedlibrarypool',K.MOLECULAR);g.react('pool_libraries','Sourcephysicalpoolthreeactuallibraryclasses; noassumedequimolarratio',S.POOL,tuple((x,R.TEMPLATE)for x in chosen),(pool,))
    lib=g.product('library','SourceinCITE-seqsequencingpool',K.LIBRARY);g.react('loading','SourceNextSeq500loading',S.GENERIC,((pool,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='AssociatedAMNextSeq500source75cyclekitdescriptorandR1:28/I1:8/R2:55retainedwithoutnormalization. WholeTotalSeq-A/bead/primer/index/readlandingsunassembled. AntibodyDNAreportercopiesarenotendogenousRNA.'))
    g.add_source_item(SourceItem('unassembled','AssociatedAMrecipe;publishedreportingsummaryidentitychecked. Fullfinalbody,exacthistorical3primev3/CITE-A/nucleushashingcomponentorderingsordinaryunreviewed.','No full-sourceabsenceorchemistryvalidationclaim.'));return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph();c('ActualGEX10of40input',g.steps['gex_library'].inputs[0].material,'wta_ten_ul');c('WTA30ulremainder','wta_thirty_ul'in g.terminals,True);c('Tag12ulremainder','tag_twelve_ul'in g.terminals,True);c('Two6ultagaliquots',len(g.steps['tag_aliquots'].outputs),3);c('HTOusesitsownaliquot',g.steps['hto_index_pcr'].inputs[0].material,'hto_six_ul');c('ADTusesitsownaliquot',g.steps['adt_index_pcr'].inputs[0].material,'adt_six_ul');c('Sourcebothnuclearpopulationscontribute',len(g.steps['pool_nuclei'].inputs),2);c('Threeactuallibraryclassespooled',len(g.steps['pool_libraries'].inputs),3);c('Sequencingstockremaindersretained',all(n+'_pool_input_remainder'in g.terminals for n in ('gex','hto','adt')),True);c('NoFACSselectionborrowed',all(s.kind!=S.SELECTION for s in g.steps.values()),True);c('KnownBSAinstocknotjustcondition',g.materials['block_stock'].kind,K.MOLECULAR);c('Privatepartialgraph',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
