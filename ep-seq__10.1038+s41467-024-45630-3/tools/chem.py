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
def opaque(g,key,label,x,ins=(),kind=S.MAGIC,state=None):
    y=g.product(key+'_product',label+'; whole molecular identities/junctions unassembled',K.MOLECULAR,state=state)
    kw={'reason':'Source-positive operation; whole kit composition, substrate context or product geometry unassembled, not publication absence.'} if kind==S.MAGIC else {}
    g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(ins),(y,),**kw);return y
def cleanup(g,key,label,x,carrier=None):
    y=g.product(key+'_kept','Source retained preparation; recovery/carryover unassigned',K.MOLECULAR);w=g.product(key+'_removed','Source workup fraction',K.WASTE)
    ins=((x,R.TEMPLATE),condition(g,key+'_condition',label))
    if carrier is not None:ins+=(add(g,key+'_carrier',carrier),)
    g.react(key,label,S.CLEANUP,ins,(y,w));g.discard(w,'Physical workup, not perfect molecular selectivity.');return y
def aliquot(g,key,label,x,remainder_label='Untransferred source stock; later fate unassigned'):
    y=g.product(key+'_chosen',label,K.MOLECULAR);r=g.product(key+'_remainder',remainder_label,K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,remainder_label));return y
def qc(g,key,label,x):
    a=g.product(key+'_analytical',label,K.MOLECULAR);k=g.product(key+'_kept','Retained production stock',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(a,k));g.mark_terminal(Terminal(a.id,Disposition.QC,'Analytical aliquot not returned; amount unassigned.'));return k
from chemdraw import feature
from nextera import ADAPTOR_S5,ADAPTOR_S7
from pcr import PCRGeometry
TITLE='EP-Seq: published selected activity-label/FACS-bin sequencing preparation'
AUTHORITY='PublishedPMC10902396main176–206 andSI1Table3originalp15. Twooriginals only; upstreammutagenesis/lookup/PacBio/Zenodoandparallel expressionbranchunassembled, notinherited.'
F_TARGET='GGCCTTTTGCCCTGCAGGC';R_TARGET='GGAGGAGAGTCTTCCTTCGGAGGG'
FORWARDS=[ADAPTOR_S5+x+F_TARGET for x in ('','T','AT','CAT')]
REVERSES=[ADAPTOR_S7+x+R_TARGET for x in ('','T','AT','GAT')]
SOURCE_PRIMERS={**{'F'+str(i+3):v for i,v in enumerate(FORWARDS)},**{'R'+str(i+3):v for i,v in enumerate(REVERSES)}}
def primer_pool(which):
    s=Scene();arr=FORWARDS if which=='F'else REVERSES
    for i,v in enumerate(arr):s.strand(which+str(i+3),[Segment(which+'_source_'+str(i),v,note='Sourcefixed0–3ntread-diversitystagger, NOTamolecularidentifier')])
    return s
def clone_tag_class():
    s=Scene();s.strand('plasmid_region',[Segment('left_flank','[plasmid flank unassembled]',placeholder=True),Segment('source_clone_id','N'*15,placeholder=True,feature=feature('ep.variant_tag','feature_barcode','random',note='Sourcecalls15N UMI; installedinupstreamlibraryfactoryandpropagatedthroughcellgrowth. Clone/variantidentity,notper-NGSmoleculededuplication.')),Segment('right_flank','[plasmid flank unassembled]',placeholder=True)]);return s

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('induced_library','SourcebarcodedDAOxyeastlibraryafter48h20Csurfaceexpression; exactmutant/barcodefactoryandpoolsizeallocationunassembled',K.BIOLOGICAL,primary=True,chemical_state={'identifier_meaning':'Upstream15Nclone/variantbarcodecalledUMIbyauthors; notfreshNGSmoleculeUMI','display':'SourcegeneticallyencodedAga2DAOxenzyme,notafreshenzymestockaddition'})
    x=opaque(g,'primary_stain','Sourceanti6His1:50030minRT inPBS.1percentBSA; productiontubeallocation20Mcellspertube/6tubebatchunassembled',x,(add(g,'anti_his','Sourceanti6HisMA1-21315antibodyphysicalstock'),add(g,'primary_bsa','Actual.1percentBSAcopresentinstainingbuffer; notanotherunrecordedaddition'),condition(g,'primary_conditions','SourcePBS/30minRT')),S.BINDING)
    x=cleanup(g,'primary_wash','Source2xicecoldPBS.1percentBSA; molecularcarrieroftheactualwashstock',x,carrier='Actualfresh0.1percentBSAinprimarywashbuffer')
    x=opaque(g,'secondary_stain','SourceAlexa594goatantimouse1:50030min4C',x,(add(g,'secondary_af594','SourceactualAlexa594-antibodyconjugateA11005,notfreefluorophore'),add(g,'secondary_bsa','Actual.1percentBSAof freshstainingbuffer'),condition(g,'secondary_conditions','Source4C30min')),S.BINDING)
    x=cleanup(g,'secondary_wash','Source2xPBS.1percentBSA/resuspend20Mcellsperml; carriedBSAcompositionknown, recoveryunassigned',x,carrier='Actualfresh0.1percentBSAinsecondarywashbuffer')
    chosen=g.product('activity_25m','Source25millionselectedstainedcells foractivitylabeling',K.MOLECULAR);rest=g.product('stained_stock_remainder','Untransferredstainedyeaststock; parallel expression/negativecontrol/replicate allocationunassembled',K.MOLECULAR);g.react('activity_input','Source25Mactivityassayaliquotfrompreparedstainedstock',S.SPLIT,((x,R.TEMPLATE),),(chosen,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Otherstocknotautomaticallysequencedinselectedactivitybranch.'));x=chosen
    x=opaque(g,'activity_label','Source500cells/ul,20min25Ctyramide/peroxidasecascade',x,(add(g,'tyramide_488','SourceAlexa488tyramidereagent1:200; actualmolecularfluorescentlabel'),add(g,'hrp','Source56.8uMHRPcatalyst',K.CATALYST,R.CATALYST),add(g,'d_alanine','Source130mMD-alaninechemicalsubstrate'),add(g,'alginate','Source.75percentsodiumalginatephysicalpolymermatrixstock'),condition(g,'activity_conditions','SourcePBS25C20min; noquantitativeradicalradius/yield/activityproof')),S.MAGIC)
    x=cleanup(g,'activity_wash','Source13000g3minpellet/2xPBS.1percentBSA+.05percentTweenwash',x,carrier='Actualfresh0.1percentBSAinactivitywashbuffer')
    bins=[g.product('bin'+str(i),'SourceactualactivityFACSgate'+str(i)+'cellpopulation; notpureactivityclass',K.MOLECULAR)for i in range(1,5)];reject=g.product('facs_uncollected','Sourceuncollected/nonsinglet/lostevents; count/fateunassigned',K.MOLECULAR)
    g.react('facs','SourceactivityFACS4bins; samegatedpopulationmergedacrosssortingcycles, cyclecountunassembled',S.SPLIT,((x,R.TEMPLATE),add(g,'facs_receiving','Actualreceivingtube medium percycle:each.5ml2x-Trp/4percentglucose/1percentBSA, totalpreparedtubeallocationunassigned'),condition(g,'gates','SourceFACSMelody488/561nm; firstbin99percentnegativecontrolrange,other3equalfractionsofhigherfluorescence; thresholdvaluesunassigned')),tuple(bins)+(reject,));g.mark_terminal(Terminal(reject.id,Disposition.UNRESOLVED,'PhysicalFACSnotcomputationalreadfilter.'))
    outputs=[]
    for i,part in enumerate(bins,1):
        k='bin'+str(i)
        part=cleanup(g,k+'_pellet','Source4000g10minpellet/supernatantdiscard',part)
        part=opaque(g,k+'_regrow','Source10ml-Trp2percentglucose/100ugmlPenstrep48h30C200rpmregrowth',part,(add(g,k+'_growth_medium','Actualfreshsourcegrowthmedium/nutrientstock; DNAreplication/metabolites/cellcountsnotmodeled'),))
        part=aliquot(g,k+'_50m','Source50millionselectedcells storedin25percentglycerolforDNAprep',part)
        part=cleanup(g,k+'_thaw','SourceRT5min/13000g1minsupernatantdiscard',part)
        part=opaque(g,k+'_wall','Source250ulGeneJETresuspension+20UZymolyase37C2h900rpm',part,(add(g,k+'_zymolyase','Source20UZymolyasewallhydrolysis',K.CATALYST,R.CATALYST),condition(g,k+'_wall_buffer','Source250ulresuspensionbuffer')),S.PERMEABILIZATION)
        part=opaque(g,k+'_plasmid','SourceGeneJETplasmidminiprep/elute15ulwater; wholeplasmid/kitfractionsunassembled',part,(add(g,k+'_miniprep_stock','ActualGeneJETDNAextractionpreparationmolecularstock; fullformulationunassembled'),),state=clone_tag_class())
        bulk=g.input(k+'_pcr_stock','ActualPCRreagentpreparation:25ulNEBNextUltraIIQ5MasterMix+5uleach1uMequimolarF3–6/R3–6pool; sourceDNA15ul separatelycontributes',K.MOLECULAR)
        cat=g.product(k+'_pcr_enzyme','SamepreparedPCRcatalystcomponent',K.CATALYST);nt=g.product(k+'_pcr_nt','SamepreparedPCRnucleotidecomponent',K.MOLECULAR);f=g.product(k+'_forward_pool','Samepreparedsource4-memberforwardpool',K.OLIGO,state=primer_pool('F'));r=g.product(k+'_reverse_pool','Samepreparedsource4-memberreversepool',K.OLIGO,state=primer_pool('R'));buf=g.product(k+'_pcr_buffer','SamepreparedPCRbuffer',K.CONDITION)
        g.react(k+'_pcr_components','Copresentcomponentsofonepreparedreagentstock, notfreshduplicateadditions',S.GENERIC,((bulk,R.TEMPLATE),),(cat,nt,f,r,buf),changes_molecule=False)
        amplified=g.product(k+'_amplicon','Sourceclone-tagampliconmixture; fullprimerlandings/unreactedproductsunassembled',K.MOLECULAR,chemical_state={'identifier_meaning':'Source15NclonevariantID preserved; fixedstaggerbasesnotUMI'})
        g.react(k+'_pcr','SourceENTIRE15ulDNAelute+35ulPCRreagents,18cycles98C10s/72C30s/72C2minthenfinal72C2min, literalprogramme',S.PCR,((part,R.TEMPLATE),(cat,R.CATALYST),(nt,R.STRUCTURAL),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(buf,R.CONDITION)),(amplified,),geometry=PCRGeometry.unavailable(unknown_scene(),reason='Sourcewholeprimersknown; plasmidbarcodeflanks/completeamplicongeometryunassembled.'));part=amplified
        part=cleanup(g,k+'_gel','Source1percentagaroseelectrophoresis/gelDNApurification',part);part=qc(g,k+'_quant','SourceDNAconcentrationanalyticalaliquot; method/volumeunassigned',part);part=cleanup(g,k+'_dna_clean','SourceDNACleanandConcentrator5secondpurification',part)
        part=opaque(g,k+'_index','SourceuniqueNexteraindexaddition; notinventedTn5tagmentation',part,(add(g,k+'_index_preparation','ActualsourceNexteraindexingpreparedmolecularstock; enzyme/wholeprimers/indexvaluesunassembled'),))
        outputs.append(aliquot(g,k+'_pool_input','Selectedindexedbinlibraryaliquot foractualmultiplexpool; amountsunassigned',part))
    pool=g.product('indexed_pool','Sourceindexedactivitybinlibrarypool; otherstudypoolmembershipunassembled',K.MOLECULAR);g.react('pool_bins','Sourceindexed samples pooled; selectedfouractivitybincontributionsretained, noassumedequimolarity',S.POOL,tuple((p,R.TEMPLATE)for p in outputs),(pool,))
    lib=g.product('library','Sourceactivity-binEP-SeqIlluminalibrarypool',K.LIBRARY);g.react('load','SourceNovaSeq6000SE100readout',S.GENERIC,((pool,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='SourceResults82NovaSeq6000SE100; wholeindexprimers/amplicon/readlandingunassembled. Source15Nisupstreamcloneidentity,notper-NGSmoleculeUMI. PacBiolookup andparallel expression/replicatebranchesnotmerged.'))
    g.add_source_item(SourceItem('unassembled','Upstreamnickingmutagenesis/barcodefactory/PacBiolookup/Zenodo8388902ordinaryunreviewed. Selectedactivityreplicateonly; actualsampleindexvalues/allstudypoolmembersandwholemoleculesunassembled.','Proximityfluorescencenotquantitativechemicalyieldproof.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    import re
    c=Check();g=protocol_graph();c('No freshDAOxenzymestock',any('freshDAOx'in m.label for m in g.materials.values()),False);c('Physicaltyramidelabelcontributes',g.materials['tyramide_488'].kind,K.MOLECULAR);c('Actual4FACSbinsplusuncollected',len(g.steps['facs'].outputs),5);c('Cloneidentityroleisfeaturebarcode',clone_tag_class().strands['plasmid_region'].segs[1].feature.role,'feature_barcode')
    for i in range(1,5):
        k='bin'+str(i);c('Wallhydrolysisnotnucleicacidcleavage '+k,g.steps[k+'_wall'].kind,S.PERMEABILIZATION);c('PCRuseswhole15ulsourceelute '+k,g.steps[k+'_pcr'].inputs[0].material,k+'_plasmid_product');c('SourcePCRNTsamepreparedstock '+k,g.producer_of(k+'_pcr_nt').id,k+'_pcr_components')
    p=sources_root()/Path(__file__).resolve().parents[1].name/'s41467-024-45630-3_supp_41467_2024_45630_MOESM1_ESM.pdf.txt'
    if have(c,Source(p,'https://doi.org/10.1038/s41467-024-45630-3')):
        t=re.sub(r'\s+','',p.read_text())
        for k,v in SOURCE_PRIMERS.items():c('OriginalwholeTable3 '+k,v in t,True)
    c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
