"""Published selected preparation; private partial molecular model, not a complete audit."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from chemdraw import Scene,Segment
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from circular import CircularizationGeometry

TITLE='sci-L3-WGS selected cultured-cell DNA route'
AUTHORITY='AssociatedAMPMC6874760 main276–298/Fig1Bcaption764/Fig2caption770;1original. Finalpublished/SI/TableS5oligowhitelistsunreviewed/gated. RNA-DNAcoassayandtargetedbranchseparate.'
INITIAL='Sourceculturedcell suspension, selectedDNA-onlyWGSroute; notgerm-cell-specificextractionorRNA/DNAcoassay'
READOUT='SourceIlluminaWGS library with threecellbarcode rounds; third-roundUMI tags individualIVTtranscripts ratherthanoriginalgenomicmolecules. SourcepairedendadapterorientationcanreverseP5/P7duetosymmetricligation; exactrunlengths/primerlandingunassembled. Transcript/germlineanalysisneverphysicalDNAremoval.'
DEFERRED='Publishedfinal/associatedSI/TableS5original andhistoricalkitfinish;bc1widthgarbledasxrefXXinJATS, notinvented. RNAprimer/foldbackRT andone-versus-two-sidedligation remain explicit unassembledsubstrateclasses.'

def inventory_state(label):return None
def inventory_chemistry(label):return None

from nextera import ME,ME_RC
from chemdraw import feature
SP2='GGGATGCAGCTCGCTCCTG'
RNA_PRIMER=ME.replace('T','U')
def inventory_chemistry(label):
    if label=='SourceTn5short5P-MErc':return {'polymer':'DNA','five_prime':'P','sequence':ME_RC}
    if label.startswith('SourceTn5long'):return {'polymer':'DNA','five_prime':'P','bc1_width':'unassigned from garbledJATS; TableS5unread'}
    if label=='SourceexternalRNA RTprimer ME-as-RNA':return {'polymer':'RNA','sequence':RNA_PRIMER}
    return None
def inventory_state(label):
    sc=Scene()
    if label=='SourceTn5short5P-MErc':sc.strand('ME_rc',[Segment('source short mosaic end',ME_RC)],mod5='P');return sc
    if label.startswith('SourceTn5long'):sc.strand('bc1_Tn5',[Segment('ligation pad','GTCTTG'),Segment('bc1','[bc1 sourcewidthunassembled]',placeholder=True,feature=feature('sci_bc1','cell_barcode','unknown')),Segment('mosaic end',ME)],mod5='P');return sc
    if label=='Sourcebc2hairpin T7loop inventory':
        sc.strand('hairpin',[Segment('ligationoverhang','CAAGAC'),Segment('bc2 reverse','n'*7,placeholder=True,feature=feature('sci_bc2','cell_barcode','unknown')),Segment('sp2reverse','CAGGAGCGAGCTGCATCCC'),Segment('sourceT7loop','AATTTAATACGACTCACTATA'),Segment('sp2',SP2),Segment('bc2','N'*7,placeholder=True,feature=feature('sci_bc2','cell_barcode','unknown'))]);return sc
    if label=='SourceexternalRNA RTprimer ME-as-RNA':sc.strand('RNAprimer',[Segment('RNA RT primer',RNA_PRIMER,note='RNA, not DNA; source explicitrA/rU/rG/rC ordering')]);return sc
    if label=='SourceSSSprimer N4UMI N6bc3 sp2':sc.strand('SSS',[Segment('IVTtranscriptUMI','N'*4,placeholder=True,feature=feature('sci_umi','umi','random')),Segment('bc3','N'*6,placeholder=True,feature=feature('sci_bc3','cell_barcode','unknown')),Segment('sp2',SP2)]);return sc
    return None

SOURCE_ORDERINGS={
'bc2_ligate_oligo_0':{'source':'j.molcel.2019.08.002_PMC6874760.xml.txt:288','raw':'1. CAAGAC 2. Y’Y’Y’Y’Y’Y’Y’ [reverse complement of 2nd round barcode] 3. CAGGAGCGAGCTGCATCCC 4. AATTTAATACGACTCACTATA 5. GGGATGCAGCTCGCTCCTG 6. YYYYYYY [2nd round barcode]','modifications':'No extra hairpin end modification asserted; full reaction ends unassembled.'},
'sss_oligo_0':{'source':'j.molcel.2019.08.002_PMC6874760.xml.txt:296','raw':'NNNN [UMI] ZZZZZZ [3rd round barcode] GGGATGCAGCTCGCTCCTG','modifications':'No extra modification stated in this ordering; N4 tags IVTtranscripts, Z6 represents cellbarcode. Concrete codes unassigned.'}}
SOURCE_OLIGO_EXCEPTIONS={'bc2_ligate_oligo_0':'Associatedmain288 prints six numbered hairpin components with annotations, not one contiguous oligo; sourcepart fixture checked locally, wholeligationgeometryunassembled.','sss_oligo_0':'Associatedmain296 prints adjacent N4UMI and Z6bc3 with bracket annotations; local fullpattern fixture preserves each distinctfeature andsp2.'}
def unknown_scene():
    sc=Scene();sc.strand('unassembled',[Segment('geometry','[source geometry unassembled]',placeholder=True)]);return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('input',INITIAL,K.BIOLOGICAL,primary=True)
    def op(key,label,kind,reagent,product,oligos=()):
        nonlocal x
        enzyme=kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.CIRCULARIZATION,S.LIGATION,S.PCR,S.DEPROTEINIZATION,S.TAGMENTATION,S.TRANSCRIPTION,S.REVERSE_TRANSCRIPTION,S.CHEMICAL_CONVERSION}
        c=g.input(key+'_stock',reagent,K.CATALYST if enzyme else K.CONDITION)
        inputs=[(x,R.TEMPLATE),(c,R.CATALYST if enzyme else R.CONDITION)]
        for i,a in enumerate(oligos):
            o=g.input(key+'_oligo_'+str(i),a,K.OLIGO,state=inventory_state(a),chemical_state=inventory_chemistry(a));inputs.append((o,R.STRUCTURAL))
        kw={'changes_molecule':False} if kind==S.GENERIC else {}
        if kind==S.MAGIC:kw['reason']='Actual published component; full source-specific molecular transition unassembled.'
        if kind==S.CIRCULARIZATION:kw['geometry']=CircularizationGeometry.unavailable('Source intramolecular closure known; full substrate and closure geometry unassembled, not a claim of source absence.')
        if kind==S.LIGATION:kw['geometry']=AdapterLigationGeometry('unavailable','unavailable',(unknown_scene(),),'Source ligation operation known; full ends and sealed junction unassembled.')
        if kind==S.PCR:kw['geometry']=PCRGeometry.unavailable(unknown_scene(),reason='Source primer inventories/whole amplicon geometry not assembled in this selected model.')
        y=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,inputs,(y,),**kw);x=y
    def clean(key,label):
        nonlocal x
        c=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','Physical source recovered preparation; carryover and recovery unquantified',K.MOLECULAR);w=g.product(key+'_waste','Actual source physically separated fractions',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(y,w));g.discard(w,'Physical workup, not a perfect molecular specificity assertion.');x=y
    def aliquot(key,label):
        nonlocal x
        y=g.product(key+'_used',label,K.MOLECULAR);r=g.product(key+'_rest','Physicalsource remainder; mass/fateunassigned',K.MOLECULAR);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,'Not implicitly discarded or counted as sequenced.'));x=y
    def qc(key,label):
        nonlocal x
        y=g.product(key+'_retained','Production preparation retained after independent analytical allocation',K.MOLECULAR);q=g.product(key+'_assay',label,K.QC);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,q));g.mark_terminal(Terminal(q.id,Disposition.QC,'Analytical aliquot; not returned to production.'));x=y
    op('fix','Source1.5percentformaldehydeRT10min',S.CROSS_LINKING,'10mlDMEM/406ul37percentformaldehyde','SourcefixedcellswithDNA/chromatin/proteins; crosslinkoccupancyunassigned')
    op('quench','800ul2.5MGlycineice5min',S.QUENCHING,'Sourceglycine','Quenchedcellpreparation')
    clean('cell_wash','Sourcepellet/wash1mllysisbuffer60mMTrisAc8.3/2mMEDTA/15mMDTT')
    op('nuclei','Source1mllysisbuffer0.1percentIGEPALice20min',S.MAGIC,'IGEPALlysis;completefractionsunassembled','Sourcenucleipreparationpluscellremnants')
    clean('nuclei_wash','Sourcepellet/washNEBuffer2.1')
    op('xsds','Source800ulNEBuffer2.1/0.3percentSDS42C30min500rpm',S.MAGIC,'SDSnucleosome-depletioncomponent','Sourcenucleosome-depletednucleiwithresidualDNA/proteins/remnants; notperfectnakedDNA')
    op('triton','180ul10percentTritonX42C30min500rpm',S.QUENCHING,'SourceTritonX','SourcequenchedSDStreatednuclei')
    clean('post_xsds','Two1mllysisbufferwashes/resuspend20000nucleiul')
    op('bc1','1.5ulnuclei/6.5ulwater/0.7ul50mMMgCl2/1.2ultransposome55C20min',S.TAGMENTATION,'LucigenTn50.2uMstock; source24first-roundwells','Sourcebc1-tagmentedDNAwithinnucleiclasses plusunprocessedsubstrates',('SourceTn5short5P-MErc','SourceTn5long5P-GTCTTG-bc1-ME'))
    op('stop1','Source5ulstop40mMEDTA/1mMspermidine',S.QUENCHING,'EDTA/spermidine','Stoppedbc1nucleipreparation')
    op('pool1','Sourcephysicallypoolwellsintrough',S.GENERIC,'Wholewellclasspool; exactindividualinputsaggregated','Sourcefirst-roundphysicalnucleipool')
    clean('pool1_wash','Add1mllysisbuffer/pellet/resuspend312uldNTP/MgTrisbuffer/water/lysisbuffer')
    op('redistribute','Source4.7ulnucleimixpernewwell',S.MAGIC,'Source64second-roundwells;actualpartitionoccupanciesunassembled','Sourceredistributednucleiclasses')
    op('bc2_ligate','0.8ul1.5uMhairpin/1.18ulligasemix20C30min',S.LIGATION,'Source0.6ul10xligasebuffer/0.48ulPEG4000/0.1ulT4ligase','Sourcebc2-hairpinligatedDNAclasses;one-sided/two-sided/unligatedclassesremain',('Sourcebc2hairpin T7loop inventory',))
    op('stop2','Add4ulstopsolution',S.QUENCHING,'SourceEDTA/spermidine','Stoppedligationpreparation')
    op('pool2_sort','Sourcepool~630ul/DAPI5ugml;sort100–300nucleiperwellinto3ullysisbuffer',S.MAGIC,'SourceFACS;wholeunselectedsortingfractionnotreconstructed','Sourcechosen100–300nucleiperwellplusunassignedsortedwellclasses;notasingle-nucleuswell')
    aliquot('chosen_sorted_well','Sourcechosenwell100–300nuclei;otherwells/uncollectednucleiremainunassigned')
    op('lysis','Source75C45min/4C;Qiagenprotease2mgml55C8h/75C30min',S.DEPROTEINIZATION,'SourceQiagenprotease','SourceDNAlysateandprotein/RNAremnants; optional-80Cstorage')
    op('gap','SourceBstWarmStart2.0/68C5min',S.EXTENSION,'2ulwater/0.7ul10xtagmentbuffer/0.35ul10mMdNTP/0.35ulBst','Sourcegap-extendedpromoterduplexclasses;one/two/noT7ligationclassesremainunassembled')
    op('ivt','Source20ulHiScribeT7QuickHighYield37C10–16h',S.TRANSCRIPTION,'2ulwater/2ulT7Polmix/10ulrNMPmix','SourceRNAampliconsandtemplateDNA/reagents;noexactlinearcopyyieldguarantee')
    op('stop_ivt','Add2.2ul0.5MEDTA',S.QUENCHING,'EDTA','StoppedRNA/DNAreactionmixture')
    clean('rna_clean','SourceRCC5/elute18ul0.1xTE')
    rp=g.input('rna_rt_primer','SourceexternalRNA RTprimer ME-as-RNA',K.OLIGO,state=inventory_state('SourceexternalRNA RTprimer ME-as-RNA'),chemical_state=inventory_chemistry('SourceexternalRNA RTprimer ME-as-RNA'))
    heat=g.input('rna_primer_anneal','Source0.6ulRNAprimer/dNTP/SUPERaseIn;70C1min90C20sec/ice',K.CONDITION)
    annealed=g.product('rt_anneal_product','Source RNA template component after thermal annealing; co-retained external RNA primer tracked separately in same tube; exact substrate classes unassembled',K.MOLECULAR)
    retained=g.product('retained_rna_rt_primer','Same-tube external RNA primer component retained after annealing; not a physical fraction or new primer addition',K.OLIGO,state=inventory_state('SourceexternalRNA RTprimer ME-as-RNA'),chemical_state=inventory_chemistry('SourceexternalRNA RTprimer ME-as-RNA'))
    g.react('rt_anneal','Source RNA primer present before annealing; template and primer outputs are components of the same tube',S.THERMAL_TREATMENT,((x,R.TEMPLATE),(rp,R.STRUCTURAL),(heat,R.CONDITION)),(annealed,retained))
    x=annealed
    e=g.input('ssiv','SourceSSIV/5xbuffer/DTT/SUPERaseIn',K.CATALYST)
    rt=g.product('rt_product','Source cDNA/RNA hybrids and earlier substrates; external-primer and possible foldback substrate classes remain unassembled',K.MOLECULAR)
    g.react('rt','Source SSIV55C15min60C10min65C12min70C8min75C5min80C10min; consumes co-retained external RNA primer, not a new stock addition or asserted intramolecular primer',S.REVERSE_TRANSCRIPTION,((x,R.TEMPLATE),(retained,R.STRUCTURAL),(e,R.CATALYST)),(rt,))
    x=rt
    op('rna_digest','Source0.5ulRNaseH/0.3ulRNaseA37C30min',S.DIGESTION,'RNaseH/RNaseA','SourcecDNAandRNAtranscript/primerfragments; noinstantaneousphysicalRNAremoval')
    op('sss','SourceQ5secondstrandsynthesis; programmeunassignedinselectedmain',S.EXTENSION,'27ulwater/20ul5xQ5buffer/20ulGCenhancer/1ulQ5/1ulSSSprimer','Sourcebc3/UMI-taggeddsDNAclassesandRNA/reagentremnants',('SourceSSSprimer N4UMI N6bc3 sp2',))
    clean('sss_clean','SourceDCC5')
    op('library_kit','SourceNEBNextUltraIIWGS A-tailing/adaptorligationcomponent',S.MAGIC,'Historicalkitinternals/loadedadaptersunreviewed','Sourceadapter-bearingWGSprecursors; notthetargetedP5P7PCRorRNA-DNAcoassay')
    op('pcr','Sourceminimal3PCRcycles forWGSsequencingadaptors',S.PCR,'HistoricalNEBNextUltraII amplificationcomponent','SourceNGSlibrarycopies; cellbarcode/IVTtranscriptUMImeaningretained',('Historicalkitforwardprimeridentityunassembled','Historicalkitreverseprimeridentityunassembled'))
    library=g.product('library','Source sequencing preparation; exact full-library geometry remains unassembled',K.LIBRARY)
    g.react('loading','Physical source loading preparation',S.GENERIC,((x,R.TEMPLATE),),(library,),changes_molecule=False)
    g.sequence(library,SequencingPlan(limitation=READOUT))
    g.add_source_item(SourceItem('pending',DEFERRED,'Selected original sections read; linked-source/whole molecular geometry review incomplete. No accepted absence claim.'))
    return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('Hairpinligation followsactualredistribution',g.steps['bc2_ligate'].inputs[0].material,'redistribute_product')
    c('IVTafterBstgapextension',g.steps['ivt'].inputs[0].material,'gap_product')
    c('RNAcleavageafterRT',g.steps['rna_digest'].inputs[0].material,'rt_product')
    c('ThirdbarcodeaddedafterRNAcleavage',g.steps['sss'].inputs[0].material,'rna_digest_product')
    c('WGSkitsnotreplacedbytargetedPCR',g.steps['pcr'].inputs[0].material,'library_kit_product')
    c('RNAprimer remains RNA notDNA',inventory_chemistry('SourceexternalRNA RTprimer ME-as-RNA')['polymer'],'RNA')
    c('RNAprimer introduced beforeheating',g.steps['rt_anneal'].inputs[1].material,'rna_rt_primer')
    c('SourceSSSfixedarmtranscription',SP2,'GGGATGCAGCTCGCTCCTG')
    c('SourceSSSwholepatterntranscription',''.join(seg.top for seg in inventory_state('SourceSSSprimer N4UMI N6bc3 sp2').strands['SSS'].segs),'NNNNNNNNNNGGGATGCAGCTCGCTCCTG')
    c('Sourcehairpinwholepatterntranscription',''.join(seg.top.upper() for seg in inventory_state('Sourcebc2hairpin T7loop inventory').strands['hairpin'].segs),'CAAGACNNNNNNNCAGGAGCGAGCTGCATCCCAATTTAATACGACTCACTATAGGGATGCAGCTCGCTCCTGNNNNNNN')
    c('Selectedsortedwellleavesremainder',g.terminals['chosen_sorted_well_rest'].disposition,Disposition.UNRESOLVED)
    c('OnlyoneexternalRNAprimerstockaddition',sum(i.material=='rna_rt_primer' for st in g.steps.values() for i in st.inputs),1)
    c('RT uses produced external primer as structural component',[(i.material,i.role) for i in g.steps['rt'].inputs if i.role in {R.STRUCTURAL,R.SELF_PRIMER}],[('retained_rna_rt_primer',R.STRUCTURAL)])
    c('Retained primer produced by same-tube annealing','retained_rna_rt_primer' in g.steps['rt_anneal'].outputs,True)
    c('Private partial model cannot be publicly promoted',bool(g.publication_problems()),True)
    c.report()
if __name__=='__main__':selftest()
