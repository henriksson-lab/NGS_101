"""Private source-qualified RELACS barcoded nuclei and ChIP preparation."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'lib'))
from chemdraw import Scene,Segment,feature,revcomp
from illumina import NEBNEXT_HAIRPIN
from hairpin import Hairpin
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem
AUTHORITY='PublishedPMC6281648main145–201; publishedSI1actualFig1p1–2/Table2p22/Fig11p18. Closed2originals; detailed2018protocol/vendororiginalsunread.'
SOURCE_ORDERS=('/5Phos/ATAATATCAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTGATATTAT*T', '/5Phos/AAAGCGGCAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTGCCGCTTT*T', '/5Phos/GCGAATGCAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTGCATTCGC*T', '/5Phos/CTGGCGGCAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTGCCGCCAG*T', '/5Phos/CCAACTCGAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTCGAGTTGG*T', '/5Phos/GAAAATTTAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTAAATTTTC*T', '/5Phos/CAAACGCTAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTAGCGTTTG*T', '/5Phos/CTGACGCCAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTGGCGTCAG*T', '/5Phos/CCAGTATGAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTCATACTGG*T', '/5Phos/AAAATGCTAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTAGCATTTT*T', '/5Phos/GTGACTAGAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTCTAGTCAC*T', '/5Phos/ATGACAGAAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTTCTGTCAT*T', '/5Phos/GGAACAGCAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTGCTGTTCC*T', '/5Phos/CTGGCTCTAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTAGAGCCAG*T', '/5Phos/AGGGACTGAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTCAGTCCCT*T', '/5Phos/CAGGTAATAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTATTACCTG*T', '/5Phos/AAAACACTAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTAGTGTTTT*T', '/5Phos/AAGGGATGAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTCATCCCTT*T', '/5Phos/GCGAAGCAAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTTGCTTCGC*T', '/5Phos/GGGAGACGAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTCGTCTCCC*T', '/5Phos/AAAGCGGCNNAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNGCCGCTTT*T', '/5Phos/CCAACTCGNNAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNCGAGTTGG*T', '/5Phos/CAGGTAATNNAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNATTACCTG*T', '/5Phos/AAGGGATGNNAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNCATCCCTT*T', '/5Phos/GGGAGACGNNAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNCGTCTCCC*T', '/5Phos/AAAGCGGCNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNGCCGCTTT*T', '/5Phos/CCAACTCGNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNCGAGTTGG*T', '/5Phos/CAGGTAATNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNATTACCTG*T', '/5Phos/AAGGGATGNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNCATCCCTT*T', '/5Phos/GGGAGACGNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNCGTCTCCC*T', '/5Phos/AAAGCGGCNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNGCCGCTTT*T', '/5Phos/CCAACTCGNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNCGAGTTGG*T', '/5Phos/CAGGTAATNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNATTACCTG*T', '/5Phos/AAGGGATGNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNCATCCCTT*T', '/5Phos/GGGAGACGNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNCGTCTCCC*T')
BARCODES=tuple(x.removeprefix('/5Phos/')[:8]for x in SOURCE_ORDERS[:20])

def hairpin_scene(barcode=None):
    bc=barcode or BARCODES[0]
    if bc not in BARCODES:raise ValueError('Selectpublishedbaselinebarcode')
    left,right=('A'+NEBNEXT_HAIRPIN).split('U')
    body=bc+left+'U'+right+revcomp(bc)+'T'
    h=Hairpin('Source-local RELACS barcodehairpin',body,21,overhang_3_len=1)
    sc=Scene();sc.strand('hairpin',[Segment('nuclear_sample_barcode',bc,feature=feature('relacs.nuclear_sample','sample_index','whitelist',whitelist='PublishedSI1Table2Bar1–20',note='Initialcellpopulation, notsinglecellUMI')),Segment('canonicalhairpinleft',left),Segment('physicaldeoxyU','U',note='SourceinternalDNAuracil; USERopensloopafterChIP'),Segment('canonicalhairpinright',right),Segment('nuclear_barcode_complement',revcomp(bc),feature=feature('relacs.nuclear_sample','sample_index','whitelist',whitelist='PublishedSI1Table2Bar1–20')),Segment('terminal3T','T',note='SourcePSbondimmediatelybeforethisT; physicalnucleotide, noextraPSbase')],mod5='P')
    sc.ligation_end('hairpin','terminal3T',kind='cohesive',polarity='3-prime',text='Source3Tprotrusion/5Psinglehairpin compatiblewithsource3AdA-tailedchromatin; occupiedwholejunctionunassembled')
    return sc

def unknown_scene():
    sc=Scene();sc.strand('unassembled',[Segment('wholemoleculesunassembled','[unassembled]',placeholder=True)]);return sc

def protocol_graph():
    g=ProtocolGraph('RELACS selectedHepG2baselinebarcodednuclei/pooledH3K27acChIP',authority=AUTHORITY,representation='published_instance')
    x=g.input('fixed_cells','ActualFA-fixedHepG2sourcepellet; sourceapprox500000nuclei processinggroup, exactidentityunassigned',K.BIOLOGICAL,primary=True)
    def inp(k,label,kind=K.CONDITION,role=R.CONDITION,**kw):return g.input(k,label,kind,**kw),role
    def op(k,label,ports=(),kind=S.MAGIC,product='Actualsourceheterogeneousmixture; wholeoccupiedmolecules/products/ends/yieldsunassembled',**kw):
        nonlocal x
        y=g.product(k+'_product',product,K.MOLECULAR)
        if kind==S.MAGIC:kw.setdefault('reason','Sourcepositiveoperation; completechemicalproducts/enddistributions/recoveryunassembled, noabsenceclaim.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        g.react(k,label,kind,((x,R.TEMPLATE),)+tuple(ports),(y,),**kw);x=y
    def clean(k,label,ports=()):
        nonlocal x
        y=g.product(k+'_kept','Actualretainedphysicalfraction; recovery/carryoverunassigned',K.MOLECULAR);w=g.product(k+'_waste','Actualremovedworkupfractions',K.WASTE)
        g.react(k,label,S.CLEANUP,((x,R.TEMPLATE),inp(k+'_conditions',label))+tuple(ports),(y,w));g.discard(w,'Actualsourcephysicalworkup; noperfectpurityclaim');x=y
    def alloc(k,label,rest='Actualuntransferredstock; quantity/fateunassigned',qc=False):
        nonlocal x
        u=g.product(k+'_used',label,K.MOLECULAR);r=g.product(k+'_rest',rest,K.MOLECULAR)
        g.react(k,label,S.SPLIT,((x,R.TEMPLATE),),(u,r));g.mark_terminal(Terminal(u.id if qc else r.id,Disposition.QC if qc else Disposition.UNRESOLVED,label if qc else rest));x=r if qc else u;return u
    def components(k,label,spec,prepared=True):
        stock=g.input(k,label,K.MOLECULAR);out=[];ports=[]
        for suffix,kind,role,text in spec:
            m=g.product(k+'_'+suffix,text,kind);out.append(m);ports.append((m,role))
        g.react(k+'_components','CopresentcomponentsofONEactualpreparedstock, nofreshduplicateadditions'if prepared else'Inventorygroupofactualseparatestocks; nopremixpackagingclaim',S.GENERIC,((stock,R.TEMPLATE),),tuple(out),changes_molecule=False);return tuple(ports)
    op('nexson','Source1mllysisbuffer/PIC/CovarisE20075W2percent200cycles60–120sHepG2',ports=(inp('lysis_supplies','ActuallysisbufferIgepal/PICingredients; packagingunassigned',K.MOLECULAR,R.STRUCTURAL),inp('nexson_conditions','Actual10mMTris8/10mMNaCl/.2percentIgepal/programme')))
    initialqc=alloc('nuclei_estimate','ActualsmallnucleiQCaliquotforDNA/cellnumberestimate',qc=True)
    # The actual analytical aliquot has its own material use; it never enters production.
    initialdna=g.product('estimate_qc_dna','Actualsonicated/de-crosslinked/columnDNAforinitialestimateQC',K.MOLECULAR)
    g.react('estimate_qc_workup','Sourcealiquotshortsonication/de-crosslink/columnDNAestimate',S.MAGIC,((initialqc,R.TEMPLATE),inp('estimate_supplies','Actualanalyticalworkupinventoryunassembled',K.MOLECULAR,R.STRUCTURAL)),(initialdna,),reason='Positiveanalyticalbranch, exactvolumes/reagentstockcompositionunassembled.');g.terminals.pop(initialqc.id);g.mark_terminal(Terminal(initialdna.id,Disposition.QC,'AnalyticalDNA/cellnumberestimate'))
    clean('nuclei_pellet','Sourceapprox500000nuclei1000g5min')
    op('swell','Source50ul.5percentSDSRT10min/25ul10percentTriton+145water+25CutSmart+2.5PIC',(inp('sds_triton','ActualSDS/Tritonchemicalsupplies',K.MOLECULAR,R.STRUCTURAL),inp('cut_buffer','Actual25ulCutSmartand2.5ul100xPIC/water; sourcevolumesunmodified')))
    op('digest','SourceCviKI1 5Uper100000nuclei/20C16h800rpm',(inp('cviki','ActualCviKI1catalyst; sourcealternative20Uwide-inputrangekeptseparate',K.CATALYST,R.CATALYST),inp('digest_conditions','Actualenzymeprogramme; RG|CYbluntcutsourceSI1')))
    clean('nuclei_wash','Source1000g5min/200ulnucleiwashTrisTritonBSA',(inp('wash_bsa','Actual.1mgmlBSAmolecularcarrierinwash',K.MOLECULAR,R.STRUCTURAL),))
    digestqc=alloc('digest_qc','Actual10–20ulnucleidigestionQCaliquot',qc=True)
    qc_mix=g.input('digest_qc_supplies','Actualanalyticalbuffer/4ul5MNaCl/2ul10mgmlRNaseA/2ul20mgmlPK inventory; packagingunassigned',K.MOLECULAR)
    qd=g.product('digest_qc_dna','Actualde-crosslinkedpurifiedDNAforQubitandcapillarydigestionQC',K.MOLECULAR)
    g.react('digest_qc_workup','SourceQCto100ulChIPelutionbuffer/37C30min65C>=2h/column',S.MAGIC,((digestqc,R.TEMPLATE),(qc_mix,R.STRUCTURAL)),(qd,),reason='ActualseparateQCbranch; fullanalyticalmolecularstates/recoveryunassembled.');g.terminals.pop(digestqc.id);g.mark_terminal(Terminal(qd.id,Disposition.QC,'ActualQubitDNAsizeQC'))
    clean('barcode_pellet','Sourceremainingnucleipellet/resuspend25ulper100–100000nuclei')
    alloc('barcode_well','Actual25ulnucleialiquotforonebarcodingwell')
    prep=components('endprep_mix','Actual1.5ulNEBNextUltraIIendprepenzymemix',[
        ('enzyme',K.CATALYST,R.CATALYST,'Actualend-repair/A-tailing enzymeclasses, identities/dosesunassembled'),('other',K.MOLECULAR,R.STRUCTURAL,'Otheractualendprepmixcomponentsunassembled')])
    op('endprep','Source25ulnuclei+1.5enzyme+3.5buffer20C30min65C5min',prep+(inp('endprep_buffer','Actual3.5ulendprepreactionbuffer; nucleotideco-components/stockformulationunassembled',K.MOLECULAR,R.STRUCTURAL),),product='Sourceend-repaired3A-tailedchromatinclassesinsidenuclei plusunprocessedclasses; wholemoleculesunassembled')
    hairpin=g.input('hairpin_stock','Actualsourcehairpinbarcodeoligo; Bar1localreference, experimentalsampleassignmentunassigned',K.OLIGO,state=hairpin_scene(),chemical_state={'raw_order':SOURCE_ORDERS[0],'five_prime':'P','internal_dU':8+len(('A'+NEBNEXT_HAIRPIN).split('U')[0]),'phosphorothioate':'Bondimmediatelybeforeterminal3T; noextraresidue','source_concentrations':'Mainliteral100mMstock/15mManneal retained, interpretationunassigned','baseline':'20nuclearbarcodes; 15UMIvariantsnotselected'})
    annealed=g.product('hairpin_annealed','Actualsourcefoldedbarcodehairpinstock; partialfold/occupancyunassigned',K.ADAPTER,state=hairpin_scene())
    g.react('hairpin_anneal','Source100mMstock→15mManneal/95C2min/slowcoolRT',S.MAGIC,((hairpin,R.TEMPLATE),inp('hairpin_buffer','Actual10mMTris8/50mMNaCl/1mMEDTA')),(annealed,),reason='Sourcehairpinannealingpositive; exactwholefold/stockamounts unassembled, sourceunitsnotnormalized.')
    hu=g.product('hairpin_used','Actual1.2ulhairpinbarcodeallocation',K.ADAPTER);hr=g.product('hairpin_rest','Actualunusedannealedhairpinstock',K.ADAPTER);g.react('hairpin_allocation','Allocateactual1.2ulhairpinstock',S.SPLIT,((annealed,R.TEMPLATE),),(hu,hr));g.mark_terminal(Terminal(hr.id,Disposition.UNRESOLVED,'Unusedadapterstockamountunassigned'))
    lig=components('ligation_master','Actual15ulNEBNextUltraIIligationmastermix',[
        ('ligase',K.CATALYST,R.CATALYST,'Actualligasecomponentclassofsame15ulstock'),('other',K.MOLECULAR,R.STRUCTURAL,'Othercopresentligationmasterco-componentsunassembled')])
    op('barcode_ligate','Source1.2ulhairpin+15ulligmaster+.5enhancer30C15min20C15min',((hu,R.STRUCTURAL),)+lig+(inp('lig_enhancer','Actual.5ulligationenhancer',K.MOLECULAR,R.STRUCTURAL),),geometry=AdapterLigationGeometry.known('cohesive',hairpin_scene(),statement='SourceSI1known5P/hairpin3T/PS-protectedterminalT compatiblewithsource3Arecipientclass; fulloccupiedchromatinjunctionsunassembled.'))
    op('inhibit_ligase','Source4.7ul3MNaClperwell',ports=(inp('lig_stop','Actual4.7ul3MNaCl'),),kind=S.GENERIC)
    other=g.input('other_barcoded_nuclei','Otheractualindependentlybarcoded/inhibitednucleicontributions; samplecount/barcodes/quantitiesunassigned; upstreampreparationoutside singledetailedbranch',K.MOLECULAR)
    pool=g.product('pooled_nuclei','Actualpooledbarcodednuclei',K.MOLECULAR);g.react('pool_nuclei','Physicallypoolactualindependentlybarcodednucleicontributions',S.POOL,((x,R.TEMPLATE),(other,R.TEMPLATE)),(pool,));x=pool
    pellet1=g.product('pellet_first','Actual5000g10min20Cnucleipellet',K.MOLECULAR);super1=g.product('super_first','ActualfirstsupernatantRETAINEDforrescue',K.MOLECULAR)
    g.react('first_spin','Source5000g10min20C; retainpelletANDsupernatantforrescue',S.SPLIT,((x,R.TEMPLATE),),(pellet1,super1))
    pellet2=g.product('pellet_rescue','Actual11000g5min20CrescuepelletfromFIRSTsupernatant',K.MOLECULAR);super2=g.product('super_final','Actualpostrescuefinalsupernatant; fateunassigned',K.MOLECULAR)
    g.react('rescue_spin','Re-spinactualfirstsupernatantat11000g5min20C',S.SPLIT,((super1,R.TEMPLATE),),(pellet2,super2));g.mark_terminal(Terminal(super2.id,Disposition.UNRESOLVED,'Postrescuesupernatantnotautomaticallyreturnedtochromatin'))
    pellets=g.product('combined_pellets','Actualcombinedfirst+rescuednucleipellets',K.MOLECULAR);g.react('combine_pellets','PoolBOTHactualpellets',S.POOL,((pellet1,R.TEMPLATE),(pellet2,R.TEMPLATE)),(pellets,));x=pellets
    op('release_chromatin','Sourceapprox130ulshearingbufferper<=500000nuclei/Covaris105W2percent200cycles5min',(inp('shear_buffer','Actual10mMTris8/.1percentSDS/1mMEDTA/sourceprogramme'),))
    op('chip_dilute','Source1:2dilutioniniC1plusPIC',ports=(inp('ic1_stock','ActualiC1/PICstocks; kitformulationunassembled',K.MOLECULAR,R.STRUCTURAL),),kind=S.GENERIC)
    alloc('chip_aliquot','Actual200ulChIPaliquot','Actualotherbarcodedchromatinstock/otherChIPcontributions; amounts/fatesunassigned')
    op('antibody','SelectedsourceH3K27ac1–2ugaccordinginputcount/10h4C',(inp('h3k27ac','ActualsourceH3K27acantibody; specificinputcountbranchunassigned',K.MOLECULAR,R.STRUCTURAL),),kind=S.BINDING)
    bound=g.product('chip_bound','ActualproteinA-beadboundimmunocomplexes',K.MOLECULAR);unbound=g.product('chip_unbound','ActualunboundChIPfraction',K.WASTE)
    g.react('chip_capture','SourceproteinAconjugatedmagneticbeads3h4C',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),inp('protein_a_beads','ActualproteinAmagneticbeads',K.SELECTION,R.SELECTION)),(bound,unbound));g.discard(unbound,'ActualexcludedChIPunboundfraction');x=bound
    clean('chip_wash','Sourcefour5minwasheswithiDealbuffers1–4')
    eluate=g.product('chip_eluate','ActualmanuallyrecoveredChIPeluate',K.MOLECULAR);resid=g.product('chip_residual','Actualresidualbeadmaterial; fateunassigned',K.MOLECULAR)
    g.react('chip_elute','Sourcemanualeluatecollection; kitelutiondetailsunassembled',S.MAGIC,((x,R.TEMPLATE),inp('chip_elution_supplies','Actualsourceelutioninventoryunassembled',K.MOLECULAR,R.STRUCTURAL)),(eluate,resid),reason='Actualmanualsourceelution; fullstockcomposition/end/productsunassembled.');g.mark_terminal(Terminal(resid.id,Disposition.UNRESOLVED,'ResidualbeadmaterialnotDNAproduction'));x=eluate
    op('de_crosslink','SourceRNaseAtreatment/de-crosslink/proteolysis37C30min65C4h',(inp('rnasea','ActualRNaseA',K.CATALYST,R.CATALYST),inp('proteolysis','Actualsourceproteolyticreagentclass; dose/packagingunassembled',K.CATALYST,R.CATALYST),inp('de_crosslink_conditions','ActualChIPsourceprogramme')))
    clean('chip_dna','SourceQiagenMinElutefinal20ul')
    alloc('pcr_input','Actual18ulof20ulChIPDNA','Actualnominal2ulof20ulChIPDNAuntransferred')
    pc=components('pcr_master','ONEactual25ulNEBNextUltraIIPCRmaster',[
        ('polymerase',K.CATALYST,R.CATALYST,'ActualPCRpolymeraseclasscopresentinsame25ulstock'),('nts',K.MOLECULAR,R.STRUCTURAL,'ActualPCRnucleotideclasscopresentinsame25ulstock'),('conditions',K.CONDITION,R.CONDITION,'ActualPCRbufferco-components')])
    user=g.input('user','Actual3ulUSERenzyme',K.CATALYST)
    primers=components('primer_supplies','Actual2ul10uMuniversal+2ul10uMindexedprimerstocks, aggregateinventorynotpremix',[
        ('universal',K.OLIGO,R.STRUCTURAL,'ActualuniversalPCRprimer; wholeordernotassembledinselectedpass'),('indexed',K.OLIGO,R.STRUCTURAL,'ActualindexedChIPPCRprimer; wholebody/indexidentityunassigned')],prepared=False)
    opened=g.product('user_open_product','SourceopenedhairpinDNA pluscopresentUSER/exciseduracil/unreactedclasses; wholeends/productsunassembled',K.MOLECULAR)
    carried=[];carryports=[]
    for m,role in pc+primers:
        c=g.product('carried_'+m.id,'SAMEphysicalcomponentcarriedthrough37CUSERincubation; nofreshaddition: '+m.label,m.kind)
        carried.append(c);carryports.append((c,role))
    g.react('user_open','Source18DNA+25PCRmaster+3USER+2universal+2indexed/37C15min; retainPCRcomponents',S.MAGIC,((x,R.TEMPLATE),)+pc+((user,R.CATALYST),)+primers,(opened,)+tuple(carried),reason='SourceUSERdUopeningpositive; fullopenedDNA/endchemistry/coproduceduracilunassembled. PCRcomponentportsinventorysamephysicalmixture, notseparationorfreshadditions.');x=opened
    op('pcr','Source98C30s/12cycles98C10s65C75s/65C5min; primerenzymeNTcarriedfromsameUSERmixture',ports=tuple(carryports),kind=S.PCR,geometry=PCRGeometry.unavailable(unknown_scene(),reason='Actualprimerclassesalreadycopresentupstream; wholeoccupiedtemplates/fullPCRproducts/primerlandingunassembled.'))
    clean('ampure_first','Source.9xAMPure')
    clean('ampure_second','Source1xAMPure')
    alloc('qubit_qc','ActualQubitDNAHSaliquot; volumeunassigned',qc=True)
    alloc('fragment_qc','ActualFragmentAnalyzeraliquot; volumeunassigned',qc=True)
    alloc('pool_allocation','Actualmolar-normalizedlibrarycontribution; amountunassigned')
    others=g.input('other_finished_libraries','OtheractualQCpassednormalizedlibrarycontributions; exactcount/identities/amountsunassigned',K.MOLECULAR)
    phix=g.input('phix','Actual10percentPhiXspike-inlibrary; wholemoleculesunassembled',K.LIBRARY)
    final=g.product('library_pool','ActualnormalizedRELACSlibrarypoolwith10percentPhiX',K.MOLECULAR);g.react('final_pool','Poolactualnormalizedcontributionsand10percentPhiX',S.POOL,((x,R.TEMPLATE),(others,R.TEMPLATE),(phix,R.TEMPLATE)),(final,));x=final
    op('denature','Sourcedenaturation/instrumentloadingprep; exactIlluminarecipeunassembled',ports=(inp('denaturation_conditions','Actualsourceinstrumentloadingprogrammeunassembled'),),kind=S.GENERIC)
    alloc('loading','ActualHiSeqloadingallocation; volumeunassigned')
    lib=g.product('library','ActualRELACSH3K27acHiSeqDNAsequencinglibrary; wholeproductsunassembled',K.LIBRARY);g.react('load','SourceHiSeq3000PE75',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(read_lengths={'Read1':75,'Read2':75},limitation='SourceHiSeq3000PE75/10percentPhiX. SI1bothreadsbeginnuclearbarcode/ChIPindexseparate; fullprimers/readlanding/indexcyclesunassembled. No baselineUMIvariantinheritance.'))
    g.add_source_item(SourceItem('source_literals','Main145prints100mMstock/15mManneal; Table2physical5P/dU/PSbondbeforeterminalT.','RetainliteralunitswithoutuMrepair/publicationerrorclaim.20baselinebarcodes+15UMIvariants; representativeBar1doesnotassignactualexperimentalsample.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();p=sources_root()/HERE.parent.name/'s42003-018-0219-z_supp_42003_2018_219_MOESM1_ESM.pdf.txt'
    if have(c,Source(p,'https://doi.org/10.1038/s42003-018-0219-z')):
        t=''.join(p.read_text().split())
        for n,raw in enumerate(SOURCE_ORDERS,1):c('WholemodifiedTable2order '+str(n),raw in t,True)
        c('SourcePSprotectsterminalT','phosphorothioatedbond'in t,True)
    main=p.parent/'s42003-018-0219-z_PMC6281648.xml.txt'
    if have(c,Source(main,'https://doi.org/10.1038/s42003-018-0219-z')):
        t=''.join(main.read_text().split());c('Source100mMstockliteral','final100mMconcentration'in t,True);c('Source15mMannealliteral','dilutedto15mM'in t,True)
    actual=''.join(s.top for st in hairpin_scene().strands.values()for s in st.segs)
    expected=SOURCE_ORDERS[0].removeprefix('/5Phos/').replace('/ideoxyU/','U').replace('*','')
    c('Wholebaselinepolymermatchesmodifiedorder',actual,expected)
    c('Nuclearbarcodemeanspopulation',hairpin_scene().strands['hairpin'].segs[0].feature.role,'sample_index')
    c('Only20baselinechoices',len(BARCODES),20)
    try:hairpin_scene('AAAAAAAA');bad=False
    except ValueError:bad=True
    c('Rejectunlistedbarcode',bad,True)
    c('DigestionQCusesactualaliquot',g.steps['digest_qc_workup'].inputs[0].material,'digest_qc_used')
    c('DigestionQCterminalisanalyticalDNA',g.terminals['digest_qc_dna'].disposition,Disposition.QC)
    c('Rescuespinsactualfirstsupernatant',g.steps['rescue_spin'].inputs[0].material,'super_first')
    c('Combinebothphysicalpellets',set(p.material for p in g.steps['combine_pellets'].inputs),{'pellet_first','pellet_rescue'})
    c('Finalsupernatantnotreturned',g.terminals['super_final'].disposition,Disposition.UNRESOLVED)
    c('ChIPuses200ulallocation',g.steps['antibody'].inputs[0].material,'chip_aliquot_used')
    c('PCRpreparationuses18of20',g.steps['user_open'].inputs[0].material,'pcr_input_used')
    c('2ulDNAstockremainderretained',g.terminals['pcr_input_rest'].disposition,Disposition.UNRESOLVED)
    c('PCRNTscarriedthroughUSERnotnewdose',g.producer_of('carried_pcr_master_nts').id,'user_open')
    c('PCRprimerscarriedthroughUSER',g.producer_of('carried_primer_supplies_universal').id,'user_open')
    c('PCRusesactualopenedDNA',g.steps['pcr'].inputs[0].material,'user_open_product')
    c('QubitQCnotproduction',g.terminals['qubit_qc_used'].disposition,Disposition.QC)
    c('PhiXisactualfinalpoolinput','phix'in{p.material for p in g.steps['final_pool'].inputs},True)
    c('WholePCRproductsunassembled',g.materials['pcr_product'].state,None)
    c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
