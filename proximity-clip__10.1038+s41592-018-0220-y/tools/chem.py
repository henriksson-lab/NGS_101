"""Protocol-owned source-qualified partial graph; no finished-publication marker."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem,SequencingPlan
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry

def clean(g,key,x,procedure):
    reagent=g.input(key+'_condition',procedure,K.CONDITION)
    out=g.product(key+'_kept','Recovered source preparation; exact molecular recovery/carryover unquantified',K.MOLECULAR)
    waste=g.product(key+'_separated','Actual separated source workup fractions; no ideal molecular partition assumed',K.WASTE)
    ports=[(x,R.TEMPLATE),(reagent,R.CONDITION)]
    if 'GlycoBlue' in procedure:
        carrier=g.input(key+'_carrier','Actual source GlycoBlue molecular carrier; carryover unassigned',K.MOLECULAR);ports.append((carrier,R.STRUCTURAL))
    g.react(key,procedure,S.CLEANUP,ports,(out,waste));g.discard(waste,'Actual source physical workup only.');return out

def op(g,key,x,label,kind,*,reagent=None,catalyst=False,product=None,reason='',changes_molecule=True,molecular=False):
    ports=[(x,R.TEMPLATE)]
    if reagent:
        stock=g.input(key+'_reagent',reagent,K.CATALYST if catalyst else K.MOLECULAR if molecular else K.CONDITION);ports.append((stock,R.CATALYST if catalyst else R.STRUCTURAL if molecular else R.CONDITION))
    if molecular and kind in {S.CHEMICAL_CONVERSION,S.CROSS_LINKING}:
        cond=g.input(key+'_conditions',label+'; sourceconditions, notanotherphysicaldose',K.CONDITION);ports.append((cond,R.CONDITION))
    out=g.product(key+'_product',product or label+'; reaction remnants/carryover remain unassigned',K.MOLECULAR)
    g.react(key,label,kind,ports,(out,),reason=reason,changes_molecule=changes_molecule);return out

def unknown_adapter(label):
    sc=Scene();sc.strand('adapter',[Segment(label,'[source end geometry unreviewed]',placeholder=True)])
    return AdapterLigationGeometry('unavailable','unavailable',(sc,),'Source operation known; exact original oligo/end/junction assembly unreviewed, no absence claim.')
SOURCE='Benhalevy et al.2018 DOI10.1038/s41592-018-0220-y; active published protocolMOESM3 steps33–109; associatedAMPMC6289640 identity/overview and HiSeq2500/3000 report remain version-qualified.'
@dataclass(frozen=True)
class SelectedCrosslinkedRNA:
    sequence:str
    four_su:frozenset[int]
    protein_adducts:frozenset[int]
    def __post_init__(self):
        if not self.sequence or set(self.sequence)-set('ACGU') or any(type(i)is not int or not 0<=i<len(self.sequence) or self.sequence[i]!='U' for i in self.four_su) or not self.protein_adducts<=self.four_su:raise ValueError('selected photoreactive class needs actual U-origin4SU and explicitly selected crosslink sites')
    def first_dna(self,diagnostic_sites=frozenset()):
        if not diagnostic_sites<=self.protein_adducts:raise ValueError('selected diagnostic copying sites must be supplied crosslink positions; no universal mismatch prediction')
        bases=[{'A':'T','C':'G','G':'C','U':'A'}[b] if i not in diagnostic_sites else 'G' for i,b in enumerate(self.sequence)]
        return ''.join(reversed(bases))

def protocol_graph():
    g=ProtocolGraph('Proximity-CLIP — published intact-RNA and footprint routes',authority=SOURCE,representation='published_instance')
    x=g.input('cells','Source targeted APEX2-expressing cells with4SU metabolic labeling; exact recording sites unassigned',K.MOLECULAR,primary=True)
    x=op(g,'bp',x,'Source BP500µM30min followed by H2O2final1mM1min',S.CHEMICAL_CONVERSION,reagent='Actual BP/H2O2 molecularcontributors; APEX2 alreadyexpressedincells',molecular=True,product='Local protein-biotin labeling plus RNA/cell material; no automatic covalent RNA-biotin assignment')
    x=op(g,'uv_crosslink',x,'Source three quencher washes then UV0.4J/cm2 at312–365nm during quenching',S.CROSS_LINKING,reagent='Actual Trolox/ascorbate/azide molecularquenchers; UV iscondition',molecular=True,product='Selected4SU–protein crosslinked RNPs with protein-biotin and other/background material')
    x=clean(g,'harvest',x,'Source cell harvest/pellet; microscopy fixation excluded')
    x=op(g,'lysis',x,'Source preparative RIPA lysis',S.EXTRACTION,reagent='source RIPA/protease inhibitors/quenchers')
    x=clean(g,'clarify',x,'Source15,000g10min4°C clarification')
    bead=g.input('beads','Source60µl streptavidin beads,1hRT capture, specified RIPA/KCl/carbonate/urea/NP40 washes',K.SELECTION)
    bound=g.product('bound','Captured protein-biotin RNP preparation plus possible background',K.MOLECULAR);ft=g.product('saved_ft','Source saved capture flowthrough for analysis; exact downstream disposition unassembled',K.MOLECULAR)
    g.react('capture','Source preparative streptavidin capture',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(bead,R.SELECTION)),(bound,ft));g.mark_terminal(Terminal(ft.id,Disposition.UNRESOLVED,'Source saves flowthrough; not silently discarded.'))
    bound=clean(g,'bound_washes',bound,'Source preparative bead wash sequence')
    intact=g.product('intact','Source intact-RNA aliquot; no RNaseT1 treatment',K.MOLECULAR);foot=g.product('foot','Source footprint aliquot',K.MOLECULAR);ms=g.product('ms','Separate source mass-spectrometry aliquot, not a proteinNGS library',K.MOLECULAR)
    g.react('aliquot','Source final-wash bead split: intactRNAseq/footprints/MS',S.SPLIT,((bound,R.TEMPLATE),),(intact,foot,ms));g.mark_terminal(Terminal(ms.id,Disposition.UNRESOLVED,'Proteomics route outside RNA NGS model scope.'))
    foot=op(g,'rnase_t1',foot,'RNaseT1 1U/µl15min22°C; retain reaction products before washing',S.DIGESTION,reagent='source RNaseT1',catalyst=True)
    foot=clean(g,'t1_wash',foot,'Source RNaseT1 and dephosphorylation-buffer washes')
    foot=op(g,'cip',foot,'Source CIP10min37°C',S.DEPHOSPHORYLATION,reagent='source CIP/end-processing reagent; full RNA terminal products unassembled',catalyst=True)
    foot=clean(g,'cip_wash',foot,'Source CIP washes and PNK-buffer exchange')
    pnk=g.input('pnk_enzyme','Actual source T4PNK forhotlabeling; nosecondfreshaddition',K.CATALYST)
    hot=g.input('hot_atp','Actual source gamma32P-ATP molecularphosphatedonor; tracer notsequenceidentifier',K.MOLECULAR)
    hot_rna=g.product('hot_rna','HotlabelingRNA mixture; copresentenzymecomponent representedseparately',K.MOLECULAR);retained_pnk=g.product('retained_pnk','Same-tube retainedPNK catalyst, notisolatedphysicalfraction',K.CATALYST)
    g.react('pnk_hot','Source hotPNKlabeling30min37C',S.PHOSPHORYLATION,((foot,R.TEMPLATE),(pnk,R.CATALYST),(hot,R.STRUCTURAL)),(hot_rna,retained_pnk))
    cold=g.input('cold_atp','Actual freshcoldATP to100uMfinal; source5min37C',K.MOLECULAR)
    out=g.product('pnk_product','SourcephosphorylatedRNA plusretainedreagents; isotope/terminalconversionfractionsunassigned',K.MOLECULAR)
    g.react('pnk','SourcecoldATPaddition withretainedPNK5min37C',S.PHOSPHORYLATION,((hot_rna,R.TEMPLATE),(retained_pnk,R.CATALYST),(cold,R.STRUCTURAL)),(out,));foot=out
    wash_condition=g.input('pnk_wash_condition','Published steps46–47: magnetize beads and wash five times with1ml PNK buffer without DTT',K.CONDITION)
    washed=g.product('pnk_wash_kept','Recovered source bead preparation after PNK washes; molecular carryover unquantified',K.MOLECULAR)
    marker=g.product('pnk_wash_qc','Source saved50µl radioactive first-wash aliquot for marking subsequent Urea-PAGEs; composition unassigned',K.QC)
    wash_rest=g.product('pnk_wash_separated','Remaining radioactive wash fractions after saved QC aliquot; exact material partition unquantified',K.WASTE)
    g.react('pnk_wash','Published steps46–47: preserve50µl first-wash marker separately from remaining washes',S.CLEANUP,((foot,R.TEMPLATE),(wash_condition,R.CONDITION)),(washed,marker,wash_rest))
    g.mark_terminal(Terminal(marker.id,Disposition.QC,'Published MOESM3 step46/source533–534: saved radioactive wash marks later Urea-PAGEs; marker composition and later fate unassembled.'))
    g.discard(wash_rest,'Remaining source wash material only; saved50µl QC marker excluded.')
    foot=washed
    released={}
    for key,current in (('intact',intact),('foot',foot)):
        current=op(g,key+'_pk',current,'Source three sequential ProteinaseK50°C30min incubations',S.DEPROTEINIZATION,reagent='source sequential ProteinaseK additions',catalyst=True,product='RNA-containing preparation plus peptide/protein/bead remnants; no automatic crosslink reversal or label transfer')
        current=clean(g,key+'_magnet',current,'Source magnetize beads, transfer RNA-containing supernatant and discard beads')
        current=clean(g,key+'_organic',current,'Source acidic phenol/chloroform then chloroform aqueous recovery')
        current=clean(g,key+'_ethanol',current,'Source GlycoBlue/ethanol RNA recovery')
        released[key]=current
    kit=g.input('nebnext','Source NEBNextUltraRNA kit, fragmentation10min; no rRNA depletion assumed for captured intactRNA',K.MOLECULAR)
    lib=g.product('intact_library','Source intact-RNAseq cDNA library; full kit/end/read geometry unassembled',K.LIBRARY)
    g.react('intact_library_prep','Source historical NEBNext library handoff',S.MAGIC,((released['intact'],R.TEMPLATE),(kit,R.STRUCTURAL)),(lib,),reason='Internal historical kit operations/oligos not borrowed; rRNA depletion applies only to separate total-cell controls.')
    g.sequence(lib,SequencingPlan(acquisition='unspecified',limitation='Source reports RNA sequencing; associatedAMHiSeq2500/3000, exact intact-library acquired roles/cycles/primer geometry unassembled.'))
    x=clean(g,'foot_size',released['foot'],'Source denaturing gel footprint-size selection and salt/filter/ethanol extraction;20–40 or41–~80nt branches')
    for key,label,enzyme in (('ligate3','Source adenylated3adapter,90°C1min denature then Rnl2(1–249)K227Q overnight on ice','source truncated Rnl2'),('ligate5','Source5adapter/ATP,90°C1min then Rnl1 1h37°C','source Rnl1')):
        adapter=g.input(key+'_adapter','Source adapter original body/ends/modifications unreviewed; no generic blocker substituted',K.ADAPTER)
        mix=g.input(key+'_mix','Source3ligation DMSO/bufferwithoutATP and32P19/35marker aliquot' if key=='ligate3' else 'Source5ligation DMSO/bufferwithATP; ATPpartofactualbuffersupply',K.MOLECULAR)
        heat=g.input(key+'_heat_condition','SourceRNA+adapter+mix90C1minthentoice BEFOREaddingligase',K.CONDITION)
        warmed=g.product(key+'_heated_rna','SameRNA/mix material afterheating; marker/reagentcarryoverunassigned',K.MOLECULAR);warm_adapter=g.product(key+'_heated_adapter','Samecopresentadapter afterheat, notfreshaddition',K.ADAPTER)
        g.react(key+'_heat','SourcecopresentRNA/adapterdenaturation',S.THERMAL_TREATMENT,((x,R.TEMPLATE),(adapter,R.STRUCTURAL),(mix,R.STRUCTURAL),(heat,R.CONDITION)),(warmed,warm_adapter))
        cat=g.input(key+'_enzyme',enzyme,K.CATALYST);out=g.product(key+'_product','Source adapter-ligated footprint population plus unligated/carryover molecules',K.MOLECULAR)
        g.react(key,label,S.LIGATION,((warmed,R.TEMPLATE),(warm_adapter,R.STRUCTURAL),(cat,R.CATALYST)),(out,),geometry=unknown_adapter(label));x=clean(g,key+'_gel',out,'Source ligation gel excision and salt/filter/ethanol RNA extraction')
    primer=g.input('rt_primer','Source reverse PCR primer; original barcoded member/sequence unreviewed',K.OLIGO);rt=g.input('ssiii','Source SuperScriptIII1h50°C',K.CATALYST)
    cdna=g.product('cdna','Source first-strand products with original RNA/peptide remnants retained; selected T-to-C copying signature is not RNA mutation',K.MOLECULAR)
    nts=g.input('rt_nt','Actual source RT dNTP substrates; supply packaging unassigned',K.MOLECULAR)
    g.react('rt','Source SSIII reverse transcription',S.REVERSE_TRANSCRIPTION,((x,R.TEMPLATE),(primer,R.STRUCTURAL),(rt,R.CATALYST),(nts,R.STRUCTURAL)),(cdna,))
    diagnostic=g.product('diagnostic_cdna','Separate source diagnostic cDNA aliquot',K.MOLECULAR);prep=g.product('preparative_cdna','Separate source preparative cDNA aliquot; not amplified diagnostic material',K.MOLECULAR)
    g.react('pcr_allocation','Source diagnostic versus preparativePCR aliquots',S.SPLIT,((cdna,R.TEMPLATE),),(diagnostic,prep))
    def pcr_array(key,template,count,total,cycles):
        supplies=g.input(key+'_supplies','AggregateactualPCRsupplies includingtwo primerclasses/PlatinumTaq/dNTP/buffer; stockcountandpackagingunassigned, notassumedpremix',K.MOLECULAR)
        prepared=g.product(key+'_prepared',str(total)+'ulactualsourcecompletePCRpreparation includingtemplate andaddedprimers; wholeoligosunassembled',K.MOLECULAR)
        g.react(key+'_prepare','Prepare sourcePCRarraybeforephysicalsplit',S.GENERIC,((template,R.TEMPLATE),(supplies,R.STRUCTURAL)),(prepared,),changes_molecule=False)
        tubes=[g.product(key+'_tube_'+str(i),str(total//count)+'ulsourcePCRreaction '+str(i+1),K.MOLECULAR) for i in range(count)]
        g.react(key+'_split','Sourceactual'+str(count)+'parallelPCRtubes',S.SPLIT,((prepared,R.TEMPLATE),),tubes)
        outputs=[]
        for i,tube in enumerate(tubes):
            f=g.product(key+'_f_'+str(i),'CopresentforwardPCRprimerclass; wholeorderingunreviewed',K.OLIGO)
            r=g.product(key+'_r_'+str(i),'CopresentreversePCRprimerclass; wholeorderingunreviewed',K.OLIGO)
            enzyme=g.product(key+'_taq_'+str(i),'CopresentPlatinumTaq fromsameactualpreparedreaction',K.CATALYST)
            nts=g.product(key+'_nt_'+str(i),'CopresentdNTPs fromsameactualpreparedreaction',K.MOLECULAR)
            target=g.product(key+'_template_'+str(i),'Templateandotherreactioncomponents ofsameactualtube; notphysicalfractionation',K.MOLECULAR)
            g.react(key+'_components_'+str(i),'Reactioncomponent bookkeeping, nofreshreagentaddition',S.GENERIC,((tube,R.TEMPLATE),),(target,f,r,enzyme,nts),changes_molecule=False)
            product=g.product(key+'_amp_'+str(i),'SourcePCRDNAproducts plusoriginaltemplates/reagents',K.MOLECULAR)
            g.react(key+'_pcr_'+str(i),'SourcePlatinumTaqPCR '+str(cycles[i])+'cycles' if key=='diagnostic' else 'SourcePlatinumTaqPCR atdiagnosticallyselectedcyclecount',S.PCR,((target,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(enzyme,R.CATALYST),(nts,R.STRUCTURAL)),(product,),geometry=PCRGeometry.unavailable(Scene(),reason='ActualparallelPCR/primerclassesknown; originaloligos andwholeproductgeometryunassembled.'))
            outputs.append(product)
        return outputs
    diag_products=pcr_array('diagnostic',diagnostic,6,60,(9,11,13,15,17,19))
    for i,product in enumerate(diag_products):
        qa=g.product('diagnostic_qa_'+str(i),'Entire10uldiagnosticPCR usedfor2.5percentagarosecyclechoice, notreusedaspreparativeDNA',K.QC)
        g.react('diagnostic_gel_'+str(i),'Sourcewhole-tubeanalyticalgelinspection',S.MAGIC,((product,R.TEMPLATE),),(qa,),reason='Actualsourceanalyticalgel; assayreagents/wholeDNAgeometryunassembled.');g.mark_terminal(Terminal(qa.id,Disposition.QC,'Spentdiagnosticreaction'))
    prep_products=pcr_array('preparative',prep,3,300,(None,)*3)
    zymo=g.input('zymo_supplies','ActualsourceZymoPCRpurification/concentration supplies',K.SELECTION)
    recovered=g.product('zymo_kept','SourcePCR-arraylibraryrecoveredto70ul; orderofcombination/cleanupunassembled',K.MOLECULAR)
    rejects=g.product('zymo_rejects','ActualsourceZymoworkuprejects',K.WASTE)
    g.react('zymo','SourceparallelPCRproductworkup/concentration, notanassertedone-tubeprecleanuppool',S.MAGIC,tuple((p,R.TEMPLATE) for p in prep_products)+((zymo,R.SELECTION),),(recovered,rejects),reason='Source108concentratesPCRproductto70ul; exactpool/cleanup scheduling ofthreePCRtubesunassembled.');g.discard(rejects,'Actualsourcecolumnworkupfractions')
    take=g.product('pippin_input','Source30ul of70ulrecoveredlibrary forPippin',K.MOLECULAR);remain=g.product('library_remainder','Untransferredlibrary after30of70ulselection; nominal40ul, recoveryunassigned',K.MOLECULAR)
    g.react('pippin_allocate','AllocateactualPippinlibrarycontribution',S.SPLIT,((recovered,R.TEMPLATE),),(take,remain));g.mark_terminal(Terminal(remain.id,Disposition.UNRESOLVED,'Untransferredsourcepreparation'))
    x=clean(g,'pippin',take,'SourcePippinsize-selection windows; fullproduct/adapterlengthassemblyunreviewed')
    g.sequence(x,SequencingPlan(acquisition='unspecified',limitation='Source footprintNGS; associatedAMHiSeq2500/3000, exact run/primer/index layout unassembled.'))
    g.declare_reviewed_scope(('cells','bp_reagent','bp','uv_crosslink_reagent','uv_crosslink','harvest_condition','harvest','lysis_reagent','lysis','clarify_condition','clarify','beads','capture','bound_washes_condition','bound_washes','aliquot','rnase_t1_reagent','rnase_t1','t1_wash_condition','t1_wash','cip_reagent','cip','cip_wash_condition','cip_wash','pnk_enzyme','hot_atp','pnk_hot','cold_atp','pnk','pnk_wash_condition','pnk_wash','ligate3_adapter','ligate3_mix','ligate3_heat_condition','ligate3_heat','ligate3_enzyme','ligate3','ligate5_adapter','ligate5_mix','ligate5_heat_condition','ligate5_heat','ligate5_enzyme','ligate5','pcr_allocation','diagnostic_supplies','diagnostic_prepare','diagnostic_split','preparative_supplies','preparative_prepare','preparative_split','zymo_supplies','zymo','pippin_allocate','pippin_condition','pippin'),SOURCE)
    g.add_source_item(SourceItem('pending_originals','Published main and original oligoTable1/kit/read geometry','Published MOESM3 protocol; exact Table1 and named predecessor/vendor originals not followed in cap.'))
    return g

def selftest():
    from checks import Check
    c=Check();x=SelectedCrosslinkedRNA('ACUU',frozenset({2,3}),frozenset({2}))
    first=x.first_dna(frozenset({2}))
    c('selected RNA diagnostic copy yields complementary firstDNA',first,'AGGT')
    c('RNA-oriented diagnostic C is a report, not physical U substitution',revcomp(first),'ACCT')
    c('full physical RNA and adduct positions preserved',x.sequence,'ACUU')
    try:x.first_dna(frozenset({3}));bad=False
    except ValueError:bad=True
    c('diagnostic selected class does not assume every4SU crosslinked',bad,True)
    g=protocol_graph();c('preparativePCR preparation uses separate unamplified aliquot',g.steps['preparative_prepare'].inputs[0].material,'preparative_cdna')
    c('six diagnostic reactions derive from actual source split',len(g.steps['diagnostic_split'].outputs),6)
    c('three preparative reactions derive from actual source split',len(g.steps['preparative_split'].outputs),3)
    c('compound workup preserves all distinct preparative products',tuple(p.material for p in g.steps['zymo'].inputs if p.role==R.TEMPLATE),tuple('preparative_amp_'+str(i) for i in range(3)))
    c('coldATP uses produced retainedPNK, no second fresh stock',g.steps['pnk'].inputs[1].material,'retained_pnk')
    c('ligation uses adapter already carried through heating',g.steps['ligate3'].inputs[1].material,'ligate3_heated_adapter')
    c('actual Pippin aliquot leaves library stock remainder',tuple(g.steps['pippin_allocate'].outputs),('pippin_input','library_remainder'))
    c('radioactive wash marker survives as source QC rather than waste',g.terminals['pnk_wash_qc'].disposition,Disposition.QC)
    c('PNK cleanup partitions marker and remaining wash separately',tuple(g.steps['pnk_wash'].outputs),('pnk_wash_kept','pnk_wash_qc','pnk_wash_separated'))
    c('published protocol leaves source/geometry publication boundary',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
