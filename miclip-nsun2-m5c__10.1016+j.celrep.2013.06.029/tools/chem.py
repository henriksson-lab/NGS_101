"""Source-owned NSUN2 methylation-iCLIP partial graph; no predecessor recipe import."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'lib'))
from chemdraw import Scene,Segment,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem,SequencingPlan
from adapter_ligation import AdapterLigationGeometry
UMI_LENGTH=5
AUTHORITY='PublishedPMC3730056 main86–90/134/158; consolidatedmmc9extendedmethods595–611/FigS1actualp11; citedpublishedKönig2010SI2098–2121 barcodecomparison only. Closed3originals.'
@dataclass(frozen=True)
class TrappedClass:
    attacking_cysteine:int=321
    release_mutant:str='C271A'
    peptide_remnant:bool=False
    def __post_init__(self):
        if self.attacking_cysteine!=321 or self.release_mutant!='C271A':
            raise ValueError('Selected source trapped class requires Cys321 and release-mutant C271A')
    def proteolyze(self):return replace(self,peptide_remnant=True)
    def descriptor(self):
        return dict(attacking_cysteine=self.attacking_cysteine,release_mutant=self.release_mutant,
                    protein_is_nucleotide=False,covalent_connection_retained=True,
                    peptide_remnant_class=self.peptide_remnant,uv_crosslinking=False,
                    scope='Selected qualitative protein/peptide–cytosine-ring class, not complete occupied molecules or atom mapping; exact linkage atom/methylation occupancy/remnant unassembled')
def rt_scene():
    s=Scene();s.strand('RT primer',[Segment('unknown5flank','[unassembled5flank]',placeholder=True),
        Segment('source5randomUMI','N'*UMI_LENGTH,feature=feature('nsun2_miclip.umi','umi','random',note='2013explicitfive-random cDNA tracing; completeprimerposition/featurecyclesunassigned')),
        Segment('inverseadapter1','[unassembledadapter1]',placeholder=True),
        Segment('BamHI region','[unassembledBamHIregion]',placeholder=True),
        Segment('inverseadapter2','[unassembledadapter2]',placeholder=True)])
    return s

def protocol_graph():
    g=ProtocolGraph('NSUN2 C271A miCLIP selected HEK293/25cycle route',authority=AUTHORITY,representation='published_instance')
    x=g.input('cells','ActualsourceHEK293cells; wholetranscriptome/complexoccupancyunassembled',K.BIOLOGICAL,primary=True)
    def inp(k,label,kind=K.CONDITION,role=R.CONDITION,**kw):return g.input(k,label,kind,**kw),role
    def op(k,label,ports=(),product='Actualsourcepreparation; wholeoccupiedmolecules/endclasses/yields/offproductsunassembled',kind=S.MAGIC,chemical_state=None,**kw):
        nonlocal x
        y=g.product(k+'_product',product,K.MOLECULAR,**({'chemical_state':chemical_state}if chemical_state is not None else{}))
        if kind==S.MAGIC:kw.setdefault('reason','Positive sourceoperation; complete occupied molecular substrates/products and end distributions unassembled, not a source absence claim.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        g.react(k,label,kind,((x,R.TEMPLATE),)+tuple(ports),(y,),**kw);x=y
    def clean(k,label):
        nonlocal x
        y=g.product(k+'_kept','Actualretainedsourcefraction; recovery/carryoverunassigned',K.MOLECULAR);w=g.product(k+'_waste','Actualexcludedworkupfraction',K.WASTE)
        g.react(k,label,S.CLEANUP,((x,R.TEMPLATE),inp(k+'_conditions',label)),(y,w));g.discard(w,'Physicalsourceworkup, no perfect purity claim');x=y
    def split(k,label,qc=False):
        nonlocal x
        y=g.product(k+'_used',label,K.MOLECULAR);r=g.product(k+'_rest','Otheractualstockfraction; amount/fateunassigned',K.MOLECULAR)
        g.react(k,label,S.SPLIT,((x,R.TEMPLATE),),(y,r))
        g.mark_terminal(Terminal((y if qc else r).id,Disposition.QC if qc else Disposition.UNRESOLVED,label if qc else 'Untransferredstock; amount/fateunassigned'));x=r if qc else y
    op('express_trap','SourceMycNSUN2C271A/Lipofectamine2000;24hharvest, noUV',(
        inp('expression_construct','ActualMycNSUN2C271Aconstruct; wholeplasmid/expressionunassembled',K.MOLECULAR,R.STRUCTURAL),
        inp('transfection','ActualLipofectamine2000/24h')),
        product='Sourcecells containing covalently trappedNSUN2/RNAclass plusfree/non-target/untransfectedclasses',chemical_state=TrappedClass().descriptor())
    op('lyse','Source50mMTris7.4/100mMNaCl/1percentNP40/.1percentSDS/.5percentdeoxycholate',(
        inp('lysis_buffer','Actualsourcelysisbuffer/stockamountsunassigned'),))
    op('nuclease','SourcehighDNase/lowRNaseI;FigS1selected1:1000RNaseI',(
        inp('dnase','ActualDNase; identity/doseunassigned',K.CATALYST,R.CATALYST),
        inp('rnase','ActualRNaseI1:1000dilution; stockactivityunassigned',K.CATALYST,R.CATALYST)))
    clean('clear','Source13,000rpm15min4C; RETAINsupernatant, removepellet')
    bound=g.product('bound','ActualProteinG/antiMyc-boundfraction; purityunassigned',K.MOLECULAR);flow=g.product('flow','Actualnonretainedcapturefraction',K.MOLECULAR)
    g.react('ip','SourceProteinGDynabeads/antiMyc9E10capture',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),inp('beads','ActualProteinGDynabeads',K.SELECTION,R.SELECTION),inp('myc_antibody','ActualantiMyc9E10antibody',K.MOLECULAR,R.STRUCTURAL)),(bound,flow));g.mark_terminal(Terminal(flow.id,Disposition.UNRESOLVED,'Nonretainedfractionlaterfateunassigned'));x=bound
    clean('wash','Sourcestringentwashing/retainboundfraction; formulationsunassembled')
    op('dephos3','SourceT4PNK3primeenddephosphorylation',(
        inp('pnk3','ActualT4PNKfor3endprocessing',K.CATALYST,R.CATALYST),inp('pnk3_conditions','Actualsourceprocessingconditions/dosesunassigned')),
        product='Source3OH-processedRNAcomplexclass plusunprocessed/copresentclasses; wholeendsunassembled')
    op('linker_ligate','Sourcepreadenylatedlinker/RNAligase after3endprocessing',(
        inp('linker','Actualpreadenylatedlinkeroligo; body/blocker/endgeometryunassembled',K.ADAPTER,R.STRUCTURAL,chemical_state={'source_activation':'preadenylated','scope':'Exact donor geometry and whole ordering unassembled; no borrowed2010 linker'}),
        inp('rna_ligase','ActualNEBRNAligase; isoform/doseunassigned',K.CATALYST,R.CATALYST)),
        geometry=AdapterLigationGeometry.unavailable('Source preadenylated donor activation is explicit; whole donor body, blocker, duplex/splint geometry and RNA junction are unassembled. No generic iCLIP adapter substituted.'))
    op('label5','SourceSEPARATET4PNK/32P-ATP5endlabeling',(
        inp('pnk5','ActualT4PNKfor5endlabeling; distinctsourceoperation',K.CATALYST,R.CATALYST),
        inp('radio_atp','Actual32P-ATPphosphatedonor; isotopomer/doseunassigned',K.MOLECULAR,R.STRUCTURAL)),
        product='Source5endradiolabeledcomplexclass plusunlabeled/copresentclasses; isotopeposition/yieldunassigned')
    op('elute','SourceeluteproteinRNAcomplexes; detailedreagentsunassembled')
    clean('gel_membrane','Sourcedenaturinggel/nitrocellulosetransfer/selecthigherMWradioactivesignal; actualdiscardedgel/membraneregions')
    op('proteinase','SourceProteinaseK/3.5MureaRNArecovery',(
        inp('proteinase_k','ActualProteinaseK',K.CATALYST,R.CATALYST),inp('urea','Actual3.5Murearecoverybuffer')),
        product='RNAwithresidualcovalentpeptideclass plusotherrecoverycomponents; not asserted unmodifiedRNA',chemical_state=TrappedClass().proteolyze().descriptor())
    rt=g.input('rt_primer','Actual2013RToligos; inverseadapters/BamHIregion/five-randombarcode; wholeordersunassembled',K.OLIGO,state=rt_scene(),chemical_state={'barcode_random_bases':UMI_LENGTH,'five_prime_phosphate':'unassigned','source_version':'2013five-random; not2010three-random/two-experiment'})
    supplies=g.input('rt_supplies','ActualRTreagentinventory; enzyme/NTidentities/doses/packagingunassembled',K.MOLECULAR)
    enzyme=g.product('rt_enzyme','PositiveRTenzymeclassfromactualinventory',K.CATALYST);nts=g.product('rt_nts','PositiveRTnucleotideclassfromactualinventory',K.MOLECULAR)
    g.react('rt_inventory','Componentinventoryonly; notclaimedpremixedstockornewadditions',S.GENERIC,((supplies,R.TEMPLATE),),(enzyme,nts),changes_molecule=False)
    op('rt','SourceRTstopsatpeptideRNAadduct; actualRTprimerusedonce',((rt,R.STRUCTURAL),(enzyme,R.CATALYST),(nts,R.STRUCTURAL)),product='SourcecDNA/recoverymixture; adduct-truncated classes andoffproducts, wholeRTjunctions/read-startrelationunassembled')
    clean('cdna_size','SourceTBEUreacDNAsizepurification; FigS1H/M/Lselectedidentityunassigned, not invented fixed sizecut')
    op('circle','SourceCircLigaseIIcircularization',(
        inp('circligase','ActualCircLigaseII',K.CATALYST,R.CATALYST),inp('circle_conditions','Actualsourcebuffer/cofactors/dosesunassembled')),
        product='SourcecDNAcirclepreparation pluslinear/offproductclasses; 5P/3OH/wholecirclesunassembled')
    op('bamhi_anneal','SourceannealBamHI-region-complementaryoligo',(
        inp('bamhi_oligo','ActualBamHI-region-complementaryoligo; whole2013bodyunassembled',K.OLIGO,R.STRUCTURAL),))
    op('bamhi_cut','SourceBamHIlinearization',(
        inp('bamhi','ActualBamHInuclease',K.CATALYST,R.CATALYST),inp('bamhi_conditions','Actualreactionconditions/dosesunassembled')))
    pcrsup=g.input('pcr_inventory','ActualPCRreagents; packaging andidentities/dosesunassembled',K.MOLECULAR)
    f=g.product('pcr_f','Sourceadapter-complementaryforwardprimerclass',K.OLIGO);r=g.product('pcr_r','Sourceadapter-complementaryreverseprimerclass',K.OLIGO);pol=g.product('pcr_pol','ActualDNApolymeraseclass',K.CATALYST);nt=g.product('pcr_nt','ActualPCRnucleotideclass',K.MOLECULAR)
    g.react('pcr_components','ActualPCRinventorycomponents; no inventedkit/premix/freeprimerextraaddition',S.GENERIC,((pcrsup,R.TEMPLATE),),(f,r,pol,nt),changes_molecule=False)
    op('pcr','Source25cyclebranch;35cycleconditionseparate',((f,R.STRUCTURAL),(r,R.STRUCTURAL),(pol,R.CATALYST),(nt,R.STRUCTURAL)),product='SourceamplifiedmiCLIPcDNApreparation; completeamplicons/indexesunassembled')
    split('gel_qc','ActualsourcepostPCRgel/EtBranalyticalallocation; amountunassigned',qc=True)
    split('loading','ActualGA2sequencingallocation; volumeunassigned')
    lib=g.product('library','SourceNSUN2miCLIPNGSlibrary; wholeconstruct/readlayoutunassembled',K.LIBRARY)
    g.react('load','SourceIlluminaGA2; cycles/loadingdetailsunassembled',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='2013five-randomUMI reported, completeprimerlanding/readcycles/sourcewholeRTsequenceunassembled; no2010barcode/GA2runlengthorTruSeqvalidationkitinheritance.'))
    g.add_source_item(SourceItem('version_boundary','2013Methods158/mmc9five-random cDNA tracing; cited2010publishedSI2098–2121three-random+two-experiment.', 'Differentmethodversions, no publication-errorclaim or silentthree-for-five replacement; wholepredecessorprimer/catalystrecipesnotinherited.'))
    g.add_source_item(SourceItem('chemical_scope','Main86/Fig1: Cys321attacksRNAcytosinering; C271Aimpairsrelease; noUV.', 'Qualitative protein/peptide–RNAclassonly; Ala271isnotdeclarednucleophile, peptideisnotnucleotide, preciseatoms/methylationoccupancyunassembled.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    p=sources_root()/HERE.parent.name/'j.celrep.2013.06.029_PMC3730056.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1016/j.celrep.2013.06.029')):
        text=p.read_text();c('Sourcefive-randomcDNAtracing','5-nucleotide random barcode to enable tracing of individual cDNAs'in text,True);c('Sourcecatalytic321','cysteine 321 of NSun2'in text,True);c('SourcenoUV','without UV crosslinking'in text,True)
    try:TrappedClass(attacking_cysteine=271);bad=False
    except ValueError:bad=True
    c('Rejectmutantpositionasattackingcysteine',bad,True)
    c('Proteolysisretainstrappedlink',TrappedClass().proteolyze().descriptor()['covalent_connection_retained'],True)
    c('Peptideisnotnucleotide',TrappedClass().proteolyze().descriptor()['protein_is_nucleotide'],False)
    c('RTusesfive-randomdeclaration',g.materials['rt_primer'].chemical_state['barcode_random_bases'],5)
    c('Laterlabelingdoesnotreplace3processing',g.steps['label5'].inputs[0].material,'linker_ligate_product')
    c('RTusespeptiderecovery',g.steps['rt'].inputs[0].material,'proteinase_product')
    c('BamHIafterannealing',g.steps['bamhi_cut'].inputs[0].material,'bamhi_anneal_product')
    c('KnownPCRprimersfromactualinventory',g.producer_of('pcr_f').id,'pcr_components')
    c('QCalicannotsequence',g.terminals['gel_qc_used'].disposition,Disposition.QC)
    c('Wholelibraryunassembled',g.materials['library'].state,None)
    c('Incompleteprivategraph',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
