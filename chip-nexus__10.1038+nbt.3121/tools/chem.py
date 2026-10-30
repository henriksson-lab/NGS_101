"""Protocol-owned source-qualified incomplete graph; no publication marker."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry

def scene(label):
    sc=Scene();sc.strand('preparation',[Segment('population','['+label+']',placeholder=True)]);return sc

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source reagents; geometry unassembled',product=None,state=None,catalyst=False):
    stock=g.input(key+'_stock',reagent,K.CATALYST if catalyst else K.CONDITION)
    out=g.product(key+'_product',product or label+'; molecular mixture/remnants unassembled',K.MOLECULAR,state=state)
    kw={'reason':'Source operation established; exact molecular ends, companion species and geometry unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    g.react(key,label,kind,((x,R.TEMPLATE),(stock,R.CATALYST if catalyst else R.CONDITION)),(out,),**kw);return out

def cleanup(g,key,x,label):
    stock=g.input(key+'_condition',label,K.CONDITION)
    out=g.product(key+'_kept','Source recovered preparation; recovery and carryover unquantified',K.MOLECULAR)
    waste=g.product(key+'_waste','Actual source-separated workup fractions',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(stock,R.CONDITION)),(out,waste));g.discard(waste,'Source physical cleanup, no perfect partition implied.');return out

def uncertain_ligation():return AdapterLigationGeometry('unavailable','unavailable',(scene('unreviewed exact ends and junction geometry'),),'Known source ligation, complete molecular geometry not assembled.')



def qc(g,key,x,label):
    keep=g.product(key+'_retained','Sourceproductionremainder; analyticalaliquotamountnotinvented',K.MOLECULAR);test=g.product(key+'_assay','Separateanalyticalaliquotforthe sourceassay',K.QC)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(keep,test));g.mark_terminal(Terminal(test.id,Disposition.QC,'Analyticalmaterialdoesnotbecomethepreparativesequencingtemplate.'));return keep

def unknown_pcr(g,key,x,label,reagent):
    f=g.input(key+'_f','Sourceforwardindexingprimerfamily; memberandfullbodyunreviewed',K.OLIGO);r=g.input(key+'_r','Sourcereverseindexingprimerfamily; memberandfullbodyunreviewed',K.OLIGO);cat=g.input(key+'_enzyme',reagent,K.CATALYST);y=g.product(key+'_product','SourcePCRcopiesplusoriginaltemplates/unextendedprimersandreactionremnants; exactbackbonesunassembled',K.MOLECULAR)
    g.react(key,label,S.PCR,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=PCRGeometry.unavailable(scene('Sourcefulladapter/primer/endgeometryunassembled'),reason='KnownPCRoperationdoesnotestablishcompletefirstcopyorfinalstrandgeometry.'));return y

import illumina as il
from chemdraw import Construct,feature,revcomp
SOURCE='PublishedMOESM434 SupplementaryProtocol1 December2014 PartsC–E/oligotable(originalp17rendered),DOI10.1038/nbt.3121. AssociatedAMNIHMS647911 Methodsqualifiedcomparison;not2019/2021updatedprotocol.'
U='GATCGGAAGAGCACACGTCTGGATCCACGACGCTCTTCC'
B_FIXED='AGATCGGAAGAGCGTCGTGGATCCAGACGTGTGCTCTTCCGATCT'
CUT='GAAGAGCGTCGTGGATCCAGACGTG'
UMI=feature('nexus.umi','umi','random',note='Source5randomnucleotides;adjacentfixedTCAGnotUMI.')
INDEX=feature('nexus.i7','sample_index','fixed',note='PublishedB01componentCGTGATinPCRprimer;tablelabelATCACGoppositeorientation;notastudysampleassignment.')
@dataclass(frozen=True)
class PSPrimer:
    sequence:str
    bond:int
    def __post_init__(self):
        if not self.sequence or set(self.sequence)-set('ACGT') or self.bond!=len(self.sequence)-2:raise ValueError('SourceonephosphorothioatebondbeforeterminalTrequired')
        if self.sequence[-1]!='T':raise ValueError('TerminalsourceTmustremainactualnucleotide')
    @property
    def extendable_3oh(self):return True
    @property
    def written(self):return self.sequence[:self.bond+1]+'*'+self.sequence[self.bond+1:]
def read_prefix_for_selected_umi(umi):
    if len(umi)!=5 or set(umi)-set('ACGT'):raise ValueError('Explicitselected5ntadapterUMIrequired')
    return revcomp(umi)+'CTGA'
def adapter_b():
    sc=Scene();sc.strand('b',[Segment('sourcefixedanchor','TCAG'),Segment('sourceUMI','N'*5,feature=UMI),Segment('sourceBbody',B_FIXED)],mod5='P');return sc
def adapter_u():
    sc=Scene();sc.strand('u',[Segment('sourceU',U)],mod5='P');return sc
def adapter_duplex():
    sc=Scene();sc.strand('u',[Segment('sourceU',U)],mod5='P')
    sc.anneal('b',[Segment('sourcefixedanchor','TCAG'),Segment('sourceUMI','N'*5,feature=UMI),Segment('source5primeextension','AGATC'),Segment('sourcepairedcore',revcomp(U)),Segment('source3primeT','T')],to='u',pair=('sourcepairedcore','sourceU'),mod5='P')
    sc.ligation_end('b','source3primeT',kind='cohesive',polarity='3-prime',owner='adapter');return sc
def reverse_primer():return Construct([Segment('sourceP7',il.P7),Segment('sourceB01physicalindex','CGTGAT',feature=INDEX),Segment('sourceannealingarm',il.TRUSEQ_READ2)])
def protocol_graph():
    g=ProtocolGraph('ChIP-nexus — published2014detailedprotocol/B01component',authority=SOURCE,representation='published_instance')
    x=g.input('chromatin','SourcecrosslinkedchromatinimmunoprecipitatedonDynabeads;PartA/Bcell/antibodycontextunassigned',K.MOLECULAR,primary=True)
    x=cleanup(g,'initial_wash',x,'SourcewashA/B/C/D/TrispH7.5;retainedbeads')
    x=op(g,'end_repair',x,'SourceT4DNApolymerase/Klenow/T4PNK/dNTPs12C30min',S.END_REPAIR,reagent='PublishedSIPartCend-repairenzymes',catalyst=True);x=cleanup(g,'wash_repair',x,'SourcewashA–D/TrispH8afterendrepair')
    x=op(g,'a_tail',x,'SourceKlenowexo-/dATP37C30min',S.TAILING,reagent='PublishedSIdATP,notAMprintedATP;noordinaryrevisionoverridden',catalyst=True);x=cleanup(g,'wash_tail',x,'SourcewashA–D/TrispH7.5afterdAtailing')
    u=g.input('adapter_u','SourceNexadapterUBamHI5phosphate',K.OLIGO,state=adapter_u());b=g.input('adapter_b','SourceNexadapterBN5BamHI5phosphate/5randomUMI',K.OLIGO,state=adapter_b())
    duplex=g.product('adapter','SourcepairedNexadapterwith5primebarcodeoverhangand3primeT',K.MOLECULAR,state=adapter_duplex());g.react('anneal_adapter','Source95C5mincoolto25C/hold30min',S.HYBRIDIZATION,((u,R.TEMPLATE),(b,R.STRUCTURAL)),(duplex,))
    lig=g.input('ligase','SourceQuickT4ligase25C60min',K.CATALYST);y=g.product('ligated','Sourceadapter-bearingchromatinplusfreeadapters/unligatedDNA/carryover;success/yieldunpredicted',K.MOLECULAR)
    geom=AdapterLigationGeometry.known('cohesive',adapter_duplex(),scene('SourcegenomicA-tailedend/fulljunctiongeometryunassembled'),statement='Sourceadapter3primeTendknown;sealedgenomicproductnotassembled.')
    g.react('adapter_ligate','SourceNexadapterligationtoA-tailedDNA',S.LIGATION,((x,R.TEMPLATE),(duplex,R.STRUCTURAL),(lig,R.CATALYST)),(y,),geometry=geom);x=cleanup(g,'wash_ligate',y,'SourcewashA–D/TrispH8afterligation')
    for key,label,reagent in [('fill','SourceKlenowexo-/dNTP37C30minadapteroverhangfill','Klenowexo-/dNTP'),('trim','SourceT4DNApolymerase12C5minendtrim','T4DNApolymerase/dNTP'),('lambda','Sourcelambda37C60min5primeexonucleaseresection','Lambda/DMSO/Triton'),('recj','SourceRecJf37C60minssDNAexonuclease','RecJf/DMSO/Triton')]:
        x=op(g,key,x,label,S.EXTENSION if key=='fill' else S.DIGESTION,reagent=reagent,catalyst=True,product='SourceprocessedDNAclassesplusoriginalDNA,removedDNA/adapterproducts,proteins/beads/reagents;noinstantperfectdepletion')
        x=cleanup(g,'wash_'+key,x,'Source'+('threeRIPAwashes' if key=='recj' else 'A–D/Triswash')+'after'+key)
    x=op(g,'elution',x,'Sourceelutionbuffer65C20min',product='SourceelutedDNA/protein/RNA/reagentsplusbeadsuntilsupernatantcollection')
    x=cleanup(g,'collect_supernatant',x,'Sourcespin1minandtransfer100ulsupernatantinto100ulTE')
    x=op(g,'reverse_crosslink',x,'Source65C6hcrosslinkreversal',product='DNA/RNA/protein/remnantmixture;noRNAorproteinclearanceasserted')
    x=op(g,'rnase',x,'SourceRNaseA37C2h',S.DIGESTION,reagent='SourceRNaseA',catalyst=True,product='DNAplusdigestedRNAproducts/proteins/reagentsuntilphenolworkup')
    x=op(g,'pk',x,'SourceproteinaseK55C2h',S.DEPROTEINIZATION,reagent='SourcePK',catalyst=True)
    x=cleanup(g,'phenol',x,'Sourcephenolthenphenol/chloroform/IAAextraction')
    x=cleanup(g,'ethanol',x,'Sourceglycogen/NaOAc/EtOHminus80C2h,pellet/70%wash,resuspend12ul')
    assay=g.product('qpcr_aliquot','Source1ulqPCRanalyticalallocation;inputcontrolreactionseparate',K.QC);keep=g.product('library_11ul','Source11ulretainedDNAlibrarypreparation',K.MOLECULAR);g.react('allocate_qpcr','Source12ulsplit1ulqPCR/11ullibrary',S.SPLIT,((x,R.TEMPLATE),),(assay,keep));g.mark_terminal(Terminal(assay.id,Disposition.QC,'SourceqPCRmaterialdoesnotbecomelibrarytemplate.'))
    x=op(g,'denature',keep,'Source95C5mindenaturethenice',S.THERMAL_TREATMENT)
    x=op(g,'circle',x,'SourceCircLigase/ATP/Mn60C1hssDNAself-circularization',reagent='SourceCircLigase;selectedcircleendgeometryunassembled',catalyst=True,product='SelectedssDNAcircleclassesplusuncircularizedDNA/enzymes/reagents;notallmoleculesassertedcircular')
    cut=g.input('cut_oligo','SourceunmodifiedNexcutBamHI',K.OLIGO,state=Construct([Segment('sourcecutoligo',CUT)]));y=g.product('cut_hybrid','Sourcecut-oligo/circleduplexwindowplusunboundcutoligo/linearDNA/remnants',K.MOLECULAR);g.react('anneal_cut','PublishedSI95C1min,coolto25C/hold30min;notAM5minsubstitution',S.HYBRIDIZATION,((x,R.TEMPLATE),(cut,R.STRUCTURAL)),(y,))
    x=op(g,'bamhi',y,'SourceFastDigestBamHI37C30min',S.DIGESTION,reagent='SourceBamHI/buffer',catalyst=True,product='SelectedrelinearizedlibraryDNAplusdigestedcutoligoproducts,uncutcircles,linearDNA/reagents;nodigestioncleanupassumed')
    x=cleanup(g,'postcut_ethanol',x,'SourceTE/NaOAc/ethanolprecipitation,pellet/cold70%wash,dry/resuspend36ul')
    f=g.input('forward','SourceNexprimerUterminalPSbond',K.OLIGO,state=Construct([Segment('sourceprimerU',il.NEBNEXT_UNIVERSAL_PRIMER)]),chemical_state=PSPrimer(il.NEBNEXT_UNIVERSAL_PRIMER,len(il.NEBNEXT_UNIVERSAL_PRIMER)-2));r=g.input('reverse','SourceB01primercomponentterminalPSbond;notstudysampleassignment',K.OLIGO,state=reverse_primer(),chemical_state=PSPrimer(il.P7+'CGTGAT'+il.TRUSEQ_READ2,len(il.P7+'CGTGAT'+il.TRUSEQ_READ2)-2));cat=g.input('phusion','SourcePhusion/HFbuffer/dNTP',K.CATALYST);xout=g.product('pcr_product','Source15–18cyclePCRcopiesplusoriginalDNA/cutoligoremnants/unextendedprimers,reagents;fullproductgeometryunassembled',K.MOLECULAR)
    g.react('pcr','PublishedSIPhusion15–18cycles;notuniversallyAM18',S.PCR,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST)),(xout,),geometry=PCRGeometry.unavailable(scene('Completecircle-cut/PCRtemplategeometryunassembled'),reason='Knownsourceoligosdonotestablishwholegenomiclibrarycopies.'))
    x=cleanup(g,'gel',xout,'Source2%gel/MinEluteRTdissolution/elute12ul;excludeprimer/adapterdimers')
    x=qc(g,'bioanalyzer',x,'SourceBioanalyzeranalyticalallocationandretainedproductionlibrary')
    lib=g.product('library','Sourcephysicalsequencinglibrary;fullmolecular/readgeometryunassembled',K.LIBRARY);g.react('library_identity','Sourcefinalpreparationidentity',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='AssociatedmainreportsIlluminaHiSeqSE50;publishedSIactivepreparation,fullreadprimerlanding/indexgeometryunassembled. Selectedreadprefix5UMI+CTGA9nt;computationalfilter/trim/dedupdoesnotdeleteDNAorconstituteperfectmoleculeclassification.'))
    g.add_source_item(SourceItem('pending','Finalmainapplicability,unreviewedPartA/B/fulllibrarygeometryanddirectsourceprotocols','Published2014SupplementaryProtocolactive;AMonlyqualifiedcomparison. No2019/2021updatedrecipeorSOLiDsourceimport.'))
    return g
def selftest():
    from checks import Check
    c=Check();g=protocol_graph();c('sourcebarcodepairedcorebuildsvalidduplex',len(adapter_duplex().strands),2)
    c('selectedUMIreadorientationreversecomplementsactualadapter',read_prefix_for_selected_umi('AACGT'),'ACGTTCTGA')
    c('onlyfivevariablebasesareninebaseprefixUMI',len(read_prefix_for_selected_umi('AACGT'))-len('CTGA'),5)
    p=PSPrimer(il.NEBNEXT_UNIVERSAL_PRIMER,len(il.NEBNEXT_UNIVERSAL_PRIMER)-2);c('terminalPSretainsTandextendability',(p.written[-3:],p.extendable_3oh),('C*T',True))
    try:PSPrimer('ACGT',1);bad=False
    except ValueError:bad=True
    c('wrongPSbondplacementrejected',bad,True)
    c('sourceqPCRdoesnotfeedcircle',g.steps['denature'].inputs[0].material,'library_11ul')
    c('cutdebrisreachessourceethanolcleanup',g.steps['postcut_ethanol'].inputs[0].material,'bamhi_product')
    c('noRNA-templateassumptionfromiCLIPancestry','reverse_transcription' not in g.steps,True)
    c('geometrylimitsremainprivate',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
