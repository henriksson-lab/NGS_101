from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
class B:
    def __init__(self,title,authority,start,kind=K.BIOLOGICAL):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',start,kind,primary=True)
    def stock(self,key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind,inputs=(),geometry=None,product='Source mixture including incomplete/unreacted/off-target classes; occupancy/recovery unassigned',reason=None):
        y=self.g.product(key+'_product',product,K.MOLECULAR);kw={}
        if kind==S.MAGIC:kw['reason']=reason or 'Positive source operation; complete occupied molecular inputs/products and historical formulation unassembled, not a missing-information conclusion.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label,extra=()):
        y=self.g.product(key+'_kept','Actual source retained physical fraction; yield/carryover unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual separated source workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label))+tuple(extra),(y,w));self.g.discard(w,'Actual workup fraction excluded; no perfect specificity or molecular yield claim.');self.x=y
    def split(self,key,label,selected,remainder):
        y=self.g.product(key+'_selected',selected,K.MOLECULAR);r=self.g.product(key+'_remainder',remainder,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r));self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,remainder));self.x=y;return y

def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('source products unassembled','[unassembled]',placeholder=True)]);return s

AUTHORITY='AssociatedAMNIHMS791916/PMC4904710main35–59/Fig1Acaption209/Fig4Acaption227. Oneoriginal, actualMaterialsMethodsDOCXgated; no finalbody/SBSvendor/predecessor recipe inherited.'
def adapter_scene():
    s=Scene();s.strand('biotinylated',[Segment('Source biotinylated adapter strand, wholebody/endpositions unassembled','[biotinylated adapter strand]',placeholder=True)]);s.note('biotinylated','Source biotin and Alexa555 adaptor fluorophore present; exact covalent positions/strand assignment of dye unassembled')
    s.strand('non_biotinylated',[Segment('Source non-biotinylated adapter strand containing uracils','[dU-containing adapter strand]',placeholder=True)],at=0);s.note('non_biotinylated','Source dU resides in this non-biotinylated strand; positions and pairing context unassembled')
    return s

def protocol_graph():
    b=B('Single-molecule nucleosome decoding followed by same-surface DNA sequencing',AUTHORITY,'Actual ESC preparation selectedforH3K27me3/H3K4me3 surfaceDNAsequencing; othercell/drug/peptidebranchesoutsidegraph');g=b.g
    b.op('mnase','SourceMNase digestion/isolationofmononucleosomes',S.MAGIC,(b.stock('mnase_stock','ActualsourceMNasecatalyst; dose/endworkupunassembled'),b.cond('mnase_conditions','Sourcecell/nucleiisolation/digestprogrammeunassembled')),reason='Positive sourceMNasenucleosomepreparation; sourceenddistribution/endprocessing/stoichiometryunassembled, no assumedligationreadyMNaseends.')
    adapter=g.input('adapter_stock','Actual fluorescent biotinylated adaptor preparation with dU-containing non-biotinylated strand; wholeendgeometryunassembled',K.MOLECULAR,state=adapter_scene(),chemical_state={'biotin':'Present on biotinylated adapterstrand, exactpositionunassigned','uracil':'Present on NON-biotinylated adapterstrand, exactpositionsunassigned','fluorophore':'SourceAlexa555 adapterlabel; exactcovalentstrand/position/linkageunassigned','whole_pairing':'Unassembled; two sourcecomponentclasses drawnwithoutinventedduplexregistration'})
    b.op('ligate','Source ligate fluorescent/biotin adaptors to free nucleosomal DNA ends',S.MAGIC,((adapter,R.STRUCTURAL),b.stock('ligation_mix','Actualsourceadapterligationpreparation; catalyst/nucleotide/buffercompositionunassembled',K.MOLECULAR,R.STRUCTURAL)),geometry=AdapterLigationGeometry('unavailable','unavailable',(adapter_scene(),),'Sourcebiotin/fluorophore/dUcomponentidentity known, complete nativeDNAends/adapterstrands/junctiongeometryunassembled.'),reason='Actual nucleosomalendligation positive; exactDNAends/adapterjunctions and productoccupancyunassembled.')
    keep=g.product('gradient_mono','Actualglycerol-gradient selected adaptor-ligated mononucleosome fraction',K.MOLECULAR);rest=g.product('gradient_other','Actualotherglycerolgradientfractions; fatesunassigned',K.MOLECULAR)
    g.react('gradient','Sourceglycerolgradientpurifiesadapter-ligatedmononucleosomes',S.SPLIT,((b.x,R.TEMPLATE),b.cond('gradient_conditions','Sourceglycerolgradient; composition/rotor/fractionnumbersunassembled')),(keep,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Actualunselectedgradientfractions notsame-surfaceproduction.'));b.x=keep
    surface=g.product('captured_surface','ActualPEG-streptavidinboundnucleosomes, captureoccupancy/purityunassigned',K.MOLECULAR);flow=g.product('capture_other','Actualnoncapturedloadingfraction',K.MOLECULAR)
    g.react('surface_capture','SourcePEG-streptavidin-coatedslidecapture',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),b.stock('surface_stock','ActualPEG-streptavidinslide, notmagneticbeads',K.SELECTION,R.SELECTION)),(surface,flow));g.mark_terminal(Terminal(flow.id,Disposition.UNRESOLVED,'Noncapturedfractionnotimplicitlyloaded/returned.'));b.x=surface
    b.op('position_image','SourceTIRF recordsadaptor-labellednucleosome positions',S.GENERIC,(b.cond('position_imaging','Sourceinitialadapterfluorophorelocalization; noDNAsequencebarcodecreated'),))
    b.op('cleave_adapter_dye','Sourcechemicalcleavageofadapterfluorophore AFTERpositionimaging',S.MAGIC,(b.stock('dye_cleavage','Actualsourcefluorophorecleavagepreparation; identity/linkageconditionsunassembled',K.MOLECULAR,R.STRUCTURAL),),reason='Source dye removedbeforehistone-antibodydecoding; exactcleavageproducts/bond and washfatesunassembled. Biotinanchor notsilentlyremovedwithdye.')
    b.op('antibodies','ActualfluorescentH3K27me3/H3K4me3 probesbind/dissociaterepeatedly',S.BINDING,(b.stock('antibody_stocks','Actualsourcefluorescenthistoneantibodypreparations; Alexa647/Alexa555 probes/preciselabelpositionsunassembled',K.MOLECULAR,R.STRUCTURAL),b.cond('antibody_conditions','Sourceprobes/timecourseconditionsunassembled; no permanentcovalentDNAbarcodeclaim')))
    b.op('modification_image','SourceTIRFrecordsbindingeventsatcapturednucleosomepositions',S.GENERIC,(b.cond('modification_imaging','Modificationcalls integrate repeatedbindingevents computationally; no idealizedperfectspecificity'),))
    dna=g.product('retained_surface_dna','Actualsurface-retained double-strandednucleosomalDNA afterhistoneeviction',K.MOLECULAR);histones=g.product('evicted_histones','Actualsalt-evictedhistonecomponents; laterfatesunassigned',K.MOLECULAR)
    g.react('histone_evict','Sourceincreasedsalt evictshistones WHILERETAININGsurfaceDNA',S.MAGIC,((b.x,R.TEMPLATE),b.cond('eviction_conditions','Sourceincreaseinsalt, exactcomposition/programmeunassembled')),(dna,histones),reason='SourceDNAretentionversushistoneeviction explicit; no wholeduplexdenaturation/biotintethercleavageinferred.');g.mark_terminal(Terminal(histones.id,Disposition.UNRESOLVED,'ActualevictedhistonesnotsequencingDNA; solublewashfateunassigned.'));b.x=dna
    b.op('user','SourceUSERexcisesdUinNON-BIOTINYLATEDadapterstrand exposesknownprimingsite',S.MAGIC,(b.stock('user_stock','ActualsourceUSERenzymepreparation; dose/timeunassembled'),),reason='Source strand-specific USER exposure; whole dUpositions/retainedDNAintervals andexcisedproductsunassembled. Biotinylatedanchoringstrandnotremovedbyanalogy.')
    primer=g.input('sequencing_primer','Actualsourcecomplementarysequencingprimer; wholebases/endlabelpositionsunassembled',K.MOLECULAR,chemical_state={'source_label':'Fig4A primerAlexa555 inbindingimage, covalentpositionunassigned'})
    b.op('primer_hybridize','SourceprimerannealstotheUSER-exposedadaptor site',S.HYBRIDIZATION,((primer,R.STRUCTURAL),b.cond('primer_conditions','Wholebindinggeometry/thermalprogrammeunassembled')))
    library=g.product('surface_library','Actualsame-surfaceDNApreparation fordirectsinglemoleculeSBS; nosolubleIlluminalibraryinvented',K.LIBRARY)
    g.react('sequencing_ready','Source directsame-surface singlemolecule sequencing substrate',S.GENERIC,((b.x,R.TEMPLATE),),(library,),changes_molecule=False)
    g.sequence(library,SequencingPlan(limitation='DirectsinglemoleculeSBS/TIRFcycles nucleotideaddition→fluorophore/terminatorcleavage; actualcyclecount/completeprimerends/instrumentreagentinventoryunassembled. No PCR/readpair/indexmaskinvented. Sameflowcellx-ycoordinatesjoinimagingandDNAreadscomputationally, notsequenceUMIs.'))
    g.add_source_item(SourceItem('strand_fates','Fig4A USERacts specifically non-biotinylated dU-containing adapterstrand afterhistoneeviction; sourceDNAdoublestrandssurface-retained.','MainFig1A chemicallycleavesadaptorfluorophore BEFOREantibodyimaging; no assumedbiotinanchorremoval. Wholemodifiedoligos/SBSreactionsunassembled.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();c('Gradient selects ligatednucleosomes',g.steps['gradient'].inputs[0].material,'ligate_product');c('Dye cleavage afterpositionimage',g.steps['cleave_adapter_dye'].inputs[0].material,'position_image_product');c('Antibodies afterdye removal',g.steps['antibodies'].inputs[0].material,'cleave_adapter_dye_product');c('USERconsumesretainedDNA',g.steps['user'].inputs[0].material,'retained_surface_dna');c('Histonesexcludedfromsequencing',g.terminals['evicted_histones'].disposition,Disposition.UNRESOLVED);c('SequencingprimerafterUSER',g.steps['primer_hybridize'].inputs[0].material,'user_product');c('OneactualsurfaceNGSclass',set(g.sequencing),{'surface_library'});c('NoinventedPCR',all(x.kind!=S.PCR for x in g.steps.values()),True)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'science.aad7701_PMC4904710.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1126/science.aad7701')):
        t=p.read_text();c('OriginalUSERnonbiotinstrand','non-biotinylated adaptor strand'in t,True);c('OriginalhistoneevictionretainsDNA','leaving behind double-stranded nucleosomal DNA'in t,True)
    c('Partialnotpublicready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
