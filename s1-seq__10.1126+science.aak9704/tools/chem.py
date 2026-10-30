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


SOURCE='AssociatedacceptedmanuscriptNIHMS841208/PMC5234563 main30–38/Fig1A440; DOI10.1126/science.aak9704. ActualchemistrySIgated; no2018detailedprotocolsubstitution.'
@dataclass(frozen=True)
class SelectedResectedEnd:
    duplex_part:str
    tail:str
    def __post_init__(self):
        if not self.duplex_part or not self.tail or set(self.duplex_part+self.tail)-set('ACGT'):raise ValueError('explicitduplex-associatedstrandand3primeDNAtailexpected')
    def conditional_tail_removal(self):
        # A declared boundary class, not a prediction of exactS1cleavagesites/ends.
        return self.duplex_part,self.tail

def protocol_graph():
    g=ProtocolGraph('S1-seq — sourcecoarseoperationalroute',authority=SOURCE,representation='published_instance')
    x=g.input('meiotic_dna','SourceyeastDNAwithresectedDSB3primetails/recombinationintermediates/backgroundDNA; exactpopulationunassigned',K.MOLECULAR,primary=True)
    x=op(g,'embed',x,'SourceagaroseplugimmobilizationprotectsDNAfromshearing',S.MAGIC,reagent='sourceagaroseplugpreparation/lysis originalSIunreviewed')
    x=op(g,'digest',x,'Source ssDNAnucleaseprocessingof3primetails',S.DIGESTION,reagent='SourceS1Seqnuclease/endprocessinginventoryandconditionsunreviewed',catalyst=True,product='SelectedprocessedssDNA-dsDNAboundaryDNA pluscleavedtailresidues/otherDNA/reagentmaterial; nophysicalremovalorrecreatedtail')
    a=g.input('first_adapter','Source5primebiotinylatedfirstadapter; exactbody/endchemistryunreviewed',K.ADAPTER);lig=g.input('first_ligase','Sourcefirstadapterligationreagentsunreviewed',K.CATALYST);y=g.product('first_ligated','SourcefirstadaptercapturedDNAclass plusunligatedDNA/adapters and retainednucleasedebris; sourceintermediateworkupunreviewed',K.MOLECULAR)
    g.react('first_ligation','Sourcein-plugbiotinadapterligationBEFOREsonication',S.LIGATION,((x,R.TEMPLATE),(a,R.STRUCTURAL),(lig,R.CATALYST)),(y,),geometry=uncertain_ligation())
    x=op(g,'extract',y,'Sourceextractionfromagaroseplugs',S.MAGIC,reagent='Actualplugextraction/recoverygeometrybehindSI',product='SourceextractedDNApreparation; partitionofdigestiondebris/freeadapterunreviewed')
    x=op(g,'sonicate',x,'Sourcepostlabelsonication',S.FRAGMENTATION,product='SourcefragmentedfirstadapterbearingDNA plusunlabelledDNA/fragment/reagentmaterial; genericnewendsnotmisidentifiedasoriginalresectionboundary')
    x=op(g,'biotin_recover',x,'Sourcepurificationofbiotinylatedfragments',S.MAGIC,reagent='Actualsourceaffinitycapture/recoverybehindgatedSI; noMyOneelutionorDTTdefault',product='SourcebiotinenrichedDNApreparation; exactbeads/recovery/retainedandremovedfatesunassembled')
    a=g.input('second_adapter','Sourceoppositeendsecondadapter; fullbody/endgeometryunreviewed',K.ADAPTER);lig=g.input('second_ligase','Sourcesecondadapterligationreagentsunreviewed',K.CATALYST);y=g.product('second_ligated','SourcetwoadapterDNAclassesplusunligatedDNA/adapters/remnants',K.MOLECULAR)
    g.react('second_ligation','SourceoppositeendadapterligationAFTERfragmentpurification',S.LIGATION,((x,R.TEMPLATE),(a,R.STRUCTURAL),(lig,R.CATALYST)),(y,),geometry=uncertain_ligation())
    f=g.input('forward','Sourcesequencinglibraryforwardprimerunreviewed',K.OLIGO);r=g.input('reverse','Sourcesequencinglibraryreverseprimerunreviewed',K.OLIGO);poly=g.input('polymerase','SourcePCRpolymeraseunreviewed',K.CATALYST);out=g.product('amp','SourcelowcyclePCRproducts plus residualoriginalDNA/primers/reagents; noquantitativerecoveryclaim',K.MOLECULAR)
    g.react('pcr','SourcelowcyclelibraryPCR',S.PCR,((y,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(poly,R.CATALYST)),(out,),geometry=PCRGeometry.unavailable(scene('Exactsourceprimer/end/backbonegeometrybehindactualgatedSI'),reason='NoEND-seq/BLESSfamilyrecipeimportedfromsimilarsubstrateorreadout.'))
    g.sequence(out,SequencingPlan(limitation='Actualsource deepsequencingandresectionendpointmapping; exactplatform/cycles/primersunreviewed. Genomicreadcoordinatesareinferredprocessed-boundarypositions,notnativeDSBcleavagepositionsoruniversalpreciseS1cuts.'))
    g.add_source_item(SourceItem('gated_si','Actual2017Supplemental_material.pdf','PMCdownloadgate; figures/MaterialsMethods/TablesS1–2originalunread. Movies/Exo1simulationdonotreplacechemistrySI. Detailed2018methodnotpartofthisoriginalscope.'))
    return g

def selftest():
    from checks import Check
    c=Check();r=SelectedResectedEnd('ACGT','GGA');retained,debris=r.conditional_tail_removal()
    c('selectedtailremovalretainsexplicitreleasednucleotides',retained+debris,r.duplex_part+r.tail)
    c('selectedreleasedtailisnotapredictednativeDSBcoordinate',debris,'GGA')
    try:SelectedResectedEnd('ACGT','');bad=False
    except ValueError:bad=True
    c('tailclassrequiresactualtail',bad,True)
    g=protocol_graph();c('firstadapterprecedessonication',g.steps['sonicate'].inputs[0].material,'extract_product')
    c('oppositeadapterfollowsbiotincapture',g.steps['second_ligation'].inputs[0].material,'biotin_recover_product')
    c('actualgatedoriginalpreventscompletepublication',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
