"""Private source-qualified partial preparation; no complete molecular/absence claim."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature,revcomp
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
class B:
    def __init__(self,title,authority,label,kind=K.MOLECULAR):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',label,kind,primary=True)
    def stock(self,key,label,kind=K.MOLECULAR,role=R.STRUCTURAL,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cat(self,key,label):return self.stock(key,label,K.CATALYST,R.CATALYST)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind=S.MAGIC,inputs=(),product='Actual source preparation; whole occupied substrates/products/side products and conversion unassembled',**kw):
        y=self.g.product(key+'_product',product,K.MOLECULAR)
        if kind==S.MAGIC:kw.setdefault('reason','Positive source operation; whole occupied molecules, end distributions, detailed molecular products and conversion unassembled, not a source-absence claim.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label,inputs=()):
        y=self.g.product(key+'_kept','Actual retained physical source fraction, carryover/recovery unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual excluded workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label))+tuple(inputs),(y,w));self.g.discard(w,'Actual source workup fraction excluded; no quantitative purification claim.');self.x=y;return y
    def aliquot(self,key,label,kept='Selected source aliquot',rest='Untransferred stock, amount/fate unassigned',qc=False):
        y=self.g.product(key+'_selected',kept,K.MOLECULAR);r=self.g.product(key+'_remainder',rest,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r))
        if qc:self.g.mark_terminal(Terminal(y.id,Disposition.QC,kept));self.x=r;return r
        self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,rest));self.x=y;return y



AUTHORITY='AssociatedAMNIHMS1809916/PMC9361928main95–101/231/keyresources1081–1085 plusactualassociatedcorrectionNIHMS1822996PMC9365200. Finalpublishedmain403; noSTAR/vendor/predecessorrecipeequivalence.'
def protocol_graph():
    b=B('RedaC:T-seq selectedWT/NAT10-null reduction preparations',AUTHORITY,'SourceWTHeLatotalRNAafterTRIzol; preparedbiologicalstock, upstreamcultureunassembled');g=b.g
    ko=g.input('nat10_null','SourceNAT10-nullHeLatotalRNAafterseparateTRIzolpreparation; notaliquotofWTstock',K.BIOLOGICAL,primary=True)
    seeds=(('wt',b.x),('null',ko));contributions=[]
    g.add_source_item(SourceItem('correction','Actualassociatedcorrection10.1016/j.molcel.2022.06.022 changescryo-EMFig7Fnonacetylatedcomplexhydrogenbond2.1to2.1or2.7angstrom afterrefinement.','InspectedbodydoesnotstateNGS/reductionrecipechange; notcompletefinal/AMrecipeequivalence. No cryo-EMatomgeometry inheritedintoRNAreductionmodel.'))
    g.add_source_item(SourceItem('kit_catalogue','Associatedmain231rRNAE7850/libraryE7760 versusKeyResourcesrRNAE6310L/libraryE7770S.','Both retained; actualusedkitvariant/wholekitcomponentsunassembled, no silentcatalogueorcurrentvendorequivalence.'))
    g.add_source_item(SourceItem('specificity','Associatedmain95–101NaBH4→tetrahydro-ac4C linkedtoC:Treadmismatches; NaBH4reactswithothermodifiedRNAs andconversionisincomplete.','SelectedWTandNAT10-nullsourcepreparations areindependentbiologicalstocks. Untreatedinputvariantfiltering/genotypecomparison arecomputational, notphysicalpurificationofac4C orabsolutestoichiometryproof. AcidicNaCNBH3ac4C-seq/TGIRTrecipe not imported.'))
    for key,start in seeds:
        b.x=start
        b.op(key+'_dnase','SourceTurboDNaseItreatmentbefore1ugselection',inputs=(b.cat(key+'_dnase_stock','ActualTurboDNaseI, dose/formulationunassembled'),),product='SourceDNase-treatedtotalRNApreparation, residualDNA/recoveryunassigned')
        b.aliquot(key+'_one_ug','Actual1ugDNase-treatedRNAforribodepletion','Actual1ugRNAselected','ActualuntransferredDNase-treatedRNAstock')
        b.op(key+'_ribodeplete','SourceNEBNextrRNAdepletion, mainE7850/keyresourceE6310L unresolved',inputs=(b.stock(key+'_ribodeplete_inventory','Actualhistoricalribodepletionreagentinventory; completeenzyme/probe/formulationandproductfractionsunassembled'),),product='SourceribodepletedRNApreparation; rRNA/offtarget/recoveryfractionsunassembled')
        converted=g.product(key+'_reduced','SourceNaBH4-treatedfragmentedRNAclasses; fullmodification/yield/distributionunassembled',K.MOLECULAR,chemical_state={'source_reduction_class':'ac4C→tetrahydro-ac4C inreactedsubset','base_identity_scope':'ModifiedCresidueclass, notliteralRNA C→U substitution','conversion':'Incomplete/unreacted/otherNaBH4-reactiveRNAclassesunassembled','fragment_class':'Sourceapprox100–200bp wording, notvalidatedsizepurity'})
        g.react(key+'_reduce','Source100mMNaBH4/20ul/55Cdark1h; openlidsevery10min',S.MAGIC,((b.x,R.TEMPLATE),b.stock(key+'_nabh4','ActualNaBH4chemicalcontribution100mM'),b.cond(key+'_reduce_conditions','Actual20ul55Cdark1h/10minventingprogramme')),(converted,),reason='SourcechemicalreductionplusRNAfragmentationcompoundboundary; occupiedmodificationpositions/stereochemistry/atommap/yields/offproductsunassembled. NoNaCNBH3substitution.');b.x=converted
        b.clean(key+'_precipitate','Source.3MNaOAc5.5/2.5xethanolprecipitation',inputs=(b.stock(key+'_linear_acrylamide','Actual15ugmllinearacrylamidemolecularcarrier'),))
        b.aliquot(key+'_library_input','ActualsourcekitRNAinputallocation; quantityunassigned','SelectedtreatedRNAforlibrarypreparation','AnyuntransferredtreatedRNAstock; amount/fateunassigned')
        inventory=g.input(key+'_kit_inventory','ActualhistoricalNEBNextUltraIIDirectionalRNAkitinventory; mainE7760/keyresourceE7770S, packaging/stockcounts/formulationunassembled',K.MOLECULAR)
        enzyme=g.product(key+'_kit_rt','ActualRTcatalystclasswithinlibraryreagentinventory; exactidentityunassembled',K.CATALYST)
        nt=g.product(key+'_kit_nt','Actualcopied-DNAnucleotideinventoryofkit; preciseformulationunassembled',K.MOLECULAR)
        oligos=g.product(key+'_kit_oligos','Actualkitoligoinventory; counts/polymer/wholebodies/endsandindexassignmentunassembled',K.OLIGO,chemical_state={'sample_identity':feature('redac.sample_'+key,'sample_index','unknown',note='Sourcemultiplexedsampleidentity; exactindexarchitecture/sequence/coordinateunassembled')})
        other=g.product(key+'_kit_other','Otheractuallibrarykitreagents/conditions; fullreactionorderingunassembled',K.MOLECULAR)
        g.react(key+'_kit_components','Inventorydecompositiononly, notclaimedonepremixoradditionalfreshsupplies',S.GENERIC,((inventory,R.TEMPLATE),),(enzyme,nt,oligos,other),changes_molecule=False)
        b.op(key+'_library_build','SourceUltraIIDirectionallibraryconstructionwithRT48C',inputs=((enzyme,R.CATALYST),(nt,R.STRUCTURAL),(oligos,R.STRUCTURAL),(other,R.STRUCTURAL),b.cond(key+'_rt48','Actualsource48CRTstep; remainingkitprogrammesunassembled')),product='SourceDNAsequencinglibrarypreparationfromtreatedRNA; reactionordering/wholeduplexes/indices/strandmarkingunassembled')
        b.aliquot(key+'_pool_input','Actualselectedmultiplexingcontribution; volume/proportionunassigned','Actualsourcegenotypelibrarycontribution','Anyuntransferredpreparedlibrarystock')
        contributions.append(b.x)
    others=g.input('other_biological_library_contributions','Otheractualbiologicalreplicatelibrarycontributions; source2biologicalreplicatesmentioned, exactpools/stocks/volumesunassembled',K.MOLECULAR)
    pool=g.product('multiplexed_library_pool','ActualHiSeqmultiplexedlibrarypool; ratiosandreplicateallocationsunassigned',K.MOLECULAR)
    g.react('multiplex','Sourcephysicalmultiplexingofpreparedlibraries, noequimolarassumption',S.POOL,tuple((x,R.TEMPLATE)for x in contributions)+((others,R.TEMPLATE),),(pool,));b.x=pool
    b.aliquot('loading','ActualHiSeqloadingallocation; amountunassigned')
    lib=g.product('library','SourceRedaC:T-seqRNA-derivedDNApool; wholemolecules/readlandingunassembled',K.LIBRARY)
    g.react('load','SourceHiSeq2500/TruSeqV4paired-end126cycles',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='Associatedmain231HiSeq2500TruSeqV4paired-end126cyclesassourcewording; exactindexmasks/primerlanding/wholelibrarystrandorientationunassembled. C:TmappingobservationnotliteralRNA C→U conversion; untreated/genotypefilteringcomputational.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('IndependentbiologicalRNAstockinputs',set(m.id for m in g.materials.values()if m.primary),{'source','nat10_null'})
    c('WTNaBH4issourcechemicalcontribution',g.steps['wt_reduce'].inputs[1].role,R.STRUCTURAL)
    c('NullNaBH4comesafteritsownribodepletion',g.steps['null_reduce'].inputs[0].material,'null_ribodeplete_product')
    c('Actualcarrierincorporatedduringprecipitation',g.steps['wt_precipitate'].inputs[2].material,'wt_linear_acrylamide')
    c('LibrarycopiesselectedreducedRNA',g.steps['wt_library_build'].inputs[0].material,'wt_library_input_selected')
    c('Knownunassembledkitoligoclassesexposed',tuple(g.materials[k+'_kit_oligos'].kind for k in ('wt','null')),(K.OLIGO,K.OLIGO))
    c('Bothgenotypelibrariescontributetopool',tuple(x.material for x in g.steps['multiplex'].inputs[:2]),('wt_pool_input_selected','null_pool_input_selected'))
    c('NoacidicNaCNBH3reagentimported',all('nacnbh3'not in m.id for m in g.materials.values()),True)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'j.molcel.2022.05.016_PMC9361928.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1016/j.molcel.2022.05.016')):
        t=''.join(p.read_text().split());c('Literalreductionconditions','100mMNaBH4ina20μlvolumeat55°Cinthedarkfor1hr'in t,True);c('Actuallinearacrylamidecarrier','15μg/mllinearacrylamide(carrier)'in t,True);c('LiteralRT48','reversetranscriptionstepperformedat48°C'in t,True)
    c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
