"""Source-qualified Chem-seq page; historical library finishing remains unassembled."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem,SourceCoverage
AUTHORITY='AssociatedAMPMC4189815/NIHMS616754 titledGenome-wide determination of drug localization, Methods114–120/Results60–64/actualFig1PDFp14; publishedSI5actualp3ligandscheme. FinalpublishedmainGenome-wide localization of small molecules fullbodyunread, nofullversionequivalence. Source-qualifiedJQ1routesonly.'
PROBE={'source_name':'bio-JQ1(S)','ligand_class':'active-enantiomerBETbromodomainligand','capture_handle':'biotin','probe_covalent_linkage':'sourcePEG-linkedamideclass','probe_target_binding':'sourceaffinity/noncovalentclass','not_an_identifier':'NoDNAbarcode, UMI orgenomicsequenceencodedbybiotin','full_atommap':'unassembled'}
IN_VIVO={'cells_per_sample':200000000,'probe_uM':5,'formaldehyde_percent':1,'minutes':20,'simultaneous':True}
IN_VITRO={'formaldehyde_percent':1,'fix_minutes':20,'probe_on_beads_uM':200,'bead_loading_hours':6,'bead_loading_washes':4}
RECOVERY={'elution_SDS_percent':10,'elution_C':65,'elution_minutes':15,'crosslink_reversal':'overnight65C'}

class Trace:
    def __init__(self,g,x):self.g=g;self.x=x
    def inp(self,key,label,kind=K.CONDITION,role=R.CONDITION,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cat(self,key,label):return self.inp(key,label,K.CATALYST,R.CATALYST)
    def op(self,key,label,ports=(),kind=S.MAGIC,product=None,chemical_state=None):
        y=self.g.product(key+'_product',product or 'Sourcepreparation; wholeoccupiedmolecules/yields/atommaps/offproductsunassembled',K.MOLECULAR,**({'chemical_state':chemical_state}if chemical_state else{}))
        kw={'reason':'Sourcepositiveoperation; fulloccupiedsubstrates/products/atommaps/conversionunassembled, notsourceabsence.'}if kind==S.MAGIC else({'changes_molecule':False}if kind==S.GENERIC else{})
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(ports),(y,),**kw);self.x=y
    def cleanup(self,key,label,ports=()):
        y=self.g.product(key+'_kept','Actualretainedsourcefraction; exactrecovery/carryoverunassigned',K.MOLECULAR);w=self.g.product(key+'_excluded','Actualexcludedsourceworkupfraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),)+tuple(ports)+(self.inp(key+'_conditions',label),),(y,w));self.g.discard(w,'Actualsourceexcludedworkupfraction; noquantitativepurityclaim');self.x=y
    def allocation(self,key,label):
        a=self.g.product(key+'_selected','ActualsourceDNAcontributionforsequencing; amountunassigned',K.MOLECULAR);r=self.g.product(key+'_remainder','UntransferredsourceDNAstock; amount/fateunassigned',K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(a,r));self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,'No inventedstockquantityorconsumptionofentireeluate'));self.x=a


# Explicitly reviewed source inventory: fixed target groups, NOT derived from graph keys.
REVIEWED_COMMON = ('cells', 'bio_jq1', 'quench', 'tris_quench', 'quench_conditions', 'cell_pellet', 'cell_pellet_conditions', 'cell_lysis', 'cell_lysis_conditions', 'nuclei', 'nuclei_conditions', 'sonication', 'sonication_conditions', 'clarification', 'clarify_conditions', 't1_beads', 'bead_bsa', 'bead_preparation_conditions', 'prepare_beads', 'capture', 'capture_washes', 'capture_washes_conditions', 'elute', 'elution_conditions', 'reverse_crosslinks', 'reverse_conditions', 'rnase', 'rnase_enzyme', 'pk', 'proteinase_k', 'dna_purification', 'dna_purification_conditions', 'dna_workup_supplies', 'sequencing_allocation', 'ngs_finishing_inventory', 'ngs_finishing', 'library')
REVIEWED_IN_VIVO = ('drug_fix', 'fa', 'drug_fix_conditions')
REVIEWED_IN_VITRO = ('fix', 'fa', 'fix_conditions', 'bead_loading_wash', 'loading_wash_bsa', 'loading_wash_conditions')
def reviewed_inventory(g, route):
    g.declare_reviewed_scope(REVIEWED_COMMON, 'AssociatedAMPMC4189815Methods114–116 JQ1capture/recovery and namedpreparedligand; publishedSI5p3 supportsbio-JQ1componentclass. Finishing33/34unreviewed, explicitlyunassembled actual boundary.')
    g.declare_reviewed_scope(REVIEWED_IN_VIVO if route=='in_vivo' else REVIEWED_IN_VITRO, 'AssociatedAMMethods114 versus116 independently reviewed selectedroute: simultaneouscellprobe/fix versusfixedchromatin/preloadedwashedbeads; recipes notmerged.')
    g.cover_source_item(SourceCoverage('source_version','modelled','library'))
    g.cover_source_item(SourceCoverage('molecular_relations','modelled','capture'))
    g.cover_source_item(SourceCoverage('finishing_scope','modelled','ngs_finishing'))

def protocol_graph():return preparation_graph('in_vivo')

def preparation_graph(route='in_vivo'):
    if route not in ('in_vivo','in_vitro'):raise ValueError('Select a sourceJQ1route, notAT7519/psoralenrecipeinheritance')
    g=ProtocolGraph('Chem-seq sourcequalifiedMM1S/bio-JQ1 '+route,authority=AUTHORITY,representation='published_instance')
    cells=g.input('cells','ActualsourceexponentiallygrowingMM1.S cells; in-vivo2x10^8/sample, upstreamcultureoutsidegraph',K.BIOLOGICAL,primary=True);b=Trace(g,cells)
    g.add_source_item(SourceItem('source_version','PublishedfinaltitleGenome-wide localization of small molecules; actualassociatedAMtitleGenome-wide determination of drug localization.','FinalmainPDFreturnedHTML. RecipequalifiedtoAMMethods; publishedSIligandsynthesisreadasnamedcomponentauthority, no fullfinal/AMequivalence orpublicationerrorclaim.'))
    g.add_source_item(SourceItem('molecular_relations','Bio-JQ1hascovalentPEG-linkedbiotin; BETbindingretrieveschromatinassociatedviaprotein/DNAformaldehydefixation; streptavidincaptureisaffinity.','Drug–DNAcovalentadductnotinferred; occupancy/bindingstoichiometry/crosslinkatommap andallwholeDNAfragmentsunassembled. No perfectcompound-specificcapture.'))
    g.add_source_item(SourceItem('finishing_scope','AM114–120reportsDNApurification34andmassivelyparallelsequencing/analysis33; selectedsourcegraphdoesnotfollowthoseoriginals.','ActualNGSpositive, but exactadapters/ends/kitstocks/PCRprogramme/platform/readlayout unassembledforselectedJQ1branch. No defaultChIP-seq orAT7519Illumina recipeinherited; noabsenceclaim.'))
    probe=g.input('bio_jq1','Actualsourcepreparedbio-JQ1(S) contribution; synthesisoutsidecapturegraph',K.MOLECULAR,chemical_state=PROBE.copy())
    if route=='in_vivo':
        b.op('drug_fix','AMsource5uMbio-JQ1AND1percentFA20minSIMULTANEOUSLY',ports=((probe,R.STRUCTURAL),b.inp('fa','Actual1percentformaldehyde',K.MOLECULAR,R.STRUCTURAL),b.inp('drug_fix_conditions','Actualcellmedium20minprogramme')),product='Sourcecell/ligand/protein/chromatinfixationclasses; occupiedcomplexesandunreactedmaterialunassembled',chemical_state={'probe':PROBE.copy(),'protein_DNA':'Sourceformaldehyde-crosslinkedclasses','drug_target':'Noncovalentaffinityclass, notdrugDNAbond'})
    else:b.op('fix','AMsourceuntreatedcells1percentFA20min',ports=(b.inp('fa','Actual1percentformaldehyde',K.MOLECULAR,R.STRUCTURAL),b.inp('fix_conditions','Actualsource20minmediumfixation')),product='Sourcefixedcell/chromatinclassesbeforeex-vivoprobeaddition; noinvivoprobeaddition')
    b.op('quench','AMsourceTRISpH7.5final300mM',ports=(b.inp('tris_quench','ActualTRISchemicalquenchcontribution',K.MOLECULAR,R.STRUCTURAL),b.inp('quench_conditions','Actualfinal300mMTRISpH7.5')),kind=S.QUENCHING)
    b.cleanup('cell_pellet','Sourceactualscrape/centrifugepellet/3xPBSwash')
    b.op('cell_lysis','AMsourceHEPES7.5/140mMNaCl/EDTA/glycerol/NP40/Triton/PIC',ports=(b.inp('cell_lysis_conditions','Actualsourcecelllysisbuffer'),))
    b.cleanup('nuclei','Sourceactualretainednuclei/washonceTris8/200mMNaCl/EDTA/EGTA/PIC')
    b.op('sonication','AMsource18W/10x30sice/30sintervals',ports=(b.inp('sonication_conditions','ActualHEPES7.5/140mMNaCl/EDTA/EGTA/1percentTriton/.1percentdeoxycholate/.1percentSDS/PIC'),),product='Sourceshearedfixedchromatin/nucleicacidclasses; fragmentsequences/enddistributionunassembled')
    supernatant=g.product('clarified_sonicate','Actualclearednuclearsonicatesupernatant',K.MOLECULAR);pellet=g.product('clarification_pellet','Actualexcludedinsolublefraction',K.WASTE)
    g.react('clarification','Sourcecentrifugeclearedsonicate/retainsupernatant',S.CLEANUP,((b.x,R.TEMPLATE),b.inp('clarify_conditions','Sourcecentrifugationparametersunassigned')),(supernatant,pellet));g.discard(pellet,'Actualexcludedsourcefraction');b.x=supernatant
    raw=g.input('t1_beads','ActualMyOneStreptavidinT1beadpreparation; amountunassigned',K.MOLECULAR);bsa=g.input('bead_bsa','Actual0.5percentBSAmolecularproteinblocker inPBS',K.MOLECULAR)
    prepared=g.product('capture_beads','ActualsourcepreparedT1selectionbeadstock; occupiedbiotin/probe/blockingclassesunassembled',K.SELECTION,chemical_state={'support':'StreptavidinT1','BSA':'Actual0.5percentblocker','probe':PROBE.copy()if route=='in_vitro'else'Bindsbio-JQ1fromsourcechromatinpreparation, notfreshbeadloadingdose'})
    prep_ports=((raw,R.TEMPLATE),(bsa,R.STRUCTURAL),b.inp('bead_preparation_conditions','ActualPBS/0.5percentBSA sourcebeadpreparation'))+(((probe,R.STRUCTURAL),)if route=='in_vitro'else())
    if route=='in_vitro':
        loaded=g.product('probe_loaded_beads','Actual6h200uMbio-JQ1/0.5percentBSA/T1beadloadingpreparation',K.MOLECULAR)
        g.react('prepare_beads','AMsourceinvitro6h200uMprobepreloading',S.MAGIC,prep_ports,(loaded,),reason='Sourcepreparedprobe-beadstock; exactoccupancy/probestoichiometry/unboundclassesunassembled.')
        unbound=g.product('unbound_loading_probe','Actualexcludedfreeprobebeadloadingwashfractions',K.WASTE)
        g.react('bead_loading_wash','Source4xPBS/0.5percentBSAwashBEFOREchromatincapture',S.CLEANUP,((loaded,R.TEMPLATE),b.inp('loading_wash_bsa','ActualBSAproteininthe4sourcewashbuffers',K.MOLECULAR,R.STRUCTURAL),b.inp('loading_wash_conditions','ActualPBSsourcewashconditions')),(prepared,unbound));g.discard(unbound,'Actualunbounddrug/loadingwashfractionexcluded, no quantitativebindingclaim')
    else:g.react('prepare_beads','AMsourcePBS/0.5percentBSApreblock; nopreblockwashasserted',S.MAGIC,prep_ports,(prepared,),reason='Actualsourceblockedbeadpreparation; no source-unstatedpreblockwashingorfreshprobedoseinvented.')
    bound=g.product('captured_chromatin','ActualselectedT1/biotinligand/protein/chromatinheterogeneouscomplexpreparation',K.MOLECULAR,chemical_state={'biotin_streptavidin':'Affinitycaptureclass','probe_target':'NoncovalentBETbindingclass','protein_DNA':'Sourcefixedchromatinclasses','exact_occupancy':'unassembled'});unbound=g.product('unbound_chromatin','Actualexcludedunboundsonicatefraction',K.WASTE)
    g.react('capture','AMsource16to20h4Cchemicalaffinitychromatincapture',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),(prepared,R.SELECTION)),(bound,unbound));g.discard(unbound,'Sourceunboundcapturefractionexcluded, notperfectspecificity');b.x=bound
    b.cleanup('capture_washes','Source2xsonicationbuffer/1x500mMNaClsonication/1xLiClbuffer/1xTRIS.1mMEDTA')
    eluate=g.product('protein_dna_eluate','Actualsourceprotein-DNAeluate; fullprobe/beadcarryoverunassembled',K.MOLECULAR);residual=g.product('residual_beads','ActualresidualT1beads/retainedmaterial',K.MOLECULAR)
    g.react('elute','AMsourceTris8/10mMEDTA/10percentSDS65C15min',S.MAGIC,((b.x,R.TEMPLATE),b.inp('elution_conditions','Actualsource10percentSDSnotAT7519onepercent')),(eluate,residual),reason='Sourcephysicaleluate; exactligand/streptavidin/protein/crosslinkreleasechemistryandyieldunassembled, notselectivebiotinbondcleavage.');g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'Residualbead/ligandmaterialfateunassigned'));b.x=eluate
    b.op('reverse_crosslinks','AMsourceeluateovernight65Ccrosslinkreversal',ports=(b.inp('reverse_conditions','Actualsourceovernight65C'),),product='Sourceintendedcrosslink-reversedDNA/protein/ligandclasses; exactbondchanges/completenessunassembled')
    b.op('rnase','AMsourceRNasedigestion; exactdose/subtypeunassigned',ports=(b.cat('rnase_enzyme','ActualsourceRNase'),),kind=S.DIGESTION)
    b.op('pk','AMsourceProteinaseKdigestion; doseunassigned',ports=(b.cat('proteinase_k','ActualsourceProteinaseK'),),kind=S.DEPROTEINIZATION)
    b.cleanup('dna_purification','AMsourceDNApurificationdelegatedref34; detailedkit/reagent/endgeometryunassembled',ports=(b.inp('dna_workup_supplies','ActualsourceDNApurificationinventory, chemicalidentities/stockpackagingunassembled',K.MOLECULAR,R.STRUCTURAL),))
    b.allocation('sequencing_allocation','ActualpurifiedDNAallocationforNGSlibrarypreparation; amountunassigned')
    finishing=g.input('ngs_finishing_inventory','ActualreportedNGSfinishing/loadinginventory; source33detailsunread, noadapters/PCR/platform/kitdefaultasserted',K.MOLECULAR)
    library=g.product('library','ActualsourceChem-seqDNAlibraryclass; completeends/molecules/readlandingunassembled',K.LIBRARY)
    g.react('ngs_finishing','AMsourcepurifiedDNAfragmentsmassivelyparallelsequenced; finishingdelegatedsource33',S.MAGIC,((b.x,R.TEMPLATE),(finishing,R.STRUCTURAL)),(library,),reason='ActualNGSendpointreported; detailedphysicalfinishingoperations/reagents/primersequences/reads unassembledinselectedsourcepass, notabsenceevidence.')
    g.sequence(library,SequencingPlan(limitation='SourceJQ1massivelyparallelDNAsequencing/analysisref33; exactplatform/runlength/index/adapter/primerlandingunassembled. AT7519explicitIlluminabranchnotautomaticallyassignedtoselectedJQ1.'))
    reviewed_inventory(g,route)
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();a=protocol_graph();v=preparation_graph('in_vitro')
    c.raises('Rejectsourcebranchmerging',lambda:preparation_graph('psoralen'))
    c('In-vivodrugentersactualsimultaneousfix',a.steps['drug_fix'].inputs[1].material,'bio_jq1')
    c('Noexvivodrugloadingininvivo',any(x.material=='bio_jq1'for x in a.steps['prepare_beads'].inputs),False)
    c('Invitrodrugentersbeadloading',any(x.material=='bio_jq1'for x in v.steps['prepare_beads'].inputs),True)
    c('Invitronofreshdrugincellfix',any(x.material=='bio_jq1'for x in v.steps['fix'].inputs),False)
    c('Actualpreloadedbeadswashedbeforecapture',v.steps['bead_loading_wash'].inputs[0].material,'probe_loaded_beads')
    c('Noimaginedinvivopreblockwash','bead_loading_wash'in a.steps,False)
    c('BSAismolecularnotcondition',a.materials['bead_bsa'].kind,K.MOLECULAR)
    c('Captureusesproducedselectionstock',a.steps['capture'].inputs[1].material,'capture_beads')
    c('Crosslinkreversalusesactualeluate',a.steps['reverse_crosslinks'].inputs[0].material,'protein_dna_eluate')
    c('Residualbeadfateexplicit','residual_beads'in a.terminals,True)
    c('NGSusesactualDNAallocation',a.steps['ngs_finishing'].inputs[0].material,'sequencing_allocation_selected')
    c('Noimaginedpool',any(s.kind==S.POOL for s in a.steps.values()),False)
    c('Noimpliedoligo/PCRrecipe',any(s.kind==S.PCR for s in a.steps.values()),False)
    p=sources_root()/HERE.parent.name/'nbt.2776_PMC4189815.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1038/nbt.2776')):
        t=p.read_text()
        c('ActualAMtitlediffers','Genome-wide determination of drug localization'in t,True)
        c('Source10percentSDSliteral','10% SDS at 65'in t,True)
        c('Sourceactual2x10power8sample','2×108 cells per sample'in t,True)
        c('Source200uMbeadloading','200 μM biotinylated drug'in t,True)
    c('Reviewedin-vivopublicationinventory',a.publication_problems(),[])
    c('Reviewedin-vitropublicationinventory',v.publication_problems(),[])
    a.input('unreviewed_added_stock','Regression-only unreviewed input',K.CONDITION)
    c('Unreviewedinventorycannotautopublish',any('unreviewed_added_stock' in issue and 'source-inventory' in issue for issue in a.publication_problems()),True)
    c.report()

PAGE_FILE = 'chem_seq.html'
def render_page():
    from html import escape
    from page import head,info
    from protocol_views import assemble_page,WorkflowSection,PlanLimitationSection
    g=protocol_graph().require_publishable()
    return assemble_page(g,[head('Chem-seq — bio-JQ1 chromatin capture'),'<div class="wrap" data-protocol-complete="true"><h1>Chem-seq</h1>',
      '<p>Map chromatin associated with a biotinylated small-molecule ligand. The selected in-vivo bio-JQ1 route combines ligand treatment with formaldehyde fixation, captures chromatin on streptavidin beads, and recovers DNA for sequencing.</p>',
      info('Source scope: <a href="https://doi.org/10.1038/nbt.2776">Anders et al., 2014</a>. Operational recipe: associated manuscript PMC4189815/NIHMS616754, titled <em>Genome-wide determination of drug localization</em>; published title: <em>Genome-wide localization of small molecules</em>. Published supplement supports the ligand component. Full final/manuscript equivalence remains unverified.'),
      '<h2>Preparation</h2><p>MM1.S cells receive 5 µM bio-JQ1 and 1% formaldehyde together for 20 min. Nuclei are isolated and sonicated; BSA-blocked streptavidin T1 beads capture the biotin-linked ligand and associated chromatin. Washing precedes eluate recovery, crosslink reversal, RNase/proteinase K treatment and DNA purification.</p>',
      '<h2>Molecular relationships</h2><p>Biotin is covalently linked to the prepared ligand through a PEG-linked amide class. Ligand–BET binding and biotin–streptavidin capture are affinity interactions. The model does not assign a covalent drug–DNA bond, binding occupancy or perfect capture specificity. No DNA barcode or UMI is encoded by biotin.</p>',
      info('The selected manuscript specifies <strong>10% SDS, 65°C for 15 min</strong> for JQ1 eluate recovery. The separate AT7519 branch uses different conditions and is not substituted. In vitro Chem-seq instead fixes untreated cells and preloads beads with 200 µM ligand for 6 h, washing four times before capture; it is not a second dose added to the in-vivo sample.'),
      '<h2>Library and reads</h2><p>Actual DNA sequencing is reported, with purification and sequencing delegated to references 34 and 33. Their detailed finishing recipe, adapter oligos, PCR programme, sequencing platform and run lengths remain unresolved in this selected source scope. No generic ChIP-seq library is supplied to fill those boundaries.</p>',PlanLimitationSection('library'),
      '<details><summary>Complete preparation graph, fractions and chemistry warnings</summary>',WorkflowSection('In-vivo bio-JQ1 preparation'),'</details>',
      '<p><a href="01_chem-seq.html">Source notes and remaining questions</a></p></div>'])

if __name__=='__main__':selftest()
