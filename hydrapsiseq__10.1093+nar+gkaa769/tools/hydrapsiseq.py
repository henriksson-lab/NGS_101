"""Source-scoped HydraPsiSeq cleavage, end repair and actual library finishing."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene, Segment
from protocol_graph import ProtocolGraph, MaterialKind, InputRole, StepKind, Terminal, Disposition, SequencingPlan
SOURCE='Marchand et al. 2020 doi:10.1093/nar/gkaa769, published Methods and Supplementary Figure S1'
BOUNDARY=('SourceNEBNextSmallRNAkit identity/revision and exact adapter/RT/PCR/read geometry remain unreviewed; actualsource libraryconstruction is represented as an opaque interface, with no modern kit substitution.')
PAGE_FILE='hydrapsiseq.html'

@dataclass(frozen=True)
class CleavageEvent:
    """One source-described U scission event, not a deterministic population predictor."""
    upstream_three_prime: str='phosphate'
    downstream_five_prime: str='phosphate'
    removed_residue: str='U'
    debris: str='released residue-derived molecular products present'
    def __post_init__(self):
        if self.removed_residue!='U':
            raise ValueError('only a selected U-cleavage event is modelled; Ψ is relatively protected')
        if self.upstream_three_prime not in {'phosphate','OH'} or self.downstream_five_prime!='phosphate':
            raise ValueError('event must retain the source-specified downstream 5-prime phosphate')
    def dephosphorylate(self):
        return replace(self,upstream_three_prime='OH')
    def scene(self):
        sc=Scene()
        sc.strand('upstream',[Segment('RNA through N−1','[upstream RNA]',placeholder=True)],
                  label='N−1 fragment',mod3=self.upstream_three_prime)
        sc.strand('downstream',[Segment('RNA from N+1','[downstream RNA]',placeholder=True)],
                  label='N+1 fragment',mod5=self.downstream_five_prime)
        sc.note('downstream','The excised U is absent from both RNA fragments; '+self.debris+'.')
        return sc

def protocol_graph():
    g=ProtocolGraph('HydraPsiSeq — hydrazine/aniline cleavage and end repair',authority=SOURCE,
                    representation='published_instance')
    rna=g.input('rna','Source50–300ngRNA population: Ucleavage stochastic; Ψ relatively protected',
                MaterialKind.MOLECULAR,primary=True)
    hydrazine=g.input('hydrazine','Actualhydrazine supply, source50percentfinal',MaterialKind.MOLECULAR)
    hydro_program=g.input('hydrazine_program','Sourceice30–60min; selectedtime/sampleorigin unassigned',MaterialKind.CONDITION)
    abasic=g.product('abasic_mixture','RNA abasic sites plus released base-derived molecular products',MaterialKind.MOLECULAR)
    g.react('hydrazinolysis','Hydrazine ring opening/abasic-site formation at susceptible U sites',
            StepKind.CHEMICAL_CONVERSION,((rna,InputRole.TEMPLATE),(hydrazine,InputRole.STRUCTURAL),(hydro_program,InputRole.CONDITION)),(abasic,))
    precip=g.product('hydrazine_pellet','RNA-containing ethanol precipitate; exact debris partition unquantified',MaterialKind.MOLECULAR)
    waste=g.product('hydrazine_supernatant','ethanol supernatant and washes',MaterialKind.WASTE)
    glycoblue=g.input('glycoblue','ActualGlycoblue carrier supply for sourcehydrazine/aniline precipitations; allocationsunassigned',MaterialKind.MOLECULAR)
    precipitation=g.input('precipitation_program','Source0.3MNaOAc,pH5.2/ethanol;−80C30min;centrifuge/two80percentEtOHwashes',MaterialKind.CONDITION)
    g.react('precipitate_hydrazine','NaOAc/Glycoblue ethanol precipitation and washes',StepKind.CLEANUP,
            ((abasic,InputRole.TEMPLATE),(glycoblue,InputRole.STRUCTURAL),(precipitation,InputRole.CONDITION)),(precip,waste))
    g.discard(waste,'Physical supernatant/wash separation; no exact composition or quantitative recovery claimed')
    aniline=g.input('aniline','Actualsource1Maniline supply,pH4.5',MaterialKind.MOLECULAR)
    aniline_program=g.input('aniline_program','Source60C15min in dark',MaterialKind.CONDITION)
    event=CleavageEvent(debris='released residue-derived products; prior hydrazine debris partition unquantified')
    cleaved=g.product('cleavage_mixture','N−1 RNA3′p + N+1 RNA5′p + released nonpolymeric products',
                      MaterialKind.MOLECULAR,state=event.scene(),chemical_state=event)
    g.react('aniline_scission','Aniline decomposes abasic sites by β/δ elimination',
            StepKind.CHEMICAL_CONVERSION,((precip,InputRole.TEMPLATE),(aniline,InputRole.STRUCTURAL),(aniline_program,InputRole.CONDITION)),(cleaved,))
    cleaned=replace(event,debris='partition of released products after precipitation unresolved')
    fragments=g.product('recovered_fragments','ethanol-recovered RNA-fragment population',
                        MaterialKind.MOLECULAR,state=cleaned.scene(),chemical_state=cleaned)
    rejected=g.product('aniline_supernatant','aniline-stage precipitation supernatant/washes',MaterialKind.WASTE)
    g.react('precipitate_aniline','Ethanol precipitation after aniline; retain RNA fragments',StepKind.CLEANUP,
            ((cleaved,InputRole.TEMPLATE),(glycoblue,InputRole.STRUCTURAL),(precipitation,InputRole.CONDITION)),(fragments,rejected))
    g.discard(rejected,'Physical supernatant/wash fraction; released products were molecular material before separation')
    pnk=g.input('pnk','Actualsource10UT4PNK3primephosphatase,ATP-free',MaterialKind.CATALYST)
    pnk_buffer=g.input('pnk_buffer','Source100mMTrisHCl,pH6.5/100mMMgOAc;37C6h then65C20min',MaterialKind.CONDITION)
    beta_me=g.input('beta_me','Actualsource5mMbeta-mercaptoethanol contribution',MaterialKind.MOLECULAR)
    repaired=cleaned.dephosphorylate()
    ends=g.product('repaired_fragments','RNA fragment3′OH; downstream5′p retained',MaterialKind.MOLECULAR,
                   state=repaired.scene(),chemical_state=repaired)
    g.react('dephosphorylate','ATP-free T4PNK3′dephosphorylation, then heat inactivation',
            StepKind.DEPHOSPHORYLATION,((fragments,InputRole.TEMPLATE),(pnk,InputRole.CATALYST),(beta_me,InputRole.STRUCTURAL),(pnk_buffer,InputRole.CONDITION)),(ends,))
    extract=g.product('kit_input','phenol/chloroform-extracted and ethanol-recovered kit input RNA',MaterialKind.MOLECULAR,
                      state=ends.state,chemical_state=repaired)
    extraction_waste=g.product('extraction_waste','extraction phases and precipitation supernatant',MaterialKind.WASTE)
    extraction=g.input('extraction_supply','Actualsourcephenol/chloroform thenethanol recovery; solventamounts/ratiosunassigned',MaterialKind.CONDITION)
    g.react('recover_ends','Phenol/chloroform and ethanol recovery after PNK',StepKind.CLEANUP,
            ((ends,InputRole.TEMPLATE),(extraction,InputRole.CONDITION)),(extract,extraction_waste))
    g.discard(extraction_waste,'Physical extraction and precipitation reject')
    kit=g.input('kit_supply','ActualsourceNEBNextSmallRNAlibrarykit reagentcollection; historicaledition,individualstocks/oligos unassembled, notapremix',MaterialKind.MOLECULAR)
    library=g.product('library','DNA-library material from actualsourcekitconversion; wholepolymer/purityunassigned',MaterialKind.LIBRARY)
    g.react('kit_conversion','Convert recovered RNA to source sequencing library',StepKind.MAGIC,
        ((extract,InputRole.TEMPLATE),(kit,InputRole.STRUCTURAL)),(library,),
        reason='PublishedMethods explicitly converts RNA to library with the namedkit. This is an actualsource-finishing interface, not sequencingRNA. Adapter/RT/PCR/workup reactions and individualstock packaging remain opaque; no modernvendorrecipe or finalDNAsequence invented.')
    qc=g.product('library_qc','ActualsourceBioAnalyzerHighSensitivityDNA/Qubit material; amountsunassigned',MaterialKind.QC)
    retained=g.product('retained_library','Source libraryretained after analyticalselection',MaterialKind.LIBRARY)
    g.react('select_qc','Reserve source final-library analyticalmaterial',StepKind.SPLIT,((library,InputRole.TEMPLATE),),(qc,retained),changes_molecule=False)
    g.mark_terminal(Terminal(qc.id,Disposition.QC,'Detailed analyticalassays outside this preparation; analytical aliquots are not returned to the library.'))
    contribution=g.product('pool_contribution','Selected librarycontribution for source multiplexing; amountunassigned',MaterialKind.LIBRARY)
    stock_rest=g.product('library_remainder','Untransferred library stock; amountunassigned',MaterialKind.LIBRARY)
    g.react('select_contribution','Select source library contribution for multiplexing',StepKind.SPLIT,((retained,InputRole.TEMPLATE),),(contribution,stock_rest),changes_molecule=False)
    g.mark_terminal(Terminal(stock_rest.id,Disposition.UNRESOLVED,'No whole-library consumption assumed.'))
    other=g.input('other_libraries','Other source libraries for actualmultiplexedrun; identities,indexmembers,countsunassigned',MaterialKind.LIBRARY)
    pool=g.product('multiplex_pool','Actualsource multiplexed librarypool; compositionunassigned',MaterialKind.LIBRARY)
    g.react('multiplex','Combine distinct source library contributions',StepKind.POOL,((contribution,InputRole.TEMPLATE),(other,InputRole.TEMPLATE)),(pool,),changes_molecule=False)
    aliquot=g.product('loading_aliquot','Selected sequencing loadingaliquot',MaterialKind.LIBRARY)
    pool_rest=g.product('pool_remainder','Untransferred multiplexpool stock; amountunassigned',MaterialKind.LIBRARY)
    g.react('select_loading','Select actualsource pool loadingaliquot',StepKind.SPLIT,((pool,InputRole.TEMPLATE),),(aliquot,pool_rest),changes_molecule=False)
    g.mark_terminal(Terminal(pool_rest.id,Disposition.UNRESOLVED,'Sourceflowcell loading doesnot imply consumption of the entire multiplexpool.'))
    concentration=g.input('loading_conditions','PublishedsourceHiSeq1000,flowcell10–12pM; exactdilutionvolumes/runprimersunassigned',MaterialKind.CONDITION)
    loading=g.product('sequencing_input','Actualselected DNA-library flowcell-preparation interface',MaterialKind.LIBRARY)
    g.react('flowcell_prep','Prepare selected source libraryloading concentration',StepKind.MAGIC,((aliquot,InputRole.TEMPLATE),(concentration,InputRole.CONDITION)),(loading,),changes_molecule=False,reason='Source10–12pMflowcellpreparation; detailedhybridization/runprimerreagents remain unassembled.')
    g.sequence(loading,SequencingPlan(acquired_roles=('Read 1',),read_lengths={'Read 1':50},limitation='PublishedHiSeq1000single-read50bp and10–12pM. Exact historicalkitoligos, indexmask/loadedprimers and wholelibrary geometry unassembled; read-end scoring doesnot validate a primerlanding.'))
    g.declare_reviewed_scope(('rna','hydrazine','hydrazine_program','hydrazinolysis','glycoblue','precipitation_program','precipitate_hydrazine','aniline','aniline_program','aniline_scission','precipitate_aniline','pnk','pnk_buffer','beta_me','dephosphorylate','extraction_supply','recover_ends','kit_supply','kit_conversion','select_qc','library_qc','select_contribution','library_remainder','other_libraries','multiplex','select_loading','pool_remainder','loading_conditions','flowcell_prep','sequencing_input'),SOURCE)
    return g.require_publishable()


def render_page():
    from page import head,info
    from chemdraw import panel
    from protocol_views import assemble_page,WorkflowSection
    return assemble_page(protocol_graph(),(head('HydraPsiSeq — pseudouridine protection from cleavage'),'<div class="wrap" data-protocol-complete="true"><h1>HydraPsiSeq</h1>',
      '<p class="research-notes"><a href="01_hydrapsiseq.html">Research notes</a></p>',info('<a href="https://doi.org/10.1093/nar/gkaa769">Published HydraPsiSeq Methods and SupplementaryFigureS1</a>'),
      '<h2>Read chemical cleavage protection</h2>',
      '<p>Hydrazine acts on susceptible U residues; subsequent aniline treatment breaks the RNA chain at abasic sites. Pseudouridine is relatively protected. The selected cleavage event removes the central U and produces an upstream3′phosphate and downstream5′phosphate. ATP-freePNK removes the3′phosphate while retaining the downstream5′phosphate.</p>',
      panel(CleavageEvent().scene().rows(),caption='One source-described U-cleavage event, not a deterministic RNA population or a Ψ cleavage prediction. Genomic/transcript sequence and conversion fractions remain unassigned.'),
      '<h2>Keep molecular debris until physical recovery</h2>',
      '<p>The source uses50–300ngRNA,50%hydrazine on ice30–60min, NaOAc/Glycoblue ethanol precipitation, then1Maniline at60°C15min. The second precipitation retains RNA fragments. T4PNK is10U,37°C6h, followed by65°C20min inactivation and phenol/chloroform/ethanol recovery. Released products stay molecular material until actual separations; exact debris/carrier partition and recovery are not predicted.</p>',
      '<h2>Actual library finishing with unknown kit geometry</h2>',
      '<p>PublishedMethods explicitly converts recovered RNA to a DNA library with NEBNextSmallRNAkit, then assesses it by BioAnalyzerHighSensitivityDNA and Qubit. The graph represents that actual finishing interface; it does not attach sequencing directly to the RNA intermediate. Historicalkit edition, adapter/RT/PCR/workups and whole DNA/read-primer geometry remain unassembled. The reagentcollection is not an invented premix or imported newer kitrecipe.</p>',
      '<h2>Multiplexing, read setup and limits</h2>',
      '<p>Source libraries are multiplexed and sequenced onHiSeq1000 with50bp single-read runs at10–12pM. Analytical material, librarycontributions and loadingaliquots leave explicit retained/untransferred stocks. Pool membership, indexmask and run-primer sites remain unknown. Read5′ends are used for protection profiling; calibrated Ψ levels, exact site identity and cleavage/background efficiency are not established by the local event model.</p>',
      WorkflowSection('Full selected RNA-to-library preparation and molecular limits'),'</div>'))

def selftest():
    from checks import Check
    c=Check();event=CleavageEvent();repaired=event.dephosphorylate()
    c('ATP-free repair preserves5p while removing3p',
      (repaired.upstream_three_prime,repaired.downstream_five_prime),('OH','phosphate'))
    try: CleavageEvent(removed_residue='Ψ')
    except ValueError: rejected=True
    else: rejected=False
    c('U event is not silently reused as a Ψ cleavage prediction',rejected,True)
    g=protocol_graph()
    c('source actualkitconversion supports qualified publication',g.publication_problems(),[])
    c('the source kit produces DNA-library material',g.materials['library'].kind,MaterialKind.LIBRARY)
    c('sequencing uses selected library, not repairedRNA',tuple(g.sequencing),('sequencing_input',))
    c('source carrier is molecular material at precipitation',any(p.material=='glycoblue' and p.role is InputRole.STRUCTURAL for p in g.steps['precipitate_hydrazine'].inputs),True)
    c('analytical material terminates separately',g.terminals['library_qc'].disposition,Disposition.QC)
    c('multiplexpool remainder survives loading',g.terminals['pool_remainder'].disposition,Disposition.UNRESOLVED)
    from checks import Source,have
    from storage import sources_root
    source=Source(sources_root()/Path(__file__).resolve().parents[1].name/'nar_gkaa769_PMC7641733.xml.txt','PublishedHydraPsiSeqMethods')
    if have(c,source,label='source actuallibrary/readout'):
        text=' '.join(source.path.read_text().split())
        c('source explicitly converts RNA to library','RNA was converted to library using NEBNext Small RNA Library kit' in text,True)
        c('source single-read50 is actualsequencing','50 bp single-read runs' in text,True)
        c('source DNAchip supports DNA-library product','High Sensitivity DNA chip' in text,True)
    c.report()

if __name__=='__main__': selftest()
