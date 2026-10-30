"""CoLa-seq accepted-version preparation; real molecules stay source-qualified."""
from pathlib import Path
from dataclasses import replace
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Scene, Segment, feature, revcomp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S,
                            Disposition as D, Terminal, SourceItem, SequencingPlan)
from rna_lariat import RNALariatTopology

AUTHORITY = ('Accepted author manuscript NIHMS1854138/PMC10448999, DOI10.1016/j.molcel.2022.11.004; '
             'introduction72/workflow90–94/STARMethods216–230. Final publication equivalence unverified. '
             'TableS4/S5 identities cited, actual oligo documents unresolved; no preprint sequence inheritance.')
UMI_LENGTH = 6
READ_LENGTH = 100


def adapter_scene():
    sc = Scene()
    sc.strand('adapter', [Segment('source UMI', 'N' * UMI_LENGTH, placeholder=True,
                                 feature=feature('cola.umi', 'umi', 'random', note='Source six-base5prime UMI; whole oligo unassembled.')),
                          Segment('source adapter body', '[TableS5 body unassembled]', placeholder=True)],
              mod5='App')
    sc.ligation_end('adapter', 'source UMI', kind='single-strand', polarity='5-prime', owner='adapter',
                    text='Source preadenylated donor; whole body/3prime chemistry unassigned')
    return sc


def lariat_example():
    """Invented sequence/coordinate for topology illustration, not a K562 molecule."""
    return RNALariatTopology('GUGCGAACUUAGU', 6, 'invented RNA example')


def cola_tail_copy(reference):
    """Source-specific copied interval only; excludes unavailable RT-primer bodies.

    This pre-adapter example copies only natural RNA tail; it excludes adapter/UMI
    and unavailable RT-primer tracts, not a complete cDNA. The accepted source's
    observed read starts are BP+1. This is not a generic
    prediction that every reverse transcriptase stops at every RNA branch.
    """
    reference.validate_connectivity()
    return revcomp(reference.tail_sequence)


def lariat_example_graph():
    g = ProtocolGraph('Illustrative RNA lariat connectivity', representation='illustrative_example',
                      example_description='RNA sequence and branch index are invented; only topology is illustrated.')
    rna = g.input('example_rna', 'Invented reference, not an experimental CoLa substrate', K.MOLECULAR,
                  primary=True, chemical_state=lariat_example())
    g.mark_terminal(Terminal(rna.id, D.UNRESOLVED, 'Topology component only, not a complete library.'))
    return g


class Prep:
    def __init__(self):
        self.g = ProtocolGraph('CoLa-seq accepted-manuscript K562 preparation', authority=AUTHORITY,
                               representation='published_instance')
        self.x = self.g.input('cells', 'Source K562 library cohort: RNA from200million cells', K.BIOLOGICAL,
                              primary=True, chemical_state={'reference_batch': 'Methods fractionation recipe uses50million cells; actual array multiplicity/scaling unassigned.',
                                                           'scope': 'Aggregate source cohort inventory, not four invented reactions or one unscaled reagent dose.'})
    def supply(self, key, label, kind=K.CATALYST, role=R.CATALYST, **kw):
        return self.g.input(key, label, kind, **kw), role
    def condition(self, key, label):
        return self.supply(key, label, K.CONDITION, R.CONDITION)
    def op(self, key, label, inputs=(), *, kind=S.MAGIC, geometry=None):
        y = self.g.product(key+'_product', label+'; occupied molecules/recovery unassembled', K.MOLECULAR)
        kw = {'reason': 'Positive accepted-source operation; complete molecular input/output classes, historical formulations and yields unassembled.'} if kind==S.MAGIC else {}
        if kind==S.GENERIC: kw['changes_molecule']=False
        if geometry is not None: kw['geometry']=geometry
        self.g.react(key, label, kind, ((self.x, R.TEMPLATE),)+tuple(inputs), (y,), **kw)
        self.x = y
    def analytical(self, key, label):
        sample = self.g.product(key+'_sample', label, K.MOLECULAR)
        rest = self.g.product(key+'_rest', 'Actual stock after analytical allocation; amount unassigned', K.MOLECULAR)
        self.g.react(key, label, S.SPLIT, ((self.x, R.TEMPLATE),), (sample,rest))
        self.g.mark_terminal(Terminal(sample.id,D.QC,'Analytical fraction excluded from production.'))
        self.x=rest
    def select(self,key,label,remainder):
        sample=self.g.product(key+'_selected',label,K.MOLECULAR)
        rest=self.g.product(key+'_rest',remainder,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(sample,rest))
        self.g.mark_terminal(Terminal(rest.id,D.UNRESOLVED,remainder));self.x=sample
    def workup(self,key,label,inputs=(),*,removed=None):
        kept=self.g.product(key+'_kept',label+': actual retained fraction; no perfect specificity claim',K.MOLECULAR)
        other=self.g.product(key+'_other',removed or 'Actual removed workup fractions',K.MOLECULAR)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.condition(key+'_conditions',label))+tuple(inputs),(kept,other))
        if removed:
            qc=self.g.product(key+'_other_qc','Source roughly5percent removed biological fraction for QC',K.MOLECULAR)
            rest=self.g.product(key+'_other_stock','Other biological stock outside selected CoLa production',K.MOLECULAR)
            self.g.react(key+'_other_qc_allocate','Actual analytical allocation from removed subcellular fraction',S.SPLIT,((other,R.TEMPLATE),),(qc,rest))
            self.g.mark_terminal(Terminal(qc.id,D.QC,'Source western/RNA-fragment QC.'))
            self.g.mark_terminal(Terminal(rest.id,D.UNRESOLVED,'Outside selected chromatin RNA production; not claimed destroyed.'))
        else:self.g.discard(other,'Removed source workup fraction; composition/recovery unassigned.')
        self.x=kept


def protocol_graph():
    p=Prep();g=p.g
    p.workup('wash','Source100g3min harvest/two15mlPBS andone1mlPBS washes per50M reference',
             (p.supply('pbs','Actual aggregate ice-coldPBS wash stocks; aggregate volume/scaling unassigned',K.MOLECULAR,R.STRUCTURAL),))
    p.op('lysis','BufferA Triton lysis9min ice',
         (p.supply('buffer_a','Actual bufferA aggregate allocation; source1ml per50M reference',K.MOLECULAR,R.STRUCTURAL,
                   chemical_state={'reference_composition':'0.1percentTriton;10mMHEPESpH7.5/10mMKCl/10percentglycerol/0.116g-mlsucrose/1mMDTT/1xproteaseinhibitor/2mMEDTA; no addedMgCl2'}),))
    p.workup('nuclei','Source800g5min4C retain nuclei',removed='Actual removed cytoplasmic fraction')
    p.op('nrb_nun','Source nuclei extraction NRB then slowly addedNUN,5min ice',
         (p.supply('nrb','Actual aggregateNRB allocation;250ul per50M reference, not an unscaled one-dose claim',K.MOLECULAR,R.STRUCTURAL),
          p.supply('nun','Actual aggregateNUN allocation;250ul per50M reference',K.MOLECULAR,R.STRUCTURAL),
          p.condition('buffer_formulas','NRB20mMHEPESpH7.5/75mMNaCl/50percentglycerol/0.5mMEDTA/1mMDTT/PI;NUN50mMHEPESpH7.5/300mMNaCl/0.2mMEDTA/1Murea/1mMDTT/1percentNP40')))
    p.workup('chromatin','Source800g5min4C retain chromatin',removed='Actual removed nucleoplasmic fraction')
    p.workup('chromatin_wash','Source500ul1:1NRB/NUN wash per50M reference',
             (p.supply('wash_premix','Actual aggregate premixedNRB/NUN wash stock; no extra separate components',K.MOLECULAR,R.STRUCTURAL),))
    p.analytical('chromatin_qc','Source roughly5percent chromatin fraction for western/RNA QC; cohort count convention unassigned')
    p.op('resuspend','Chromatin pellet resuspension inNRB',
         (p.supply('pellet_nrb','Actual NRB stock; acceptedJATSprints200mL, intended unit/aggregate amount unassigned',K.MOLECULAR,R.STRUCTURAL),))
    p.workup('trizol','SourceTrizol extraction followed by phenol/chloroform/isopropanol recovery',
             (p.supply('extraction_reagents','Actual source extraction/precipitation reagent allocations; volumes unassigned',K.MOLECULAR,R.STRUCTURAL),))
    p.op('genomic_dnase','SourceTurboDNase removal of residual genomicDNA',
         (p.supply('dnase_genome','Actual firstTurboDNase allocation; historical formulation unassembled'),))
    p.workup('dnase_workup','Source phenol/chloroform extraction after genomicDNase')
    p.workup('size1','SourceSPRI1:1.2 beads:sample removes short RNA roughly200nt; some valid short lariats lost',
             (p.supply('spri1','Actual firstSPRI bead allocation',K.SELECTION,R.SELECTION),))
    rz=g.input('ribozero_kit','Actual historicalRiboZeroMagneticGold inventory; packaging/formulation unassembled',K.MOLECULAR)
    rz_probes=g.product('ribozero_probes','Known oligo-based depletion probe class from SAME actual kit; whole orders/polymer chemistry unassembled',K.OLIGO)
    rz_other=g.product('ribozero_other','Remaining same-kit depletion supplies; selection chemistry/quantities unassembled',K.MOLECULAR)
    g.react('ribozero_inventory','Copresent oligo and other kit supply classes; no new probe dose or physical separation',S.GENERIC,((rz,R.TEMPLATE),),(rz_probes,rz_other),changes_molecule=False)
    p.op('ribozero','SourceRiboZeroMagneticGold rRNA depletion',
         ((rz_probes,R.STRUCTURAL),(rz_other,R.STRUCTURAL)))
    p.op('aso_anneal','Source abundant-ncRNA probe annealing75C5min→22C0.1C-sec',
         (p.supply('asos','ActualTableS4 biotinylated depletion oligo mixture;82.5pmol total perugRNA, whole orders unassembled',K.OLIGO,R.STRUCTURAL),
          p.condition('aso_solution','Actual2xSSC/5mMEDTA; total reaction volume unassigned')))
    p.workup('c1','SourceC1 capture: retain RNA SUPERNATANT, remove probe-bound ncRNA',
             (p.supply('c1_beads','ActualMyOneStreptavidinC1 bead allocation; amount unassigned',K.SELECTION,R.SELECTION),))
    p.workup('rna_clean','SourceRNA Clean and Concentrator cleanup')
    p.op('probe_dnase','SeparateTurboDNase addition removes excess antisenseDNA',
         (p.supply('dnase_probes','Actual secondTurboDNase allocation after probe depletion'),))
    p.workup('probe_cleanup','Source purification before enzyme enrichment; exact recovery unassigned')
    p.op('rpph','SourceRppH treatment described as linearRNA decapping',
         (p.supply('rpph_stock','ActualRppH NEB allocation; doses/conditions delegated and unassembled'),))
    p.op('terminator','SourceTerminator5prime→3prime exonuclease treatment',
         (p.supply('terminator_stock','ActualTerminator Lucigen allocation; not RNaseR'),))
    p.workup('size2','SecondSPRI1:1.2 selection removes≤200nt; retained mixture includes lariat classes, not pure/exclusive lariats',
             (p.supply('spri2','Actual secondSPRI bead allocation',K.SELECTION,R.SELECTION),))
    p.x=replace(p.x,chemical_state={'known_target_topology':'RNA lariat with2prime5prime branch and3prime tail; actual RNA sequences and branch coordinates unassigned', 'specificity':'Enrichment mixture is not exclusively lariats; background classes remain unassembled.'})
    g.materials[p.x.id]=p.x
    p.select('library_rna','Actual approximately500ng enrichedRNA library allocation','Untransferred enrichedRNA stock; source recovery0.5–1percent is empirical, not imposed yield')
    geom=AdapterLigationGeometry.known('single-strand',adapter_scene(),statement='Source known5prime preadenylation/6baseUMI; wholeTableS5 body and3prime chemistry unassembled.')
    p.op('adapter_ligate','Source3prime adapter ligation22C3h',
         (p.supply('adapter','SourceTableS5 preadenylated3end adapter; whole order unknown',K.OLIGO,R.STRUCTURAL,state=adapter_scene()),
          p.supply('rnl2','ActualtruncatedK227Q T4RNA ligase2; source10U-ul'),
          p.condition('ligation_conditions','Source17.5percentPEG8000/1xligasebuffer/1uMadapter; total reaction volume unassigned')),
         kind=S.LIGATION,geometry=geom)
    p.workup('adapter_cleanup','SourceSPRI1.6:1 beads:sample removes excess adapters',
             (p.supply('spri3','Actual adapter-cleanup bead allocation',K.SELECTION,R.SELECTION),))
    rt_inventory=g.input('rt_supplies','Actual RT functional inventory; stock count/packaging unassigned, not a claimed premix',K.MOLECULAR)
    rt_enzyme=g.product('ssiii','SSIII component of actualRT inventory; source10U-ul',K.CATALYST)
    rt_nt=g.product('rt_nt','Nucleotide component of actualRT inventory; identities/quantities unassembled',K.MOLECULAR)
    rt_buffer=g.product('rt_buffer','Buffer component of actualRT inventory; historical formulation unassembled',K.CONDITION)
    g.react('rt_inventory','Copresent RT supply classes, not physical separation or fresh duplicate additions',S.GENERIC,((rt_inventory,R.TEMPLATE),),(rt_enzyme,rt_nt,rt_buffer),changes_molecule=False)
    p.op('rt','SourceSSIII reverse transcription50C1h; branchpoint stop is a local source readout claim',
         (p.supply('rt_primer','ActualTableS5 RT primer0.25uM; whole order unassembled',K.OLIGO,R.STRUCTURAL),
          (rt_enzyme,R.CATALYST),(rt_nt,R.STRUCTURAL),(rt_buffer,R.CONDITION)))
    p.op('exosap','Actual3.5ulExoSAP-IT addition,37C30min',
         (p.supply('exosap_stock','Actual3.5ulExoSAP-IT preparation; copresent components unassembled'),))
    p.op('stop','Actual10ul0.5MEDTA then10ul1NNaOH additions stop/degradeRNA',
         (p.supply('edta','Actual10ul0.5MEDTA',K.MOLECULAR,R.STRUCTURAL),p.supply('naoh','Actual10ul1NNaOH',K.MOLECULAR,R.STRUCTURAL)))
    p.workup('cdna_clean','SourceOligo Clean and Concentrator cleanup')
    p.op('circle','SourceCircLigase cDNA circularization',
         (p.supply('circligase','ActualCircLigase allocation; historical conditions/whole cDNA ends unassembled'),))
    primers=g.input('pcr_supplies','Actual PCR primer inventory; two functional classes, stock count/packaging unassigned',K.MOLECULAR)
    fw=g.product('pcr_forward','ActualTableS5 forward primer class, whole order unassembled',K.OLIGO)
    rv=g.product('pcr_reverse','ActualTableS5 reverse-index primer class, whole order/assignedindex unassembled',K.OLIGO,
                 chemical_state={'unpositioned_identifier':feature('cola.sample','sample_index','unknown',note='Source index is in reversePCRprimers; width, position and assignment unassembled.')})
    g.react('primer_classes','Copresent functional primer inventory, not an invented premix or physical separation',S.GENERIC,((primers,R.TEMPLATE),),(fw,rv),changes_molecule=False)
    pcr_inventory=g.input('pcr_reagents','Actual PCR functional inventory; stock count/packaging unassigned',K.MOLECULAR)
    pol=g.product('phusion','Phusion enzyme component from actualPCR supplies',K.CATALYST)
    pcr_nt=g.product('pcr_nt','Nucleotide component from actualPCR supplies; quantities unassembled',K.MOLECULAR)
    pcr_buffer=g.product('pcr_buffer','Buffer component from actualPCR supplies; conditions unassembled',K.CONDITION)
    g.react('pcr_inventory','Copresent PCR supply classes; no invented premix or duplicate additions',S.GENERIC,((pcr_inventory,R.TEMPLATE),),(pol,pcr_nt,pcr_buffer),changes_molecule=False)
    template_view=Scene()
    template_view.strand('template',[Segment('source circular cDNA','[whole circular cDNA unassembled]',placeholder=True)])
    template_view.note('template','Placeholder reference only; no actual linearization, circle junction or sequence inferred.')
    p.op('pcr','SourcePhusionPCR15–18cycles; reverse primers introduce sample indices',
         ((fw,R.STRUCTURAL),(rv,R.STRUCTURAL),(pol,R.CATALYST),(pcr_nt,R.STRUCTURAL),(pcr_buffer,R.CONDITION)),
         kind=S.PCR,geometry=PCRGeometry.unavailable(template_view,reason='Whole circular template and TableS5 primers are unassembled; source-positive PCR does not establish binding geometry.'))
    p.select('gel_load','Actual selectedPCR gel-loading fraction, amount unassigned','Unloaded PCR stock')
    p.workup('gel','Source6percentTBE gel selects175–750bp DNA')
    p.select('load','Actual3nM library loading allocation; volume unassigned','Untransferred final library stock')
    lib=g.product('library','SourceCoLa-seq library; whole occupied molecules and full primer landing unassembled',K.LIBRARY)
    g.react('sequence','AcceptedsourceHiSeq4000 paired-end100',S.GENERIC,((p.x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(read_lengths={'R1':READ_LENGTH,'R2':READ_LENGTH},limitation='Accepted-version sourcePE100. Full TableS5 oligos, index assignment/acquisition and whole products unassembled.'))
    g.add_source_item(SourceItem('accepted_main',AUTHORITY,'Actual source-positive preparation, not proof of final-version equivalence.'))
    g.add_source_item(SourceItem('lariat_topology','Acceptedintroduction72/workflow90–94','Known2prime5prime branch topology; actual genomic sequences/branch coordinates unknown. Invented example is separate and not a reaction substrate.'))
    return g


def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    from chemistry_validation import check_graph,Status
    c=Check();g=protocol_graph();r=lariat_example()
    c('Local source-readout interval excludes branch, not ordinaryRNA circular traversal',cola_tail_copy(r),'ACTAAG')
    c('Invented RNA is not a real preparation input','example_rna' in g.materials,False)
    c('Topology reference remains labelled as illustrative',lariat_example_graph().representation,'illustrative_example')
    c('Topology checker does not establish any CoLa reaction',[x.status for x in check_graph(lariat_example_graph())],[Status.ESTABLISHED])
    c('Terminator is its own named enzyme, not RNaseR',g.steps['terminator'].inputs[1].material,'terminator_stock')
    c('Actual input is selected after two enrichment selections',g.steps['library_rna'].inputs[0].material,'size2_kept')
    c('C1 actual retained fraction has a separate removed output',tuple(g.steps['c1'].outputs),('c1_kept','c1_other'))
    c('ExoSAP addition acts on the actualRT mixture',g.steps['exosap'].inputs[0].material,'rt_product')
    c('PCR consumes both functional primer classes',tuple(i.material for i in g.steps['pcr'].inputs[1:3]),('pcr_forward','pcr_reverse'))
    c('Source PCR is typed without invented primer landing',g.steps['pcr'].kind,S.PCR)
    c('Known depletion probes belong to same actual RiboZero inventory',g.steps['ribozero_inventory'].inputs[0].material,'ribozero_kit')
    c('RNA depletion consumes actual kit probe class',g.steps['ribozero'].inputs[1].material,'ribozero_probes')
    c('Loaded library comes only from selected final allocation',g.steps['sequence'].inputs[0].material,'load_selected')
    c('Adapter N6 means UMI',adapter_scene().strands['adapter'].segs[0].feature.role,'umi')
    c('Whole experimental library stays unassembled',g.materials['library'].state,None)
    c('Partial library is not a ready schematic',bool(g.publication_problems()),True)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'j.molcel.2022.11.004_PMC10448999.xml.txt'
    if have(c,Source(p,'https://pmc.ncbi.nlm.nih.gov/articles/PMC10448999/')):
        text=p.read_text()
        for label,phrase in [('Source2prime5prime','2′−5′ phosphodiester linkage'),('SourceRNA500ng','~500 ng of RNA'),('SourceUMI6','6 nt UMI'),('Sourceacceptedunit','200 mL of NRB'),('SourceTerminator','Terminator 5′−3′ exonuclease'),('SourcePE100','100 bp paired-end mode')]:c(label,phrase in text,True)
    c.report()

if __name__=='__main__':
    out=Path(__file__).resolve().parents[1]/'lariat_example.svg'
    out.write_text(lariat_example().svg(caption='Illustrative lariat: invented RNA sequence'))
    selftest()
