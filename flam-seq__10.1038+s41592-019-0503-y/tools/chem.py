"""Protocol-owned source-qualified partial chemistry. No public schematic endpoint."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Scene, Segment, revcomp
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan, Terminal, Disposition
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry

def scene(name, text):
    sc = Scene(); sc.strand(name, [Segment(name, text, placeholder=True)]); return sc

def oligo(g, key, name, bases, **mods):
    parts = []
    for i, block in enumerate(bases.split('[BtndT]')):
        if i: parts.append(Segment('Internal biotin-dT', 'T', note='Source [BtndT]; physical modification not inherited by copies'))
        if block: parts.append(Segment(key + str(i), block, placeholder=bool(set(block)-set('ACGTU')), note='IUPAC synthesis pattern retained; full landing unassembled'))
    sc = Scene(); sc.strand(key, parts, **mods)
    return g.input(key, name, K.OLIGO, state=sc)

def op(g, key, x, label, kind=S.MAGIC, *, reagent='Source mixture; internal chemistry unassembled', product=None, structural=()):
    catalyst = kind in {S.DIGESTION, S.END_REPAIR, S.EXTENSION, S.TAILING, S.REVERSE_TRANSCRIPTION, S.PCR, S.LIGATION}
    stock = g.input(key + '_reagent', reagent, K.CATALYST if catalyst else K.CONDITION)
    inputs = [(x, R.TEMPLATE), (stock, R.CATALYST if catalyst else R.CONDITION)]
    inputs.extend((a, R.STRUCTURAL) for a in structural)
    kwargs = {}
    if kind == S.MAGIC:
        kwargs['reason'] = 'Source operation established; exact molecular transition/kit internals unassembled. No accepted absence claim.'
    if kind == S.PCR:
        kwargs['geometry'] = PCRGeometry.unavailable(Scene(), reason='Exact template/primer endpoint geometry unassembled; selected oligos do not establish whole amplicon.')
        if len(structural) < 2:
            for arm in ('f', 'r'):
                inputs.append((g.input(key + '_' + arm, 'Source PCR primer; exact identity/landing unassembled', K.OLIGO, state=scene(key + '_' + arm, '[source primer]')), R.STRUCTURAL))
    if kind == S.REVERSE_TRANSCRIPTION and not structural:
        inputs.append((g.input(key + '_primer', 'Source RT primer; geometry unassembled', K.OLIGO, state=scene(key + '_primer', '[source RT primer]')), R.STRUCTURAL))
    if kind == S.LIGATION:
        kwargs['geometry'] = AdapterLigationGeometry('unavailable', 'unavailable', (scene(key + '_ends', '[source ligation; complete substrate/end geometry unassembled]'),), 'This incomplete model does not claim publication absence.')
        if not structural:
            inputs.append((g.input(key + '_adapter', 'Source ligating adapter; full geometry unassembled', K.ADAPTER, state=scene(key + '_adapter', '[source adapter]')), R.STRUCTURAL))
    y = g.product(key + '_product', product or label + '; physical companions retained until separation', K.MOLECULAR)
    g.react(key, label, kind, inputs, (y,), **kwargs); return y

def clean(g, key, x, label):
    condition = g.input(key + '_condition', label, K.CONDITION)
    kept = g.product(key + '_kept', 'Source physically recovered fraction; carryover/yield unquantified', K.MOLECULAR)
    waste = g.product(key + '_waste', 'Actual separated cleanup fractions; no perfect classification assumed', K.WASTE)
    g.react(key, label, S.CLEANUP, ((x, R.TEMPLATE), (condition, R.CONDITION)), (kept, waste))
    g.discard(waste, 'Actual physical separation rather than computational removal.'); return kept

def qc(g, key, x, label):
    kept = g.product(key + '_kept', 'Source retained preparation', K.MOLECULAR)
    assay = g.product(key + '_aliquot', 'Separate analytical aliquot; no return to library assumed', K.QC)
    g.react(key, label, S.SPLIT, ((x, R.TEMPLATE),), (kept, assay))
    g.mark_terminal(Terminal(assay.id, Disposition.QC, 'Source analytical aliquot.')); return kept

def finish(g, x, label, limitation):
    lib = g.product('library', label, K.LIBRARY)
    g.react('library_identity', 'Physical library interface; complete molecular state unassembled', S.GENERIC, ((x, R.TEMPLATE),), (lib,), changes_molecule=False)
    g.sequence(lib, SequencingPlan(limitation=limitation))
    g.add_source_item(SourceItem('pending', 'Historical exact kit, substrate and sequencing geometry', 'Private partial graph; no full source audit or accepted absence.'))
    return g

SOURCE_OLIGO_EXCEPTIONS = {'iso_tso': 'Independent DNA implementation reader verified published MOESM3 original p2 literal iCiGiCAAGCAGTGGTATCAACGCAGAGTACATrGrGrG on2026-10-29. Guard macro-normalization cannot match physical noncanonical iC/iG/iC placeholders; exact DNA body ACAT and terminal RNA retained. Narrow input-only text-match exception, no ordinary CGC substitution.'}
SOURCE = 'Published2019 supplementary protocol MOESM3 all preparative/setup sections; published SI MOESM1 pp1–2 selected tail-anchoring limits. Main final text inaccessible; no preprint recipe.'
RT1 = 'GGTAATACGACTCACTATAGCGAGA' + 'N'*10 + 'C'*9 + 'TTT'
RT2 = 'TGAGTCGGCAGAGAACTGGCGAA' + 'N'*10 + 'C'*9 + 'TTT'
PCR1 = 'GGTAATACGACTCACTATAGCGAG'
PCR2 = 'TGAGTCGGCAGAGAACTGGCGAA'
from rt import SMART_HANDLE
TSO_BODY = SMART_HANDLE + 'ACAT'
TSO_WRITTEN = 'iCiGiC' + TSO_BODY + 'rGrGrG'

@dataclass(frozen=True)
class AnchoredCopy:
    body: str
    tail_length: int
    def __post_init__(self):
        if not self.body or set(self.body)-set('ACGU') or type(self.tail_length) is not int or self.tail_length < 3:
            raise ValueError('Explicit selected RNA body and at least three terminal adenosines required')
    def first_strand(self, primer):
        if not primer or set(primer)-set('ACGT') or not primer.endswith('C'*9+'TTT'):
            raise ValueError('Supply an actual resolved source-pattern DNA primer register')
        return primer + 'T'*(self.tail_length-3) + revcomp(self.body.replace('U','T'))

def tso_scene():
    sc=Scene();sc.strand('tso',[Segment('5prime iso-dC iso-dG iso-dC','[iso-dC][iso-dG][iso-dC]',placeholder=True,note='Physical noncanonical bases, not ordinary CGC'),Segment('TSO DNA body',TSO_BODY),Segment('terminal RNA GGG','GGG',note='Physical riboguanosines rGrGrG')]);return sc

def protocol_graph():
    g=ProtocolGraph('FLAM-seq — published2019 supplementary preparation',authority=SOURCE,representation='published_instance')
    x=g.input('total_rna','Source totalRNA2–10ug; extraction/application branch unassembled',K.MOLECULAR,primary=True)
    beads=g.input('purification_beads','SourceTruSeqRNApurificationbeads',K.SELECTION)
    y=g.product('first_bound','Source bead-bound RNA fraction with carryover',K.MOLECULAR);w=g.product('first_supernatant','Physically discarded first capture supernatant',K.WASTE)
    g.react('capture1','Source65C/ice/roomtemperature then magnetic capture',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(y,w));g.discard(w,'Protocolexplicitdiscard')
    x=clean(g,'wash1',y,'Sourcebeadwash and supernatantdiscard')
    x=op(g,'elute1',x,'Source80C2min elution before rebinding',product='SourceRNA/elutionbufferandbeads together; no magneticdiscardatthisstep')
    x=op(g,'rebind',x,'Source bindingbuffer/roomtemperature5min rebinding',product='ReboundRNA/beads plus unboundRNAandbuffer until magnetic separation')
    x=clean(g,'capture2wash',x,'Source magnetic supernatantdiscard and second beadwash')
    x=op(g,'elute2',x,'Sourcewater70C2min/ice1min elution',product='ElutedRNAplusbeadsbeforemagneticcollection')
    x=clean(g,'rna_elute',x,'Source magnetically collect16ul RNA supernatant; bead fraction separated')
    x=qc(g,'rna_qc',x,'Source optional1ul Bioanalyzer aliquot; allocated RNA scope')
    x=op(g,'gi_tail',x,'Source GI tailing37C60min followed by stopsolution',S.TAILING,reagent='USBpolyAtail-lengthassay5xtailbuffer/10xtailenzymemix; exactenzyme/nucleotidecompositionunreviewed',product='SourceGI-taggedRNA classes plus original/freeRNA and kit/remnant reagents; no knownGIlength/composition asserted')
    x=clean(g,'gi_clean',x,'Source1.8xRNAbeads/80percentEtOHtwice/elutesupernatant')
    rt=oligo(g,'rt1','PublishedRTprimer1; N10function unassigned in selected model',RT1)
    tso=g.input('iso_tso','PublishedisoTSO with noncanonical5primebases and terminalRNA',K.OLIGO,state=tso_scene())
    x=op(g,'rna_primer_preheat',x,'SourceRNAplusRTprimer72C3minthen coolto42C; noRTenzyme/TSOyet',reagent='SourceRNA/RTprimerpreheat',structural=(rt,),product='PreheatedRNA/RTprimerpreparation; noenzyme/TSOexposedto72C')
    x=op(g,'rt',x,'At42CaddRTmastermix/isoTSO toalreadyprimedRNA; nosecondRTprimeraddition',S.REVERSE_TRANSCRIPTION,reagent='SMARTScribe/dNTPs/DTT/RNaseinhibitor; sourceRTmastermixaddedafterRNA/primerpreheatand42Ccooling',structural=(tso,),product='Selectedfull-lengthfirst-strandDNAclassesplusoriginalRNA,unextendedRTprimer/TSOandkitremnants; fullTSincorporationgeometryunassembled')
    x=clean(g,'rt_clean',x,'Source0.6xDNAbeads/80percentEtOHtwice/elutesupernatant')
    x=qc(g,'cdna_qc',x,'Source optional1ul Bioanalyzer aliquot')
    p=oligo(g,'pcr1','PublishedPCRprimer1 pairedwithRTprimer1',PCR1)
    q=g.input('pcr_iia','Source5primePCRPrimerIIA; exactbody/landingnotassignedfromTSObody',K.OLIGO,state=scene('pcr_iia','[named source PCR primer]'))
    x=op(g,'pcr',x,'SourceAdvantage2PCR; optimize18/20/22/24cycles if new sample',S.PCR,reagent='SourceAdvantage2SA/polymerasemix/dNTPs',structural=(p,q),product='AmplifieddsDNAclassesplusoriginalRNA/DNA/primers/TSO/reagents; noRNAhydrolysisstepinvented')
    x=clean(g,'pcr_clean',x,'Source0.6xDNAbeads/80percentEtOHtwice/elutesupernatant')
    x=qc(g,'library_qc',x,'Source1ulFragmentAnalyzer and Qubit aliquots')
    x=op(g,'damage_repair',x,'SourceSMRTbell1.0SPv3DNArepair37C60min',S.END_REPAIR,reagent='NamedDNArepairmix/buffer/NAD+/ATP/dNTPs',product='Repair-treateddsDNAandoriginalmolecular/reagentcompanions')
    x=op(g,'end_repair',x,'Sourceend-repairmix25C5min',S.END_REPAIR,reagent='HistoricalSMRTbellEndRepairMix')
    x=clean(g,'repair_clean',x,'Source0.6xAMPure/80percentEtOHtwice/EBsupernatant')
    x=op(g,'smrtbell_ligation',x,'Sourceannealedbluntadapterligation25C15min/65C10min',S.LIGATION,reagent='Sourceannealedbluntadapter/TemplatePrepBuffer/ATPlow/ligase',product='SMRTbellclassesplusunligated/partlyligateddsDNA/freeadapters/remnants; exactadaptersequenceunassembled')
    x=op(g,'exo',x,'SourceExoIIIandVII37C60min',S.DIGESTION,reagent='ExonucleasesIII/VII',product='Sourceexonuclease-treatedSMRTbellpreparationplusdegradedDNAfragments/freeadapters/reagents; noinstantdiscard')
    x=clean(g,'smrt_clean1',x,'Source0.6xAMPure/80percentEtOHtwice/EBsupernatant')
    x=clean(g,'smrt_clean2',x,'Source repeatedAMPurecleanup/elute11ulEB')
    x=qc(g,'smrt_qc',x,'SourceBioanalyzer/Qubitanalyticalaliquots')
    seq_primer=g.input('sequencing_primer_stock','HistoricalSequencingPrimerv3 stock; exactsequenceunassembled',K.OLIGO,state=scene('seq_primer','[source sequencingprimer]'))
    conditioned=g.product('conditioned_sequencing_primer','SourceSequencingPrimerv3conditionedALONE80C2minthen4C',K.OLIGO,state=scene('conditioned_seq','[conditionedsourceprimer]'))
    primer_buffer=g.input('primer_condition_buffer','Source1xelutionbuffer',K.CONDITION)
    g.react('condition_sequencing_primer','ConditionSequencingPrimerv3stockonly80C2min; noDNAlibrarypresent',S.MAGIC,((seq_primer,R.TEMPLATE),(primer_buffer,R.CONDITION)),(conditioned,),reason='Sourceprimerconditioning; exactmolecularstateunassembled')
    x=op(g,'primer_bind',x,'Mixsamplewithconditionedv3primer/PrimerBufferv2;20C60minanneal',reagent='Sourceprimerannealingbuffer',structural=(conditioned,))
    x=op(g,'polymerase_bind',x,'SourceSequelPolymerase3.0binding30C60min',reagent='SequelPolymerase3.0/Sequelbindingbuffer/dNTP/DTT',product='Polymerase-boundtemplateclassesplusunboundpolymerase/primersandkitremnants')
    x=clean(g,'complex_clean',x,'SourceMagbeadbindingbuffer/AMPurecleanup; NOethanolwash; collectsupernatantafterresuspension')
    x=qc(g,'recovered_qc',x,'SourceQubitrecoveredDNAaliquotbeforeloadingmix')
    control=g.input('internal_control','SourceDILUTEDINTERNALCONTROL stock; actualidentity/amountunassigned',K.MOLECULAR)
    load_mix=g.input('load_mix','SourceDTT/SequelAdditive/SequelComplexDilutionBuffer',K.CONDITION)
    loaded=g.product('loaded_preparation','Sourcepreparedtemplateplusindependentinternalcontrolandloadingreagents',K.MOLECULAR)
    g.react('loading_pool','Source85ulperSMRTcellphysicaltemplate/controlloadingmix',S.POOL,((x,R.TEMPLATE),(control,R.TEMPLATE),(load_mix,R.CONDITION)),(loaded,))
    return finish(g,loaded,'PhysicalFLAMSMRTbell/polymerasepreparation; fullgeometryunassembled','SourceSequel/SMRTCell1Mv3/SequencingKit3.0; diffusionloading360–600minmovies and CCSminimum3passes. N10identifiermeaning/adapter/primergeometryunassigned; G/Iandnativeguanylationnotseparable. AnchoredTTTexcludesnative3primeendmodifications in source-selectedscope.')

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();p='ACGT'+'C'*9+'TTT';a=AnchoredCopy('ACGU',9)
    c('Anchored RT copies remaining tail and upstream RNA without artificialGI',a.first_strand(p),p+'T'*6+'ACGT')
    try:a.first_strand('ACGT'); rejected=False
    except ValueError:rejected=True
    c('Selected extension requires actual anchored register',rejected,True)
    g=protocol_graph();c('Tailing stop/cleanup precedes RT',g.steps['rna_primer_preheat'].inputs[0].material,'gi_clean_kept')
    c('PCR retains source RNA remnant scope', 'originalRNA' in g.materials['pcr_product'].label,True)
    c('Exonuclease digestion precedes actual separation',g.steps['smrt_clean1'].inputs[0].material,'exo_product')
    c('TwoSMRTbellcleanupsarephysical',g.steps['smrt_clean2'].inputs[0].material,'smrt_clean1_kept')
    c('Sequencingprimerconditioningexcludeslibrary',tuple(i.material for i in g.steps['condition_sequencing_primer'].inputs),('sequencing_primer_stock','primer_condition_buffer'))
    c('Independentcontroladdedaftercomplexcleanup',g.steps['loading_pool'].inputs[1].material,'internal_control')
    c('Partial graph cannot publish',bool(g.publication_problems()),True)
    src=Source(sources_root()/'flam-seq__10.1038+s41592-019-0503-y'/'s41592-019-0503-y_41592_2019_503_MOESM3_ESM.pdf.txt','get_sources --slugs flam-seq__10.1038+s41592-019-0503-y')
    if have(c,src,label='published FLAM protocol'):
        t=src.path.read_text();c('Published primer1 transcription',RT1 in t,True);c('Source says noethanol complex wash','do not wash with ethanol' in t,True)
    c.report()
if __name__=='__main__':selftest()
