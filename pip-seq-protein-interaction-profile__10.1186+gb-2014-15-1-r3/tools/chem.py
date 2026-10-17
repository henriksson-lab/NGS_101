"""Protocol-owned source, molecular states and reaction graph for protein-interaction-profile PIP-seq.

Authority: Silverman et al., Genome Biology 15:R3 (2014), PMC4053792 JATS
(Materials and methods: Cross-linking experiments, PIP-seq library preparation;
Results/Figure 1A overview).  Library construction itself is delegated to refs 29/30
(Li et al. 2012 Cell Reports; Zheng et al. 2010 PLoS Genetics); it is a warned boundary
here rather than a reconstructed adapter chemistry.
"""
from chemdraw import Construct
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)
from protocol_model import molecule, role
import seqprimers as sp


TITLE = 'PIP-seq — protein interaction profile sequencing'
NOTES = '01_protein_pip_seq.html'
SOURCE = 'Silverman et al., Genome Biology (2014), <a href="https://doi.org/10.1186/gb-2014-15-1-r3">doi:10.1186/gb-2014-15-1-r3</a>.'
SUMMARY = ('Cross-linked (formaldehyde or UV) or untreated HeLa lysate is split four ways: '
           'single- or double-strand-specific RNase, applied either with proteins present '
           '(footprint) or after SDS/proteinase K removal (digestion control). '
           'Comparing footprint with control libraries reveals protein-protected sites.')
CAVEAT = ('This is protein-interaction-profile PIP-seq, not the unrelated particle-templated '
          'single-cell RNA-seq protocol with the same acronym. The paper performs strand-specific '
          'library construction "as previously described" (refs 29, 30) and does not print '
          'adapter, RT or PCR oligos. That step is therefore shown as an unresolved boundary.')
READ_LENGTH = 50
FINAL_LIBRARY = Construct((
    role('5-prime library adapter (predecessor method)'),
    role('RNase-protected RNA-derived insert, strand retained'),
    role('3-prime library adapter (predecessor method)')),
    name='protein-interaction-profile PIP-seq library')
FINAL_CAPTION = ('One strand-specific, DSN-normalized library is made for each of the four '
                 'branches (ss/ds RNase × footprint/control) per replicate. The paper does not '
                 'describe adapter identities, indexing or multiplexing.')
SEQ_PRIMERS = ()
READ_LENGTHS = None
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = ('Each library was sequenced on an Illumina HiSeq2000 using the standard '
                     f'protocol for {READ_LENGTH}-bp single-read sequencing. The paper does not state whether libraries were '
                     'indexed or how they were assigned to lanes.')
SEQUENCING_UNAVAILABLE = ('The adapter and sequencing-primer sequences come from the cited '
                          'predecessor library method and are not given in the defining paper.')
READ_SPECS = (sp.read_spec('Read 1', READ_LENGTH,
                           'RNase-protected fragment (strand-specific); Additional file 1 reports raw, trimmed and mapped reads',
                           'include'),)


def oligos():
    return ()


def protocol_graph():
    source = SOURCE
    g = ProtocolGraph(TITLE, authority=source)
    cells = g.input("cells", "adherent HeLa cells at 90% confluence (~18 million per 15-cm dish)",
                    MaterialKind.BIOLOGICAL,
                    state=molecule("cellular RNPs", (role("RNA"), role("bound proteins"))),
                    primary=True)
    crosslinker = g.input("crosslinker",
                          "1% formaldehyde 10 min RT then 125 mM glycine 5 min; or UV-C 254 nm "
                          "400 mJ/cm² (Stratalinker 2400); or no cross-linking",
                          MaterialKind.CONDITION)
    crosslinked = g.product("crosslinked_cells", "cross-linked (or untreated) HeLa cells",
                            MaterialKind.BIOLOGICAL, state=cells.state)
    g.react("crosslink", "Cross-link RNA–protein complexes in cells", StepKind.CROSS_LINKING,
            ((cells, InputRole.TEMPLATE), (crosslinker, InputRole.CONDITION)), (crosslinked,))
    rip = g.input("rip_lysis", "RIP buffer (25 mM Tris pH 7.4, 150 mM KCl, 5 mM EDTA, 0.5% NP40, "
                  "10 µM DTT, protease inhibitors; 850 µL per 10 million cells), manual grinding",
                  MaterialKind.CONDITION)
    dnase = g.input("dnase", "RNase-free DNase (Qiagen)", MaterialKind.CATALYST)
    lysate = g.product("lysate", "DNA-depleted cross-linked lysate", MaterialKind.MOLECULAR,
                       state=molecule("crosslinked RNPs", (role("RNA"), role("bound proteins"))))
    g.react("lyse", "Lyse cells in RIP buffer and treat with DNase", StepKind.EXTRACTION,
            ((crosslinked, InputRole.TEMPLATE), (rip, InputRole.CONDITION),
             (dnase, InputRole.CATALYST)), (lysate,))
    branches = []
    for identifier, label in (
            ("bound_ss", "footprint, ssRNase (RNaseONE)"),
            ("bound_ds", "footprint, dsRNase (RNaseV1)"),
            ("removed_ss", "digestion control, ssRNase (RNaseONE)"),
            ("removed_ds", "digestion control, dsRNase (RNaseV1)")):
        branches.append(g.product(identifier, label, MaterialKind.MOLECULAR, state=lysate.state))
    g.react("split", "Split lysate into footprint and RNase-digestion-control samples for each RNase",
            StepKind.SPLIT, ((lysate, InputRole.TEMPLATE),), tuple(branches),
            changes_molecule=False)
    reviewed = ["cells", "crosslinker", "crosslink", "crosslinked_cells", "rip_lysis",
                "dnase", "lyse", "lysate", "split", *(b.id for b in branches)]
    rnases = {
        "ss": ("RNaseONE (Promega) 100 U/mL, 200 µg/mL BSA, 1× RNaseONE buffer, 1 h RT",
               "single-strand-RNase-resistant RNA"),
        "ds": ("RNaseV1 (Ambion) 2.5 U/mL in 1× RNA structure buffer, 1 h 37 °C",
               "double-strand-RNase-resistant RNA"),
    }

    def step(prefix, name, label, kind, current, extra, state):
        out = g.product(f"{prefix}_{name}_out", label, MaterialKind.MOLECULAR, state=state)
        g.react(f"{prefix}_{name}", label, kind, ((current, InputRole.TEMPLATE), *extra), (out,))
        reviewed.extend((f"{prefix}_{name}", out.id, *(m.id for m, _ in extra)))
        return out

    for branch in branches:
        p, nuclease = branch.id, branch.id[-2:]
        current = branch
        if p.startswith("removed"):
            sds = g.input(f"{p}_strip_condition", "1% SDS + 0.1 mg/mL proteinase K (Roche), then ethanol "
                          "precipitation", MaterialKind.CONDITION)
            current = step(p, "strip", "Remove proteins before RNase digestion",
                           StepKind.DEPROTEINIZATION, current, ((sds, InputRole.CONDITION),),
                           molecule("protein-depleted RNA", (role("RNA"),)))
        enzyme = g.input(f"{p}_rnase", rnases[nuclease][0], MaterialKind.CATALYST)
        protected = ("protein-protected and " if p.startswith("bound") else "") + rnases[nuclease][1]
        current = step(p, "digest", "Digest with the structure-specific RNase", StepKind.DIGESTION,
                       current, ((enzyme, InputRole.CATALYST),),
                       molecule(protected, (role(protected),)))
        reverse = g.input(f"{p}_reverse_condition", "1% SDS + 0.1 mg/mL proteinase K 15 min RT; then 65 °C 2 h "
                          "(formaldehyde) or re-treatment with 1% SDS + 1 mg/mL proteinase K 30 min (UV)",
                          MaterialKind.CONDITION)
        current = step(p, "reverse", "Denature/digest proteins and reverse cross-links",
                       StepKind.DEPROTEINIZATION, current, ((reverse, InputRole.CONDITION),),
                       molecule("RNase-resistant RNA fragments", (role("RNase-resistant RNA fragment"),)))
        kit = g.input(f"{p}_mirneasy", "Qiagen miRNeasy RNA isolation", MaterialKind.CONDITION)
        current = step(p, "isolate", "Isolate digested RNA (miRNeasy)", StepKind.EXTRACTION,
                       current, ((kit, InputRole.CONDITION),), current.state)
        method = g.input(f"{p}_libprep", "strand-specific library preparation as previously "
                         "described (refs 29, 30)", MaterialKind.CONDITION)
        raw = g.product(f"{p}_raw_library", f"{branch.label}: strand-specific library before "
                        "normalization", MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True)
        g.react(f"{p}_libprep_step", "Prepare a strand-specific sequencing library", StepKind.MAGIC,
                ((current, InputRole.TEMPLATE), (method, InputRole.CONDITION)), (raw,),
                reason="The defining paper delegates library construction to Li et al. 2012 (unavailable "
                       "here) and Zheng et al. 2010; the exact inherited adapters, ligation/RT/PCR "
                       "steps and primers are not reconstructed.")
        dsn = g.input(f"{p}_dsn", "Illumina DSN normalization: 100 ng library denatured 95 °C, "
                      "re-annealed 5 h 68 °C, 2 µL DSN enzyme (1 U/µL)", MaterialKind.CATALYST)
        library = g.product(f"{p}_library", f"{branch.label}: DSN-normalized PIP-seq library",
                            MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True)
        g.react(f"{p}_dsn_step", "Deplete rapidly re-annealing abundant duplexes with DSN",
                StepKind.DIGESTION, ((raw, InputRole.TEMPLATE), (dsn, InputRole.CATALYST)),
                (library,))
        g.sequence(library, SequencingPlan(
            limitation="Adapter and primer sequences belong to the cited predecessor method.",
            acquired_roles=("Read 1",), read_lengths={"Read 1": READ_LENGTH},
            read_uses={"Read 1": "include"}))
        reviewed.extend((f"{p}_libprep", f"{p}_raw_library", f"{p}_libprep_step", f"{p}_dsn",
                         f"{p}_library", f"{p}_dsn_step"))
    g.annotate(SourceDiscrepancy(
        "bound_ss_library", "Raw reads per library: Results range vs Additional file 1",
        (SourceClaim("Results, PIP-seq overview paragraph (PMC4053792 JATS line 71)",
                     "approximately 31 to 60 million raw reads per library (Additional file 1)"),
         SourceClaim("Additional file 1 (S1.xlsx, Supplemental Table 1, raw reads column)",
                     "per-library raw reads from 24,604,010 (UV dsRNase footprint replicate 2) "
                     "to 103,702,805 (formaldehyde dsRNase footprint replicate 2)")),
        "The page gives no per-library read depth; only the 50-bp single-read format is modelled.",
        "Both are in the published article and its own supplement.",
        "It is unclear whether the Results range refers to a subset of the libraries."))
    g.declare_reviewed_scope(reviewed, source)
    return g.require_publishable()
