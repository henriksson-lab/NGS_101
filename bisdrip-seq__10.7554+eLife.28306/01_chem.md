# bisDRIP-seq: source-qualified private model

Dumelie and Jaffrey (2017), **Defining the location of promoter-associated R-loops at near-nucleotide resolution using bisDRIP-seq**, eLife6:e28306. [Published original v2](https://doi.org/10.7554/eLife.28306).

🟢 Read the published v2 main Methods514–546 and584–592, with control identities556–582 and conceptual attribution56–74. This is one original; no preprint or unreviewed historical kit recipe is merged. The standard fetch acquired associated code, reporting form and result tables, which remain unreviewed. Mandatory scraper yielded3 candidates in one promoter-results workbook; their printed contexts resemble SCARNA gene-name expansions rather than adopted oligos. They are not used as sequences and do not establish sequence absence.

🟢 This method combines non-denaturing bisulfite during lysis with S9.6 RNA-DNA-hybrid capture. It preserves a chemical mark on exposed DNA before the hybrid is processed. This differs from ordinary DRIP-seq capture and from BisMapR's targeted nuclease release. Methods538 explicitly say the IP is modeled after Ginno2012 (bibliography930–932, DOI10.1016/j.molcel.2012.01.017); this establishes the S9.6 IP strategy scope, not equality of all precursor steps or a constructor import.

🟢 Standard source route: MCF7 dish →PBS wash/scrape/pellet →1.5ml bisulfite lysis preparation →two475µl aliquots each into25µl LDS →37°C30min dark →increased ammonium bisulfite/hydroquinone37°C2h dark →nominal200fold dialysis, then5000fold into ProteinaseK buffer →proteinaseK/SDS overnight37°C →two-tube phenol/chloroform/chloroform/ethanol workup and combination of precipitated DNA →restriction cocktail HindIII/EcoRI/BsrGI/XbaI/SspI37°Covernight.

The private graph selects one475µl tube as a scoped route. The other source475µl tube and arithmetic550µl remainder are explicitly unresolved, without inventing their final fate or rejoining them before a source-established pooling point. The internal organic-workup split/recombine remains one opaque physical component, not an extra invented purification. This graph deliberately does not claim complete source coverage.

🟢 Separate0.8%agarose and Quant-iT/SpectraMax measurements precede IP. Source samples were used only if digest mostly1–10kb; the graph does not predict a QC outcome. Source20µg digested DNA goes to20µgS9.6 IP overnight4°C, then150µl DynabeadsProteinG2h4°C, magnetic separation and three washes. Remaining digest and analytical aliquots remain separate physical fates.

🟢 On beads,150µl0.3NNaOH37°C20min completes the bisulfite reaction, followed by HCl/Tris neutralization. RNaseI50°C20min →RNaseH/Mg37°C20min →proteinaseK50°C1h precede magnetic supernatant transfer. Digested RNA/protein products physically remain until actual workup; no ideal protection or instant molecular deletion is claimed. Phenol/chloroform workup includes an actual400µl aqueous transfer, whose untransferred portion remains unresolved, then chloroform and glycogen-assisted ethanol precipitation48h.

🟢 Recovered DNA is sonicated to approximately300bp with CovarisS2. Historical Zymo Pico-MethylSeq random-priming/index/amplification is an opaque component. The source explicitly omits its bisulfite-treatment steps; the graph never adds a second kit conversion. Separate Bioanalyzer/Quant-iT QC precedes loading five or six libraries per lane. Acquisition is HiSeq2500, single-index100bp paired-end. Exact historical kit oligos, junctions, amplification and read landing remain unassembled.

Input-only, RNaseH/control and delayed-bisulfite variants are identified in the source but not fully modeled here. They are not additional counted protocols, nor silently merged into the standard IP route. The source RNaseH matched experiment physically splits one sample, whereas study replicates use independent dishes.

🟡 `SelectedDeamination` is an inferred, caller-selected DNA C-to-U witness. Its physical product retains U; a separate reported-reference view uses T. It neither converts RNA nor predicts accessibility, hybrid protection, methylation or conversion efficiency. Bismark conversion/read filtering is computational, not physical disappearance. No finished schematic or accepted missing-information claim.

Open scope: historical Pico-MethylSeq guide/revision; source-cited Ginno2012 IP, Yu2003 non-denaturing mapping and Hayatsu2004 conversion original scopes. No linked original is presumed unavailable merely because it was deferred.
