# Golden Gate animation source review — 1 November 2026

Selected protocol-owned illustrative example, not the experimental plasmids.

- Published defining article: Engler et al.2008, DOI10.1371/journal.pone.0003647,
  external `sources/golden-gate__10.1371+journal.pone.0003647/` main PDF and JATS.
  Read Results basic cloning strategy, subcloning experiment and Figure1 caption.
  Inspected the rendered original PDF page2, including Figure1 (temporary rendering
  `/tmp/golden-gate-original-2.png`). The figure shows cuts outside recognition
  sequences, complementary numbered fusion ends and loss of sites from productD.
  Methods paragraph on page2 names BsaI and T4 DNA ligase in ATP-containing buffer.
  Original methods explain restriction/ligation in one tube; separate animation
  shots distinguish molecular relations and are not separate source incubations.
- Ran the source bundle primer scraper without a known-only filter. No source
  oligo is adopted: displayed exact bases come from the explicitly invented
  linear-part example in the protocol module, not the source's donor/vector maps.
- Protocol owner `golden_gate.py`: `example_fragments()` and `precursor()` supply
  all input sequences; canonical BsaI geometry is imported. Its offsets1/5 are
  measured after the six-base recognition site, so simulator cuts are7/11 from
  motif start, for both forward and reverse orientations. Both physical cuts are
  generated on each substrate; no strand interval is hand repaired to fit product.

Simulation accounts for all128 initial nucleotides. Cleavage removes eight backbone
bonds, no atoms. Explicit end exchange unpairs cross-fragment contacts without
claiming degradation. Four complementary pairs at each selected seam anneal.
Four ATP-dependent3′OH→5′P ligations close both28nt strands. The selected circular
product loses both site orientations; all72 flanking nucleotides retain their
actual chains and remaining hybrid contacts. Pink marks require actual bonds.

Scope: invented nonfunctional two-part construct; selected successful ligation,
not alternate-product kinetics, clone validity, temperature optimization or yield.
The full native vector/insert/backbone population and functional maps remain outside
this animation. Nothing establishes confirmed absence of source information.
