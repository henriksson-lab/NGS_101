# External reference manifest

Originals and standard text twins live at `/husky/henriksson/NGS101_papers/sources/paired-damage-seq__10.1038+s41592-025-02632-3/`, outside this checkout.

Re-obtain with `python3 tools/get_sources.py --slugs paired-damage-seq__10.1038+s41592-025-02632-3`; then run `python3 tools/scrape_primers.py /husky/henriksson/NGS101_papers/sources/paired-damage-seq__10.1038+s41592-025-02632-3/`. Only declared standard routes are used; failed/subscription responses are not chemistry originals.

- Original: https://pmc-oa-opendata.s3.amazonaws.com/PMC12182908.1/PMC12182908.1.xml
  Locator: `/husky/henriksson/NGS101_papers/sources/paired-damage-seq__10.1038+s41592-025-02632-3/s41592-025-02632-3_PMC12182908.xml.txt`.
  Scope: AssociatedauthorMS main59/125–185/223/227; actualjointrecipe/sourcecomponentattribution andcomputation; noentirefinal-equivalenceclaim.
- Original: https://www.nature.com/articles/s41592-025-02632-3
  Locator: `/husky/henriksson/NGS101_papers/sources/paired-damage-seq__10.1038+s41592-025-02632-3/s41592-025-02632-3.html.txt`.
  Scope: PublishedfinalNaturepage identity/abstractandattachmentmapping; actualMethodsnotreceived(subscriptionpreview). MapsMOESM3asTables1–5,notSourceData.
- Original: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41592-025-02632-3/MediaObjects/41592_2025_2632_MOESM3_ESM.xlsx
  Locator: `/husky/henriksson/NGS101_papers/sources/paired-damage-seq__10.1038+s41592-025-02632-3/s41592-025-02632-3_41592_2025_2632_MOESM3_ESM.xlsx`.
  Scope: PublishedoriginalTable1_1B3:B29,Table1_2actualDNA/polyT/N6rows2–37,first18rowsTable1_3..6andall6Table1sheetidentities. WholepMENTsddC andselectedoligomodificationstringsverifieddirectlyinoriginalXML. Resulttabsnotchemistryauthority. Standardtexttwinconvertedexplicitlyforguard;targetedscraper39exactlocations.

The28MB published workbook is explicitly converted using `python3 tools/doctext.py FILE.xlsx` when automatic acquisition skips conversion. Targeted `scrape_primers.py DIR --find AGGCCAGAGCATTCG` reads that large original. The source guard uses its external standard text twin; no exception or shared guard change.
The failed associated-manuscript SupplementaryTables.zip is not requested as missing material: the actual published Table1 workbook has been retrieved and read.

No originals, extracted text or rendered original images are checked in. Exact source versions and capped selected sections are not complete source/linked-information/absence reviews.
