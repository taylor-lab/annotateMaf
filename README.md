# annotatemaf

A set of functions to add variant annotation to a MAF file.
Sources currently include OncoKB, BRCA Exchange and somatic hotspots from the Taylor Lab.

## Installation

Load and install the library this way:

``` r
devtools::install_github('taylorlab/annotate-maf')
library(annotatemaf)
```

## Example

Run the functions simply with your MAF (as a `data.table`, not the file path) as the input.

``` r
annotated_maf1 = brca_exchange_annotate_maf(maf1)
annotated_maf2 = oncokb_annotate_maf(maf2)
annotated_maf3 = hotspot_annotate_maf(maf3)
```

### Annotation sources:
- OncoKB: Queries latest version of [OncoKB](http://oncokb.org), version number included but currently no support for querying older versions.
- BRCA Exchange: Queries latest version of [BRCA Exchange](https://brcaexchange.org), also currently does not support versioning. 
- Somatic hotspots: List generated from PMIDs 26619011, 29247016, 28115009. Semi-manual curation was carried out to remove false-positive germline variants that were in the oldest publication. 

