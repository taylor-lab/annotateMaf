# annotateMaf
 [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Build Status](https://travis-ci.org/taylor-lab/annotateMaf.svg?branch=master)](https://travis-ci.org/taylor-lab/annotateMaf)
[![Coverage status](https://codecov.io/gh/taylor-lab/annotateMaf/branch/master/graph/badge.svg)](https://codecov.io/github/taylor-lab/annotateMaf?branch=master)

A set of functions to add variant annotation to a MAF file.
Sources currently include OncoKB, BRCA Exchange and somatic hotspots from the Taylor Lab.

## Installation

~~This package python modules `urllib3` and `ga4gh`, the latter of which only works in python (< 3.0).~~\
_2019-02-14: Disabled due to this [issue](https://github.com/BRCAChallenge/brca-exchange/issues/981). Instead of querying the API, we now use a static table._

Load and install the library this way:
``` r
devtools::install_github('taylor-lab/annotateMaf')
library(annotateMaf)
```

## Examples

Run the functions simply with your MAF (as a `data.table`, not the file path) as the input.

hotspot_annotate_maf requires a VEP-annotated MAF file.

``` r
# Note that the BRCA Exchange database is geared towards germline variants but by default the variant allele in a MAF is called Tumor_Seq_Allele2
annotated_maf = brca_exchange_annotate_maf(input_maf)

# Only retain oncogenic or likely oncogenic mutations after OncoKB annotation
maf %>% 
    oncokb_annotate_maf(input_maf) %>% 
    filter(oncogenic %like% 'Oncogenic) 
```

### Annotation sources:
- OncoKB: Queries latest version of [OncoKB](http://oncokb.org), version number included but currently no support for querying older versions.
- BRCA Exchange: Queries latest version of [BRCA Exchange](https://brcaexchange.org), also currently does not support versioning. 
- Somatic hotspots: List generated from PMIDs 26619011, 29247016, 28115009. Semi-manual curation was carried out to remove false-positive germline variants that were in the oldest publication. 
