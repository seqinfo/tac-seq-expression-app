# TAC-seq expression app
This repository contains TAC-seq expression app for TAC-seq data analysis.

## Input
* [TAC-seq data analysis](https://github.com/cchtEE/TAC-seq-data-analysis) output file with:
  * column `sample` - sample ID
  * column `locus` - locus ID
  * column `molecule_count` - molecule count
  
* Target file with:
  * column `target` - values have to match with `locus` values in [TAC-seq data analysis](https://github.com/cchtEE/TAC-seq-data-analysis) output file
  * column `type` - possible values: `biomarker`, `housekeeper` or `spike_in`

* Control file with:
  * column `sample` - sample ID
  * column `label` - class label
  * column for each `biomarker`

## Output
Normalized count table with:
* column `sample` - sample ID
* column for each `biomarker`
