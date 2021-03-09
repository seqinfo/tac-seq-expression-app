# TAC-seq expression app
This repository contains TAC-seq expression app for TAC-seq data analysis.

## Input
* [TAC-seq data analysis](https://github.com/cchtEE/TAC-seq-data-analysis) output file with:
  * column `sample` - sample ID
  * column `locus` - locus ID
  * column `molecule_count` - molecule count

* Target gene subset for modelling:
  * column `locus` - only one column with header "locus" and a list of selected gene names. Have to be a subset of `target` -> `locus` values.
  
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

Score table with:
* column `sample` - sample ID
* column `pre-rec` - sample score for the pre-receptivity
* column `rec` - sample score for the receptivity
* column `post-rec` - sample score for the post-receptivity
* column `predcted group` - the predicted class of the sample
* column `score` - general score to be used in the receptometer ranging from (pre=) 0 - (rec=) 100 - (post=) 200
