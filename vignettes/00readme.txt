
do place vignettes here. When the vignettes are built, the following files will be copied to inst/doc:
- .Rmd
- .html
- .R

Hints for building vignettes are here:
 http://stackoverflow.com/questions/19372260/how-to-get-rstudio-to-automatically-compile-r-markdown-vignettes
 
 Options:
 1. use knit button
 2. call devtools::build_vignettes()
 3. build a source package and extract the vignettes from the tar.gz
 
 If you need a pdf version of a vignette. Do specify in the .Rmd as output: pdf_document
 
 