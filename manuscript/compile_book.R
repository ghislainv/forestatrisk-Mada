#!/usr/bin/Rscript

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# license         :CC-BY-SA 4.0
# ==============================================================================

# Libraries
library(bookdown)
library(knitr)
library(kableExtra)

# Working directory
setwd("manuscript")

# Bookdown
# pdf
options(knitr.table.format="latex")
pdf_format <- bookdown::pdf_document2(citation_package="natbib", fig_caption=TRUE, keep_tex=FALSE,
                                      latex_engine="pdflatex", number_sections=FALSE, toc=FALSE,
                                      includes=list(in_header="header.tex", before_body="doc_prefix.tex"))
params <- list(title="",author="",date="")
bookdown::render_book("index.Rmd", output_format=pdf_format)

# html
# Don't indicate output_format to take into account YAML options
options(knitr.table.format="html")
# Dynamic YAML options
title_html <- "Accounting for spatial regional variability in modelling and forecasting deforestation -- the fate of Madagascar's forests"
author_html <- "Ghislain Vieilledent and Frédéric Achard"
data_html <- format(Sys.time(), "%d %B, %Y")
params <- list(title=title_html,author=author_html ,date=data_html)
bookdown::render_book("index.Rmd")
