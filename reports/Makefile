# Makefile to use knitr for package vignettes
# Original source from Yihui Xie:
# http://yihui.github.com/knitr/demo/vignette/

# For more about Makefiles and R, see:
# http://www.stat.auckland.ac.nz/~stat782/downloads/make-tutorial.pdf

.PHONY: clean

# Put all PDF targets here, separated by spaces.
PDFS= paper-rda.pdf

all: $(PDFS) 

%.pdf: %.Rnw
	Rscript -e "library(knitr); knit2pdf('$*.Rnw')"

clean:
	@echo "Cleaning Up"
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary demo* *.synctex.gz
