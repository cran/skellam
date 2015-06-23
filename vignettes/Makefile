
all: skellam.pdf 

clean:
	rm -rf *.tex *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ cache/ *.dvi *.4* *.css *.idv *.lg *.tmp *.xref *.gz

skellam.tex: skellam.rnw
	'$(R_HOME)/bin/Rscript' -e "knitr::knit('skellam.rnw')" $(Rargs)

skellam.pdf: skellam.tex
	xelatex skellam;
	xelatex skellam
