# Build the conjugate-heat-transfer design note.
#   make          -> conjugate_ibm.pdf
#   make clean    -> remove LaTeX intermediates (keeps the pdf)
#   make distclean-> also remove the pdf
# Needs pdflatex + latexmk with tikz/pgfplots (texlive-pictures).

DOC = conjugate_ibm

all: $(DOC).pdf

$(DOC).pdf: $(DOC).tex
	latexmk -pdf -interaction=nonstopmode -halt-on-error $(DOC).tex

clean:
	latexmk -c $(DOC).tex

distclean:
	latexmk -C $(DOC).tex

.PHONY: all clean distclean
