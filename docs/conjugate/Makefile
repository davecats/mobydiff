# Build the conjugate-heat-transfer design note.
#   make          -> conjugate_ibm.pdf
#   make clean    -> remove LaTeX intermediates (keeps the pdf)
#   make distclean-> also remove the pdf
# Needs pdflatex + latexmk with tikz/pgfplots (texlive-pictures).

# conjugate_ibm          the DESIGN note: alternatives, the full derivation,
#                        the corner model -- why the scheme looks as it does
# conjugate_ibm_asbuilt  the AS-BUILT note: only what the code does today
DOCS = conjugate_ibm conjugate_ibm_asbuilt

all: $(addsuffix .pdf,$(DOCS))

%.pdf: %.tex
	latexmk -pdf -interaction=nonstopmode -halt-on-error $<

clean:
	$(foreach d,$(DOCS),latexmk -c $(d).tex;)

distclean:
	$(foreach d,$(DOCS),latexmk -C $(d).tex;)

.PHONY: all clean distclean
