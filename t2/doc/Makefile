FILE=Capa

LODFIGS = $(patsubst %.odg,%.pdf,$(wildcard *.odg))

all: $(FILE).pdf
	evince $< &
	
$(FILE).pdf: $(wildcard *.tex) $(LODFIGS)
	cp -f ../mat/*.eps .
	cp -f ../mat/*.tex .
	cp -f ../sim/*.pdf .
	cp -f ../sim/*.tex .
	cp -f ../sim/*.ps .
	pdflatex $(FILE).tex
	pdflatex $(FILE).tex
	pdflatex $(FILE).tex
	
%.pdf: %.odg
	soffice --headless --convert-to pdf $<
	
clean:
	@rm -rf *.aux *.bbl *.blg *.glg *.glo *.gls *.ilg *.ist *.lof
	@rm -rf *.log *.lot *.nlo *.nls *.out *.toc *~ *.*% ./*.pdf ./*.eps
	
.PHONY: all clean
