all:
	latex Galaxy.tex
	bibtex Galaxy
	latex Galaxy.tex
	latex Galaxy.tex
	dvips -t a4 Galaxy.dvi
	ps2pdf Galaxy.ps
	rm Galaxy.ps Galaxy.dvi Galaxy.tex Galaxy.aux Galaxy.bbl Galaxy.blg Galaxy.out Galaxy.log
	latex Faithful.tex
	bibtex Faithful
	latex Faithful.tex
	latex Faithful.tex
	dvips -t a4 Faithful.dvi
	ps2pdf Faithful.ps
	rm Faithful.ps Faithful.dvi Faithful.tex Faithful.aux Faithful.bbl Faithful.blg Faithful.out Faithful.log
	latex Tandmob.tex
	bibtex Tandmob
	latex Tandmob.tex
	latex Tandmob.tex
	dvips -t a4 Tandmob.dvi
	ps2pdf Tandmob.ps
	rm Tandmob.ps Tandmob.dvi Tandmob.tex Tandmob.aux Tandmob.bbl Tandmob.blg Tandmob.out Tandmob.log
	rm ./figures/figFaithful01.ps ./figures/figFaithful02.ps ./figures/figFaithful03.ps \
	./figures/figFaithful04.ps ./figures/figFaithful05.ps ./figures/figFaithful06.ps \
	./figures/figFaithful11.ps ./figures/figFaithful12.ps ./figures/figFaithful13.ps \
	./figures/figFaithful14.ps ./figures/figFaithful15.ps ./figures/figFaithful19.ps
	rm ./figures/figGalaxy01.ps ./figures/figGalaxy02.ps ./figures/figGalaxy03.ps \
	./figures/figGalaxy04.ps ./figures/figGalaxy05.ps ./figures/figGalaxy06.ps
	rm ./figures/figTandmob01.ps ./figures/figTandmob02.ps ./figures/figTandmob03.ps \
	./figures/figTandmob04.ps ./figures/figTandmob05.ps
