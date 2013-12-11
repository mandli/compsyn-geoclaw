
MANUSCRIPT_FILES = \
	database.bib \
	commands.tex \
	paper.tex

OTHER_FILES = \
	highlights.pdf \
	cover_letter.pdf \
	response_2.pdf

FIGURES_FILES = \
	earth_quake_tsunami_source.png \

.PHONY: manuscript figures

all: clean manuscript figures

clean:
	-rm manuscript.tgz figures.tgz

manuscript: manuscript.tgz

figures: figures.tgz

# Note that COPYFILE_DISABLE=1 is included in an alias
manuscript.tgz: $(MANUSCRIPT_FILES) $(OTHER_FILES)
	tar cvzf $@ $(MANUSCRIPT_FILES)

# Note that COPYFILE_DISABLE=1 is included in an alias
figures.tgz: $(FIGURES_FILES)
	tar cvzf $@ $(FIGURES_FILES)

