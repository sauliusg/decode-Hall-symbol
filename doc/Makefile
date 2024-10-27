
# Generate CSV and TSV files from the edited ODS files.

SOURCE_ODS = $(wildcard *.ods)

TARGET_TSV = ${SOURCE_ODS:%.ods=%.tsv}

.PREVIOUS: %.csv

.PHONY: all clean distclean mostlyclean cleanAll display

all: ${TARGET_TSV}

display:
	@echo ${TARGET_TSV}

%.csv: %.ods
	ods2csv $< > $@

%.tsv: %.csv
	csv2tab $< > $@

clean:
	rm -f ${SOURCE_ODS:%.ods=%.csv}

mostlyclean: clean

distclean: mostlyclean
	rm -f ${TARGET_TSV}

cleanAll: distclean

