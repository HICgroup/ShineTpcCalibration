#############################################################
#  User definitions section:

#  Calibration targets (final output files of calibration)
CALIBRATIONTARGETS := TPCPhaseCalibration.xml MHTDCCalibration.xml
#  Postprocessing calibration analysis program names
CALIBRATIONANALYZERS := TPCPhaseCalibrationCalculator MHTDCCalibrationCalculator
#  Generated input XML files to Shine
GENERATEDXMLS := $(patsubst %.xml.in,%.xml,$(wildcard *.xml.in))

#  Symbolic targets (not real generated files)
.PHONY: tools clean

#  Generate necessary input XMLs and programs
tools: $(GENERATEDXMLS) $(CALIBRATIONANALYZERS)

#  Run Shine module sequence and generate calibration data dump (production)
produce: $(GENERATEDXMLS)
#	ShineOffline -b bootstrap.xml && touch $@
	./runModuleSeq.sh $(TESTFILE) && touch $@

#  Merge and time slice the produced calibration data dumps
merge: produce
	calibMerge -i `ls tpcPhases*.root | head -c -1 | tr '\n' ','` -o mergedTPCPhases.root && \
	calibMerge -i `ls mhtdcDiffs*.root | head -c -1 | tr '\n' ','` -o mergedMHTDCDiffs.root && \
	touch $@

#  Run calibration analysis on merged and time sliced calibration data dumps
calibrate: merge $(CALIBRATIONANALYZERS)
	./TPCPhaseCalibrationCalculator -i mergedTPCPhases.root && \
	./MHTDCCalibrationCalculator -i mergedMHTDCDiffs.root && \
	touch $@

#  Alias for 'calibrate'
run: calibrate
	touch $@


#############################################################
#  The parts below are for experts only (steering the compilation):
#  set up input XML configuration and compile analyzer programs

SHINEOFFLINECONFIG := shine-offline-config
CONFIGFILES := $(shell $(SHINEOFFLINECONFIG) --config)
DBFILES := $(shell $(SHINEOFFLINECONFIG) --db-path)
DOCPATH := $(shell $(SHINEOFFLINECONFIG) --doc)
#TESTFILE := $(DOCPATH)/SampleEvents/run-009917x054_exampleevent.raw
TESTFILE := $(DOCPATH)/SampleEvents/run-014660x000_10exampleevent.raw

ifeq "$(strip $(CXX))" ""
  CXX := g++
endif

#  For generating input XMLs
%: %.in
	@echo -n "Generating $@ file..."
	@ sed -e 's!@''CONFIGDIR@!$(CONFIGFILES)!g;s!@''SHINEDBDIR@!$(DBFILES)!g' $< >$@
	@echo "done"
	@chmod +x runModuleSeq.sh

#  For compiling analyzer programs
%: %.C
	$(CXX) -o $@ $< `root-config --cflags` `root-config --libs`
%: %.cpp
	$(CXX) -o $@ $< `root-config --cflags` `root-config --libs`
%: %.cc
	$(CXX) -o $@ $< `root-config --cflags` `root-config --libs`
%: %.cxx
	$(CXX) -o $@ $< `root-config --cflags` `root-config --libs`

#  For cleaning up
clean:
	rm -f $(GENERATEDXMLS) $(CALIBRATIONANALYZERS) $(CALIBRATIONTARGETS) *.root *.pdf produce merge calibrate run
