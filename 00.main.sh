#!/bin/bash

# ####################################
# # ******************************** #
# # * SF sequence analysis script  * #
# # ******************************** #
# ####################################

# start all star formation sequence analysis
# scripts in correct order

# run this script as
# scripts/00.main.sh 2>&1 | tee pipeline_run.log

#############################################################

scripts/01.masks.sh
scripts/02.load_into_CLASS.sh
scripts/03.list_emission_pixels.sh

for j in `seq 1 6`
do
	scripts/04.fit_NH3.sh NH3_${j}-${j}

ipython << ipythonINPUT
J=${j}
execfile('scripts/05.make_fit_list.py')
ipythonINPUT

	scripts/06.make_map.sh NH3_${j}-${j}
done

scripts/07temperature_map.sh

ipython << ipythonINPUT2
execfile('scripts/08.everything_list.py')
execfile('scripts/09.fit_results_histogram.py')
execfile('scripts/10.sequence_plots.py')
ipythonINPUT2
