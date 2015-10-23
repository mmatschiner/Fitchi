# This script assumes that 'fitchi.py' is located in the same directory, together with the example file 'example.nex'.
# Running Fitchi requires the python3, as well as the installation of graphviz, and the python modules pygraphviz, biopython, scipy, and numpy.
# See http://www.evoinformatics.eu/fitchi.htm for more information.
# Start this script with 'bash run_simulation.sh'.

fitchi.py example.nex example.html -m 3 -p pop1 pop2
