# LIGO-Signal-Analysis-for-Cosmic-String-Cusps
This is a python script that can be used to try and detect cosmic string cusp events within LIGO strain data. 
Cosmic string cusps are hypothetical events that would be indicators that the theoretical object cosmic strings would exist.
It was predicted that cosmic strings would have formed during a phase transition during the very early universe.
Based on a model from Xavier Siemens of what the gravitational wave signal from a cusp event would be characterized as, 
raw strain data from LIGO can be compared using Optimal Matched Filtering to see if there is a cusp event signal within the data. 
Currently this overall method is used to identify the merger of two black holes, otherwise known as a compact binary coalescence.
For any section of data that is selected the program can also try to find a CBC signal as well and compare the Signal to Noise Ratio (SNR),
to that found while searching for a cusp signal.  

The main file will produce two graphs for the SNR results from the two templates. After downloading data from the LIGO portal simply change the name of the data file
where the arguemnet for rl.loaddata. 
