This repository is developed to construct ultracold neutron (UCN) events from the hit-photon array of two paired PMTs. 
The code begins with reading a single TTress from the blinded replayed root files. 
The hits are stored in the form of a vector of pairs, where the pair comprises the time stamp of the photon hit and the PMT it fired in. 
The single vectors are filtered out to remove dead events using a 20 ns time window. The dead-corrected vectors from two PMTs are merged and sorted to construct UCN events.
The paired vector is supplied to a coincidence method which constructs UCN events by looking for events in two different PMTs with a 100 ns timing window. 
The outcome of the coincidence function is a struct carrying all the information of a UCN event (time, number of hits, length of the event).
For the 2022 segmented dagger, there are four independent structs with information about UCN events from PMT12, PMT34, PMT56, and PMT78. 
Each of the UCN events list is independently corrected for deadtime correction. 
The UCN events are later stored in a root TTree and also as histograms in a root file.
