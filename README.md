# Fluka-Scripts

DS_20k_fluka_updated.inp;

This is the input file which contains the information about the materials of the detector, dimensions, and densities for the DarkSide-20k experiment. 
This files is updated according to the recent TDR document from DarkSide collaboration. 

Mgdraw_ds_updated.f 
This files provides the information about the collision tape which is the output of fluka. This file is a user routine file in fluka simulation. We use this file
file to write the collision tape, i.e which region we want to record and then it provides us the text file with energies, particle type, time etc.

source_000_new_step123.f:
This files provides the information about reading the source file. The source file contains muon primary events. we use this source for the muon shower in the cavern. The source file contain primary events which we read by this source file.
