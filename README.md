# TCSGen
Generator for Timelike Compton Scattering process.

It includes Bethe-Heitler (BH) process cross section, and The interferenc (between TCS and BH) term cross section 

==================== I N S T A L L A T I O N =======================
In the base directory make d directory named lib 
> mkdir lib

Then make sure you have installed root, It will compile the code using thr $ROOTSYS enviromental variable

Just type make
> make

If make is successful, then make sure the lib directory, i.e. `pwd`/lib is in the $LD_LIBRARY_PATH

To run from an arbitrary location, $TSCGEN_DIR must be set to this repo's directory.

======================================================================================




===================== R U N N I N G ============================================
1) Make sure when running, the file CFFs_DD_Feb2012.dat is present in the directory, from
where you run the executable. It contains the grid of CFFs, and for interference term
cross section calculation, it uses that grid.


2) Also some basic variables are inside the config file named "GenOptions.dat"
- Nsim          Number of events to generate
- Eb            The beam energy in GeV
- tLim          The t_Max value in GeV
- EgMin         Photon energy minimum
- EgMax         Photon energy maximum
