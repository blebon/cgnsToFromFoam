Description
===========
This branch was forked and adapted from svn://svn.code.sf.net/p/openfoam-extend/svn/trunk/Breeder_2.0/OSIG/TurboMachinery/ - via mirror https://github.com/Unofficial-Extend-Project-Mirror/openfoam-extend-Breeder2.0-OSIG-TurboMachinery

That original repository is part of the TurboMachinery SIG, at http://openfoamwiki.net/index.php/Sig_Turbomachinery


This is currently an experimental port for OpenFOAM 5.x and it relies on Ubuntu's more recent CGNS and HDF5 libraries. To install those packages, run:

    sudo apt install libhdf5-dev libcgns-dev


This branch OF5x is provided at https://github.com/wyldckat/cgnsToFromFoam

All of the non-CGNS related files were removed, given that this repository focuses on porting CGNS to more recent OpenFOAM versions.


Structure
=========
The folder "applications" contains new applications.

The folder "src" contains new libraries.


License
=======

Part of the code contained in this repository is following the same license as OpenFOAM, namely the GNU General Public License v3.
Other parts are in GNU General Public License v2, but was indicated to be upgradable on a need basis.
