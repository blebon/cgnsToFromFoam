Description
===========

This branch was forked and adapted from svn://svn.code.sf.net/p/openfoam-extend/svn/trunk/Breeder_2.0/OSIG/TurboMachinery/
via mirror https://github.com/Unofficial-Extend-Project-Mirror/openfoam-extend-Breeder2.0-OSIG-TurboMachinery

That original repository is part of the TurboMachinery SIG, at http://openfoamwiki.net/index.php/Sig_Turbomachinery

All of the non-CGNS related files were removed, given that this repository focuses on porting CGNS to more recent OpenFOAM versions.

For any problems with this fork, please report it here: https://github.com/wyldckat/cgnsToFromFoam/issues


Structure
=========

The folder "applications" contains new applications.

The folder "src" contains new libraries.


Requirements
============

This is currently an experimental port for OpenFOAM 5.x and it relies on more
recent CGNS and HDF5 libraries, which usually are provided by your Linux
Distribution.

Therefore, beyond the usual requirements which are already needed for building
OpenFOAM from source code, you also need to install those libraries, which can
be done as follows:

  * To install those packages on Ubuntu 14.04 or newer, run:

    ```
    sudo apt install libhdf5-dev libcgns-dev autoconf git
    ```


  * To install those packages on CentOS 7, run (one line at a time):

    ```
    su -
    yum install epel-release
    yum install cgnslib-devel git
    exit
    ```


Installation and respective git branches
========================================

Currently there are the following branches:

  * [OF5x](https://github.com/wyldckat/cgnsToFromFoam/tree/OF5x) - This branch
    will build with OpenFOAM 5.x

Quick usage instructions:

  1. Clone the repository:

     ```
     git clone https://github.com/wyldckat/cgnsToFromFoam.git
     ```

  2. Checkout the right branch, e.g. for `OF5x`:

     ```
     cd cgnsToFromFoam
     git checkout OF5x
     ```

  3. Now build:

     ```
     wmake all -j src
     wmake all -j applications
     ```


License
=======

Part of the code contained in this repository is following the same license as OpenFOAM, namely the GNU General Public License v3.
Other parts are in GNU General Public License v2, but was indicated to be upgradable on a need basis.
