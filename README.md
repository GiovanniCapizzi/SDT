SOUND DESIGN TOOLKIT: INSTALLATION INSTRUCTIONS
===============================================

### Note

This **fork** adds the C implementation of the [SeamCarvingAudio](https://github.com/GiovanniCapizzi/SeamCarvingAudio) algorithm. This is still **WIP** and it may not be integrated in the future. Please check the [SDT](https://github.com/SkAT-VG/SDT) repository of this project and the commit changes.

#### Diffs

I wrote 4 new files:

```
Pd/seamcarving~-help.pd
src/Pd/seamcarving~.c
src/SDT/SDTSeamCarving.c
src/SDT/SDTSeamCarving.h
```

I also updated these:

```
src/Pd/SDT.c
build/macosx/Makefile
```

The first one simply introduced the seam carving setup in `purrdata`. The second one has new flags to compile and link [GSL](https://www.gnu.org/software/gsl/) library I used to work with matrices.

#### Todo

- [ ] Update the documentation in `SDTSeamCarving.h` and `SDTSeamCarving.c` 📝;
- [ ] Fix all build files (rethink the GSL integration)  💥;
- [ ] Use the same method of erasing a seam (same as [here](https://github.com/GiovanniCapizzi/SeamCarvingAudio/blob/master/fixed_frequencies.py)) 🔨;
- [ ] Verify and validate the algorithm ✔️; 
- [ ] Max?
- [ ] Make a coffee ☕.



PRECOMPILED BINARIES
--------------------

Precompiled binaries are available for Max and PureData, running on Mac OS X and Windows.
Simply download the universal SDT package for from the official website:

http://www.soundobject.org/SDT

then unpack it and copy the branch for your operating system and platform into the desired
destination folder.

COMPILING FROM SOURCE
---------------------

**MAC OS X**

1. In a terminal, type the following commands to compile the software in all its flavors
(Max package, Pd library, Apple framework):

        cd build/macosx
        make

2. Install one or more products. The script will install the desired product in the given
<path>, creating a SDT subfolder:

        make install_max DSTDIR=<path>
        make install_pd DSTDIR=<path>
        make install_core DSTDIR=<path>

3. To clean the source directories after compilation:

        make clean
	

**WINDOWS**

To compile the Sound Design Toolkit under Windows, you need a distribution of the
GNU C Compiler and a UNIX style shell, as provided in MinGW + MSYS 
(http://www.mingw.org, recommended) or Cygwin (http://www.cygwin.com).

1. Once the compilation environment is installed, open its shell and issue the following
commands to compile the software in all its flavors (Max package, Pd library, Shared DLL):

        cd build/win32 (or cd build/win64 for the x64 version)
        make

2. Install the desired products. The script will install the desired product in the given
<path>, creating a SDT subfolder:

        make install_max DSTDIR=<path>
        make install_pd DSTDIR=<path> (only for 32 bit)
        make install_core DSTDIR=<path>

3. To clean the source directories after compilation:

        make clean
	

**LINUX**

1. In a terminal, type the following commands to build the SDT:

        cd build/linux
        make
    
2. Install the SDT. By default, the building environment will install a shared library in
/usr/lib and a collection of PureData externals and patches in /usr/lib/pd/extras/SDT.
Root privileges may be required to access the default install path. If you want to change
the install path, provide a PREFIX argument:
        
        make install
        make install PREFIX=<path>
	
3. To clean the source directories after compilation:

        make clean
