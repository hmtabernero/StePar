How to install StePar
======================

StePar is a code used to automatically infer stellar atmospheric parameters using the EW method.For the time being this is the code and the installation instructions. A paper describing this code has been recently accepted in A&A. 

If you use this code please cite Tabernero et al. (2019), in press. This Code is under the two clause BSD-licence.

Perhaps these instructions need more beta-testing. If you find any bugs, let me now.

Requirements: GNU/Linux, Mac OS, BSD.
Dependencies: MOOG, numpy, scipy, gfortran, make

I provided a MOOG folder with everything you need. However, it is wise to check MOOGS's webpage: https://www.as.utexas.edu/~chris/moog.html 

The MOOG folder includes a pretty awesome hack written by Sergi Blanco Cuaresma to avoid supermongo (Thanks, Sergi!).

Highly recommended auxiliary tools:

	ARES: https://github.com/sousasag/ARES 
	TAME: http://astro.snu.ac.kr/~wskang/tame/

You must install numpy and scipy by yourself (pip, anaconda, apt, pacman, yum, etc), as well as development tools (make, gfortran, ...)

StePar is known to work on Mac OS, and some popular distros (Arch, Debian, Fedora, Manjaro, Solus, and Ubuntu). I assume it will run on a BSD system as well. If you try it let me now.

StePar will only derive parameters, it will not measure any EWs for you. However, we do RECOMMEND to grab ARES or TAME to do the task. 

To install dependencies, please use pip or your distro's repo.

pip:

	$ pip3 install --user numpy
	$ pip3 install --user scipy
	
	Note: Install gfortran from repo, source or whatever. If you know how to compile your own things, this step is a piece of cake, :D.

Debian or Debian-based (i.e. Ubuntu) distro:

	$ sudo apt update
	$ sudo apt upgrade
	$ sudo apt install python3-numpy python3-scipy gfortran build-essential

	Note: Debian is so stable that it hurts. Its software is incredibly reliable, but also incredibly old (according to Matt Hartley, see: https://www.youtube.com/watch?v=IC2op3AZMfM).
	
Arch or Arch-based (i.e. Manjaro):

	$ sudo pacman -Syu
	$ sudo pacman -S python-numpy python-scipy gfortran base-devel

	Note: Arch is a bleeding edge distro, be careful as it fetches updates directly from upstream. Manjaro is more or less the same.  

Mac OS:
	Use a ports tree. You are on a BSD-like system. Mac ports is a good one:

	$ sudo port -d selfupdate
        $ sudo port install py37-numpy py37-scipy gcc8

Clone the repo:

	git clone https://github.com/hmtabernero/StePar

Now it is time to compile MOOG. This is the same thing in each operating system (GNU/Linux, BSD or Mac OS)

In the StePar folder: 

	$ cd MOOG
	$ make -f Makefile.fake clean
	$ make -f Makefile.fake

If it compiled you are good to go:

	$ cp MOOGnointro ..

If you have cloned the repo you will have a folder structure already set up. However, I think it is important to know the folders:

	LOGS -> the whole iterative process is stored here
 	
	EW -> here you place your EW file (in MOOG format)
	
	STD_OUT -> The last MOOG iteration (aka out1)
	
	ABUN -> The last MOOG iteration is stored here (aka out2)

	DATA -> ABO data included with MOOG (Barklem.dat and BarklemUV.dat). It is important to not skip these files.

	TXT -> MOOG response files

	PAR -> MOOG par file
	
Just in case:
	
	$ mkdir LOGS EW STD_OUT ABUN DATA TXT PAR

Model grid:
	
	MARCS1M.bin -> A grid of MARCS atmospheric models. 
	
For the time being this version only accepts MARCS models. In the near future I will include KURUCZ models. If there is a particular need for a GRID (i.e., APOGEE or PHOENIX), please let me know.

If you did everything correctly, please place your EW file in EW/.  

There are already two example files in that particular folder: 

	Sun HARPS  -> HARPS.Archive_Sun-4Fe.l
	Sun NARVAL -> NARVAL_Sun-1Fe.l

Just in case:

	$ chmod +x MOOGnointro
	$ chmod +x *.sh

Now it is time to properly run StePar for these two examples. Check the contents of runStePar.sh:

	./runstar.sh HARPS.Archive_Sun-4
	./runstar.sh NARVAL_Sun-1

If you want to use another starfile just do this:

put a file called MyStarFe.l  in EW/
put a line in runStePar.sh called:

	./runstar.sh MyStar

Wait for a little while, then check StePar_results when it is done. There will be one line for each StePar iteration (by default 2).

Enjoy! Drop me a line if you need help.
