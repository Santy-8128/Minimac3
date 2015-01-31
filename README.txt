This is a sample repository setup for a program that uses the statgen library.
It is a starting point you can use for creating your own program.

Update Makefile, replacing any occurrances of SAMPLE_PROGRAM with your program name.
Update EXE and add your own cpp/h files to the appropriate lines of the Makefile.

--------------------------------------------------------------------------------
To use git to clone the required statgen library:
  make cloneLib
Next, to build libStatGen and this program:
  make
To install:
  make install INSTALLDIR=pathToInstall
