ifort -c  -lstdc++ -fPIC  -O3 *.f   # compile all of the .f files to produce .o files
ar rv libblas.a *.o    #  combine the .o files into a library
