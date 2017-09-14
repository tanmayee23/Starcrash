#! /usr/sbin/smake
#
MPID	  = /Users/tmg9722/include
FFTDIR    = /Users/tmg9722/Downloads/fftw-2.1.5
SFFTDIR   = $(FFTDIR)/../../lib
RFFTDIR   = $(FFTDIR)/rfftw
MPIDIR    = $(FFTDIR)/mpi
TMPIDIR   = $(FFTDIR)/mpi
INCL      = -I$(FFTDIR)/../include -I$(RFFTDIR) -I$(MPIDIR) -I$(TMPIDIR) -I$(FFTDIR)/fftw -I$(MPID)
FC        = gfortran
CC        = gcc
LD        = $(FC)
F77       = $(FC)
SWP       = swplist
RM        = /bin/rm -f
####### The -w flag gets rid of annoying compiler warnings
MP        = #-w
ABI       = #-n32  This flag slows down the code considerably!!!
ISA       = #-mips4
PROC      = ip27
ARCH      = $(MP) $(ABI) $(ISA)
#######OLEVEL    = -Ofast=$(PROC) This flag makes the code give nan's...
OLEVEL    = -O3
FOPTS     = #-OPT:IEEE_arithmetic=3
COPTS     = #-OPT:IEEE_arithmetic=3
FFLAGS    = $(ARCH) $(OLEVEL) $(FOPTS)
CFLAGS    = $(ARCH) $(OLEVEL) $(COPTS)
LIBS      = $(INCL) -L$(SFFTDIR) -lfftw -lfftw_mpi -lrfftw -lrfftw_mpi -lmpich -lfmpich
LDFLAGS   = $(ARCH) $(OLEVEL)
PROF      =

RFFT_LIBS = $(RFFTDIR)/librfftw.la -lm
LIBTOOL = $(SHELL) $(FFTDIR)/libtool
LINK = $(LIBTOOL) --mode=link

FOBJS     = advance.o balAV.o calcdivv.o claAV.o densplot.o \
                doqgrav.o dots.o getcurlv.o grav.o grav2.o gravrad.o \
                init.o kernels.o main.o nelist.o newAV.o output.o \
                poly.o ran1.o setup1ns.o setup2ns.o setup2q.o setupirr.o \
		setuphyper.o

COBJS     = fftw_f77.o
FFTOBJS   = $(MPIDIR)/fftwnd_mpi.o $(TMPIDIR)/transpose_mpi.o \
                $(TMPIDIR)/TOMS_transpose.o \
                $(RFFTDIR)/rfftwnd.o
OBJS      = $(FOBJS) $(COBJS) $(FFTOBJS)
EXEC      = sphagr

testinput : testinput.f
	$(FC) -o testinput testinput.f

b2pos	: bin2pos.f
	$(FC) -o b2pos bin2pos.f
b2pos2	: bin2pos2.f
	$(FC) -o b2pos2 bin2pos2.f

$(EXEC):  $(OBJS) 
	$(LINK) $(FC) -Wall -o $@ $(PROF) $(FFLAGS) $(QV_LOPT) $(OBJS) $(LIBS)

fftw_f77.o : fftw_f77.c spha.h
	$(CC) -c $(CFLAGS) fftw_f77.c $(INCL)

clean:
	$(RM) $(EXEC) $(OBJS)

.SUFFIXES:
.SUFFIXES: .o .c .f

.f.o: spha.h
	$(FC) -Wall -c $(FFLAGS) $(QV_OPT) $(DEFINES) $< -I$(MPID)

