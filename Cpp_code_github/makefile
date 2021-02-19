#~ CFLAGS= -std=c++11 -pg -g -c  
CFLAGS= -std=c++11 -O3 -c  

#~ CFLAGS= -std=c++0x -O3 -c  
#LDFLAGS= -pg  CFLAGS= -std=c++11 -pg -g -Wall -Wextra -c  
#LDFLAGS= -g
LDFLAGS= -O3


#~ NFFT_DWNLD = /home/jaco/nfft-3.5.1
#~ NFFT_ROOT = /home/jaco/nfft3_install
#~ FFTW_ROOT = /home/jaco/fftw3_install

NFFT_DWNLD =$${HOME}/nfft-3.5.1
NFFT_ROOT = $${HOME}/nfft3_install
FFTW_ROOT = $${HOME}/fftw3_install

libdir = ${NFFT_ROOT}/lib 
#libdir2 = ${NFFT_ROOT}/.libs 
libdir_FFTW = ${FFTW_ROOT}/lib
#libdir_FFTW_libs = ${FFTW_ROOT}/.libs 
includedir = ${NFFT_ROOT}/include 
includedir_config = ${NFFT_DWNLD}/include 
includedir_FFTW = ${FFTW_ROOT}/include

CONFIG_HEADER = $(top_builddir)/include/config.h

CPPFLAGS        = -I. \
                  -I$(NFFT_ROOT) \
                  -I$(includedir) \
                  -I$(includedir_config) \
                  -I$(includedir_FFTW)
LDFLAGS         += -L$(NFFT_ROOT) \
                  -L$(libdir) \
                  -L$(libdir_FFTW) 
                  #-L$(libdir2) \
                  #-L$(libdir_FFTW_libs) \
                  
#RFLAGS            = -Wl,-rpath,$(libdir2) -Wl,-rpath,$(libdir_FFTW_libs)
RFLAGS            = -Wl,-rpath,$(libdir) -Wl,-rpath,$(libdir_FFTW)
#-rpath,$(NFFT_ROOT),-rpath,$(libdir),-rpath,$(libdir_FFTW),
                 
                  
LDLIBS          = -lfftw3 -lm  -lnfft3 


coarse_grained_coevo: classi.o coarse_grained_coevo.o fastsum_NFFT_simple.o
	g++ $(RFLAGS)  $(LDFLAGS) classi.o coarse_grained_coevo.o fastsum_NFFT_simple.o $(LOADLIBES) $(LDLIBS) -o coarse_grained_coevo
classi.o: classi.cpp classi.h fastsum_NFFT_simple.h
	g++ $(CFLAGS) $(CPPFLAGS) classi.cpp
fastsum_NFFT_simple.o: fastsum_NFFT_simple.cpp fastsum_NFFT_simple.h 
	g++ $(CFLAGS)  $(CPPFLAGS) fastsum_NFFT_simple.cpp
coarse_grained_coevo.o: coarse_grained_coevo.cpp classi.h fastsum_NFFT_simple.h
	g++ $(CFLAGS) $(CPPFLAGS) coarse_grained_coevo.cpp 

clean: 
	@rm -f *.o coarse_grained_coevo 

