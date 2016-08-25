OBJS = delta_sigma_src/xi_nfw/xi_nfw.o delta_sigma_src/sigma_r/sigma_r.o delta_sigma_src/sigma_r/sigma_r_at_r.o delta_sigma_src/delta_sigma/delta_sigma.o delta_sigma_src/delta_sigma/delta_sigma_at_r.o delta_sigma_src/miscentered_sigma_r/miscentered_sigma_r.o delta_sigma_src/miscentered_sigma_r/miscentered_sigma_r_at_r.o delta_sigma_src/miscentered_delta_sigma/miscentered_delta_sigma.o delta_sigma_src/miscentered_delta_sigma/miscentered_delta_sigma_at_r.o delta_sigma_src/wrapper/wrapper.o

CC = gcc
ifdef SHARED
ifeq ($(SHARED),yes)
$(info Building shared library)
EXEC = Build_Delta_Sigma_Library.so
CFLAGS = -fPIC
OFLAGS = -shared 
#-W1,-soname=$(EXEC)
endif
else
$(info Building executable)
EXEC = main.exe
CFLAGS = 
OFLAGS = 
endif

#Note the paths to GSL
INCL = -I/home/tom/code/gsl/include/ -fopenmp -O3
LIBS = -lgsl -lgslcblas -L/home/tom/code/gsl/lib -lm -fopenmp -O3
.SUFFIXES : .c .o
%.o: %.c
	$(CC) $(CFLAGS) $(INCL) -c $< -o $@

$(EXEC): $(OBJS)
	$(CC) $(OFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

#$(OBJS): $(INCL)

.PHONY : clean

clean:
	rm -f $(OBJS) main.exe Build_Delta_Sigma.so
	rm -f *~
