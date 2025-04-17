PROG = ggmcalc
FC = gfortran
FCFLAGS = -Wall -O3 -Isrc -Jsrc -ffree-form -ffree-line-length-none -ffast-math -Ofast -fopenmp
LD = $(FC)
INSTALL_DIR = /opt/bin/
SRCS = src/nrtype.f90 \
       src/inout_mod.f90 \
       src/coordinates_mod.f90 \
       src/legendre_mod.f90 \
       src/c_2n_mod.f90 \
       src/gamma_mod.f90 \
       src/w_mod.f90 \
       src/u_mod.f90 \
       src/undulation_mod.f90 \
       src/height_anomaly_mod.f90 \
       src/gravity_disturbance_mod.f90 \
       src/gravity_anomaly_mod.f90 \
       src/progressbar_mod.f90 \
       src/date_sub.f90 \
       src/duration.f90 \
       src/color_mod.f90 \
       src/help_mod.f90 \
       src/ggmcalc_mod.f90
SRCP = src/main.f90
OBJS = $(SRCS:%.f90=%.o)
MODS = $(SRCS:%.f90=%.mod)
OBJP = $(SRCP:%.f90=%.o)
RM = rm -f

all: $(PROG)

$(PROG): $(OBJS) $(OBJP)
	$(LD) $(FCFLAGS) $^ -o $@

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	$(RM) $(PROG) $(OBJS) $(OBJP) $(MODS)
	
install:
	mkdir -p $(INSTALL_DIR)
	cp $(PROG) $(INSTALL_DIR)
