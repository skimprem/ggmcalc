PROG = ggmcalc
FC = gfortran
FCFLAGS = -Wall -O3 -Isrc -Jsrc -ffree-form -ffree-line-length-none -ffast-math -Ofast -fopenmp
LD = $(FC)
INSTALL_DIR = /opt/bin/
SRCS = src/nrtype.f \
       src/coordinates_mod.f \
       src/legendre_mod.f \
       src/c_2n_mod.f \
       src/gamma_mod.f \
       src/w_mod.f \
       src/u_mod.f \
       src/undulation_mod.f \
       src/height_anomaly_mod.f \
       src/gravity_disturbance_mod.f \
       src/gravity_anomaly_mod.f \
       src/progressbar_mod.f \
       src/date_sub.f \
       src/duration.f \
       src/color_mod.f \
       src/help_mod.f \
       src/ggmcalc_mod.f
SRCP = src/main.f
OBJS = $(SRCS:%.f=%.o)
MODS = $(SRCS:%.f=%.mod)
OBJP = $(SRCP:%.f=%.o)
RM = rm -f

all: $(PROG)

$(PROG): $(OBJS) $(OBJP)
	$(LD) $(FCFLAGS) $^ -o $@

.f.o:
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	$(RM) $(PROG) $(OBJS) $(OBJP) $(MODS)
	
install:
	mkdir -p $(INSTALL_DIR)
	cp $(PROG) $(INSTALL_DIR)
