DBFLAGS	=  #-g  #-fno-underscoring

ARCH = `uname -s`

# Linux, g95
# -fdefault-real-8
# -fdefault-double-8
FFLAGS	= -O2 -Wuninitialized  -fdefault-real-8 -fdefault-double-8 -Wall -ffpe-trap=zero #-fopenmp -lgomp

DEST	      = .

EXTHDRS	      = 

HDRS	      =

INSTALL	      = /etc/install

LD	      = gfortran # Linux, gcc

LDFLAGS       = -O2 -Wuninitialized  -fdefault-real-8 -fdefault-double-8 -Wall -ffpe-trap=zero #-lgomp

LIBS	      = 

MAKEFILE      = makefile

OBJS	      = \
                A_vec.o\
                Hbar_vec.o\
                beta.o\
                bicgstab.o\
                bicgstab_p.o\
                calc_b.o\
                calc_error.o\
                cg_cond.o\
                cg_solve_ns.o\
                chol.o\
                cor_vel_un.o\
                divergence.o\
                divergence_rc.o\
                dp_centers.o\
                dp_edges.o\
                driver.o\
                filename.o\
                forcing.o\
                get_kovas_pbc.o\
                get_p_diag.o\
                get_pbc_navier_stokes.o\
                inflow_bcs.o\
                initialize.o\
                int_a_tilde.o\
                kovas.o\
                l3_solve.o\
                matlabvtk.o\
                outputvtk_ns.o\
                outputvtk_bound.o\
                p_mat_vec.o\
                pres_rhs.o\
                rhie_chow.o\
                sub_ave.o\
                u3_solve.o\
                uniform_bcs.o\
                update_navier_stokes.o\
                vel_Hmat_vec.o\
                vel_rhs.o\
                vtk_scalar_add_ns.o\
                vtk_scalar_out_ns.o\
                vtk_vector_add_ns.o\

PRINT	      = pr

PROGRAM       = navier_stokes_2d

SHELL	      = /bin/sh

SRCS	      = \
                A_vec.f\
                Hbar_vec.f\
                beta.f\
                bicgstab.f\
                bicgstab_p.f\
                calc_b.f\
                calc_error.f\
                cg_cond.f\
                cg_solve_ns.f\
                chol.f\
                cor_vel_un.f\
                divergence.f\
                divergence_rc.f\
                dp_centers.f\
                dp_edges.f\
                driver.f\
                filename.f\
                forcing.f\
                get_kovas_pbc.f\
                get_p_diag.f\
                get_pbc_navier_stokes.f\
                kovas.f\
                inflow_bcs.f\
                initialize.f\
                int_a_tilde.f\
                l3_solve.f\
                matlabvtk.f\
                outputvtk_ns.f\
                outputvtk_bound.f\
                p_mat_vec.f\
                pres_rhs.f\
                rhie_chow.f\
                sub_ave.f\
                u3_solve.f\
                uniform_bcs.f\
                update_navier_stokes.f\
                vel_Hmat_vec.f\
                vel_rhs.f\
                vtk_scalar_add_ns.f\
                vtk_scalar_out_ns.f\
                vtk_vector_add_ns.f\

SYSHDRS	      =

.f.o:
	$(LD) $(FFLAGS) $(DBFLAGS) -c $<

.F.o:
	$(LD) $(FFLAGS) $(DBFLAGS) -c $<

all:		$(PROGRAM)

linux-fort77:   linux-fort77-vars $(PROGRAM)

linux-g95:      linux-g95-vars $(PROGRAM)

linux-fort77-vars:
	@echo "Building for linux-fort77"
	@FFLAGS='-a -r8 -ext -O3 -funroll-loops'
	@LD='fort77'
	@LDFLAGS='-O3 -funroll-loops' # -pg

redo:		clean $(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS)
	$(LD) $(LDFLAGS) $(DBFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
	@echo done

clean:;		@rm -f $(OBJS) *.c core *~ \#*.*\# *.o
	@/usr/bin/strip $(PROGRAM)

clobber:;	@rm -f $(OBJS) $(PROGRAM) core tags *~ \#*.*\# Makefile.bak *.o

dep:;		@$(HOME)/bin/make_depend $(SRCS)

depend:;	@makedepend -f$(MAKEFILE) -- -- $(SRCS)

echo:;		@echo $(HDRS) $(SRCS)

index:;		@ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		@echo Installing $(PROGRAM) in $(DEST)
		@-strip $(PROGRAM)
		@if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(PROGRAM); $(INSTALL) -f $(DEST) $(PROGRAM)); fi

print:;		@$(PRINT) $(HDRS) $(SRCS)

tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)

update:		$(DEST)/$(PROGRAM)

$(DEST):	$(SRCS) $(LIBS) $(HDRS) $(EXTHDRS)
		@$(MAKE) -f $(MAKEFILE) ROOT=$(ROOT) DEST=$(DEST) install

# DO NOT DELETE THIS LINE -- make  depend  depends  on  it.
# driver.o: param.h
