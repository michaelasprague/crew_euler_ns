DBFLAGS	=  #-g  #-fno-underscoring

ARCH = `uname -s`

# Linux, g95
FFLAGS	= -O2 -Wuninitialized  -fdefault-real-8 -fdefault-double-8 -Wall -ffpe-trap=zero

DEST	      = .

EXTHDRS	      = 

HDRS	      =

INSTALL	      = /etc/install

LD	      = gfortran # Linux, gcc

LDFLAGS       =  -O2 -Wuninitialized -fdefault-real-8 -fdefault-double-8 -Wall -ffpe-trap=zero

LIBS	      = 

MAKEFILE      = makefile

OBJS	      = \
                add.o\
                copy.o\
                deflux.o\
                driver.o\
                error.o\
                filename.o\
                generate_ix_iy.o\
                init_euler.o\
                max_speeds_euler.o\
                outputvtk.o\
                outputvtk_err.o\
                r_advect5.o\
                update_euler.o\
                update_pressure.o\
                vtk_scalar_add.o\
                vtk_scalar_out.o\
                vtk_vector_add.o\

PRINT	      = pr

PROGRAM       = euler2d

SHELL	      = /bin/sh

SRCS	      = \
                add.f\
                copy.f\
                deflux.f\
                driver.f\
                error.f\
                filename.f\
                generate_ix_iy.f\
                init_euler.f\
                max_speeds_euler.f\
                outputvtk.f\
                outputvtk_err.f\
                r_advect5.f\
                update_euler.f\
                update_pressure.f\
                vtk_scalar_add.f\
                vtk_scalar_out.f\
                vtk_vector_add.f\

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
