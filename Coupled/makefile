#### makefile automatically created by makenek 01/16/2012 14:41:41 ###
BINNAME=euler_ns
SRC_C:=src
SRC_E:=../Euler
SRC_NS:=../Navier_Stokes
OBJDIR=obj
F77:=gfortran
CC:=gcc
P:=-O1 -Wuninitialized -fdefault-real-8 -fdefault-double-8 -Wall -ffpe-trap=zero
USR:=
USR_LFLAGS:=

################################################################################

# CORE##########################################################################
COUPLED       = \
                add_ke.o\
                add_ke3.o\
                add_value.o\
                calc_energy_e.o\
                calc_energy_n.o\
                copy_e_to_n.o\
                driver.o\
                error_ns.o\
                error_ns_ke.o\
                get_conv_bcs.o\
                get_max_div.o\
                get_navier_bcs.o\
                grab_bound_vels.o\
                gradient_u.o\
                init_navier_stokes.o\
                proj_cg_solve.o\
                proj_mat_vec.o\
                projection.o\
                put_bound_vels.o\
                overwrite_n_to_e.o\

EULER         = \
                add.o\
                copy.o\
                deflux.o\
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
                vtk_vector_add.o\
                vtk_scalar_add.o\
                vtk_scalar_out.o

NAVIER_STOKES = \
                A_vec.o\
                Hbar_vec.o\
                beta.o\
                bicgstab.o\
                bicgstab_p.o\
                calc_b.o\
                cg_cond.o\
                cg_solve_ns.o\
                chol.o\
                cor_vel_un.o\
                divergence.o\
                divergence_rc.o\
                dp_centers.o\
                dp_edges.o\
                forcing.o\
                get_p_diag.o\
                get_pbc_navier_stokes.o\
                initialize.o\
                int_a_tilde.o\
                l3_solve.o\
                matlabvtk.o\
                outputvtk_ns.o\
                outputvtk_bound.o\
                p_mat_vec.o\
                pres_rhs.o\
                rhie_chow.o\
                sub_ave.o\
                u3_solve.o\
                update_navier_stokes.o\
                vel_Hmat_vec.o\
                vel_rhs.o\
                vtk_scalar_add_ns.o\
                vtk_scalar_out_ns.o\
                vtk_vector_add_ns.o

################################################################################
NOBJS_C = $(patsubst %,$(OBJDIR)/%,$(COUPLED))

NOBJS_E = $(patsubst %,$(OBJDIR)/%,$(EULER))

NOBJS_NS = $(patsubst %,$(OBJDIR)/%,$(NAVIER_STOKES))

NOBJS  = $(NOBJS_C) $(NOBJS_E) $(NOBJS_NS)

FL   = $(P) -I$(SRC_C) -I$(SRC_E) -I$(SRC_NS)

################################################################################
all : buildit

objdir: 
	@mkdir $(OBJDIR) 2>/dev/null; cat /dev/null 

buildit: objdir $(NOBJS)
	$(F77) -o ${BINNAME} $(NOBJS) $(lFLAGS)
	@if test -f ${BINNAME}; then \
	echo "#############################################################"; \
	echo "#                  Compilation successful!                  #"; \
	echo "#############################################################"; \
        size ${BINNAME}; \
        echo ""; \
	else \
	echo -e "\033[1;31;38m" "ERROR: Compilation failed!"; \
	echo -e "\033[0m"; \
	fi
clean:
	rm -rf ./obj ${BINNAME}

# CORE     ##################################################################

# COUPLER ROUTINES
$(OBJDIR)/add_ke.o              :$(SRC_C)/add_ke.f;               $(F77) -c $(FL) $< -o $@
$(OBJDIR)/add_ke3.o             :$(SRC_C)/add_ke3.f;              $(F77) -c $(FL) $< -o $@
$(OBJDIR)/add_value.o           :$(SRC_C)/add_value.f;            $(F77) -c $(FL) $< -o $@
$(OBJDIR)/calc_energy_e.o       :$(SRC_C)/calc_energy_e.f;        $(F77) -c $(FL) $< -o $@
$(OBJDIR)/calc_energy_n.o       :$(SRC_C)/calc_energy_n.f;        $(F77) -c $(FL) $< -o $@
$(OBJDIR)/copy_e_to_n.o         :$(SRC_C)/copy_e_to_n.f;          $(F77) -c $(FL) $< -o $@
$(OBJDIR)/driver.o              :$(SRC_C)/driver.f;               $(F77) -c $(FL) $< -o $@
$(OBJDIR)/error_ns.o            :$(SRC_C)/error_ns.f;             $(F77) -c $(FL) $< -o $@
$(OBJDIR)/error_ns_ke.o         :$(SRC_C)/error_ns_ke.f;          $(F77) -c $(FL) $< -o $@
$(OBJDIR)/get_conv_bcs.o        :$(SRC_C)/get_conv_bcs.f;          $(F77) -c $(FL) $< -o $@
$(OBJDIR)/get_max_div.o         :$(SRC_C)/get_max_div.f;          $(F77) -c $(FL) $< -o $@
$(OBJDIR)/get_navier_bcs.o      :$(SRC_C)/get_navier_bcs.f;       $(F77) -c $(FL) $< -o $@
$(OBJDIR)/grab_bound_vels.o     :$(SRC_C)/grab_bound_vels.f;      $(F77) -c $(FL) $< -o $@
$(OBJDIR)/gradient_u.o          :$(SRC_C)/gradient_u.f;           $(F77) -c $(FL) $< -o $@
$(OBJDIR)/init_navier_stokes.o  :$(SRC_C)/init_navier_stokes.f;   $(F77) -c $(FL) $< -o $@
$(OBJDIR)/proj_cg_solve.o       :$(SRC_C)/proj_cg_solve.f;        $(F77) -c $(FL) $< -o $@
$(OBJDIR)/proj_mat_vec.o        :$(SRC_C)/proj_mat_vec.f;         $(F77) -c $(FL) $< -o $@
$(OBJDIR)/projection.o          :$(SRC_C)/projection.f;           $(F77) -c $(FL) $< -o $@
$(OBJDIR)/put_bound_vels.o      :$(SRC_C)/put_bound_vels.f;       $(F77) -c $(FL) $< -o $@
$(OBJDIR)/overwrite_n_to_e.o    :$(SRC_C)/overwrite_n_to_e.f;     $(F77) -c $(FL) $< -o $@


# EULER ROUTINES
$(OBJDIR)/add.o                   :$(SRC_E)/add.f;               $(F77) -c $(FL) $< -o $@
$(OBJDIR)/copy.o                  :$(SRC_E)/copy.f;              $(F77) -c $(FL) $< -o $@
$(OBJDIR)/deflux.o                :$(SRC_E)/deflux.f;            $(F77) -c $(FL) $< -o $@
$(OBJDIR)/error.o                 :$(SRC_E)/error.f;             $(F77) -c $(FL) $< -o $@
$(OBJDIR)/filename.o              :$(SRC_E)/filename.f;          $(F77) -c $(FL) $< -o $@
$(OBJDIR)/generate_ix_iy.o        :$(SRC_E)/generate_ix_iy.f;    $(F77) -c $(FL) $< -o $@
$(OBJDIR)/init_euler.o            :$(SRC_E)/init_euler.f;        $(F77) -c $(FL) $< -o $@
$(OBJDIR)/max_speeds_euler.o      :$(SRC_E)/max_speeds_euler.f;  $(F77) -c $(FL) $< -o $@
$(OBJDIR)/outputvtk.o             :$(SRC_E)/outputvtk.f;         $(F77) -c $(FL) $< -o $@
$(OBJDIR)/outputvtk_err.o         :$(SRC_E)/outputvtk_err.f;     $(F77) -c $(FL) $< -o $@
$(OBJDIR)/r_advect5.o             :$(SRC_E)/r_advect5.f;         $(F77) -c $(FL) $< -o $@
$(OBJDIR)/update_euler.o          :$(SRC_E)/update_euler.f;      $(F77) -c $(FL) $< -o $@
$(OBJDIR)/update_pressure.o       :$(SRC_E)/update_pressure.f;   $(F77) -c $(FL) $< -o $@
$(OBJDIR)/vtk_vector_add.o        :$(SRC_E)/vtk_vector_add.f;    $(F77) -c $(FL) $< -o $@
$(OBJDIR)/vtk_scalar_add.o        :$(SRC_E)/vtk_scalar_add.f;    $(F77) -c $(FL) $< -o $@
$(OBJDIR)/vtk_scalar_out.o        :$(SRC_E)/vtk_scalar_out.f;    $(F77) -c $(FL) $< -o $@

# NAVIER-STOKES ROUTINES
$(OBJDIR)/A_vec.o                 :$(SRC_NS)/A_vec.f;                 $(F77) -c $(FL) $< -o $@
$(OBJDIR)/Hbar_vec.o              :$(SRC_NS)/Hbar_vec.f;              $(F77) -c $(FL) $< -o $@
$(OBJDIR)/beta.o                  :$(SRC_NS)/beta.f;                  $(F77) -c $(FL) $< -o $@
$(OBJDIR)/bicgstab.o              :$(SRC_NS)/bicgstab.f;              $(F77) -c $(FL) $< -o $@
$(OBJDIR)/bicgstab_p.o            :$(SRC_NS)/bicgstab_p.f;            $(F77) -c $(FL) $< -o $@
$(OBJDIR)/calc_b.o                :$(SRC_NS)/calc_b.f;                $(F77) -c $(FL) $< -o $@
$(OBJDIR)/cg_cond.o               :$(SRC_NS)/cg_cond.f;               $(F77) -c $(FL) $< -o $@
$(OBJDIR)/cg_solve_ns.o           :$(SRC_NS)/cg_solve_ns.f;           $(F77) -c $(FL) $< -o $@
$(OBJDIR)/chol.o                  :$(SRC_NS)/chol.f;                  $(F77) -c $(FL) $< -o $@
$(OBJDIR)/cor_vel_un.o            :$(SRC_NS)/cor_vel_un.f;            $(F77) -c $(FL) $< -o $@
$(OBJDIR)/divergence.o            :$(SRC_NS)/divergence.f;            $(F77) -c $(FL) $< -o $@
$(OBJDIR)/divergence_rc.o         :$(SRC_NS)/divergence_rc.f;         $(F77) -c $(FL) $< -o $@
$(OBJDIR)/dp_centers.o            :$(SRC_NS)/dp_centers.f;            $(F77) -c $(FL) $< -o $@
$(OBJDIR)/dp_edges.o              :$(SRC_NS)/dp_edges.f;              $(F77) -c $(FL) $< -o $@
$(OBJDIR)/forcing.o               :$(SRC_NS)/forcing.f;               $(F77) -c $(FL) $< -o $@
$(OBJDIR)/get_pbc_navier_stokes.o :$(SRC_NS)/get_pbc_navier_stokes.f; $(F77) -c $(FL) $< -o $@
$(OBJDIR)/get_p_diag.o            :$(SRC_NS)/get_p_diag.f;            $(F77) -c $(FL) $< -o $@
$(OBJDIR)/initialize.o            :$(SRC_NS)/initialize.f;             $(F77) -c $(FL) $< -o $@
$(OBJDIR)/int_a_tilde.o           :$(SRC_NS)/int_a_tilde.f;           $(F77) -c $(FL) $< -o $@
$(OBJDIR)/l3_solve.o              :$(SRC_NS)/l3_solve.f;              $(F77) -c $(FL) $< -o $@
$(OBJDIR)/matlabvtk.o             :$(SRC_NS)/matlabvtk.f;             $(F77) -c $(FL) $< -o $@
$(OBJDIR)/outputvtk_ns.o          :$(SRC_NS)/outputvtk_ns.f;          $(F77) -c $(FL) $< -o $@
$(OBJDIR)/outputvtk_bound.o       :$(SRC_NS)/outputvtk_bound.f;       $(F77) -c $(FL) $< -o $@
$(OBJDIR)/p_mat_vec.o             :$(SRC_NS)/p_mat_vec.f;             $(F77) -c $(FL) $< -o $@
$(OBJDIR)/pres_rhs.o              :$(SRC_NS)/pres_rhs.f;              $(F77) -c $(FL) $< -o $@
$(OBJDIR)/rhie_chow.o             :$(SRC_NS)/rhie_chow.f;             $(F77) -c $(FL) $< -o $@
$(OBJDIR)/sub_ave.o               :$(SRC_NS)/sub_ave.f;             $(F77) -c $(FL) $< -o $@
$(OBJDIR)/u3_solve.o              :$(SRC_NS)/u3_solve.f;              $(F77) -c $(FL) $< -o $@
$(OBJDIR)/update_navier_stokes.o  :$(SRC_NS)/update_navier_stokes.f;  $(F77) -c $(FL) $< -o $@
$(OBJDIR)/vel_Hmat_vec.o          :$(SRC_NS)/vel_Hmat_vec.f;          $(F77) -c $(FL) $< -o $@
$(OBJDIR)/vel_rhs.o               :$(SRC_NS)/vel_rhs.f;               $(F77) -c $(FL) $< -o $@
$(OBJDIR)/vtk_vector_add_ns.o     :$(SRC_NS)/vtk_vector_add_ns.f;     $(F77) -c $(FL) $< -o $@
$(OBJDIR)/vtk_scalar_add_ns.o     :$(SRC_NS)/vtk_scalar_add_ns.f;     $(F77) -c $(FL) $< -o $@
$(OBJDIR)/vtk_scalar_out_ns.o     :$(SRC_NS)/vtk_scalar_out_ns.f;     $(F77) -c $(FL) $< -o $@
                
