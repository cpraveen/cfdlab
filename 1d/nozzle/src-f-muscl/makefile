# Location of tapenade executable
# TAPENADE_HOME must be defined, probably in your bashrc
# Else define it first
TPN  = ${TAPENADE_HOME}/bin/tapenade

#
# fortran flags
#
 
FC = gfortran
CC = gcc
RM = /bin/rm
FFLAGS   = -O3
FFLAGS2  = -DCOMPLEX
CFLAGS   = -O3
CPPFLAGS = -O3
LFLAGS   = -O3

TARGETS =  flo adj error

#
# transformation rules
#

.F.o:     ; ${RM} -f $*.o; ${FC} ${FFLAGS} -c $*.F
.c.o:     ; ${RM} -f $*.o; ${CC} ${CFLAGS} -c $*.c

#
# program dependencies
#

noz_flo = nozzle.o common.o residu.o routines.o

noz_adj = noz_adj.o common.o routines.o     \
          adBuffer.o adStack.o              \
          flux_ausm_bq.o  flux_ausm_ba.o    \
          flux_roe_bq.o  flux_roe_ba.o    \
          flux_kfvs_bq.o  flux_kfvs_ba.o    \
          flux_lf_bq.o  flux_lf_ba.o        \
			 flux_in_bq.o flux_in_ba.o         \
			 flux_out_bq.o flux_out_ba.o       \
			 flux_out2_bq.o flux_out2_ba.o       \
			 source_term_bq.o source_term_ba.o \
			 dmax1_bq.o dmin1_bq.o errf_bq.o   \
			 costfunc_bq.o costfunc2_bq.o \
			 muscl_bq.o nozzle_x.o flux_in_steger_bq.o \
			 flux_in_steger_bmach.o

error = error.o common.o routines.o nozzle_x.o InterpolateFlux.o \
        InterpolateFlux_x.o InterpolateSource.o adStack.o adBuffer.o 

#
# programs
#

ALL:		$(TARGETS)

flo:	${noz_flo}
		${FC} ${FFLAGS} -o noz_flo ${noz_flo}

adj: ${noz_adj}
		${FC} ${FFLAGS} -o noz_adj ${noz_adj}
		${RM} -f *~ *.msg

error: ${error}
		${FC} ${FFLAGS} -o error ${error}

#
# Real version of nonlinear routines
#

InterpolateFlux.o:	   InterpolateFlux.f
		${FC} ${FFLAGS} -c InterpolateFlux.f

InterpolateSource.o:	   InterpolateSource.f
		${FC} ${FFLAGS} -c InterpolateSource.f

routines.o:	   routines.f
		${FC} ${FFLAGS} -c routines.f

noz_adj.o: noz_adj.f flux_in_steger_bmach.o
		${FC} ${FFLAGS} -c noz_adj.f

##############################################################################
# Adjoint versions of nonlinear routines
##############################################################################

dmax1_bq.o:	routines.f
		${TPN} -backward                     \
		       -head       dmax1             \
		       -vars       "a b"             \
		       -outvars    "dmax1"           \
		       -difffuncname "_bq"           \
		        routines.f;
		${FC} ${FFLAGS} -c dmax1_bq.f;
#-----------------------------------------------------------------------------
dmin1_bq.o:	routines.f
		${TPN} -backward                     \
		       -head       dmin1             \
		       -vars       "a b"             \
		       -outvars    "dmin1"           \
		       -difffuncname "_bq"           \
		        routines.f;
		${FC} ${FFLAGS} -c dmin1_bq.f;
#-----------------------------------------------------------------------------
# Error function used by KFVS
errf_bq.o:	routines.f
		${TPN} -backward                     \
		       -head       errf              \
		       -vars       "x"               \
		       -outvars    "x errf"          \
		       -difffuncname "_bq"           \
		        routines.f;
		${FC} ${FFLAGS} -c errf_bq.f;
#-----------------------------------------------------------------------------
nozzle_ba.o:	routines.f
		${TPN} -backward                     \
		       -head       nozzle            \
		       -vars       "c d"             \
		       -outvars    "nozzle"          \
		       -difffuncname "_ba"           \
		        routines.f;
		${FC} ${FFLAGS} -c nozzle_ba.f;
#-----------------------------------------------------------------------------
nozzle_x.o:	routines.f
		${TPN} -backward                     \
		       -head       nozzle            \
		       -vars       "x"               \
		       -outvars    "nozzle"          \
		       -difffuncname "_x"            \
		        routines.f;
		${FC} ${FFLAGS} -c nozzle_x.f;
#-----------------------------------------------------------------------------
# Differentiation of flux routines
# This takes care of all the flux routines
# They must be named "flux_name" where name is the name of the flux function
# for example flux_ausm, flux_kfvs, flux_lf, etc.
flux_%_bq.o:	routines.f
		${TPN} -backward                     \
		       -head       flux_$*           \
		       -vars       "ql qr resl resr" \
		       -outvars    "ql qr resl resr" \
		       -difffuncname "_bq"           \
		        routines.f;
		${FC} ${FFLAGS} -c flux_$*_bq.f

flux_%_ba.o:	routines.f
		${TPN} -backward                 \
		       -head       flux_$*       \
		       -vars       "a resl resr" \
		       -outvars    "a resl resr" \
		       -difffuncname "_ba"       \
		        routines.f;
		${FC} ${FFLAGS} -c flux_$*_ba.f;

#-----------------------------------------------------------------------------
# Flux for inflow face
flux_in_bq.o:	routines.f
		${TPN} -backward             \
		       -head       flux_in   \
		       -vars       "res"   \
		       -outvars    "res"   \
		       -difffuncname "_bq"   \
		        routines.f;
		${FC} ${FFLAGS} -c flux_in_bq.f;

flux_in_ba.o:	routines.f
		${TPN} -backward             \
		       -head       flux_in   \
		       -vars       "a res"   \
		       -outvars    "a res"   \
		       -difffuncname "_ba"   \
		        routines.f;
		${FC} ${FFLAGS} -c flux_in_ba.f;

#-----------------------------------------------------------------------------
# Flux for inflow face
flux_in_steger_bq.o:	routines.f
		${TPN} -backward             \
		       -head       flux_in_steger   \
		       -vars       "q res"   \
		       -outvars    "q res"   \
		       -difffuncname "_bq"   \
		        routines.f;
		${FC} ${FFLAGS} -c flux_in_steger_bq.f;

# Flux for inflow face
flux_in_steger_bmach.o:	routines.f
		${TPN} -backward             \
		       -head       flux_in_steger   \
		       -vars       "minf res"   \
		       -outvars    "res"   \
		       -difffuncname "_bmach"   \
		        routines.f;
		${FC} ${FFLAGS} -c flux_in_steger_bmach.f;

#-----------------------------------------------------------------------------
# Flux for outflow face
flux_out_bq.o:	routines.f
		${TPN} -backward                 \
		       -head       flux_out      \
		       -vars       "ql qll res"  \
		       -outvars    "ql qll res"  \
		       -difffuncname "_bq"       \
		        routines.f;
		${FC} ${FFLAGS} -c flux_out_bq.f;

flux_out_ba.o:	routines.f
		${TPN} -backward                 \
		       -head       flux_out      \
		       -vars       "a res"       \
		       -outvars    "a res"       \
		       -difffuncname "_ba"       \
		        routines.f;
		${FC} ${FFLAGS} -c flux_out_ba.f;

flux_out2_bq.o:	routines.f
		${TPN} -backward                 \
		       -head       flux_out2     \
		       -vars       "ql res"  \
		       -outvars    "ql res"  \
		       -difffuncname "_bq"       \
		        routines.f;
		${FC} ${FFLAGS} -c flux_out2_bq.f;

flux_out2_ba.o:	routines.f
		${TPN} -backward                 \
		       -head       flux_out2     \
		       -vars       "a res"       \
		       -outvars    "a res"       \
		       -difffuncname "_ba"       \
		        routines.f;
		${FC} ${FFLAGS} -c flux_out2_ba.f;

#-----------------------------------------------------------------------------
source_term_bq.o:	routines.f
		${TPN} -backward                 \
		       -head       source_term   \
		       -vars       "q res"       \
		       -outvars    "q res"       \
		       -difffuncname "_bq"       \
		        routines.f;
		${FC} ${FFLAGS} -c source_term_bq.f;

source_term_ba.o:	routines.f
		${TPN} -backward                 \
		       -head       source_term   \
		       -vars       "al ar res"   \
		       -outvars    "al ar res"   \
		       -difffuncname "_ba"       \
		        routines.f;
		${FC} ${FFLAGS} -c source_term_ba.f;

#-----------------------------------------------------------------------------
costfunc_bq.o:	routines.f
		${TPN} -backward             \
		       -head       costfunc  \
		       -vars       "q cost"  \
		       -outvars    "q cost"  \
		       -difffuncname "_bq"   \
		        routines.f;
		${FC} ${FFLAGS} -c costfunc_bq.f;

#-----------------------------------------------------------------------------
costfunc2_bq.o:	routines.f
		${TPN} -backward             \
				 -head       costfunc2 \
				 -vars       "q cost"  \
				 -outvars    "q cost"  \
				 -difffuncname "_bq"   \
				  routines.f;
		${FC} ${FFLAGS} -c costfunc2_bq.f;

#-----------------------------------------------------------------------------
muscl_bq.o:	routines.f
		${TPN} -backward             \
				 -head       muscl     \
				 -vars       "qjm1 qj qjp1 qjph"  \
				 -outvars    "qjm1 qj qjp1 qjph"  \
				 -difffuncname "_bq"   \
				  routines.f;
		${FC} ${FFLAGS} -c muscl_bq.f;


#-----------------------------------------------------------------------------
InterpolateFlux_x.o:	InterpolateFlux.f
		${TPN} -backward                     \
		       -head       InterpolateFlux   \
		       -vars       "xgauss"          \
				 -outvars    "flux"            \
		       -difffuncname "_x"           \
		        InterpolateFlux.f;
		${FC} ${FFLAGS} -c InterpolateFlux_x.f;
#-----------------------------------------------------------------------------
# Flux for inflow face
residu_bq.o:	residu.f
		${TPN} -backward             \
		       -head       residu   \
		       -vars       "q res"   \
		       -outvars    "q res"   \
		       -difffuncname "_bq"   \
				 -output residu \
		        residu.f;
		${FC} ${FFLAGS} -c residu_bq.f;

#-----------------------------------------------------------------------------
adBuffer.o:	$(TAPENADE_HOME)/ADFirstAidKit/adBuffer.f
		${FC} ${FFLAGS} -c $<

adStack.o:	$(TAPENADE_HOME)/ADFirstAidKit/adStack.c
		${CC} ${FFLAGS} -c $<

##############################################################################
# clean things up
##############################################################################

clean:	
	${RM} -f *.o *.msg *~ noz_flo noz_adj error

allclean:	
	@make clean
	${RM} -f *_bq.f *_ba.f *_bmach.f *_bq.inc *_ba.inc *_bmach.inc *_bx.inc
	${RM} -f *_bx.f *_x.f
