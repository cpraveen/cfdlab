#
# fortran flags
#
FC       = gfortran

TARGETS  = scalar

#
# program dependencies
#

SRCF  = $(wildcard *.f)
OBJF  = $(patsubst %.f,%.o,$(SRCF))
scalar= $(OBJF)

#
# programs
#

ALL:  $(TARGETS)

scalar:  ${scalar}
		${FC} -o scalar ${scalar}

##############################################################################
# transformation rules
##############################################################################

.f.o:     ; ${FC} -Wall -O3 -c $*.f

##############################################################################
# clean things up
##############################################################################

clean:
	${RM} *.o $(TARGETS)
