ADOLC_HOME = /opt/adolc

#
# fortran flags
#
CC       = gcc
CPP      = g++
FC       = gfortran
ADOLCINC = -I$(ADOLC_HOME)/include/adolc
CCFLAGS  = ${CFLAGS}
CPPFLAGS = ${CFLAGS} ${ADOLCINC} -DNO_KFVS_

INDENT   = indent -br -brs -npcs -npsl -nsaf -nsai -nsaw -cli3 -l78 -i3 -ts0

TARGETS = noz_flo noz_adj

#
# program dependencies
#

FLO_OBJC  = noz_flo.o \
		  		routines.o \
		  		common_cc.o
noz_flo= $(FLO_OBJC)

ADJ_OBJC  = noz_adj.o \
            routines_q.o \
            routines_x.o \
            common_cpp.o
noz_adj= $(ADJ_OBJC)

#
# programs
#

ALL:  $(TARGETS)

noz_flo:  ${noz_flo}
		${CC} ${noz_flo} -o noz_flo -lm

noz_adj:  ${noz_adj}
	${CPP} -o noz_adj ${noz_adj} -L$(ADOLC_HOME)/lib -ladolc -lm

##############################################################################
# transformation rules
##############################################################################
noz_flo.o: noz_flo.c
	${CC} ${CFLAGS} -DDIF0 -c $*.c

common_cc.o: common.c
	${CC} ${CFLAGS} -DDIF0 -c $< -o $*.o

routines.o: routines.c
	${CC} ${CFLAGS} -DDIF0 -c $*.c

##############################################################################
noz_adj.o: noz_adj.c
	${CPP} ${CPPFLAGS} -DDIF0 -c $*.c

common_cpp.o: common.c
	${CPP} ${CPPFLAGS} -DDIF0 -c $< -o $*.o

routines_q.o: routines.c
	${CPP} ${CPPFLAGS} -DDIFQ -c $< -o $*.o

routines_x.o: routines.c
	${CPP} ${CPPFLAGS} -DDIFX -c $< -o $*.o
##############################################################################
indent:
	${INDENT} *.c
	${RM} *.c~
##############################################################################
# clean things up
##############################################################################

clean:
	${RM} *.o $(TARGETS)
