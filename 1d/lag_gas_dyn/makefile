FC=gfortran

TARGETS = con noncon1 noncon2

all: $(TARGETS)

con: lag_gas_con.f
	$(FC) -Wall -fdefault-real-8 $< -o $@

noncon1: lag_gas_noncon1.f
	$(FC) -Wall -fdefault-real-8 $< -o $@

noncon2: lag_gas_noncon2.f
	$(FC) -Wall -fdefault-real-8 $< -o $@

clean:
	rm -f $(TARGETS)
