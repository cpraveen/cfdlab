#ifndef MYROUTINES_H_

#define MYROUTINES_H_ 1

void flux_in(REALX, REALA*);

void flux_ausm(REALX, REALQ*, REALQ*, REALA*, REALA*);
#ifndef NO_KFVS_
void flux_kfvs(REALX, REALQ*, REALQ*, REALA*, REALA*);
#endif
void flux_lf(REALX, REALQ*, REALQ*, REALA*, REALA*);

void flux_out(REALX, REALQ*, REALQ*, REALA*);

void source_term(REALX, REALX, REALQ*, REALA*);

#endif
