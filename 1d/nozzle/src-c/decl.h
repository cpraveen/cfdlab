#ifndef MYDECL_H_

#define MYDECL_H_ 1

void init(void);

void read_shape(int, double*, double*);

void read_flow(int, double**);

void time_step(int, double**, double*);

void save_old(int, double**, double**);

//void residu(int, int, double*, double**, double**);

void sol_residue(int, double**, double*);

//void costfunc(double*, double, double*);

void result(int, int, double*, double**);

#endif
