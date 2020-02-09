"""
To run without control
   python ./run.py nocontrol
To run with control
   python ./run.py
   python plot_opt.py
"""
import os
import sys
import numpy as np
import scipy.optimize as sop

def write_param(u):
    f = open('param.dat','w')
    f.write(str(initc) +"  initc\n")
    f.write(str(tstart)+"  tstart\n")
    f.write(str(tend)  +"  tend\n")
    f.write(str(u[0])  +"  u1\n")
    f.write(str(u[1])  +"  u2\n")
    f.close()

def cost_value(u):
    # save params to file
    write_param(u)
    # Run forward solver
    os.system("./fvm > forward.log")
    # Read cost function value
    f = open('cost.dat','r')
    c = f.read()
    c = float(c)
    return c

def cost_gradient(u):
    # save params to file
    write_param(u)
    # Run forward solver
    os.system("./adj > adjoint.log")
    # Read gradient value
    g = np.loadtxt('gradient.dat')
    return g

#------------------------------------------------------------------------------
nint  = 10
T     = 5.0
u     = np.array([0.0, 0.0])  # Initial control
initc = 0                     # Set initial condition internally

# Run without optimization
if len(sys.argv) > 1:
    if sys.argv[1]=="nocontrol":
        print("Running uncontrolled case ...")
        tstart, tend = 0.0, T
        c     = cost_value(u)
        os.system("cp sol.dat sol_nocontrol.dat")
        print("Cost = ", c)
        exit()

print("Running optimization ...")
t = np.linspace(0.0, T, nint+1, True)

u_all = np.zeros((nint,2))
c_all = np.zeros(nint)
for i in range(nint):
    tstart, tend = t[i], t[i+1]
    print('tstart, tend =', tstart, tend)
    res = sop.minimize(cost_value, u, \
                       method='BFGS', jac=cost_gradient, \
                       options={'disp': True})
    u = res.x; c = cost_value(u)
    print("   Control variable = ", u)
    print("   Objective func   = ", c)
    u_all[i,:] = u; c_all[i] = c;
    # Final solution becomes initial condition for next interval
    os.system("cp sol.dat init.dat")
    initc = 1 # From now on, read initial cond from file
    print("----------------------------------------------")

# Save final solution
os.system("cp sol.dat sol_control.dat")

# Save history to file
f = open('opt.dat','w')
for i in range(nint):
    print(u_all[i,:], c_all[i])
    st = str(u_all[i,0])+" "+str(u_all[i,1])+" "+str(c_all[i])+"\n"
    f.write(str(t[i])  +" "+st)
    f.write(str(t[i+1])+" "+st)
f.close()
