import plotUtil.simplePlot;

// setup simulation constants
config const nx = 41,           // Number of Finite Difference points
             nt = 20,           // Number of time steps
             nu = 0.3,          // The fluid's viscosity
             sigma = 0.2;       // A simulation stability parameter

const dx = 2.0 / (nx - 1),      // Spatial step-size (Delta-x)
      dt = sigma * dx**2 / nu;  // Temporal step-size (Delta-t)

// define initial conditions
var u : [0..<nx] real = 1;
u[(0.5 / dx):int..(1.0 / dx):int] = 2;

// plot initial conditions
writeln("'u' initial conditions:");
simplePlot(u, minY = 1, maxY = 2, minX=0, maxX=2);

// run FD simulation
var un = u;       // temporary copy of 'u'
for n in 1..nt {  // simulate 'nt' time steps
    u <=> un;     // swap 'u' and 'un'
    // apply FD equation
    forall i in 1..<(nx-1) do
        u[i] = un[i] + nu * dt / dx**2 * (un[i+1] - 2 * un[i] + un[i-1]);
}

// plot results
writeln("'u' after ", nt, " iterations:");
simplePlot(u, minY = 1, maxY = 2, minX=0, maxX=2);
