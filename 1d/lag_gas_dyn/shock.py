from pylab import *

gam = 1.4

rhor = 0.125
pr = 0.1
pl = 1.0

# Stationary shock relations
Ml = sqrt( (gam+1)/(2*gam)*(pr/pl - 1) + 1 )

rhol = rhor * (2 + (gam-1)*Ml**2) / ((gam+1)*Ml**2)
cl = sqrt(gam * pl / rhol)
ul = -cl * Ml

ur = rhol * ul / rhor

# Make right velocity zero
# Now the shock is moving

ul = ul - ur
ur = 0.0

print("rhol, rhor = ", rhol, rhor)
print("ul  , ur   = ", ul, ur)
print("pl  , pr   = ", pl, pr)

# Shock speed from RH condition
s1 = (rhor*ur - rhol*ul)/(rhor - rhol)
s2 = ((pr+rhor*ur**2) - (pl+rhol*ul**2))/(rhor*ur - rhol*ul)
El = pl/(gam-1) + 0.5*rhol*ul**2
Er = pr/(gam-1) + 0.5*rhor*ur**2
s3 = ((Er+pr)*ur - (El+pl)*ul)/(Er - El)
print("Shock speed = ",s1,s2,s3)
