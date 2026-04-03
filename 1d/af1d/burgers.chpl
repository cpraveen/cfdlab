// Semi-discrete active flux
module PDE
{

proc flux(u)
{
   return 0.5 * u**2;
}

proc jacobian(u)
{
   return u;
}

}
