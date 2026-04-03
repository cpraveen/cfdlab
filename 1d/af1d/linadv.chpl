// Semi-discrete active flux
module PDE
{

proc flux(u)
{
   return jacobian(u) * u;
}

proc jacobian(u)
{
   return 1.0;
}

}
