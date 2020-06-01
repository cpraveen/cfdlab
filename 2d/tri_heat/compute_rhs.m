function compute_rhs(t)

global n_ei n_eb eb rhs gradu v_ei v_eb nbr_ei nbr_eb norm_ei norm_eb
global mu cpen u

compute_gradient(t);

rhs(:) = 0;

for j=1:n_ei
   vl = v_ei(j,1);
   tl = nbr_ei(j,1);
   rhs(vl) = rhs(vl) + 0.5 * mu(tl) * dot(gradu(tl,:), norm_ei(j,:));

   vr = v_ei(j,2);
   tr = nbr_ei(j,2);
   rhs(vr) = rhs(vr) - 0.5 * mu(tr) * dot(gradu(tr,:), norm_ei(j,:));
end

% boundary flux
for j=1:n_eb
   vl = v_eb(j);
   tl = nbr_eb(j);
   rhs(vl) = rhs(vl) + 0.5 * mu(tl) * dot(gradu(tl,:), norm_eb(j,:));

   % boundary flux to two points on edge
   ve = eb(j,:);
   rhs(ve) = rhs(ve) + 0.5 * mu(tl) * dot(gradu(tl,:), norm_eb(j,:));

   % nitsche penalty
   u1 = boundary_value(t, j, 1);
   u2 = boundary_value(t, j, 2);
   ds = norm(norm_eb(j,:));
   he = ds;
   rhs(ve(1)) = rhs(ve(1)) + 0.5*(cpen/he)*mu(tl)*(u1 - u(ve(1)))*ds;
   rhs(ve(2)) = rhs(ve(2)) + 0.5*(cpen/he)*mu(tl)*(u2 - u(ve(2)))*ds;
end
