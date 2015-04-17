function res = residual(u,w)

np = length(u);
nf = np + 1;
res = zeros(size(u));

for f=2:nf-1
   wf = 0.5*(w(f-1) + w(f));
   flux = god_flux(u(f-1), u(f), wf);
   res(f-1) = res(f-1) + flux;
   res(f)   = res(f)   - flux;
end

res(1) = 0;
res(np) = 0;
