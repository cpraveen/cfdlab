function minmax(ubar)
globals;

N = length(ubar);

uu = [ubar(1:N-2); ubar(2:N-1); ubar(3:N)];

umin(2:N-1) = min(uu);
umax(2:N-1) = max(uu);

if periodic==yes
   umin(1) = min([ubar(N), ubar(1), ubar(2)]);
   umin(N) = min([ubar(N-1), ubar(N), ubar(1)]);
   umax(1) = max([ubar(N), ubar(1), ubar(2)]);
   umax(N) = max([ubar(N-1), ubar(N), ubar(1)]);
else
   umin(1) = min([ubar(1), ubar(2)]);
   umin(N) = min([ubar(N-1), ubar(N)]);
   umax(1) = max([ubar(1), ubar(2)]);
   umax(N) = max([ubar(N-1), ubar(N)]);
end
