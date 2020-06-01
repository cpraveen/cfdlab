function uavg = exact_cell_average(t)

global np nt tri X Y tarea carea

uavg = zeros(np,1);

for j=1:nt
   v = tri(j,:);
   xc = sum(X(v))/3;
   yc = sum(Y(v))/3;
   ue = exact_solution(t,xc,yc);
   uavg(v) = uavg(v) + ue * tarea(j) / 3;
end

uavg = uavg ./ carea;
