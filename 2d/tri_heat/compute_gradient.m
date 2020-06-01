function [] = compute_gradient(time)

global gradu u n_ei n_eb ei eb nbr_ei nbr_eb norm_ei norm_eb tarea

gradu(:,:) = 0;

% interior edges
for j=1:n_ei
   v1 = ei(j,1);
   v2 = ei(j,2);
   uavg = 0.5 * (u(v1) + u(v2));
   tl = nbr_ei(j,1);
   tr = nbr_ei(j,2);
   gradu(tl,:) = gradu(tl,:) + uavg * norm_ei(j,:);
   gradu(tr,:) = gradu(tr,:) - uavg * norm_ei(j,:);
end

% boundary edges
for j=1:n_eb
   u1 = boundary_value(time, j, 1);
   u2 = boundary_value(time, j, 2);
   uavg = 0.5*(u1+u2);
   tl = nbr_eb(j);
   gradu(tl,:) = gradu(tl,:) + uavg * norm_eb(j,:);
end

gradu(:,1) = gradu(:,1) ./ tarea;
gradu(:,2) = gradu(:,2) ./ tarea;

%min(gradu(:,1))
%max(gradu(:,1))
%min(gradu(:,2))
%max(gradu(:,2))
%pause
