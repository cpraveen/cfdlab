function [A,b] = compute_matrix_weak()

global np n_ei n_eb ei eb nbr_ei nbr_eb norm_ei norm_eb tarea
global v_ei v_eb tri tri_normal mu cpen e_dir e_neu

fprintf(1,'Computing matrix for weak bc\n');

A = sparse(np,np);
b = zeros(np,1);

% interior edges
for e=1:n_ei
   vj = ei(e,1);
   vk = ei(e,2);
   tria = nbr_ei(e,:);

   for n=1:2
      vi = v_ei(e,n);
      i = find(tri(tria(n),:)==vi);
      j = find(tri(tria(n),:)==vj);
      k = find(tri(tria(n),:)==vk);

      ni= tri_normal(tria(n),i,:);
      nj= tri_normal(tria(n),j,:);
      nk= tri_normal(tria(n),k,:);

      af = mu(tria(n)) * 0.25 / tarea(tria(n));

      A(vi,vi) = A(vi,vi) - af*dot(ni, ni);
      A(vi,vj) = A(vi,vj) - af*dot(nj, ni);
      A(vi,vk) = A(vi,vk) - af*dot(nk, ni);
   end
end

% boundary edges: ignore bc
for e=1:n_eb
   vi = v_eb(e);
   vj = eb(e,1);
   vk = eb(e,2);

   tl = nbr_eb(e);

   i = find(tri(tl,:)==vi);
   j = find(tri(tl,:)==vj);
   k = find(tri(tl,:)==vk);

   ni= tri_normal(tl,i,:);
   nj= tri_normal(tl,j,:);
   nk= tri_normal(tl,k,:);

   af = mu(tl) * 0.25 / tarea(tl);

   % vertex opposite edge
   A(vi,vi) = A(vi,vi) - af*dot(ni, ni);
   A(vi,vj) = A(vi,vj) - af*dot(nj, ni);
   A(vi,vk) = A(vi,vk) - af*dot(nk, ni);

   % vertex on edge: j
   A(vj,vj) = A(vj,vj) - af*dot(nj,ni);
   A(vj,vi) = A(vj,vi) - af*dot(ni,ni);
   A(vj,vk) = A(vj,vk) - af*dot(nk,ni);

   % vertex on edge: k
   A(vk,vk) = A(vk,vk) - af*dot(nk,ni);
   A(vk,vi) = A(vk,vi) - af*dot(ni,ni);
   A(vk,vj) = A(vk,vj) - af*dot(nj,ni);
end

% boundary edges: correct for dirichlet bc
for ee=1:length(e_dir)
   e  = e_dir(ee);
   vi = v_eb(e);
   vj = eb(e,1);
   vk = eb(e,2);

   tl = nbr_eb(e);

   i = find(tri(tl,:)==vi);
   j = find(tri(tl,:)==vj);
   k = find(tri(tl,:)==vk);

   ni= tri_normal(tl,i,:);
   nj= tri_normal(tl,j,:);
   nk= tri_normal(tl,k,:);

   af = mu(tl) * 0.25 / tarea(tl);

   % vertex opposite edge
   A(vi,vj) = A(vi,vj) - af*dot(ni, ni);
   A(vi,vk) = A(vi,vk) - af*dot(ni, ni);

   % vertex on edge: j
   A(vj,vj) = A(vj,vj) - af*dot(nj,ni);
   A(vj,vk) = A(vj,vk) - af*dot(nj,ni);

   A(vj,vj) = A(vj,vj) - af*dot(ni,ni);
   A(vj,vk) = A(vj,vk) - af*dot(ni,ni);

   % vertex on edge: k
   A(vk,vk) = A(vk,vk) - af*dot(nk,ni);
   A(vk,vj) = A(vk,vj) - af*dot(nk,ni);

   A(vk,vk) = A(vk,vk) - af*dot(ni,ni);
   A(vk,vj) = A(vk,vj) - af*dot(ni,ni);

   % penalty term
   ap = 0.5*cpen*mu(tl);
   A(vj,vj) = A(vj,vj) - ap;
   A(vk,vk) = A(vk,vk) - ap;

   % right hand side vector
   fj = boundary_value(0, e, 1);
   fk = boundary_value(0, e, 2);
   fa = 0.5*(fj + fk);
   b(vi) = b(vi) + 0.5*mu(tl)*fa*dot(ni,ni)/tarea(tl);
   b(vj) = b(vj) + 0.5*mu(tl)*fa*dot(nj,ni)/tarea(tl);
   b(vk) = b(vk) + 0.5*mu(tl)*fa*dot(nk,ni)/tarea(tl);

   b(vj) = b(vj) + 0.5*mu(tl)*fa*dot(ni,ni)/tarea(tl);
   b(vk) = b(vk) + 0.5*mu(tl)*fa*dot(ni,ni)/tarea(tl);

   b(vj) = b(vj) + ap*fj;
   b(vk) = b(vk) + ap*fk;
end

% boundary edges: correct for neumann bc
for ee=1:length(e_neu)
   e  = e_neu(ee);
   vi = v_eb(e);
   vj = eb(e,1);
   vk = eb(e,2);

   tl = nbr_eb(e);

   i = find(tri(tl,:)==vi);
   j = find(tri(tl,:)==vj);
   k = find(tri(tl,:)==vk);

   ni= tri_normal(tl,i,:);
   nj= tri_normal(tl,j,:);
   nk= tri_normal(tl,k,:);

   af = mu(tl) * 0.25 / tarea(tl);

   % vertex on edge: j
   A(vj,vj) = A(vj,vj) + af*dot(nj,ni);
   A(vj,vi) = A(vj,vi) + af*dot(ni,ni);
   A(vj,vk) = A(vj,vk) + af*dot(nk,ni);

   % vertex on edge: k
   A(vk,vk) = A(vk,vk) + af*dot(nk,ni);
   A(vk,vi) = A(vk,vi) + af*dot(ni,ni);
   A(vk,vj) = A(vk,vj) + af*dot(nj,ni);

   % right hand side vector
   % we assume zero neumann condition
end
