function [A,b] = compute_matrix_strong()

global np n_ei n_eb ei eb nbr_ei nbr_eb norm_ei norm_eb tarea
global v_ei v_eb tri tri_normal mu cpen

fprintf(1,'Computing matrix for strong bc\n');

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

   % vertex on edge: j
   % nothing to add

   % vertex on edge: k
   % nothing to add

   % right hand side
   fj = boundary_value(0, e, 1);
   fk = boundary_value(0, e, 2);
   b(vi) = b(vi) - af*fj*dot(nj,ni) - af*fk*dot(nk,ni);
end

% boundary edges: correct for bc
for e=1:n_eb
   vj = eb(e,1);
   vk = eb(e,2);
   A(vj,:) = 0;
   A(vk,:) = 0;
   b(vj)   = 0;
   b(vk)   = 0;
end
