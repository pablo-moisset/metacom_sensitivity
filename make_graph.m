function [G]= make_graph(L)
   A = L > 0 ;
   A=L
   assert(all(diag(A)==0)) ;
   D = diag(sum(A)) ;
%   Laplacian = D-A
%   [V,E] = eig(Laplacian) ;
%   [s,idx] = sort(abs(diag(E))) ;
   W = A + inv(D);
   [V,E] = eigs(W', 3, 'sm');
   [s,idx] = sort(abs(diag(E))) ;
   X = V(:,idx(1)) ;
   Y = V(:,idx(2)) ;
   G = graph(A) ;
   plot(G, 'XData',X,'YData',Y) ;
end
   