function TPsp=TP_shortestpath(A)
% Trophic positions computed as shortests paths to a basal
% Pablo Moisset and Rodrigo Ramos 2021
%
% A is a 0-1 n x n matrix (n: number of species). A(i,j)==1 iff j preys on i.
% TPsp is an n x 1 vector.
% TPsp(i) is the shortest length (hop distance) to species i from a basal.
% 

if isempty(A)
    TPsp = [] ;
else
   %Deleting self effects (Loops)0
   M = A - diag(diag(A));  % M = A after setting the main diagonal to 0.

   basals = find( sum(M)==0 );

   G =  digraph(M) ;
   D = distances(G, basals) ;
   TPsp = min(D,[],1)' ;
end

end