function  [A, TP, q]  = gppm_coherence( S, B, L, T )
% Generalized preferential prey model
% Creates a trophic network using the algorithm found in 
% "From neurons to epidemics: How trophic coherence affects spreading
%  processes", Janis Klaise and Samuel Johnson, Chaos, 2016 ;
% 
% The output network may contain cycles but not self loops (canibalism)
%
% A: (0-1) Adjacency matrix. A(i,j)=1 iff j preys on i.
% TP: 1xlength(A) vector. Tropic position of each species (Basals are at
% level 0)
% q: positive real. Throphic incoherence (see Klaise et al 2016).
% S: Species richness
% B: Number of basal species
% L: Number of links
% T: temperature. The higher the temp, the lower the coherence.
%    do not use very small temperatures (say T<0.01) to avoid 
%    floating point underflows.

MAX_ITERS = 20 ;

TD = zeros(S) ;

iters=0 ;
while sum(TD(:)>0) <= (L-S+B) && iters < MAX_ITERS

A = zeros(S) ;

for j=(B+1):S
    i  = randsample((j-1),1) ;
    A(i,j) = 1 ;
end

TP = tropos(A)' ;

TP = repmat(TP,length(TP),1);

TD = exp(-(TP-TP'-1).^2/(2*T^2)) ; %Favor tropic interactions between species
                                   %whose trophic positions differ by one.
TD(:,1:B) = 0 ; % Basals cannot prey on anything
TD = TD-diag(diag(TD)) ; % Eliminate the possibility of self loops
TD = TD.*(1-A) ; %Prob of adding an existing arc is set to zero
TD = TD/sum(TD(:)) ; %Normalize probabilities
iters = iters + 1;
end


if sum(TD(:)>0) >= (L-S+B)
for l=1:(L-S+B)
   idx_candidate = find(A(:), 1, 'first') ; %find an existing arc
   assert(~isempty(idx_candidate)) ;
   while A(idx_candidate) == 1 % while candidate is already there, 
       TD(idx_candidate) = 0 ; % set its entry in the probability matrix TD  to
                               % zero
       assert(sum(TD(:)>0)>0) ; %There must be at least one possible edge to add
       TD = TD/sum(TD(:)) ;   %Re-normalize probabilities
       idx_candidate = find(mnrnd(1,TD(:))) ; %Choose a new candidate
   end
   A(idx_candidate) = 1 ;
end

TL = tropos(A) ; % Recompute trophic leves after all additions
TL = repmat(TL,1,length(TP));

D = (TL'-TL);
tmp = D.*A;
%mean(tmp(tmp~=0))
q  = sqrt(mean(tmp(tmp~=0).^2)-1) ;
else
A=[] ;
TP=[] ;
q=[] ;
assert(false);
end
