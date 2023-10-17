function  [F, X, R] = is_feasible(M, max_r, min_mort)
% Returns true iff there is an r such that the system dx/dt = x.*(r + Mx) has positive equilibrium
% M[i,j] is a real number describing the per-capita effect of species j over species i
% max_r is the maximum allowed (absolute) value for r.
% min_mort is a positive number

sign_of_M = sign(M) ;
tmp = sign_of_M + sign_of_M' ;  % tmp[i,j]~=0 iff there is a non-trophic relation between i and j
Trophic_M = sign_of_M ;
Trophic_M(tmp~=0) = 0 ; % qualitative community matrix with only the signs of trophic interactions

%basals = find(sum(Trophic_M==-1)==0) ;  % detects basal species
non_basals = find(sum(Trophic_M==-1)>0) ; % detects non-basal species

n=length(M) ;
Aeq=[eye(n) M zeros(n,1)] ;
beq = zeros(n,1) ;

% tmp = ones(n,1) ;
% tmp(basals) = -1 ;
% A = [ zeros(n) -eye(n) ones(n,1) ; diag(tmp) zeros(n) ones(n,1)] ;
% b=zeros(2*n,1) ;

A = [ zeros(n) -eye(n) ones(n,1)] ;
b=zeros(n,1) ;

ub = [max_r*ones(n,1) ; inf(n+1,1)] ;
lb = -ub ;
ub(non_basals) = -min_mort ;

f = [zeros(2*n,1) ; -1] ;
%[x,~,exit_flag]=linprog(f,A,b,Aeq,beq,lb,ub,optimoptions('linprog','Display','off'));
options = optimoptions('linprog','Algorithm', 'interior-point', 'Display','none');
[x,~,exit_flag]=linprog(f,A,b,Aeq,beq,lb,ub,[],options);
F = exit_flag >0 && x(end) > 0 ;

if  F>0
    X=x(n+1:end-1);
    R=x(1:n);
else
    X=zeros(n,1) ;
    R=X ;
end
end