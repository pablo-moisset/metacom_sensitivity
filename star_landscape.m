function [L, sources, non_intermitent, xP, yP] = star_landscape(land_par)
% Star topology. Seed is always 1

land_seed = land_par.seed ;
patch_qty = land_par.land_size ;
ni_qty = land_par.non_intermitent ;

if isfinite(land_seed)
    s=rng  ;
    rng(land_seed) ;
end

L = zeros(patch_qty) ;

L(1,2:patch_qty) = 1 ;

L = L + L' ;

sources = 1 ;

tmp = setdiff(1:patch_qty, sources) ;
idx = randperm(patch_qty-1) ;
non_intermitent = tmp(idx(1:ni_qty)) ;
%   non_intermitent = [] ;

xP = zeros(patch_qty,1) ;
yP = zeros(patch_qty,1) ;
for i = 2:patch_qty
    xP(i) = 300*cos(2*pi/(patch_qty-1)*i ) ;
    yP(i) = 300*sin(2*pi/(patch_qty-1)*i ) ;
end
if isfinite(land_seed)
    rng(s)
end
end
