function [L, sources, non_intermitent] = chain_landscape(land_par)
% chain landscape. Edges from i to i+1 and from i+1 to i are present
% 
   land_seed = land_par(LAND_SEED_IDX) ;
   patch_qty = land_par(LAND_SIZE_IDX) ;
   ni_qty = land_par(NI_QTY_IDX) ;

   if ~isempty(land_seed)
      s=rng  ;
      rng(land_seed) ;
   end
   
   directed=false ;

   L = diag(ones(patch_qty-1,1),1) + (1-directed)*diag(ones(patch_qty-1,1), -1);
   sources = 1 ;
   
   tmp = setdiff(1:patch_qty, sources) ;
   idx = randperm(patch_qty-1) ;
   non_intermitent = tmp(idx(1:ni_qty)) ;

%   non_intermitent = [] ;
   if ~isempty(land_seed)
      rng(s)
   end

end
