function [L, sources, non_intermitent] = sgrid_center_landscape(land_par)
% square grid p=patch_qty x patch_qty with source at center

   land_seed = land_par(LAND_SEED_IDX) ;
   patch_qty = land_par(LAND_SIZE_IDX) ;
   ni_qty = land_par(NI_QTY_IDX) ;
   
   if isfinite(land_seed)
       s=rng  ;
       rng(land_seed) ;
   end
   
   p = patch_qty * patch_qty ;
   
   two2lin = zeros(patch_qty) ;
   two2lin(1:p) = 1:p ;
   
   L  = zeros(p) ;
   
   for i=1:(patch_qty-1)
       for j=1:(patch_qty-1)
           center = two2lin(i,j) ;
           right = two2lin(i,j+1) ;
           down = two2lin(i+1,j) ;
           L(center,right) = 1 ;
           L(center,down) = 1 ;
       end
   end
   
   for i=1:(patch_qty-1)
       center = two2lin(i,patch_qty) ;
       down = two2lin(i+1,patch_qty) ;
       L(center,down) = 1 ;
   end
   
   for j=1:(patch_qty-1)
       center = two2lin(patch_qty,j) ;
       right = two2lin(patch_qty,j+1) ;
       L(center,right) = 1 ;
   end
   
   L = L + L' ;
   
   sources = floor(patch_qty/2)*patch_qty + floor(patch_qty/2) +1 ;
   
   tmp = setdiff(1:p, sources) ;
   idx = randperm(p-1) ;
   non_intermitent = tmp(idx(1:ni_qty)) ;
   
   if isfinite(land_seed)
       rng(s)
   end
end
