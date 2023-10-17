function  [L, sources, non_intermitent, xP, yP] = modular_landscape(land_par)

   land_seed = land_par.seed ;
   patch_qty = land_par.land_size ;
   desired_edges = land_par.edges ;
   E = land_par.excess ; %Excess factor
   ni_qty = land_par.non_intermitent ;
   M = land_par.centers ;

   assert( ni_qty<patch_qty );
   assert( E>=1.0 );

   if isfinite(land_seed)
       s=rng  ;
       rng(land_seed) ;
   end

   point_qty = ceil(patch_qty*E) ;
   xy_puntos = rand(point_qty,2)*512 ;
   centros = 1:M ;
   DParches = squareform(pdist(xy_puntos,'euclidean')) ; 
   Dclosest = min(DParches(:,centros),[],2);
   t = Dclosest/sum(Dclosest(:)) ;
   mt = t(:) ;

   bitmap_landscape = true(1,point_qty) ;
   chosen = numel(bitmap_landscape) ;

   while chosen > patch_qty 
       R = mnrnd(chosen - patch_qty, mt) > 0;
       bitmap_landscape = bitmap_landscape & (~R) ;
       chosen = sum(bitmap_landscape) ; 
   end

   xP = xy_puntos(bitmap_landscape,1);
   yP = xy_puntos(bitmap_landscape,2);
   L = DParches(bitmap_landscape,bitmap_landscape) ; 

   G = graph(L) ;
   TR = minspantree(G) ;
   sorted_edges = sortrows(G.Edges,2) ;
   current_edges = TR.adjacency ;

   candidate_idx = 1 ;
   while nnz(current_edges) < 2*desired_edges
      current_edges(sorted_edges.EndNodes(candidate_idx,1), sorted_edges.EndNodes(candidate_idx,2)) = 1 ;
      current_edges(sorted_edges.EndNodes(candidate_idx,2), sorted_edges.EndNodes(candidate_idx,1)) = 1 ;
      candidate_idx = candidate_idx + 1 ;
   end
   L = full(L.*current_edges) ;

   sources = 1 ;
   non_intermitent = 2:(ni_qty+1) ;

   if isfinite(land_seed)
       rng(s)
   end
end

