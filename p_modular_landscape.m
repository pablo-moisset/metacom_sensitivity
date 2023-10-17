function  [L, sources, non_intermitent, xP, yP] = p_modular_landscape(land_par)

   land_seed = land_par.seed ;
   patch_qty = land_par.land_size ;
   desired_edges = land_par.edges ;
   E = land_par.excess ; %Excess factor
   ni_qty = land_par.non_intermitent ;
   M = land_par.centers ;
   p = land_par.p_exponent ;
  
   
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
   l_edges = G.Edges.Weight ;
   inv_edges = min(l_edges)./l_edges ; % normalized inverse of edges lengths.
                                       % the largest number in inv_edges
                                       % should be 1.0
   p_edge = inv_edges.^p ;
   
   TR = minspantree(G) ;
   current_edges = TR.adjacency ;
   
   qty_to_add = desired_edges - nnz(current_edges)/2 ; 
   while qty_to_add > 0
       p_edge = p_edge / sum(p_edge) ;
       to_add = mnrnd(qty_to_add, p_edge) > 0 ;
       idx = sub2ind(size(current_edges),G.Edges.EndNodes(to_add,1),G.Edges.EndNodes(to_add,2));
       current_edges(idx) = 1 ;
       idx = sub2ind(size(current_edges),G.Edges.EndNodes(to_add,2),G.Edges.EndNodes(to_add,1));
       current_edges(idx) = 1 ;
       qty_to_add = desired_edges - nnz(current_edges)/2 ;
       p_edge(to_add) = 0 ;
   end
   L = full(L.*current_edges) ;

   sources = 1 ;
   non_intermitent = 2:(ni_qty+1) ;

   if isfinite(land_seed)
       rng(s)
   end
end

