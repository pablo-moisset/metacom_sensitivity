function [L, sources, non_intermitent, xP, yP] = fractal_landscape(land_par)
% Fractal landscape

   global LAND_SIZE_IDX NI_QTY_IDX LAND_SEED_IDX LAND_EDGES_IDX

   land_seed = land_par(LAND_SEED_IDX) ;
   patch_qty = land_par(LAND_SIZE_IDX) ;
   desired_edges = land_par(LAND_EDGES_IDX) ;
   ni_qty = land_par(NI_QTY_IDX) ;
   
   if isfinite(land_seed)
       s=rng  ;
       rng(land_seed) ;
   end
   
   [~, L, xP, yP] = gen_landscape(patch_qty, 0.0, 512) ;

   G = graph(L) ;
   TR = minspantree(G) ;
   sorted_edges = sortrows(G.Edges,2) ;
   %plot(TR, 'XData', xP, 'YData', yP)
   current_edges = TR.adjacency ;

   candidate_idx = 1 ;
   while nnz(current_edges) < 2*desired_edges
      current_edges(sorted_edges.EndNodes(candidate_idx,1), sorted_edges.EndNodes(candidate_idx,2)) = 1 ;
      current_edges(sorted_edges.EndNodes(candidate_idx,2), sorted_edges.EndNodes(candidate_idx,1)) = 1 ;
      candidate_idx = candidate_idx + 1 ;
   end
   L = full(L.*current_edges) ;

   nodes = randperm(patch_qty) ;
   sources = nodes(1) ;
   non_intermitent = nodes(2:(ni_qty+1)) ;
   
   if isfinite(land_seed)
       rng(s)
   end
end
