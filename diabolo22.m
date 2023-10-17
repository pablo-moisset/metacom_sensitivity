function [L, sources, non_intermitent, xP, yP] = diabolo22(land_par)
% diabolo landscape with 22 nodes. Just for testing 
   sources = 1 ;
   non_intermitent = [] ;
   
   L=zeros(22) ;
   
   L(1,2:11) = 1 ;
   L(2:11,12) = 1 ;
   L(12,13:22) = 1 ;
   
   L = ((L+L') > 0)*200 ;
   xP = [0, ones(1,10), 2, 3*ones(1,10) ]'*150; 
   yP = [1.5, linspace(0,3,10), 1.5, linspace(0,3,10) ]'*150; 
end