function  A = gen_niche_model(connectance , species_qty, contiguity)
% Version creada 14 Julio 2020. Rescritura de version Leslie
% Based on Williams & Martínez 2000; stouffer et al 2006
% 
% A is a species_qty x species_qty 0-1 matrix. A(i,j)=1 iff j preys on i

assert( connectance<=1 & connectance>=0 ) ;
assert( contiguity<=1 & contiguity>=0 ) ;
assert( species_qty >= 1 ) ;

%NICHE VALUE
niche_value = sortrows( rand(species_qty, 1) ) ;

%RANGE VALUE
%assign the range value to each i species (beta no puede ser <0)
Beta = (1/(2*connectance))-1;
X = betarnd(1,Beta,species_qty,1);
range_NM = niche_value.*X;  %ranges for the Niche Model 
                            %length of "prey selection interval"

%GENERALIZED NICHE MODEL
range_GNM = contiguity * range_NM; % ranges adjusted by contiguity param.

range_GNM(1) = 0 ; % This forces species 1 to be basal

%CENTRE OF "PREDATION INTERVALS"
center = rand(species_qty,1).*(niche_value-range_GNM/2) + range_GNM/2;

left_boundary  = center - range_GNM/2;
right_boundary = center + range_GNM/2;

A = zeros(species_qty);

for j = 1:species_qty
   A(:,j)=left_boundary(j)<niche_value & niche_value<right_boundary(j);
end

deltak = round( (1-contiguity)*range_NM*species_qty ) ;

predators = find(deltak>0) ;
for predator_idx = 1:numel(predators) %add deltak(j) preys to predator j
   j = predators(predator_idx) ;
   candidate_preys =  find(niche_value < niche_value(j)) ;
   candidate_preys =  setdiff(candidate_preys, A(:,j)) ;
   n = numel(candidate_preys) ;
   candidate_preys = candidate_preys(randperm(n)) ;
   A(candidate_preys(1:min(deltak(j),n)),j) = 1 ; %add extra preys
end

end
