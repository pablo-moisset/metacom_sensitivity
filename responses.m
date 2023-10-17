function [persistence_alpha, persistence_beta, persistence_gamma,...
          biomass_beta, biomass_gamma,...
          shannon_alpha, shannon_beta, shannon_gamma, t]...
          = responses(T, p, include_sources)
% T is a 1d cell array of metacom structures
% This version redefines shannon alpha index
% p is the value for the generalized mean
% include_sources is a boolean defining if permanent sites should be considered
%   typically, it makes sense to use false here
% outputs are scalars except t, which is an array of time stamps (one for event)

richness_master = length(T{1}.A) ; 

T  = T(1:end-1) ; %remove the last element to avoid inf timestamp
CM = T{1}.CM ;
R  = T{1}.R ;

L = length(T) ;
patch_qty = size(T{1}.present,1) ;

consider_patches = 1:patch_qty ;
if ~include_sources
   consider_patches = setdiff(consider_patches, T{1,1}.sources) ;
end

L

block_size = 100 ; %somewhat arbitrary

norm_p = zeros(1, length(consider_patches)) ;
sys_present = zeros(L,1);

p_alpha = 0 ; % power mean of instantaneous persistence alpha
p_gamma = 0 ; % power mean of instantaneous persistence gamma
cv = 0 ;      % power mean of coefficient of variation
biomass_gamma = 0 ;
shannon_alpha = 0 ;
shannon_gamma = 0 ;

tspan = T{end}.timestamp - T{1}.timestamp ;

for block_idx_start = 1:block_size:(L-1)
   block_idx_end = min(block_idx_start + block_size, L ) ; %Overlap between end of a block and beginning of new one is intentional
   block_actual_size =  block_idx_end - block_idx_start + 1 ;

%   patch_present = zeros(block_actual_size, patch_qty);
   p_alpha_inst = zeros(block_actual_size,1) ; %p_alpha(k) is p_alpha at time t(block_idx_start+k-1)
   p_gamma_inst = zeros(block_actual_size,1) ; %p_gamma(k) is p_gamma at time t(block_idx_start+k-1)
   biomass = zeros(block_actual_size, patch_qty, richness_master, 'single');
   t = zeros(block_actual_size,1);

   for k = 1:block_actual_size %k is a time index
      t(k) = T{block_idx_start+k-1}.timestamp ;
      P = T{block_idx_start+k-1}.present ; %P(i,j)==1 iff species j is in patch i at time t(k)
      patch_present = sum(P(consider_patches,:),2)' ; % Richnes at time t(k), for every patch
      p_alpha_inst(k) = mean(patch_present/richness_master) ;
      p_gamma_inst(k) = sum(sum(P(consider_patches,:))>0)/richness_master ; % Gamma persistence at time t(k) for this block
      for i = 1:patch_qty
         present_ki = P(i,:) ; %present at time t(k), patch i
         CM_present = CM(present_ki,present_ki) ;
         R_present = R(present_ki) ;
         biomass(k,i,present_ki) = single(-CM_present\R_present) ; %time t(k), patch i, species
      end
   end

   % Should be a scalar. It is a PARTIAL computation of norm p
   p_alpha = p_alpha + cont_p_norm(p_alpha_inst, t, p) ;

   % Should be a scalar. It is a PARTIAL computation of norm p
   p_gamma = p_gamma + cont_p_norm(p_gamma_inst, t, p) ;

   biomass_for_r = biomass(:,consider_patches,:) ; % Biomasses to compute responses
   %clear biomass ;

   norm_p = cont_p_norm(sum(sum(biomass_for_r,3),2), t, p) ; %norm_p is 1 x length(consider_patches). It is a PARTIAL computation of norm p
   %should be a scalar
   biomass_gamma = biomass_gamma + norm_p ; %accumulated

   cv_block = std(sum(biomass_for_r,3),0,2)./(mean(sum(biomass_for_r,3),2)+realmin('single')) ;
   norm_p = cont_p_norm(cv_block,t,p) ; 

   %Should be a scalar
   cv = cv + norm_p ;

   %tmp1 = reshape(norm_p, richness_master, []) ; %tmp1(i,j) is the p-avg (over time) biomass value of species i patch j
   %tmp2 = diag(1./(sum(tmp1,2)+realmin('single')))*tmp1  ;
   tmp1 = sum(biomass_for_r,3) ; assert(all(size(tmp1) == [block_actual_size, numel(consider_patches)])) ;

   % n_biomass(k,i,s), normalized biomass of species s, patch k, time t(k)
   n_biomass = biomass_for_r./(repmat(tmp1,1,1,richness_master)+realmin('single')) ;
   shannon_alpha_block = mean(exp(-sum(n_biomass.*log(n_biomass+realmin('single')), 3))/richness_master,2) ;
   shannon_alpha = shannon_alpha + cont_p_norm(shannon_alpha_block,t,p) ;

   tmp1 = squeeze(sum(biomass_for_r,2)) ; % tmp1(k,i) is the biomass of species i at time t(k) (entire metacom)
   n_biomass = tmp1./(repmat(sum(tmp1,2),1,richness_master)+realmin('single')) ; %Normalize dividing by total species biomass
   shannon_gamma_block = exp(-sum(n_biomass.*log(n_biomass+realmin('single')), 2))/richness_master ;
   shannon_gamma = shannon_gamma + cont_p_norm(shannon_gamma_block,t,p) ;
end

tspan = T{end}.timestamp - T{1}.timestamp ;

persistence_alpha = (p_alpha/tspan)^(1/p) ;
persistence_gamma = (p_gamma/tspan)^(1/p) ;
%P_\beta is one if P_\alpha=P_\gamma=0. Hence the addition ofrealmins
persistence_beta =  (persistence_gamma+realmin)/(persistence_alpha+realmin) ;

biomass_gamma = (biomass_gamma/tspan)^(1/p) ;
biomass_beta = (cv/tspan)^(1/p) ;

shannon_alpha = (shannon_alpha/tspan)^(1/p) ;
shannon_gamma = (shannon_gamma/tspan)^(1/p) ;
shannon_beta = shannon_gamma - shannon_alpha ;
end
