function main_cluster_v16(proc_id,...
                         rng_seed,... 
                         ext_thr_par,... 
                         mig_rate_par,... 
                         variabilities,... 
                         land_par,... 
                         exp_name,... 
                         packet,...
                         foodweb_seed,...
                         event_seed,...
                         check_feasibility,...
                         len_on_par,...
                         foodweb_richness,...
                         foodweb_connectance,...
                         foodweb_temperature,...
                         foodweb_diagonal,...
                         foodweb_basals_frac,...
                         foodweb_std,...
                         foodweb_bcomp)

proc_id
%rng_seed      : main random seed for this simulation
%ext_thr_par   : 1x2 matrix with the extintion thresholds for the smallest 
%              : (shortest on period) and largest patches.

%mig_rate_par  : 1D cell array of parameters to compute migration rates.
%                For an element of this list, call it e, e(1)=k_m*m_0, 
%                e(2)=b_m^alpha_m, and e(3) is the threshold for rounding
%                the migration rate to zero. Finilly mr = e(1)/D*e(2)^T,
%                where T is the trophic level and D is the distance between
%                patches. mr is rounded to 0 if mr<e(3)

%variabilities : List of scalars Variability of ON (wet) period beginnings.

%land_par: 1D vector with landscape related parameters
%             land_par(1): Alias land_f. Landscape topology 
%                1:chain, 2:ER, 3:tree, 4:grid, 5:grid, single central
%                source, 6:Star, 7:fractal, 8:modular, 9:p-modular. 
%                For 1, the single source is at an  endpoint (patch 1). 
%                For 5, and 6 there is a single central source. 
%                For 2,3 and 4 the source is at a random node. For 8 and 9,
%                the source is at node 1. 
%                In all cases (except 8 and 9) the non intermitent patches are 
%                chosen with uniform probability among the non-source patches.
%                For 8 and 9, the non intermitent patches are 2..M+1. This
%                coincides as much as possible con cluster centers.
%             land_par(2): Alias land_size. For land_f in {4,5} this 
%                will produce a square grid with land_size^2 patches.
%                For land_f = 3 this will produce a balanced tree with 
%                land_size levels (hence with 2^land_size-1 patches) 
%                For other values, land_size is the number of patches.
%             land_par(3): Alias ni_qty. Number of non intermitent patches
%             land_par(4): Alias land_seed. Either an integer or Inf. 
%                An integer fixes the landscape seed. This allows to use 
%                the same landscape for all samples in this packet. An 
%                empty array does not force a seed and therefore all 
%                samples use independently generated landscapes.
%             land_par(5): alias land_edges. Number of total edges
%                 It is ignored for trees, grids and chains.
%             land_par(6): Only used for (p-)modular landscapes. It is the
%                number of cluster centers.
%             land_par(7): nodes excess factor for (p-)modular landscapes
%             land_par(8): p_exponent for p-modular landscapes

%exp_name      : String. Some short name describing the esperiment to be
%                added to the output filenames. 'output_'+exp_name is the
%                output dir.
%packet        : integer. This is a hack to divide a treatment into many 
%                packets (subsets of a treatment), each one with 50
%                samples. This helps with memory issues, otherwise save 
%                of timeseries crashes.
% foodweb_seed  : Same idea but for generating the trophic network.
% event_seed    : Same idea but for generating the migration events.
% check_feasibility: boolean. Should check for network feasibility after
%                   trying a species migration. If false, all migrations
%                   are successful and cause no extintion in the
%                   destination patch. Warning, setting it to false may
%                   produce equilibria with negative abundances.
% len_on_par   : 1D vector of parameters for the ON period generation
%    [0 min max] generates lengths using a uniform random distribution
%    [1 min max a] generates lengths using a truncated power_law
%    in both cases, 0.0<=min<=max<=1.0

global LAND_EXCESS_IDX
global MIG_R_KM0_IDX MIG_R_BA_IDX MIG_R_TRS_IDX


LAND_P.F = land_par(1) ;
LAND_P.land_size = land_par(2) ;
LAND_P.non_intermitent = land_par(3) ;
LAND_P.seed = land_par(4) ;
LAND_P.edges = land_par(5) ;
LAND_P.centers = land_par(6) ;
LAND_P.excess = land_par(7) ;
LAND_P.p_exponent = land_par(8) ;

MIG_R_KM0_IDX = 1 ;
MIG_R_BA_IDX = 2 ;
MIG_R_TRS_IDX = 3 ;
rng(rng_seed) ;
MASTER_P.max_r = 1.0 ;
MASTER_P.min_mort = 0.01 ;
MASTER_P.connectance_min = foodweb_connectance ;
MASTER_P.connectance_max = foodweb_connectance ;
%MASTER_P.contiguity = 1.0 ;
MASTER_P.richness_min = foodweb_richness ;
MASTER_P.richness_max = foodweb_richness ;
MASTER_P.sigma_min = foodweb_std ;
MASTER_P.sigma_max = foodweb_std ;
MASTER_P.temperature =  foodweb_temperature ;
MASTER_P.fraction_basals = foodweb_basals_frac ;
MASTER_P.feeding_efficiency = 0.4 ;
MASTER_P.basals_comp = foodweb_bcomp ;
MASTER_P.algorithm = 1 ; % 0: generalized niche
                         % 1: generalized preferential prey
MASTER_P.seed = foodweb_seed ;
MASTER_P.diagonal = foodweb_diagonal ;
MASTER_P.basals_diag = MASTER_P.diagonal ;

% They all return [L, sources, non_intermitent]
% L is the"Landscape matrix". L(i,j)>0 iff from i we can reach 
% j in one hop. A positive L(i,j) is the distance from i to j.

land_f = LAND_P.F ;
if land_f==1
   LAND_P.creation_function = @(land_p)chain_landscape(land_p) ;
elseif land_f==2
   LAND_P.creation_function = @(land_p)fixed_landscape(land_p) ;
elseif land_f==3
   LAND_P.creation_function = @(land_p)ftree_landscape(land_p) ;
%   patch_qty = 2^(landscape_size-1) ;
elseif land_f==4
   LAND_P.creation_function = @(land_p)sgrid_landscape(land_p) ;
%   patch_qty = landscape_size^2 ;
elseif land_f==5
   LAND_P.creation_function = @(land_p)sgrid_center_landscape(land_p) ;
%   patch_qty = landscape_size^2 ;
elseif land_f==6
   LAND_P.creation_function = @(land_p)star_landscape(land_p) ;
elseif land_f==7
   LAND_P.creation_function = @(land_p)fractal_landscape(land_p) ;
elseif land_f==8
   LAND_P.creation_function = @(land_p)modular_landscape(land_p) ;
elseif land_f==9
   LAND_P.creation_function = @(land_p)p_modular_landscape(land_p) ;
else
   assert(false) ;
end

LAND_P.cf_parameter = [] ;
LAND_P.len_on_period = len_on_par ;

warning('Sample size changed to 1 for star experiment') ;
%sample_size = 10;%% Samples per packet. Magic constant
sample_size = 1;%% Samples per packet. Magic constant
METACOM_SIM_P.check_feasibility = check_feasibility ;
METACOM_SIM_P.ext_thr_par = ext_thr_par ;
METACOM_SIM_P.len_year = 1.0 ; % length of year
METACOM_SIM_P.num_years = 5 ;  % Run sim for this many years
METACOM_SIM_P.event_seed = event_seed ;

results = [] ;
for i = 1:numel(mig_rate_par)
    mig_rate_p = mig_rate_par{i} ;
    for variability = variabilities
       METACOM_SIM_P.mig_rate_p = mig_rate_p ;
       METACOM_SIM_P.variability = variability ;
       [mig_rate_p, variability]
       filename_prefix = ['outputs_' exp_name '/tseries_' proc_id] ;

       filename_suffix = assemble_suffix(ext_thr_par, land_par,...
                                         mig_rate_p, len_on_par,...
                                         variability, packet, exp_name) ;

       run_sims(MASTER_P, LAND_P, METACOM_SIM_P, sample_size,...
                filename_prefix, filename_suffix) ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [L, sources, non_intermitent] = fixed_landscape(land_par)
        % ER networks.
        % Iterates until all patches are reachable from node 1
        % an integer value for seed fixes the landscape, an []
        % value for seed does not set the random seed.
        land_seed = land_par.seed ;
        land_size = land_par.size ;
        land_conn = land_par.edges ;
        ni_qty = land_par.non_intermitent ;
        
        if isfinite(land_seed)
            s=rng  ;
            rng(land_seed) ;
        end
        
        is_connected = false ;
        while ~is_connected
            L = er_landscape(land_size, land_conn) ;
            G =  digraph(L) ;
            D = distances(G, 1) ;
            is_connected = all(isfinite(D)) ;
        end
        
        sources = nodes(1) ;
        
        non_intermitent = 2:(1+ni_qty) ;
        
        if isfinite(land_seed)
            rng(s)
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [L, sources, non_intermitent] = ftree_landscape(land_par)
        % Full and complete binary tree generator
        land_seed = land_par(LAND_SEED_IDX) ;
        levels = land_par(LAND_SIZE_IDX) ;
        ni_qty = land_par(NI_QTY_IDX) ;
        
        if isfinite(land_seed)
            s=rng  ;
            rng(land_seed) ;
        end
        
        node_qty = 2^levels - 1 ;
        L = zeros(node_qty) ;
        for i = 1:((node_qty-1)/2)
            L(i,2*i)   = 1 ;
            L(i,2*i+1) = 1 ;
        end
        L = L + L' ;
        
        nodes = randperm(node_qty) ;
        
        sources = nodes(1) ;
        non_intermitent = nodes(2:(ni_qty+1)) ;
        
        if isfinite(land_seed)
            rng(s)
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function s = assemble_suffix(ext_thr_par, land_par, mig_rate_p,...
                                 len_on_par, variability, packet, exp_name)
        
%{
        s = sprintf(['_x%1.3g,%1.3g',...
                    '_l%d,%d,%d,%d,%1.3g'...
                    '_m%1.3g,%1.3g,%1.3g'...
                    '_w%d,%1.3g,%1.3g'...
                    '_v%1.3g'...
                    '_P%02d_%s.mat']Niklaus66
,...
                     ext_thr_par,...
                     land_par,...
                     mig_rate_p,...
                     len_on_par,...
                     variability,...
                     packet, exp_name) ;
%}
       s = tempname('.') ;
       s = ['_' s(3:end) sprintf('_P%02d_%s', packet, exp_name)] ;
    end
end
