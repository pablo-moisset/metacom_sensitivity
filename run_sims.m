function run_sims(MASTER_P, LAND_P, METACOM_SIM_P, sample_size,...
                  filename_prefix, filename_suffix)
%run sims V16
%March 24, 2022
%after code cleanup and better parameterization

% MASTER_P: Structure with parameters for master matrix creation
% MASTER_P.max_r
% MASTER_P.min_mort
% MASTER_P.connectance_min
% MASTER_P.connectance_max
% MASTER_P.contiguity
% MASTER_P.richness_min
% MASTER_P.richness_max
% MASTER_P.sigma_min
% MASTER_P.sigma_max
% MASTER_P.basals_diag
%
%MASTER_P.algorithm = 1 ; % 0: generalized niche
                          % 1: generalized preferential prey
% For the GPPM we also need
% MASTER_P.temperature
% MASTER_P.fraction_basals
% MASTER_P.feeding_efficiency


% LAND_P: Structure with parameters for landscape creation
% LAND_P.creation_function
% LAND_P.cf_parameter
% LAND_P.len_on_pars
%
% METACOM_SIM_P: Structure with parameters for metacommunity dynamics
% METACOM_SIM_P.mig_rate
% METACOM_SIM_P.variability
% METACOM_SIM_P.ext_threshold
% METACOM_SIM_P.len_year % length of year
% METACOM_SIM_P.num_years  % Run sim for this many years

resilience = zeros(sample_size,1) ;
C = resilience ;
S = resilience ;
sigma = resilience ;

stability = zeros(sample_size,1) ;
total_biomass = zeros(sample_size,1) ;

event_series = cell(sample_size,1) ;

warning('sample size changed to 1 for the star experiment. Change back!!!')
for sample = 1:1
%for sample = 1:sample_size
   fprintf('Sample %d\n', sample) ;

warning('Foodweb forced for the star experiment. Change back!!!')
%{
   
   if isfinite(MASTER_P.seed)
       s = rng ; 
       rng(MASTER_P.seed) ;
   end

   success = false ; 
   while ~success
       connectance = unifrnd(MASTER_P.connectance_min, MASTER_P.connectance_max) ;
       richness_master = round(unifrnd(MASTER_P.richness_min, MASTER_P.richness_max)) ;
       std_off_d = unifrnd(MASTER_P.sigma_min, MASTER_P.sigma_max) ;
       
       S(sample) = richness_master ;
       C(sample) = connectance ;
       sigma(sample) = std_off_d ;
       [success, A, CM, R] = feas_and_stable_GLV(connectance,...
                                                 MASTER_P,...
                                                 richness_master,...
                                                 std_off_d,...
                                                 10) ; %max iters
   end
   if isfinite(MASTER_P.seed)
      rng(s)
   end
%}
   richness_master = 45 %forced
   S(sample) = richness_master ;
   C(sample) = 0.2 ;
   sigma(sample) = 0.25 ;
   tmp_matrices;

   [landscape, sources, non_intermitent,xP,yP] = LAND_P.creation_function(LAND_P) ;
   assert(all(diag(landscape)==0)) ;
   M{sample} = create_metacom(landscape, sources, non_intermitent, A, METACOM_SIM_P.mig_rate_p) ;
   M{sample}.xP = xP ;
   M{sample}.yP = yP ;

   if isfinite(METACOM_SIM_P.event_seed)
       s = rng ; 
       rng(METACOM_SIM_P.event_seed) ;
   end

   if LAND_P.len_on_period(1)==0
       M{sample}.len_on_period = unifrnd(LAND_P.len_on_period(2),LAND_P.len_on_period(3),...
                                         length(M{sample}.landscape) ,1) ;
   elseif LAND_P.len_on_period(1)==1
       M{sample}.len_on_period = rnd_trunc_power(LAND_P.len_on_period(4),...
                                                 LAND_P.len_on_period(2),...
                                                 LAND_P.len_on_period(3),...
                                                 length(M{sample}.landscape), 1) ;
   else
       assert(0) ;
   end
   Q = patch_events(length(landscape), union(M{sample}.sources,M{sample}.non_intermitent),...
                    METACOM_SIM_P.len_year, METACOM_SIM_P.num_years,...
                    M{sample}.len_on_period, METACOM_SIM_P.variability) ;

   M{sample}.CM = CM ;
   M{sample}.R = R ;
   M{sample}.ext_thr = extintion_thresholds(LAND_P.len_on_period(2:3),...
                                            METACOM_SIM_P.ext_thr_par,...
                                            M{sample}.len_on_period) ;

   event_series{sample} = metacom_sim(M{sample}, Q, METACOM_SIM_P) ;

   if isfinite(METACOM_SIM_P.event_seed)
       rng(s) ;
   end
   LAND_P.seed = LAND_P.seed+1 ;
   MASTER_P.seed = MASTER_P.seed+1 ;
   METACOM_SIM_P.event_seed = METACOM_SIM_P.event_seed + 1;
end
%ts_bin_filename = strcat(filename_prefix, 'tseries', filename_suffix, '.bin')
%serial_data = hlp_serialize(event_series) ;
%f = fopen(ts_bin_filename,'w') ;
%fwrite(f,serial_data);
%fclose(f);

ts_bin_filename = strcat(filename_prefix, filename_suffix,'.hlp')
results = struct() ;
results.M = M ;
results.METACOM_SIM_P = METACOM_SIM_P ;
results.LAND_P = LAND_P ;
results.C = C ;
results.S = S ;
results.MASTER_P = MASTER_P ;
results.sigma = sigma ;
results.event_series = event_series ;
serial_data = hlp_serialize(results) ;
f = fopen(ts_bin_filename,'w') ;
fwrite(f,serial_data);
fclose(f);
[status,cmdout] = system(sprintf('gzip -f %s', ts_bin_filename));
assert(status==0) ;
end

function ext_thr = extintion_thresholds(period_range, ext_thr_pars, len_on_period)
   shortest = min(period_range) ; %Lm
   longest =  max(period_range) ; %LM
   b = log(ext_thr_pars(1)/ext_thr_pars(2))/log(longest/shortest) ;
   a = ext_thr_pars(2)*longest^b ;
   ext_thr = a*len_on_period.^(-b) ;
end
