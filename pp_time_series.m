function pp_time_series(filename_in, filename_out, p)
% Post processes time series file (filename_in) leaving response variables in an output file 
% (filename_out)
% Input file has a time series per sample
% Output file has 11 one-dimensional arrays with the response variables
% and the experimental factors. Each array has 1 item per sample

fprintf('Input: %s\nOutput: %s\np: %f\n',filename_in, filename_out, p) ;

all_C = [] ;
all_C_measured = [] ;
all_C_land_measured = [] ;
all_S = [] ;
all_umig_rate = [] ;
all_sigma = [] ;
all_stability = [] ;
all_total_biomass = [] ;
all_p_alpha = [] ;
all_p_beta_sub = [] ;
all_p_beta_div = [] ;
all_p_gamma = [] ;
all_b_alpha = [] ;
all_b_beta_sub = [] ;
all_b_beta_div = [] ;
all_b_gamma = [] ;
all_ext_tsh = [] ;
all_variability = [] ;
all_shannon_alpha = [] ;
all_shannon_beta = [] ;
all_shannon_gamma = [] ;

r = load_serialized(filename_in) ;
   
LAND_P = r.LAND_P ;
METACOM_SIM_P = r.METACOM_SIM_P ;
M = r.M ;
S = r.S
MASTER_P = r.MASTER_P ;
sigma = r.sigma ;

all_C =  r.C ; %Desired foodweb connectance. Not always exactly achieved
all_S =  r.S ;
all_sigma = sigma ;
   
all_C_land_measured = zeros(numel(M),1) ;
for j = 1:numel(M)
    patch_qty = length(M{j}.landscape) ;
    all_C_land_measured(j) = sum(M{j}.landscape(:))/(patch_qty^2-patch_qty) ;
end

all_C_measured = zeros(numel(M),1) ;
for j = 1:numel(M)
    all_C_measured(j) = sum(M{j}.A(:))/(S(j)^2-S(j)) ;
end

%warning('8 should be 1. This value is for debugging only')
for j = 1:length(r.event_series) 
   'interpreting'
   ts = interpret_events(r.event_series{j}, M{j}) ;
   [p_alpha, p_beta_sub, p_beta_div, p_gamma, b_alpha, b_beta_sub, b_beta_div, b_gamma, s_alpha, s_beta, s_gamma] = responses(ts, p, false) ;
   all_p_alpha = [all_p_alpha; p_alpha] ;
   all_p_beta_sub = [all_p_beta_sub; p_beta_sub] ;
   all_p_beta_div = [all_p_beta_div; p_beta_div] ;
   all_p_gamma = [all_p_gamma; p_gamma] ;
   all_b_alpha = [all_b_alpha; b_alpha] ;
   all_b_beta_sub = [all_b_beta_sub; b_beta_sub] ;
   all_b_beta_div = [all_b_beta_div; b_beta_div] ;
   all_b_gamma = [all_b_gamma; b_gamma] ;
   all_shannon_beta = [all_shannon_beta; s_beta] ;
   all_shannon_gamma = [all_shannon_gamma; s_gamma] ;
   all_shannon_alpha = [all_shannon_alpha; s_alpha] ;
   
   %all_ext_tsh(i,j) is the extintion threshold for patch j, sample i 
   all_ext_tsh = [all_ext_tsh ; M{j}.ext_thr'] ; 
   
   %all_umig_rate(i,j) is the "unit migration rate" for species j, sample i
   % migration rate is the unit migration rate / distance between parches
   all_umig_rate = [all_umig_rate ; M{j}.unit_mig_rates'] ;

end

all_variability = [METACOM_SIM_P.variability*ones(length(r.event_series),1)] ;

for i = 1:length(M)
   M{i} = rmfield(M{i},'present') ;
   M{i} = rmfield(M{i},'active') ;
   M{i} = rmfield(M{i},'basals') ;
   M{i} = rmfield(M{i},'last_event_type') ;
   M{i} = rmfield(M{i},'last_touched_patch') ;
   M{i} = rmfield(M{i},'unit_mig_rates') ;
   M{i} = rmfield(M{i},'CM') ;
   M{i} = rmfield(M{i},'R') ;
end

save(filename_out, 'all_p_alpha', 'all_p_beta_sub', 'all_p_beta_div', 'all_p_gamma', ...
                   'all_b_alpha', 'all_b_beta_sub', 'all_b_beta_div', 'all_b_gamma',...
                   'all_shannon_beta', 'all_shannon_gamma', 'all_shannon_alpha',...
                   'all_ext_tsh', 'all_variability', 'all_umig_rate', ...
                   'all_C_measured', 'all_C', 'all_C_land_measured', 'all_S', 'all_sigma',...
                   'LAND_P', 'METACOM_SIM_P', 'M', 'MASTER_P') ;     
end
