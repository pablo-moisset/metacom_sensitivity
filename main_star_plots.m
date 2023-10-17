% Makes graphs to show the mechanism by which, for unstable foodwebs low
% alpha diversity comes with high gamma
% run animate_metacom first to load data



richness_master=45
%fname_out = '3_star_stable' ;


filename_list = {...
'outputs_star_experiment/tseries_star_test_tp1b5bb97d_f540_4c15_8ac2_5762be3da813_P00_star_experiment.hlp.gz',...   
'outputs_star_experiment/uns_tseries_star_test_tpd991a2be_b50c_4a11_907a_a5dd9124e7b8_P00_star_experiment.hlp.gz',...
'outputs_star_experiment/tseries_star_test_tp830546d8_967a_4882_9b36_f46bff78cfcd_P00_star_experiment.hlp.gz',...
'outputs_star_experiment/uns_tseries_star_test_tp626681ec_8d78_48d9_96b0_76a5d4c81f1b_P00_star_experiment.hlp.gz'
} ;

% Star, three points, unstable foodweb
%r = load_serialized('outputs_star_experiment/uns_tseries_star_test_tp626681ec_8d78_48d9_96b0_76a5d4c81f1b_P00_star_experiment.hlp.gz');

% Star, ten points, unstable foodweb
%r = load_serialized('outputs_star_experiment/uns_tseries_star_test_tpd991a2be_b50c_4a11_907a_a5dd9124e7b8_P00_star_experiment.hlp.gz') ;

% Star, three points, stable foodweb
%r = load_serialized('outputs_star_experiment/tseries_star_test_tp830546d8_967a_4882_9b36_f46bff78cfcd_P00_star_experiment.hlp.gz');

% Star, ten points, stable foodweb
%r = load_serialized('outputs_star_experiment/tseries_star_test_tp1b5bb97d_f540_4c15_8ac2_5762be3da813_P00_star_experiment.hlp.gz');

fh={} ;
ph={} ;
figure ;
for idx=1:numel(filename_list)
    subplot(2,2,idx)
    r = load_serialized(filename_list{idx});
    M = r.M ;
    % METACOM_SIM_P = r.METACOM_SIM_P ;
    % LAND_P = r.LAND_P ;
    % C = r.C ;
    % S = r.S ;
    % sigma = r.sigma ;
    event_series = r.event_series ;
    
    one_M = M{1} ;
    one_E = event_series{1} ;
    
    biomass_gamma =[];
    persistence_alpha = [] ;
    
    [M_TS,col,ext] = interpret_events(one_E, one_M) ;
    
    [persistence_alpha(1), persistence_beta, persistence_gamma,...
        biomass_beta, biomass_gamma(1),...
        shannon_alpha, shannon_beta, shannon_gamma, t]=responses(M_TS,4,false);
    
    persistence_gamma
    persistence_alpha(1)*3/10
    persistence_beta*10/3
    
    
    [R, persistence_gamma, t] = responses_low_mem(M_TS,false);
    R = R';
    
    RR = R(:,2:end)';
    p_alpha = sum(RR)/richness_master/10 ;  % P_alpha for a 10 point star regardless of the active ones
    stairs(t,p_alpha);
    
    hold on;
    p_gamma = persistence_gamma'/richness_master ; %persistence_gamma is actually species richness
    stairs(t,p_gamma) ;
    %stairs(t,p_gamma./p_alpha) ;
    
    xlim([0.97 1.3]);
    ylim([0 1.1]);
    xticklabels([])
    box off
    hold off
    set(gcf,'position',[10,10,500,300])
end


%set(gcf,'position',[10,10,500,300])
%savefig(gcf,[fname_out '.fig']);
%saveas(gcf,[fname_out '.pdf']);
