MR_a = logspace(log10(30),log10(3000),5) %[50,100,200,300,500,1000]  %Migration rates
MR_b = [1.0]
MR_t = 0.1 %linspace(0.2, 3, 3)
ET = {[1e-3,1e-3]}     %Extinction thresholds

PQ = [5,25,50,100]       %Patch quantities (is 2^PQ-1 patches if it is a tree, PQ^2 for grids)
%LE = [2*PQ]       %Land edges, ignored for trees, grids, chains, fractal
topology = 9    %1:chain landscape, 2:ER, 3:tree, 4:grid, 5:centered grid, 6:star, 7: fractal, 8:modular, 9:p-modular. Source patches assumed [1]
NI = 0         % Non intermitent nodes
centers = 5
%excess = 50
p_exponent = 2

%variabilities = 0.5 %linspace(0,0.8,5) %[0.0]
len_on_par = [0,0.2,0.3]
%len_on_par = [1, 0.08, 2, 1.06]

check_feasibility=true
foodweb_richness=45
foodweb_connectance=0.2
foodweb_basals_frac=0.2
foodweb_std=0.25
foodweb_bcomp=0
foodweb_temperatures = linspace(1e-3, 0.7, 5) ;
foodweb_diagonals =    linspace(-0.2, -1, 5) ;


for patch_qty = PQ
    for land_edges = 2*patch_qty %LE
        base_name = sprintf('p%dle%dS%dc%dx%du_var%g',patch_qty,land_edges,foodweb_richness,centers,excess,variabilities) ;
        
        for foodweb_temperature = foodweb_temperatures
            for foodweb_diagonal = foodweb_diagonals
                
                name = sprintf( '%s_T%1.3g_D%1.3g_px%1.3g_CON%g', base_name, foodweb_temperature, foodweb_diagonal, p_exponent, foodweb_connectance ) ;
                
                f=fopen(['condor_jobs_' name], 'w') ;
                fprintf(f,'getenv = true\n') ;
                fprintf(f,'universe = vanilla\n') ;
                fprintf(f,'executable = ./run_sims.sh\n');
                fprintf(f,'output = ./log/$(Cluster).$(Process).out\n');
                fprintf(f,'error = ./log/$(Cluster).$(Process).err\n');
                fprintf(f,'log = ./log/$(Cluster).$(Process).log\n');
                %fprintf(f,'requirements = (Machine =!= "node06.econets.org")\n') ;
                                
                rng_seed = 1 ;
                landscape_seed = 1000 ;
                foodweb_seed = 2000 ;
                event_seed = 3000 ;
                
                mig_rate_pars = {} ;
                for a = MR_a
                    for b = MR_b
                        for t = MR_t
                            mig_rate_pars = {mig_rate_pars{:},[a,b,t]} ;
                        end
                    end
                end
                
                for i = 1:length(mig_rate_pars)
                    MR_pars = mig_rate_pars{i};
                    for ext_thr = ET
                        
                        for packet = 0:4
                            format_s = 'arguments = ''$(Cluster).$(Process)'' %d [%f,%f] {[%f,%f,%f]} ';
                            format_s = [format_s, '['  repmat('%f,', 1, length(variabilities)-1) '%f] '] ;
                            format_s = [format_s, '[%f,%f,%f,%f,%d,%d,%f,%f] ''%s'' %d %d %d %s '] ;
                            format_s = [format_s, '[%d,' repmat('%f,', 1, length(len_on_par)-2) '%f] %d'] ;
                            format_s = [format_s, ' %d %f %f %f %f %f %f'] ;
                            fprintf(f, format_s,...
                                rng_seed,...
                                ext_thr{1},...
                                MR_pars,...
                                variabilities,...
                                [topology,patch_qty,NI,landscape_seed,land_edges,centers,excess,p_exponent],...
                                name,...
                                packet,...
                                foodweb_seed,...
                                event_seed,...
                                mat2str(check_feasibility),...
                                len_on_par,...
                                foodweb_richness,...
                                foodweb_connectance,...
                                foodweb_temperature,...
                                foodweb_diagonal,...
                                foodweb_basals_frac,...
                                foodweb_std,...
                                foodweb_bcomp...
                                );
                            fprintf(f,'\nqueue\n');
                            rng_seed = rng_seed + 1 ;
                            sample_size = 10 ; %Should be at least as large as sample_size as defined in main
                            landscape_seed = landscape_seed + sample_size ;
                            foodweb_seed = foodweb_seed + sample_size ;
                            event_seed = event_seed + sample_size ;
                        end
                    end
                end
            end
        end
        fclose(f);
    end
end
%{
function main_cluster_v12(rng_seed,...
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
%}
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
%                source, 6:Star, 7:fractal, 8:modular.
%                For 1, the single source is at an  endpoint (patch 1).
%                For 5, and 6 there is a single central source.
%                For 2,3 and 4 the source is at a random node. For 8, the
%                source is at node 1.
%                In all cases (except 8) the non intermitent patches are
%                chosen with uniform probability among the non-source patches.
%                For 8, the non intermitent patches are 2..M+1. This
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
%             land_par(6): Only used for modular landscapes. It is the
%                number of cluster centers.
%             land_par(7): nodes excess factor for modular landscapes

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
