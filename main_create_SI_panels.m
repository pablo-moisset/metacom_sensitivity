L_b_beta.minval  = -0.8; L_b_beta.maxval = 1; 
L_b_gamma.minval = -300; L_b_gamma.maxval = 50; 
L_p_alpha.minval = -0.35; L_p_alpha.maxval = 0.05; 
L_p_beta.minval  = -18; L_p_beta.maxval = 5;
L_p_gamma.minval = -0.12; L_p_gamma.maxval = 0.22; 
L_shannon_alpha.minval = -0.22; L_shannon_alpha.maxval = 0.06; 
L_shannon_beta.minval  = -0.1; L_shannon_beta.maxval = 0.25; 
L_shannon_gamma.minval  = -0.04; L_shannon_gamma.maxval = 0.1; 

subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.1 0.01], [0.1 0.01]);

if false
filename = 'landscape_panel_data.csv' ;

T = csvimport(filename,'noHeader',true) ;

response_vars  = unique(T(:,1)) ;
patch_qtys     = {'p10','p25','p50','p100'} ;
excess_factors = {'x5','x10','x50'} ;
for i = 1:numel(response_vars)
    response_var = response_vars{i}
    %figure starts here
    figure('Name', response_var) ;

    plot_number = 1 ;
    for j = 1:numel(patch_qtys)
        patch_qty = patch_qtys{j}
        for k = 1:numel(excess_factors)
            excess_factor = excess_factors{k}
            idx = find(strcmp(T(:,1), response_var));
            Tf = T(idx,:) ;
            idx = find(strcmp(Tf(:,2),patch_qty));
            Tf = Tf(idx,:) ;
            idx = find(strcmp(Tf(:,3),excess_factor));
            Tf = Tf(idx,:) ;
            Tf;
            subplot(4,3,plot_number)
            bar_height = cell2mat(Tf(1,4:12)) ;
            bar_se     = cell2mat(Tf(1,13:end)) ;
            bh=bar(bar_height, 'FaceColor',[100,170,166]/255) ;
            bh.EdgeColor = 'none' ;
            ylim([eval(['L_' response_var '.minval']) eval(['L_' response_var '.maxval'])]) ;
            hold on ;
            h=errorbar(1:numel(bar_height), bar_height, bar_se,bar_se);
            h.LineStyle = 'none' ;
            h.Color = 'Black' ;
            if j==4
               xlabel('Predictor')
            else
                set(gca,'XTickLabel',[]);
            end
            if k==1
               ylabel('Effect size')
            else
                set(gca,'YTickLabel',[]);
            end
            plot_number = plot_number+1 ;
        end
    end
    set(gcf,'position',[10,10,650,650])
    
    export_filename = sprintf('Lpanel_%s',response_var);
    saveas(gcf, [export_filename '.eps'],'epsc');
 %   saveas(gcf, [export_filename '.tif'],'tiff');
end

end

F_b_beta.minval  = -0.46; F_b_beta.maxval = 0.86; 
F_b_gamma.minval = -160; F_b_gamma.maxval = 45; 
F_p_alpha.minval = -0.3; F_p_alpha.maxval = 0.1; 
F_p_beta.minval  = -4.5; F_p_beta.maxval = 2;
F_p_gamma.minval = -0.15; F_p_gamma.maxval = 0.2; 
F_shannon_alpha.minval = -0.15; F_shannon_alpha.maxval = 0.05; 
F_shannon_beta.minval  = -0.1; F_shannon_beta.maxval = 0.25; 
F_shannon_gamma.minval  = -0.07; F_shannon_gamma.maxval = 0.12; 


if true
filename = 'foodweb_panel_data.csv' ;

%T = readtable(filename) ;
T = csvimport(filename,'noHeader',true) ;

response_vars  = unique(T(:,1)) ;
richness_list  = {'S30','S45','S60','S75'} ;
connectances = {'CON0.15','CON0.20','CON0.25'} ;

for i = 1:numel(response_vars)
    response_var = response_vars{i}
    %figure starts here
    figure('Name', response_var) ;
    plot_number = 1 ;
    for j = 1:numel(richness_list)
        richness = richness_list{j};
        for k = 1:numel(connectances)
            connectance = connectances{k};
            idx = find(strcmp(T(:,1), response_var));
            Tf = T(idx,:) ;
            idx = find(strcmp(Tf(:,2),richness));
            Tf = Tf(idx,:) ;
            idx = find(strcmp(Tf(:,3),connectance));
            Tf = Tf(idx,:) ;
            Tf;
            subplot(4,3,plot_number)
            bar_height = cell2mat(Tf(1,4:12));
            bar_se     = cell2mat(Tf(1,13:end)) ;
            bh=bar(bar_height, 'FaceColor',[100,170,166]/255) ;
            bh.EdgeColor = 'none' ;
            ylim([eval(['F_' response_var '.minval']) eval(['F_' response_var '.maxval'])]) ;
            hold on ;
            h=errorbar(1:numel(bar_height), bar_height, bar_se,bar_se);
            h.LineStyle = 'none' ;
            h.Color = 'Black' ;
            if j==4
               xlabel('Predictor')
            else
                set(gca,'XTickLabel',[]);
            end
            if k==1
               ylabel('Effect size')
            else
                set(gca,'YTickLabel',[]);
            end
            plot_number = plot_number+1 ;
        end
    end
    set(gcf,'position',[10,10,650,650])
    
    export_filename = sprintf('Fpanel_%s',response_var);
    saveas(gcf, [export_filename '.eps'],'epsc');
 %   saveas(gcf, [export_filename '.tif'],'tiff');
end
end