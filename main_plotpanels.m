close all
h11 = openfig('p_alpha_MOD_p50S45c5x50px2_CON0.2.fig') ;
h12 = openfig('p_beta_MOD_p50S45c5x50px2_CON0.2.fig') ;
h13 = openfig('p_gamma_MOD_p50S45c5x50px2_CON0.2.fig') ;
h21 = openfig('p_alpha_MOD_p50S75c5x50px2_CON0.25.fig') ;
h22 = openfig('p_beta_MOD_p50S75c5x50px2_CON0.25.fig') ;
h23 = openfig('p_gamma_MOD_p50S75c5x50px2_CON0.25.fig') ;
h31 = openfig('p_alpha_MOD_p10S45c5x5px2_CON0.2.fig') ;
h32 = openfig('p_beta_MOD_p10S45c5x5px2_CON0.2.fig') ;
h33 = openfig('p_gamma_MOD_p10S45c5x5px2_CON0.2.fig') ;

figure
h(1)=subplot(3,3,1);xticks(1:9)
ylim([-0.3 0.101])
h(2)=subplot(3,3,2);xticks(1:9)
ylim([-4 2])
h(3)=subplot(3,3,3);xticks(1:9)
ylim([-0.12 0.18])
h(4)=subplot(3,3,4);xticks(1:9)
ylim([-0.3 0.101])
h(5)=subplot(3,3,5);xticks(1:9)
ylim([-4 2])
h(6)=subplot(3,3,6);xticks(1:9)
ylim([-0.12 0.18])
h(7)=subplot(3,3,7);xticks(1:9)
ylim([-0.3 0.101])
h(8)=subplot(3,3,8);xticks(1:9)
h(9)=subplot(3,3,9);xticks(1:9)
ylim([-0.12 0.18])

% Paste figures on the subplots
copyobj(allchild(get(h11,'CurrentAxes')),h(1));
copyobj(allchild(get(h12,'CurrentAxes')),h(2));
copyobj(allchild(get(h13,'CurrentAxes')),h(3));
copyobj(allchild(get(h21,'CurrentAxes')),h(4));
copyobj(allchild(get(h22,'CurrentAxes')),h(5));
copyobj(allchild(get(h23,'CurrentAxes')),h(6));
copyobj(allchild(get(h31,'CurrentAxes')),h(7));
copyobj(allchild(get(h32,'CurrentAxes')),h(8));
copyobj(allchild(get(h33,'CurrentAxes')),h(9));
