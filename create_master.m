function [A, CM, R, X_eq, stability] = create_master(connectance, MASTER_P, richness, ...
                                          diagonal, std_off_d)
%
if MASTER_P.algorithm == 0 %Generalized niche model
    A = gen_niche_model(connectance, richness, MASTER_P.contiguity) ;
    f = tril(A+A')==2 ;
    
    A(f) = 0 ; % Avoid mutual predation leaving only one interaction
    basals = sum(A)==0 ;
    
    CM = -A ;
    CM = CM - CM' ;
    %    CM(basals,basals) = -1 ; %Add interspecies competition for basals
    warning('Competition for basals is OFF') ;
    %M = M.*rand(length(M)) ; % TODO: Maybe we should set a minimum interaction strength
    CM = CM*std_off_d.*abs(normrnd(0, 1, length(CM))) ; % TODO: Maybe we should set a minimum interaction strength
    new_diag = basals*MASTER_P.basals_diag + (1-basals)*diagonal ;
    CM = CM - diag(diag(CM)) + diag(new_diag) ; %replace diagonal
elseif MASTER_P.algorithm == 1 %GPPM
    basals_qty = round(MASTER_P.fraction_basals*richness) ;
    assert ( basals_qty < richness) ;
    
    A = gppm_coherence(richness,basals_qty,...
        connectance*(richness^2-richness*basals_qty-richness+basals_qty),...
        MASTER_P.temperature) ;
    
    f = tril(A+A')==2 ;
    
    A(f) = 0 ; % Avoid mutual predation leaving only one interaction
    basals = sum(A)==0 ;
    
    CM = -A ;
    mu_ef = 0.4 ; % Magic number. Efficiency
    %CM = A.*(CM+std_off_d*normrnd(0, 1, length(CM) )); %GPPM style
    %CM = CM*std_off_d.*abs(normrnd(0, 1, length(CM))) ; % May style

    ln_mu = -0.5*log(1+std_off_d^2) ;
    ln_sigma = sqrt(log(std_off_d^2+1)) ;
    LN = lognrnd(ln_mu, ln_sigma,length(CM)) ; %lognormal mean=1, std=std_off_d
    CM = CM.*LN ;

    CM = CM - mu_ef*CM' ;

    CM(basals,basals) = -MASTER_P.bcomp*LN(basals,basals) ; %Add interspecies competition for basals
else
    assert(false) ;
end

new_diag = basals*MASTER_P.basals_diag + (1-basals)*diagonal ;
CM = CM - diag(diag(CM)) + diag(new_diag) ; %replace diagonal

[feasible, X_eq, R] = is_feasible(CM, MASTER_P.max_r, MASTER_P.min_mort) ;


stability = -max(real(eig(diag(X_eq)*CM))) ;
end

