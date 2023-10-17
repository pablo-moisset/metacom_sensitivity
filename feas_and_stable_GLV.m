function [success, A, CM, R, X_eq] = ...
       feas_and_stable_GLV(connectance, MASTER_P, richness, ...
                           std_off_d, max_iters)
% Feasible and stable generalized LV system
% A is a 0-1 richness x richness adj. matrix. A(i,j)==1 iff j preys on i
% CM is a the community matrix (real richness x richness)
% R is the intrinsic growth rate (real richness x 1)
% X_eq are the abundances at equilibrium
% connectance is a real between 0 and 1
% CON is the diet contiguity (real between 0 and 1)
% richness (positive integer) is the number of species
% std_off_d (positive real) is the standard deviation non-zero off-diagonal elements of CM
% max_r (positive real) is an upper bound for the intrinsic growth rate.
% min_mort (positive real) is an upper bound for the mortality rate of
%    non-basal specie.
% max_iters positive integer, stating the maximum number of attempts at
% creating a stable and feasible matrix. After that many tries,
% the function exits unsuccesfully

tries = 0 ;

success = false ;

while tries < max_iters && ~success
    feasible = false ; is_stable = false ;
    
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
        CM = CM.*abs(normrnd(0, std_off_d, length(CM))) ; % TODO: Maybe we should set a minimum interaction strength
        new_diag = basals*MASTER_P.basals_diag + (1-basals)*MASTER_P.diagonal ;
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
        
        %CM = A.*(CM+std_off_d*normrnd(0, 1, length(CM) )); %GPPM style
        %CM = CM*std_off_d.*abs(normrnd(0, 1, length(CM))) ; % May style
        
        ln_mu = -0.5*log(1+std_off_d^2) ;
        ln_sigma = sqrt(log(std_off_d^2+1)) ;
        LN = lognrnd(ln_mu, ln_sigma,length(CM)) ; %lognormal mean=1, std=std_off_d
        CM = CM.*LN ;
        
        CM = CM - MASTER_P.feeding_efficiency*CM' ;
        
        CM(basals,basals) = -MASTER_P.basals_comp*LN(basals,basals) ; %Add interspecies competition for basals
    else
        assert(false) ;
    end
    
    new_diag = basals*MASTER_P.basals_diag + (1-basals)*MASTER_P.diagonal ;
    CM = CM - diag(diag(CM)) + diag(new_diag) ; %replace diagonal
    
    [feasible, X_eq, R] = is_feasible(CM, MASTER_P.max_r, MASTER_P.min_mort) ;
    
    
    if feasible
        is_stable = -max(real(eig(diag(X_eq)*CM))) > 0 ;
    end
    
    tries = tries + 1 ;
    success = feasible && is_stable ;
end
end