function [to_add, to_remove] = eq_and_try_add(to_idx, species_idx, M,...
                                              check_feasibility)
% [to_add, to_remove] logical arrays 
% present: logical array
% species_idx: integer
% M: metacommunity structure
% r: community intrinsic growth rates.
% check_feasibility: shoud feasibility be checked for attempted additions 
% because of migrations?
% ext_thr (scalar) is the extintion threshold for the to_idx patch

present = M.present(to_idx,:) ;
CM = M.CM ;
r = M.R ;

assert(size(present,2)==size(CM,1) && size(CM,1)==size(CM,2)) ;
assert(size(CM,1)==size(r,1) && size(r,2)==1) ;
assert(~present(species_idx)) ;

if ~check_feasibility %skip fesibility test. Migration attempt is successful
    to_add = present & false ;
    to_remove = to_add ;
    to_add(species_idx) = true ;
else
    keep_testing = true ;
    
    new_present = present ;
    assert(~new_present(species_idx)) ;
    new_present(species_idx) = 1 ;
    
    while keep_testing
        sub_CM = CM(new_present, new_present) ;
        sub_r = r(new_present) ;
        %    if rank(sub_CM) < length(sub_CM)
        %        warning('sub_CM is not a full rank matrix') ;
        %    end
        X_eq = -sub_CM\sub_r ;
        non_positive = find(X_eq <= M.ext_thr(to_idx)) ;
        if ~isempty(non_positive)
            idx = find(new_present) ;
            %       fprintf('Removing %d\n', idx(non_positive));
            new_present(idx(non_positive)) = false ;
            keep_testing = ~all(~new_present) ; %there is at least a species left
        else
            keep_testing = false ;
        end
    end
    
    to_add    = (new_present - present) > 0 ;
    to_remove = (new_present - present) < 0 ;
end
end
