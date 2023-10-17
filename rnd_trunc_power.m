function X = rnd_trunc_power(a,lo,hi,rows,cols)
%    [a,lo,hi,rows,cols]
    X = randraw('pareto', [lo a] , rows, cols) ;

    too_big = find(X>hi) ;
    while ~isempty(too_big)
       fprintf('too big = %d\n', sum(too_big)) ;
       fprintf('retrying\n')
       X(too_big) = randraw('pareto', [lo a] , numel(too_big), 1) ;
       too_big = find(X>hi) ;
    end
end