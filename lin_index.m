% Conceptually, the inverse of linspace function
% i is the index at which x should be in linspace(xmin, xmax, N)
% or as close as possible

function i = lin_index(x, xmin, xmax, N)
   i = round((x-xmin+realmin)/(xmax-xmin+realmin)*(N-1))+1 ;
end
