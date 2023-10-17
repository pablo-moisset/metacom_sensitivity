function N = cont_p_norm(X, t, p)
% Continuous time p norm of a time series
% t is a 1d matrix with timestamps
% X is a matrix with numel(t) rows. Size(X,2) is the number of univariated time series
% For the time series, we assume variable i takes value X(k,i) at time t(k)
% and does not change it until t(k+1)
% N is a 1 x size(X,2) matrix. N(i) is the p-norm of time series i.
% This version misses the division by (t(end)-t(1)) and the raising to 1/p
% it is meant to process time series by blocks to save memory

   X = abs(X(1:end-1,:)) ;
   N = zeros(1, size(X,2)) ;
   for i = 1:size(X,2)
      N(i) = sum(X(:,i).^p.*(t(2:end)-t(1:end-1)));
   end
end
