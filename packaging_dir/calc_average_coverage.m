function avg = calc_average_coverage(cov, func, include_zeros, is_log)
% calculate the average of the coverage (cov)
% func can be @mean, @median, %mode
% is_log indivates if to aveeage the log or the converage itself
% include_zero=true -> to average also on zeros (only when is_log=false)

if nargin<2
    func = @median;
end

if nargin<3 
    include_zeros = false;
end

if nargin<4
    is_log = false;
end

if ~include_zeros
    cov = cov(cov>0);
end

if is_log
    avg = exp(func(log(cov)));
else
    avg = func(cov);
end

end
