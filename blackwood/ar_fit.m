function [params, noise] = ar_fit(x, order, varargin)
% ar_fit: fit an AR model of order (order) to the given data (x).
% Pads start of data if necessary to avoid an error (but warns if doing so).
%
% Syntax: params = ar_fit(x, order, [method]='burg', [mean_models]=false)
%
% If mean_models is set to true, averages the AR model estimates of each
% column of the matrix x. Otherwise, if x is a matrix, returns a separate
% model for each column in the rows of 'params'.
%
% Choices for method are 'burg' and 'yule' (Yule-Walker).
% Ani: edited to return estimate of noise
if isvector(x)
    x = x(:);
end
[lenx, nx] = size(x);

noptarg = length(varargin);

arFunc = @arburg;
if noptarg > 0 && ~isempty(varargin{1})
    method = varargin{1};
    if strcmp(method, 'yule')
        arFunc = @aryule;
    else
        assert(strcmp(method, 'burg'), 'Invalid method argument "%s"', method);
    end
end

mean_models = false;
if noptarg > 1 && ~isempty(varargin{2})
    mean_models = varargin{2};
end

npad = max(0, order + 1 - lenx);
if npad ~= 0
    warning('Padding data by %d samples to fit AR model', npad);
end
x = [zeros(npad, nx); x];

[params,noise] = arFunc(x, order);
if mean_models
    params = mean(params, 1);
end

end
