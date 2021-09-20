function [x_ext, params, noise,predError] = ar_extend(x, nsamp, varargin)
% ar_extend: extend x by nsamp samples using an AR model.
%
% Syntax: [x_ext, params] = ar_extend(x, nsamp, [nsamp_ar]=<size(x,1)>, [order]=20, ...
%                                    [mean_models]=false, [method]='burg')
%
% If x is a vector, a column vector is returned. Otherwise, operates along
% the columns of x (fitting a separate model to each column).
%
% If nsamp_ar is provided, uses only the last nsamp_ar samples of x to
% calculate the AR model. nsamp_ar (or length(x) if not specified) must be
% at least order+1.
%
% If order is provided, uses an AR(order) model.
%
% If mean_models is set to true, averages the AR model estimates of each
% column of the matrix x and uses this mean model to predict forward for
% each column. (Experimental feature!)
%
% 'method', if provided, specifies the AR estimation method - either 'burg'
% (the default) or 'yule' (for Yule-Walker).
%
% If requested, 'params' contains the AR model(s) used for prediction
% (formatted as in the output of arburg). If mean_models is true, this is
% the single average model.
% 
%Ani:  Now generating a prediction error as well to have confidence bounds on
% estimated phase


if isvector(x)
    x = x(:);
end
lenx = size(x, 1);

noptarg = length(varargin);

nsamp_ar = lenx;
if noptarg > 0 && ~isempty(varargin{1})
    nsamp_ar = varargin{1};
end

order = 20;
if noptarg > 1 && ~isempty(varargin{2})
    order = varargin{2};
end

mean_models = false;
if noptarg > 2 && ~isempty(varargin{3})
    mean_models = varargin{3};
end

method = 'burg';
if noptarg > 3 && ~isempty(varargin{4})
    method = varargin{4};
end

% use only nsamp_ar points to fit, if requested
x_clipped = x(lenx - (nsamp_ar - 1):lenx, :);

% fit model
[params, noise] = ar_fit(x_clipped, order, method, mean_models);

% extend signal using model
[x_ext,predError] = ar_predict(x, nsamp, params, noise);
x_ext = [x; x_ext];
end

