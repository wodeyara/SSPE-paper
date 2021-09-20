function [phase, lowerPhase, upperPhase, analytic] = hilbert_transformer_causal(data, buffer_len, band, varargin)
% hilbert_transformer_phase: estimate real-time phase using a hilbert transformer,
% using AR prediction to offset the group delay.
%
% Syntax: [phase, estimate_mask, analytic] = hilbert_transformer_phase(...
%   data, buffer_len, [ht_b], [band], [Fs], [upsample]);
%
% Input:
%    data:       raw data at 30 kHz
%    buffer_len: number of 30 kHz samples per buffer
%    [ht_b]:     (designed using firpm) numerator of (FIR) Hilbert
%                transformer to use, at 500 Hz sample rate.
%    [band]:     ([4 8]) target frequency band, in Hz
%    [Fs]:       (30000) sample rate in Hz, must be a multiple of 500.
%    [upsample]: (false) whether to upsample using linear interpolation
%
% Output:
%    phase:         estimated phase in radians
%    lower/upperPhase: estimated confidence bounds for phase
%    analytic:      analytic signal (same size as "phase" output)

assert(isvector(data), 'Data must be a vector');

% sample rates (#have changed them from 30000 and 500)
Fs = 1000;
Fs_ds = 1000;

% if length(varargin) >= 3 && ~isempty(varargin{3})
%     Fs = varargin{3};
%     assert(isscalar(Fs) && isreal(Fs) && Fs > 0, 'Invalid sample rate');
%     assert(mod(Fs, Fs_ds) == 0, 'Sample rate must be a multiple of 500');
% end

ds_factor = Fs / Fs_ds;

% bandpass filter
if length(varargin) >= 2 && ~isempty(varargin{2})
    band = varargin{2};
    assert(isvector(band) && isnumeric(band) && isreal(band) && length(band) == 2 ...
        && band(1) > 0 && band(2) > band(1) && band(1) < Fs_ds/2, ...
        'Invalid frequency band.');
    band = band(:);
end
[b, a] = butter(2, band/(Fs/2));
data_filt = filter(b, a, data); 

% also there's no correction for phase alteration inherent in this
% filtering step 

% hilbert transformer
if ~isempty(varargin) && ~isempty(varargin{1})
    ht_b = varargin{1};
    ht_order = length(ht_b) - 1;
else
    ht_order = 18; % changed it to a lower value since this improved performance
    ht_freq = [band(1) (Fs_ds/2)-band(1)] / (Fs_ds/2);
    ht_b = firpm(ht_order, ht_freq, [1 1], 'hilbert');
end

ht_offset_samp = ht_order / 2;
f_resp = freqz(ht_b, 1, (band(1):0.1:band(2)) * 2 * pi / Fs_ds);
mag_resp = abs(f_resp);
% normalize by geometric mean of magnitude response range
scale_factor = 1/sqrt(min(mag_resp) * max(mag_resp));

% ar prediction (want to save ~1 second of data for training)
ar_order = 5; % ALTERED FROM 20 -- note that thinking in the oscillator fromework an order of 5 allows two potential oscillators.
ar_buf_len = max(Fs_ds, 2 * ar_order);
ar_buf = zeros(ar_buf_len, 1);

% state setup for filter
ht_state = zeros(ht_order, 1);
ht_state_upper = zeros(ht_order,1);
ht_state_lower = zeros(ht_order,1);

%outputs
analytic_buf = zeros(size(data_filt)-buffer_len);
analytic_buf_upper = zeros(size(data_filt)-buffer_len);
analytic_buf_lower = zeros(size(data_filt)-buffer_len);

max_n_samp = ceil(buffer_len / ds_factor);

ar_out = zeros(ar_buf_len + ht_offset_samp + 1, 1);
ht_buf = zeros(max_n_samp + ht_offset_samp + 1, 1);
ht_out = ht_buf;

for kBuf = 1+buffer_len:length(data_filt) % runs across all data points
    % get downsampled samples in this buffer
    buf_ds = data_filt(kBuf-buffer_len+1:kBuf);
    n_samp = length(buf_ds);
    ht_buf_offset = max_n_samp - n_samp;
    
    % push new samples to ar buffer
    ht_buf = buf_ds;
    ar_buf = buf_ds;
    
    % extend by the hilbert transformer offset
    % only fit the AR model every second:
    if mod(kBuf-1,1000) ==0
        [ar_out(:),params, noise,predError] = ar_extend(ar_buf, ht_offset_samp + 1, [], ar_order);
        ht_buf = ar_out;
        x_ext = ar_out(end-ht_offset_samp:end);
    else
        % else predict out the AR model at every new sample
        [x_ext,predError] = ar_predict(ht_buf, ht_offset_samp+1, params, noise);
        ht_buf = [ht_buf; x_ext];
    end

    ht_buf_upper =[ht_buf(1:end-ht_offset_samp-1); x_ext + sqrt(predError)*1.96];
    ht_buf_lower =[ht_buf(1:end-ht_offset_samp-1); x_ext - sqrt(predError)*1.96];
    
    % apply hilbert transformer (2 steps)
    ht_out(:) = 0;
    % save state at t0
    [ht_out(ht_buf_offset + 1:max_n_samp), ht_state] = filter(ht_b, 1, ...
        ht_buf(ht_buf_offset + 1:max_n_samp),ht_state);
    
    % now do this for the upper and lower bounds of the estimate separately
    [ht_out_upper(ht_buf_offset + 1:max_n_samp),ht_state_upper] = filter(ht_b, ...
        1, ht_buf_upper(ht_buf_offset + 1:max_n_samp),ht_state_upper);
    
    [ht_out_lower(ht_buf_offset + 1:max_n_samp),ht_state_lower] = filter(ht_b, ...
        1, ht_buf_lower(ht_buf_offset + 1:max_n_samp),ht_state_lower);
    
    % apply to extension w/o saving state
    [ht_out(max_n_samp+1:max_n_samp+1+ht_offset_samp)] = filter(ht_b, 1, ...
        ht_buf(max_n_samp+1:max_n_samp+1+ht_offset_samp),ht_state);
    
    [ht_out_upper(max_n_samp+1:max_n_samp+1+ht_offset_samp)] = filter(ht_b, ...
        1, ht_buf_upper(max_n_samp+1:max_n_samp+1+ht_offset_samp),ht_state_upper);
    
    [ht_out_lower(max_n_samp+1:max_n_samp+1+ht_offset_samp)] = filter(ht_b, ...
        1, ht_buf_lower(max_n_samp+1:max_n_samp+1+ht_offset_samp),ht_state_lower);

    
    % get analytic signal and confidence bounds
    % note that analytic signal is coming only from the AR forecasted
    % samples alone
    analytic_buf(kBuf-1000+1) = complex(ht_buf(max_n_samp+1), scale_factor * ht_out(max_n_samp+1+ht_offset_samp));
    % the following values seem to severely underestimate the confidence limits so I don't trust them
    analytic_buf_upper(kBuf-1000+1) = complex(ht_buf_upper(max_n_samp+1), scale_factor * ht_out_upper(max_n_samp+1+ht_offset_samp));
    analytic_buf_lower(kBuf-1000+1) = complex(ht_buf_lower(max_n_samp+1), scale_factor * ht_out_lower(max_n_samp+1+ht_offset_samp));
    
end

% reshape and get phase estimates
analytic  =[zeros(1,length(data_filt) - length(analytic_buf)),analytic_buf];
phase = angle(analytic);
lowerPhase =[zeros(1,length(data_filt) - length(analytic_buf)),angle(analytic_buf_lower)];
upperPhase = [zeros(1,length(data_filt) - length(analytic_buf)),angle(analytic_buf_upper)];
end

