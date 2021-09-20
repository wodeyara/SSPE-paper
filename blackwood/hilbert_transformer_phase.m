function [phase, estimate_mask, analytic,analytic_buf] = hilbert_transformer_phase(data, buffer_len, varargin)
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
%    estimate_mask: true for indices of source data with corresponding phase estimates.
%                   (should be all true if 'upsample' is true)
%    analytic:      analytic signal (same size as "phase" output)

assert(isvector(data), 'Data must be a vector');

% sample rates (#have changed them from 30000 and 500)
Fs = 1000;
Fs_ds = 1000;

if length(varargin) >= 3 && ~isempty(varargin{3})
    Fs = varargin{3};
    assert(isscalar(Fs) && isreal(Fs) && Fs > 0, 'Invalid sample rate');
    assert(mod(Fs, Fs_ds) == 0, 'Sample rate must be a multiple of 500');
end

ds_factor = Fs / Fs_ds;

% bandpass filter
band = [4 8];
if length(varargin) >= 2 && ~isempty(varargin{2})
    band = varargin{2};
    assert(isvector(band) && isnumeric(band) && isreal(band) && length(band) == 2 ...
        && band(1) > 0 && band(2) > band(1) && band(1) < Fs_ds/2, ...
        'Invalid frequency band.');
    band = band(:);
end
[b, a] = butter(2, band/(Fs/2));
data_filt = filter(b, a, data); %
% also there's no correction for phase alteration inherent in this
% filtering step.  

% hilbert transformer
if ~isempty(varargin) && ~isempty(varargin{1})
    ht_b = varargin{1};
    ht_order = length(ht_b) - 1;
else
    ht_order = 18; % farr too low. The response in theta range is very weak.
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

% output setup
estimate_mask = false(size(data));
estimate_mask(1:ds_factor:end) = true;
estimate_mask_buf = buffer(estimate_mask, buffer_len);

ht_state = zeros(ht_order, 1);
data_buf = buffer(data_filt, buffer_len);

analytic_buf = zeros(size(data_buf));

max_n_samp = ceil(buffer_len / ds_factor);

ar_out = zeros(ar_buf_len + ht_offset_samp + 1, 1);
ht_buf = zeros(max_n_samp + ht_offset_samp + 1, 1);
ht_out = ht_buf;

% for interpolation
do_upsamp = false;
if length(varargin) >= 4 && varargin{4}
    do_upsamp = true;
    
    lastcs = 0;
    lastcs_offset = 1 - ds_factor;
    
    estimate_mask = true(length(data), 1);
    valid_samp = buffer(estimate_mask, buffer_len);
end

for kBuf = 1:size(data_buf, 2)
    % get downsampled samples in this buffer
    curr_mask = estimate_mask_buf(:, kBuf);
    buf_ds = data_buf(curr_mask, kBuf);
    n_samp = length(buf_ds);
    ht_buf_offset = max_n_samp - n_samp;
    
    % push new samples to ar buffer
    ar_buf(1:end-n_samp) = ar_buf(1+n_samp:end);
    ar_buf(end-n_samp+1:end) = buf_ds;
    
    % extend by the hilbert transformer offset
    [ar_out(:)] = ar_extend(ar_buf, ht_offset_samp + 1, [], ar_order);
    ht_buf(ht_buf_offset + 1:end) = ar_out(end-(n_samp+ht_offset_samp+1)+1:end);

    
    % apply hilbert transformer (2 steps)
    ht_out(:) = 0;
    % save state at t0
    [ht_out(ht_buf_offset + 1:max_n_samp), ht_state] = filter(ht_b, 1, ht_buf(ht_buf_offset + 1:max_n_samp), ht_state);
    % apply to extension w/o saving state
    ht_out(max_n_samp+1:end) = filter(ht_b, 1, ht_buf(max_n_samp+1:end), ht_state);
    
    % get analytic signal
    % problem I see with this is that we are including the AR samples when the way I now understand this
    % code is that it is intended to offset the fact that a filter
    % introduces a group delay. We shouldn't be including the AR predicted
    % samples. Not to mention this is not placing the timing correctly.
    analytic_buf(curr_mask, kBuf) = complex(buf_ds, scale_factor * ht_out(end-n_samp:end-1)); 
    
    if do_upsamp        
        x_interp = (lastcs_offset : ds_factor : lastcs_offset + (n_samp+1)*ds_factor)';
        lastcs_offset = x_interp(end-1) - buffer_len;
                
        buf_interp = [lastcs;
                      analytic_buf(curr_mask, kBuf);
                      complex(ar_out(ar_buf_len + 1), ht_out(end))];                  
        lastcs = buf_interp(end-1);
        
        xq = find(valid_samp(:, kBuf));
        analytic_buf(xq, kBuf) = ...
            interp1(x_interp, abs(buf_interp), xq) ...
            .* exp(1j * interp1(x_interp, unwrap(angle(buf_interp)), xq));                    
    end
end

% reshape and mask phase
analytic = analytic_buf(:);
analytic = analytic(estimate_mask);
phase = angle(analytic);

end

