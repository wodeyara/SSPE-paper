function [phi, Q, M] = genParametersSoulatMdl_sspp(freqs, Fs,ampVector, sigmaFreqs)

    rotMat = zeros(length(freqs), 2,2);
    varMat = zeros(length(freqs),2,2);
    cnt = 1;

    for freqRunThrough = freqs
        rotMat(cnt,:,:) = createRotMat(freqRunThrough,Fs);
        rotMat(cnt,:,:) = rotMat(cnt,:,:) * ampVector(cnt);
        varMat(cnt, :,:) = [sigmaFreqs(cnt), 0; 0 , sigmaFreqs(cnt)];
        cnt = cnt + 1;
    end
    phi = stackBlockMat(rotMat);
    Q = stackBlockMat(varMat);

    M = reshape([ones(length(freqs),1) , zeros(length(freqs),1)]', length(freqs)*2,1); % for extracting real part of signal from analytic 'x'

     function rotMat = createRotMat(freq, Fs)
        % create the appropriate rotation matrix for the frequency
        rotMat = [cos(2*pi*freq/Fs), -sin(2*pi*freq/Fs); ...
                  sin(2*pi*freq/Fs), cos(2*pi*freq/Fs)];
    end        

    function fullBlockMat = stackBlockMat(allMat)
        % pull all the rotation matrices together into a block diagonal
        % structure -> improves efficiency in computation
        for i  = 1:size(allMat,1)
           fullBlockMat(size(allMat,2) * (i-1) + 1:size(allMat,2) * (i), ...
               size(allMat,2) * (i-1) + 1:size(allMat,2) * (i)) = squeeze(allMat(i,:,:)); %don't do squeeze, do a reshape and extract points we want.
        end
    end

end