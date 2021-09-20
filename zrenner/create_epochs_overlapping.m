function epochs = create_epochs_overlapping(data,Fs)

    N = length(data);
    cnt =1;
    for i = 1:(length(data) - (Fs*2))+1
        epochs(:,cnt) = data(i: (i-1) + Fs*2);
        cnt = cnt + 1;
    end
%     epochs = reshape(data,Fs*2,N/(2*Fs));
    
end