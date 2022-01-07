function s=ifty(fs);
% ========================================================================
% Chirp Scaling Algorithm - Inverse Fourier Transform in Azimuth

% ========================================================================
s=ifftshift(ifft(fftshift(fs,2),[],2),2);
% s=fftshift(ifft(fftshift(fs.'))).';

% = end ==================================================================
