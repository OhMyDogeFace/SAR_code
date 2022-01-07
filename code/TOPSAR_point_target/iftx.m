% ========================================================================
% Chirp Scaling Algorithm - Inverse Fourier Transform in Range

% ========================================================================

function s=iftx(fs)

s=ifftshift(ifft(fftshift(fs,1),[],1),1);

% = end ==================================================================