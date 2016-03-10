function [ signal_filt ] = frequencyFilterSignal( signal_nofilt, cutFreq, sensorSamplingFrequency)
% FREQUENCYFILTERSIGNAL computes the filtering of a signal in the frequency domain.
%
%  - signal_nofilt :            no filtered signal in teh form (lengthSignal,3);
%  - cutFreq       :            cutOff frequency of the filter;
%  - sensorSamplingFrequency :  sampling frequency.

signal_filt = zeros(length(signal_nofilt),3);

    for i = 1:3
        signal = signal_nofilt (:,i);
        len = length(signal);
        Fs = sensorSamplingFrequency;            % sampling frequency of the sensor acquisition
        Ts = 1/Fs;                          % sampling time
        time = 0:Ts:(Ts * (len-1));

        % FFT of signal
        signal_fft = fft(signal);
        freq = Fs * (0:(len/2))/len; 

        % shift signal_fft as the FFT is symmetric (fftshift for plotting signal)
        signal_fft = fftshift(signal_fft);
        
        % build a gaussian low pass filter
        HLPF_gauss=(fftshift(lpfilter('gaussian', len, 1, cutFreq)));

        % apply gaussian filter to signal_fft
        signal_fftFilt = HLPF_gauss.* signal_fft;
        
        % compute the inverse FFT of the unshifted signal
        signal_iFFT = ifft(ifftshift(signal_fftFilt));
        
        % fill signal_filt
        signal_filt(:,i) = signal_iFFT;

    end
end
