clear; close all; clc;

%% testOptions
plotFrequencyFilter = true;
plotComparisonData = true;

%% selected subjects and trials
subjectList = 1;
trialList = 1;  

for subjectID = subjectList
    fprintf('\n---------\nSubject : %d',subjectID);
    for trialID = trialList
        fprintf('\nTrial : %d\n ',trialID);
         
        
        %% load raw data
        load('./experiments/humanFixedBase/data/VICONsaveDataGen16.mat')
        load('./experiments/humanFixedBase/data/imuExtractedDataGen16.mat');
        
        % ====force plate
        f_raw = subjectData(subjectID,trialID).analogsFOR;
        mom_raw = subjectData(subjectID,trialID).analogsMOM;
        
        % ====IMU 
        accl_raw = imuData(subjectID,trialID).accln;
        omega_raw = imuData(subjectID,trialID).gyro;

        % ====MARKERS  --> TO DO EVENTUALLY!
        
        
        %% compute FFT
        
        for i = 1:3
            data = omega_raw(:,i);
            len = length(data);
            Fs = 1e+3;              % sampling frequency
            Ts = 1/Fs;              % sampling time
            time = 0:Ts:(Ts)*(len-1);

            % fft of data
            data_fft = fft(data);
            freq = Fs*(0:(len/2))/len; 

            % shift data_fft as the FFT is symmetric (fftshift for plotting data)
            data_fft = fftshift(data_fft);

        %       data_fft = data_fft(2:((len/2+1)));
        %       data_fft(1:end-1) = 2*data_fft(1:end-1); %throw away first peak



            %% filter in frequency domain 

            % ====build a gaussian low pass filter
            cutOffFreq = 20; %20Hz
            HLPF_gauss=(fftshift(lpfilter('gaussian', len, 1, cutOffFreq)));

            % ====apply gaussian filter to data_fft
            data_fftFilt = HLPF_gauss.* data_fft;

            plot(HLPF_gauss);
            axis tight;
            grid on;

            if (plotFrequencyFilter)
                fig = figure();
                axes1 = axes('Parent',fig,'FontSize',16);
                box(axes1,'on');
                hold(axes1,'on');
                grid on;

                subplot(211);
                plot(abs(data_fft)/len, 'lineWidth',2.0);hold; % plotting the abs value/magnitude of data_fft
                plot(HLPF_gauss, 'lineWidth',2.0)
                leg = legend('dataFFT','gaussian LPF','Location','southeast');
                set(leg,'FontSize',12);
                xlabel('Frequeny [Hz]','FontSize',12);
                axis tight;
                grid on;
                title(sprintf('Subject : %d, Trial : %d',subjectID,trialID));

                subplot(212);
                plot((abs(data_fftFilt)/len), 'lineWidth',2.0); % plotting the abs value/magnitude of data_fftFilt
                leg = legend('filtered data','Location','southeast');
                set(leg,'FontSize',12);
                xlabel('Frequeny [Hz]','FontSize',12);
                axis tight;
                grid on;
            end

            %% compute IFFT

            %the inverse of FFT needs to be done on unshifted data
            dataAfterFFT = ifft(ifftshift(data_fftFilt));

            if (plotComparisonData)
                fig = figure();
                axes1 = axes('Parent',fig,'FontSize',16);
                box(axes1,'on');
                hold(axes1,'on');
                grid on;

                plot(time,data,time,dataAfterFFT, 'lineWidth',2.0);hold;
                leg = legend('original data','anti-transformed data','Location','southeast');
                set(leg,'FontSize',12);
                xlabel('Time [s]','FontSize',12);
                ylabel('Force [N]','FontSize',12);
                axis tight;
                grid on;
                title(sprintf('Subject : %d, Trial : %d',subjectID,trialID));
            end
        end 
     end
     fprintf('\n');
end
