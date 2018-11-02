%% Script_modulation_spectra
% 
% This script computes the amplitude- and frequency-modulation spectra for
% a sound or group of sounds. Four different spectra are computed (see [1]
% for a complete description):
% - AMa spectrum: classic Fourier AM spectrum (see e.g. [2])
% - AMi spectrum: modulation index spectrum (see e.g. [3])
% - FM spectrum: Fourier FM spectrum (see e.g. [2]) 
% - f0M spectrum: Fourier spectrum of the modulations in f0 trajectory
% The computation of f0 modulation spectra requires the YIN toolbox
% developped by de Cheveigné & Kawahara (http://audition.ens.fr/adc/) [4]
%
% refs:
% [1] A cross-linguistic study of speech modulation spectra. Varnet L,
% Ortiz-Barajas MC, Erra RG, Gervain J, Lorenzi C. J Acoust Soc Am. 2017
% Oct;142(4):1976. doi: 10.1121/1.5006179.
% [2] Effects of age and hearing loss on the relationship between
% discrimination of stochastic frequency modulation and speech perception.
% Sheft S, Shafiro V, Lorenzi C, McMullen R, Farrell C. Ear Hear. 2012
% Nov-Dec;33(6):709-20. doi: 10.1097/AUD.0b013e31825aab15. 
% [3] A review of the MTF concept in room acoustics and its use for
% estimating speech intelligibility in auditoria. T. Houtgast and H. J. M.
% Steeneken. J. Acoust. Soc. Am. 1985; 77:3, 1069-1077
% [4] YIN, a fundamental frequency estimator for speech and music. de
% Cheveigné A, Kawahara H. J Acoust Soc Am. 2002 Apr;111(4):1917-30.
%
% Leo Varnet - 2017

close all
clear all

D = dir('*.wav'); % name of the audiofile (or group of audiofiles) to be processed

flim_gammabank = [70 6700]; % gammatone range (Hz)
flim_spectra = [0.1 500]; % modulation rate range for the AMa, FM and f0M spectra (Hz)
flim_Eoct = [0.1 500]; % modulation rate range for the AMi spectrum (Hz)
NthOct = 3; % width of modulation filters for the AMi spectrum (1/X octave filters) - determines the resolution of the AMi spectrum
N_fsamples = 250; % number of (log-spaced) frequency samples for the AMa, FM and f0M spectra

N_wav = length(D);
fc = ERBlinspace( flim_gammabank(1), flim_gammabank(end), 1 );
Nchan = length(fc); % number of gammatones for the modulation spectra
f_spectra_intervals = logspace(log10(flim_spectra(1)), log10(flim_spectra(2)), N_fsamples+1);
f_spectra = logspace(log10(sqrt(f_spectra_intervals(1)*f_spectra_intervals(2))), log10(sqrt(f_spectra_intervals(end)*f_spectra_intervals(end-1))), N_fsamples);

f_oct = 10.^(log10(flim_Eoct(1)):log10(2^(1/NthOct)):log10(flim_Eoct(2))); 
cutoff_oct = 10.^(log10(flim_Eoct(1))-log10(sqrt(2^(1/NthOct))):log10(2^(1/NthOct)):log10(flim_Eoct(2))+log10(sqrt(2^(1/NthOct)))); 
N_fsamples_oct=length(f_oct);

E_spectrum = nan(Nchan, N_fsamples, N_wav);
Eoct_spectrum = nan(Nchan, N_fsamples_oct, N_wav);
m_spectrum = nan(Nchan, N_fsamples_oct, N_wav); 
FM_spectrum = nan(Nchan, N_fsamples, N_wav);
f0_spectrum = nan(N_fsamples, N_wav);

for i_wav=1:N_wav
    
    fprintf(['\nprocessing stim ' num2str(i_wav) ' of ' num2str(N_wav) '\n']);
    
    %% loading sound
    
    fprintf('sound loading\n');
    NameWav = D(i_wav).name;
    [son, fs] = audioread(NameWav);
    %[son,fs] = wavread(NameWav);
    son = son/sqrt(mean(son.^2));
    Nsamples = length(son);
    
    %% gammatone filtering
    
    [gamma_responses, ~, fc] = gammatone_filtering( son, flim_gammabank(1), flim_gammabank(end), 1, [], fs);
    
    %% AM and FM extraction
  
    fprintf('E/FM extraction\n');
    [E, FM] = hilbert_extraction( gamma_responses', fs, 0.05);
    
    %% AM spectra
    
    fprintf('calculating envelope spectra\n');
    for ichan=1:Nchan
        clear Efft
        [Efft, f] = periodogram(E(ichan,:),[],f_spectra,fs); Efft=sqrt(Efft);
        E_spectrum(ichan,:,i_wav) = Efft;
    end
    
    clear Efft f Nfft
    
    %% Nth octave-band spectra
    
    fprintf('calculating Nth octave-band spectra\n');
    
%     if ~exist('Bmod')
        for i=1:N_fsamples_oct
            [Bmod(i,:), Amod(i,:)] = butter(1,[cutoff_oct(i) cutoff_oct(i+1)]/(fs/2));
        end
%     end
    
    for ichan=1:Nchan
        DC = mean(E(ichan,:));

        for i=1:N_fsamples_oct
            clear F
            F = filter(Bmod(i,:), Amod(i,:), E(ichan,:));
            Eoct_spectrum(ichan, i, i_wav) = rms(F)*sqrt(2);
            m_spectrum(ichan, i, i_wav) = rms(F)*sqrt(2)/DC;
        end
    end
    
    clear E F ichan DC
     
    %% FM spectra
    
    fprintf('calculating FM spectra\n');
    
    t=(1:Nsamples)/fs;
    undersample = 1;
    for ichan=1:Nchan
        clear f_periodo FM_P
        FMwithoutnan=FM(ichan,1:undersample:end);FMwithoutnan(isnan(FMwithoutnan))=[];
        twithoutnan=t(1:undersample:end);twithoutnan(isnan(FM(ichan,1:undersample:end)))=[];
        if length(FMwithoutnan)<100
            FM_spectrum(ichan,:,i_wav) = nan(1,N_fsamples,1);
        else
            [FM_P,f_periodo] = plomb(FMwithoutnan,twithoutnan,100,16); 
            FM_P = sqrt(FM_P);
            %[FM_P,f_periodo] = fastlomb(FMwithoutnan,twithoutnan,0,1,20); FM_P=sqrt(FM_P/length(twithoutnan));%
            FM_spectrum(ichan,:,i_wav) = interpmean( f_periodo, FM_P, f_spectra_intervals );
            clear FMwithoutnan twithoutnan f_periodo FM_P
        end
    end
    
    clear FM
    
    %% f0 extraction
    
    fprintf('f0 extraction\n');
    
    % Parameters for the YIN algorithm
    ap0_thres = 0.6;
    P=[]; P.hop = 1;P.sr = fs/undersample;
    
    % Parameters for artifact removing in the f0 trajectory
    maxjump = 10;
    minduration = 0.02;
    minf = 50;
    maxf = 600;
    
    R = yin(son(1:undersample:end), P);
    f0 = 440*2.^(R.f0);
    f0_withnan = f0; 
    f0_withnan(R.ap0>ap0_thres)=NaN;
    
    f0_withnan = remove_artifact_FM( f0_withnan, fs, maxjump, minduration, minf, maxf );
    
    clear P R f0
    
    %% f0 modulation spectrum
    
    fprintf('calculating f0 modulation spectrum\n');
    
    f0withoutnan=f0_withnan;f0withoutnan(isnan(f0_withnan))=[];
    twithoutnan=t(1:undersample:end);twithoutnan(isnan(f0_withnan))=[];
    
    clear t f0_withnan
    
    %[f0_P,f_periodo] = fastlomb(f0withoutnan,twithoutnan,0,1,20); f0_P = sqrt(f0_P/length(twithoutnan));
    [f0_P,f_periodo] = plomb(f0withoutnan,twithoutnan,100,16); 
    f0_P = sqrt(f0_P);
    
    f0_spectrum(:,i_wav) = interpmean( f_periodo, f0_P, f_spectra_intervals );
    
    clear f0withoutnan twithoutnan f0_P f_periodo
end

%% plotting

figure

subplot(4,1,1)
semilogx(f_spectra, 10*log10(squeeze(nanmean(nanmean(E_spectrum,1),3))));hold on
ylabel('amplitude (dB)'); xlabel('rate (Hz)'); title('AM spectra')

subplot(4,1,2)
semilogx(f_oct, 10*log10(squeeze(nanmean(nanmean(m_spectrum,1),3))));hold on
ylabel('modulation index'); xlabel('rate (Hz)'); title('modulation index spectra')

subplot(4,1,3)
semilogx(f_spectra, 10*log10(squeeze(nanmean(nanmean(FM_spectrum,1),3))));hold on
ylabel('amplitude (dB)'); xlabel('rate (Hz)'); title('FM spectra')

subplot(4,1,4)
semilogx(f_spectra, 10*log10(squeeze(nanmean(f0_spectrum,2))));hold on
ylabel('amplitude (dB)'); xlabel('rate (Hz)'); title('f0 modulation spectra')
