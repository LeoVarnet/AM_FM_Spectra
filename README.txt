This MATLAB toolbox runs through audio files in a designated directory and calculate both the amplitude- and frequency-modulation spectra, as described in Varnet et al. 2017 [1]. It allows the characterization of sounds in terms of their temporal envelope and temporal fine structure content. 
Four different spectra are computed (see [1] for a complete description):
- AMa spectrum: classic Fourier AM spectrum (see e.g. [2])
- AMi spectrum: modulation index spectrum (see e.g. [3])
- FM spectrum: Fourier FM spectrum (see e.g. [2]) 
- f0M spectrum: Fourier spectrum of the modulations in f0 trajectory
The computation of f0 modulation spectra requires the YIN toolbox developped by de Cheveigné & Kawahara (http://audition.ens.fr/adc/) [4].

refs:
[1] A cross-linguistic study of speech modulation spectra. Varnet L, Ortiz-Barajas MC, Erra RG, Gervain J, Lorenzi C. J Acoust Soc Am. 2017 Oct;142(4):1976. doi: 10.1121/1.5006179.
[2] Effects of age and hearing loss on the relationship between discrimination of stochastic frequency modulation and speech perception. Sheft S, Shafiro V, Lorenzi C, McMullen R, Farrell C. Ear Hear. 2012 Nov-Dec;33(6):709-20. doi: 10.1097/AUD.0b013e31825aab15. 
[3] A review of the MTF concept in room acoustics and its use for estimating speech intelligibility in auditoria. T. Houtgast and H. J. M. Steeneken. J. Acoust. Soc. Am. 1985; 77:3, 1069-1077
[4] YIN, a fundamental frequency estimator for speech and music. de Cheveigné A, Kawahara H. J Acoust Soc Am. 2002 Apr;111(4):1917-30.
