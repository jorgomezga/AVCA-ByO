# AVCA-ByO toolbox

Matlab *Automatic Voice Condition Analysis* toolbox, used in the ByO lab (www.byo.ics.upm.es/BYO). 
It contains the code used in the a series dedicated to the Automatic Voice Condition Analysis. If you use this sofware please cite the following sources:

3. Gómez-García, J.A. Moro-Velázquez, L. Arias-Londoño, J.D. Godino-Llorente. J.I. "On the design of automatic voice condition analysis systems. Part III: review of acoustic modelling strategies." Biomedical Signal Processing and Control. 

2. Gómez-García, J. A., Moro-Velázquez, L., & Godino-Llorente, J. I. (2019). "On the design of automatic voice condition analysis systems. Part II: Review of speaker recognition techniques and study on the effects of different variability factors." Biomedical Signal Processing and Control, 48, 128-143. [full paper](https://doi.org/10.1016/j.bspc.2018.09.003) - [preprint](https://zenodo.org/record/2624815)

1. Gómez-García, J. A. Moro-Velázquez, L., & Godino-Llorente, J. I. (2019). "On the design of automatic voice condition analysis systems. Part I: Review of concepts and an insight to the state of the art." Biomedical Signal Processing and Control, 51, 181-199. [full paper](https://doi.org/10.1016/j.bspc.2018.12.024) - [preprint](https://zenodo.org/record/2624638)

# Contents

# Prerrequisites

Some external toolboxes were used for the computation of certain features in the paper, including:
- The VoiceBox toolbox for the MFCC computation and use of some auxiliary functions (http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html) [M Brookes. Voicebox: Speech processing toolbox for matlab. software,[may 2020], 2020.]

- The Covarep toolbox for the computation of CPP (https://covarep.github.io/covarep/) [Gilles Degottex, John Kane, Thomas Drugman, Tuomo Raitio, and Stefan Scherer. Covarep A collaborative voice analysis repository for speech technologies. In 2014 ieee international conference on acoustics, speech and signal processing (icassp), pages 960–964. IEEE, 2014.] 

- The implementation in (http://www.maxlittle.net/software/fastdfa.zip) to compute DFA [Max A Little, P. E. McSharry, Irene M Moroz, and Stephen J. Roberts. Nonlinear, biophysically890 informed speech pathology detection. In IEEE International Conference on Acoustics, Speech and Signal Processing, 2006. ICASSP 2006 Proceedings., volume 2, pages II–II. IEEE, 2006.]

- The code in (http://www.maxlittle.net/software/rpde.zip) to compute RPDE [Max A Little, Patrick E. McSharry, Stephen J. Roberts, Declan Ae Costello, and Irene M Moroz. Ex893 ploiting Nonlinear Recurrence and Fractal Scaling Properties for Voice Disorder Detection. BioMedical Engineering OnLine, 6(1):23, 2007.]

- The implementation in (https://www.mathworks.com/matlabcentral/fileexchange/19148-hurst-parameter-estimate) to compute the Hurst exponent

- the functions in (https://github.com/jdariasl/ME) for the computation of the Markovian entropies [J.D. Arias-Londoño and J.I Godino-Llorente. Entropies from Markov Models as Complexity Measures of Embedded Attractors. Entropy, 17(6):3595–3620, 2015.]

- The HCTSA toolbox for the computation of D2 and LLE (https://github.com/benfulcher/hctsa) [Ben D Fulcher, Max A Little, and Nick S Jones. Highly comparative time-series analysis: the empirical structure of time series and their methods. Journal of the Royal Society Interface, 10(83):20130048, 2013]
- the Modulation Toolbox library ver 2.1 to calculate MS (https://sites.google.com/a/uw.edu/isdl/projects/modulation-toolbox) [Les Atlas, Pascal Clark and Steven Schimmel, Modulation Toolbox Version 2.1 for MATLAB, https://sites.google.com/a/uw.edu/isdl/projects/modulation-toolbox, University of Washington, September 2010.]

- The implementation in (https://www.ee.columbia.edu/~dpwe/resources/matlab/rastamat/) to calculate PLP [Daniel P. W. Ellis. PLP and RASTA (and MFCC, and inversion) in Matlab, 2005. online web resource.]

- The Modulation spectrum features in [Moro-Velázquez, L., Gómez-García, J. A., Godino-Llorente, J. I., & Andrade-Miranda, G. (2015). Modulation spectra morphological parameters: a new method to assess voice pathologies according to the grbas scale. BioMed research international, 2015.]

<!--
[a] M Brookes. Voicebox: Speech processing toolbox for matlab. software,[may 2020], 2020.
[b] Gilles Degottex, John Kane, Thomas Drugman, Tuomo Raitio, and Stefan Scherer. Covarep A collaborative voice analysis repository for speech technologies. In 2014 ieee international conference on acoustics, speech and signal processing (icassp), pages 960–964. IEEE, 2014.
[c] Max A Little, P. E. McSharry, Irene M Moroz, and Stephen J. Roberts. Nonlinear, biophysically890 informed speech pathology detection. In IEEE International Conference on Acoustics, Speech and Signal Processing, 2006. ICASSP 2006 Proceedings., volume 2, pages II–II. IEEE, 2006.
[d] Max A Little, Patrick E. McSharry, Stephen J. Roberts, Declan Ae Costello, and Irene M Moroz. Ex893 ploiting Nonlinear Recurrence and Fractal Scaling Properties for Voice Disorder Detection. BioMedical Engineering OnLine, 6(1):23, 2007.
[e] J.D. Arias-Londoño and J.I Godino-Llorente. Entropies from Markov Models as Complexity Measures of Embedded Attractors. Entropy, 17(6):3595–3620, 2015.
[f] Ben D Fulcher, Max A Little, and Nick S Jones. Highly comparative time-series analysis: the empirical structure of time series and their methods. Journal of the Royal Society Interface, 10(83):20130048, 2013
[g] Les Atlas, Pascal Clark and Steven Schimmel, Modulation Toolbox Version 2.1 for MATLAB, https://sites.google.com/a/uw.edu/isdl/projects/modulation-toolbox, University of Washington, September 2010.
[h] Daniel P. W. Ellis. PLP and RASTA (and MFCC, and inversion) in Matlab, 2005. online web resource.
[i] Moro-Velázquez, L., Gómez-García, J. A., Godino-Llorente, J. I., & Andrade-Miranda, G. (2015). Modulation spectra morphological parameters: a new method to assess voice pathologies according to the grbas scale. BioMed research international, 2015. 
-->
