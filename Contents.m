% SRC
% 
%   Contents file for SRC and its subfolders.
%   
%   SRC/AUX
%   /src/Aux.CalcularEmbDim                                                 - Function to obtain the embedding dimension. First it tries to use the
%   /src/Aux.CalcularTao                                                    - Function to obtain the delay time (tau). First it tries to use the automutual
%   /src/Aux.cao_dim                                                        - Tstoolbox/@signal/cao
%   /src/Aux.DeFrame                                                        - Reconstructs the vSenal signal from the windows frames created by the
%   /src/Aux.LogEnergy                                                      - Calculates the energy of a frame in logarithmic and linear scale,
%   /src/Aux.LPC2                                                           - Computes the coeffiecients vA and the gain of the filter using a LPC
%   /src/Aux.QuitaExtremos                                                  - This function removes the 0's (or low amplitude values ) at the beginning or at the end
%   /src/Aux.SilentDetectorThreshold                                        - Calculates if the vFrame is silent or not according to the energy of
%   /src/Aux.ZeroCrossingRate                                               - Calculate the zero-crossing rate of the frame
%   
%   SRC/COMPLEXITY
%   /src/Complexity.CalculateRegularity                                     - Calculates a series of regularity features from the inputted attractor
%   /src/Complexity.embeb                                                   - Reconstruct state space vectors in a iDim dimensional space, using iTao delays
%   
%   SRC/PERTURBATIONFLUCTUATION
%   /src/PerturbationFluctuation.PerturbationFluctuation                    - Calculation of the perturbation and fluctuation parameters
%   
%   SRC/PERTURBATIONFLUCTUATION/PERTURBATION/CHNR
%   /src/PerturbationFluctuation/Perturbation/CHNR.CHNR                     - Calculates the harmonic-to-noise ratio (HNR) of a voice signal, using
%   /src/PerturbationFluctuation/Perturbation/CHNR.CHNRi                    - Calculates the harmonic-to-noise ratio (HNR) of a voice segment, using
%   
%   SRC/PERTURBATIONFLUCTUATION/PERTURBATION/GNE
%   /src/PerturbationFluctuation/Perturbation/GNE.CalculaEnvolventesHilbert - Calculates the Hilbert envelopes in different frequency bands of each of the input voice
%   /src/PerturbationFluctuation/Perturbation/GNE.CalculaParesHilbert       - Calculates the center frequencies for which the cross correlations of Hilbert envelopes are to be computed
%   /src/PerturbationFluctuation/Perturbation/GNE.EstructuraFina            - Calculate the fine structure or the glottic excitation signal (FineStructure)
%   /src/PerturbationFluctuation/Perturbation/GNE.GNE                       - Calculates the GNE -Glottal-to-Noise Excitation ratio- according to
%   
%   SRC/PERTURBATIONFLUCTUATION/PERTURBATION/HNR
%   /src/PerturbationFluctuation/Perturbation/HNR.HNRiYum                   - Calculates the value of the HNR ratio for a speech segment corresponding to
%   /src/PerturbationFluctuation/Perturbation/HNR.HNRiYumBoy                - Calculates the value of the HNR ratio for a speech segment corresponding to
%   /src/PerturbationFluctuation/Perturbation/HNR.HNRYum                    - Calculates the harmonic ratio to HNR noise in dB according to the Yumoto method
%   
%   SRC/PERTURBATIONFLUCTUATION/PERTURBATION/JITTER_SHIMMER
%   /src/PerturbationFluctuation/Perturbation/Jitter_Shimmer.apq            - Calculates the amplitude disturbance ratio in % (APQ)
%   /src/PerturbationFluctuation/Perturbation/Jitter_Shimmer.feijoo_pert    - Calculate the disturbance ratio proposed by Feijoo as a percentage, which is similar to Kasuya's PQ perturbation ratio
%   /src/PerturbationFluctuation/Perturbation/Jitter_Shimmer.jitter         - Calculates the jitter of the voice signal. Calculate its absolute value
%   /src/PerturbationFluctuation/Perturbation/Jitter_Shimmer.jitter_env     - Computes iNumPoints jitter values (absolute jitter in microseconds, and relative
%   /src/PerturbationFluctuation/Perturbation/Jitter_Shimmer.ParamPitch     - Calculates the pitch period and pitch amplitude sequences of the input voice
%   /src/PerturbationFluctuation/Perturbation/Jitter_Shimmer.ppq            - Calculates the disturbance quotient of the pitch period in % (RAP)
%   /src/PerturbationFluctuation/Perturbation/Jitter_Shimmer.pq             - Calculates the disturbance ratio as a percentage of the signal vSequence,
%   /src/PerturbationFluctuation/Perturbation/Jitter_Shimmer.rap            - Calculates the relative average disturbance of the pitch period in % (RAP)
%   /src/PerturbationFluctuation/Perturbation/Jitter_Shimmer.sapq           - Calculates the smoothed amplitude disturbance ratio in % (sAPQ)
%   /src/PerturbationFluctuation/Perturbation/Jitter_Shimmer.SecuenPitch    - Calculates the pitch or amplitude sequences given a vSignal
%   /src/PerturbationFluctuation/Perturbation/Jitter_Shimmer.SeparaCiclos   - Separates the cycles of a known voice signal given the moments in which
%   /src/PerturbationFluctuation/Perturbation/Jitter_Shimmer.shimmer        - Calculates the shimmer of the signal. Calculate its absolute value
%   /src/PerturbationFluctuation/Perturbation/Jitter_Shimmer.shimmer_env    - Calculates iNumPoints shimmer values (absolute value in dB, and relative value
%   /src/PerturbationFluctuation/Perturbation/Jitter_Shimmer.sppq           - Calculates the disturbance ratio of the smoothed amplitude in % (sPPQ)
%   
%   SRC/PERTURBATIONFLUCTUATION/PERTURBATION/NNE
%   /src/PerturbationFluctuation/Perturbation/nne.NNE                       - Calculates the normalized noise energy (NNE) of a speech segment (according to
%   /src/PerturbationFluctuation/Perturbation/nne.NNEi                      - Calculates the normalized noise energy (NNE) of a voice segment according to
%   /src/PerturbationFluctuation/Perturbation/nne.NNEs                      - Calculates the normalized noise energy (NNE) of a voice signal,
%   
%   SRC/PERTURBATIONFLUCTUATION/TREMOR
%   /src/PerturbationFluctuation/tremor.Tremor                              - Calculates the tremor (parameters that modulate the pitch period and
%   
%   SRC/PITCHDETERMINATION
%   /src/PitchDetermination.CalculatePitch                                  - Calculates the pitch of a given frame using different pitch
%   
%   SRC/PITCHDETERMINATION/BOYANOV
%   /src/PitchDetermination/boyanov.DynamicClip                             - Dynamic central clipping with threshold adaptation
%   /src/PitchDetermination/boyanov.Picos                                   - Computes the amplitudes P and the i_P positions of the positive or negative peaks
%   /src/PitchDetermination/boyanov.PicosMayores                            - Calculate the highest energy peaks, cycle by cycle of the input voice
%   /src/PitchDetermination/boyanov.PitchBoyanov                            - Calculates the pitch of a given frame using the method described in:
%   /src/PitchDetermination/boyanov.PitchDetecBoy                           - Calculates the pitch using the Boyanov algorithm
%   
%   SRC/PITCHDETERMINATION/KASUYA_FEIJOO
%   /src/PitchDetermination/kasuya_feijoo.PicosKas                          - Computes the amplitudes P and the i_P positions of the positive or negative peaks
%   /src/PitchDetermination/kasuya_feijoo.PitchKas                          - Calculate pitch using the Kasuya pitch determination algorithm
%   /src/PitchDetermination/kasuya_feijoo.PitchMedio                        - Calculates the average pitch in Hz of the input voice
%   /src/PitchDetermination/kasuya_feijoo.PitchSegKas                       - Calculates the pitch of a given frame using the method described by
%   
%   SRC/PITCHDETERMINATION/PICOS
%   /src/PitchDetermination/Picos.PicosNeg                                  - Obtains the amplitudes P and the positions i_P of the negative peaks,
%   /src/PitchDetermination/Picos.PicosNegKas                               - Computes the amplitudes P and the positions i_P of the negative peaks,
%   /src/PitchDetermination/Picos.PicosPos                                  - Gets the amplitudes P and the i_P positions of the positive peaks,
%   /src/PitchDetermination/Picos.PicosPosKas                               - Computes the amplitudes P and the positions i_P of the positive peaks,
%   
%   SRC/PITCHDETERMINATION/RABINER
%   /src/PitchDetermination/rabiner.PitchCeps                               - Pitch estimation using the LPC method of Rabiner
%   /src/PitchDetermination/rabiner.PitchClip                               - Pitch estimation using static clipping
%   /src/PitchDetermination/rabiner.PitchCorr                               - Pitch estimation using the LPC method of Rabiner
%   /src/PitchDetermination/rabiner.StaticClip                              - Static clipping
%   
%   SRC/SPECTRALCEPSTRAL
%   /src/SpectralCepstral.MFCCi                                             - Calculate the mel cepstrum of a signal frame
%    
%   This file was generated by updateContents.m on 11 May 2020 at 00:34:24.
