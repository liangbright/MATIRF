function [F J Ib]=GaussMix_FitErr(Param, ParticleNum, PSFSigma, I_Data, X, Y, ...
                                  Ib_Prior, Ib_SigmaPrior, A_Prior, A_SigmaPrior)
                              
if nargin == 6
    [F J Ib]=GaussMix_FitErr_NoPrior(Param, ParticleNum, PSFSigma, I_Data, X, Y);
elseif nargin == 8
    [F J Ib]=GaussMix_FitErr_with_b_Prior(Param, ParticleNum, PSFSigma, I_Data, X, Y, ...
                                          Ib_Prior, Ib_SigmaPrior);
elseif nargin == 10
    [F J Ib]=GaussMix_FitErr_with_b_Prior(Param, ParticleNum, PSFSigma, I_Data, X, Y, ...
                                          Ib_Prior, Ib_SigmaPrior, A_Prior, A_SigmaPrior);
else
    disp('error in GaussMix_FitErr !')
end                              