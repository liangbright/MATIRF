function [BgDensity, SNR_output, SNRError]=SetParticleSNR(SNR, C, z, R, Range_BgDensity, MATIRF_Param)

Options=optimset('TolX',1e-3, 'Display', 'off');

BgDensity=fminbnd(@SNRDiff, Range_BgDensity(1), Range_BgDensity(2), Options, C, z, R, MATIRF_Param, SNR);

SNR_output=Cal_ParticleSNR(C, z, R, BgDensity, MATIRF_Param);

SNRError=abs(SNR-SNR_output);

end

function Value=SNRDiff(BgDensity, C, z, R, MATIRF_Param, SNR)

tempSNR=Cal_ParticleSNR(C, z, R, BgDensity, MATIRF_Param);

Value=abs(SNR-tempSNR);
end
