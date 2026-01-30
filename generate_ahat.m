function [interferencePow , a_hat,ck_wnoise,SINR] = generate_ahat(ck,interference,CutoffCoeff_value,SNRdB_value,Normalized_Coeff,P_value,VardsidashNormalized,M)
deltaTHP=0;

if CutoffCoeff_value>=1
    SNR_adjust = 0;
else
    SNR_adjust = 10*log10(1./CutoffCoeff_value);
end

interferencePow = var(interference);
SNR=db2pow(SNRdB_value +SNR_adjust);

% SNRdB_value = P_wTHP_value/Noise power; : Continuous Energy/time /sampling noise power
% SINR = symbol minimum distance / (sampling interference power + sampling noise power)

ck_wnoise_Norm = ck.*Normalized_Coeff./ sqrt(P_value) + wgn(1,numel(ck),1/SNR,'linear') ;
ck_wnoise = ck_wnoise_Norm .*sqrt(P_value)./Normalized_Coeff;

SINR = VardsidashNormalized.SignalD.^2/(VardsidashNormalized.ik + 1/SNR);

a_hat=wrapToM(ck_wnoise-deltaTHP,M)+deltaTHP;
end
