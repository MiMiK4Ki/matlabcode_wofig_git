function [xnorm,ynorm,evm]  = CalcEVM(rx,tx)

hEVM = comm.EVM;

StartUpDelay = 1;
SamplesPerSymbol = 1;

x = rx(StartUpDelay+1:SamplesPerSymbol:end);
y = tx(StartUpDelay+1:SamplesPerSymbol:end);

c = xcorr(x,y);
[XcorrPeak, XcorrMaxInd] = max(c);

DelayInd = XcorrMaxInd - (mean(numel(x), numel(y)))+1;
if DelayInd < 1
    DelayInd = 1;
end
xDelayed = x(DelayInd:1:end);
yCancelled = y(1:1:(end-DelayInd+1));

SelForNInd = find(abs(yCancelled) <= 100000);
yForRot = yCancelled(SelForNInd);
xForRot = xDelayed(SelForNInd);

angleMean = angle(XcorrPeak);

xnorm = xDelayed./(sum(abs(xForRot))).*exp(-1i*angleMean);
ynorm = yCancelled./(sum(abs(yForRot)));

evm = hEVM(ynorm,xnorm);
end
