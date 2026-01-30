function [y3,P,y1] = calculate_Power(input,TXD_N,h_butterworth,hFLpS,NoSpS,deltat,h,Tc,fc)
impulse = zeros(1,(TXD_N-1)*NoSpS+1);
i = 1;
for data = input
    impulse(i) = data;
    i = i + NoSpS;
end

y1=conv(impulse,h); %
y2=y1; %upconvert
y3=conv(y2,h_butterworth).*deltat; % butterworth

[~,maxindex] = max( conv(h,h_butterworth) );

t=0:deltat:((2*hFLpS+TXD_N-1)*NoSpS+numel(h_butterworth)-1)*deltat;


% calculte power

sampling_index_P=(maxindex:1:maxindex+(TXD_N-1)*NoSpS+1 ) ;

P = sum(y3(sampling_index_P).^2).*deltat ./((TXD_N-1)*Tc) ;

% figure
% hold on
% plot(t(sampling_index_P),y3(sampling_index_P))
% 
% plot(t,y3,'--')

end
