
function [ak,bk,ck,dk,bkfromD,bk_wowrap,ck_wowrap,ck_woTHP]=IIRoutput(K,Input,FilterCoeff,ChannelCoeff,M)
%  M:bmax
%  ak:Input of System (Input of THP)
%  bk:Output of IIR (Input of Channel)
%  ck:Output of Channel (Output of System)
%  dk:Output of THP (Input of IIR)
%  Input(ak) -> THP(dk) -> IIR(bk) -> Channel (ck) 

ak=zeros(1,K);
bk=zeros(1,K);
ck=zeros(1,K);
dk=zeros(1,K);
bkfromD=zeros(1,K);
ck_woTHP=zeros(1,K); %蛻晄悄蛹?

ak=Input;
% ak = randi([0 1],1,K);
Kindex = 1;

bk(Kindex) = ak(1)/(1+FilterCoeff(1)); % bk:output
bkfromD(Kindex) =bk(Kindex);
dk(Kindex) =ak(Kindex);

bk_wowrap(Kindex) = bk(Kindex);
bk(Kindex) = wrapToM(bk(Kindex),M);


ck(Kindex) = filtering(ChannelCoeff,bk(1:Kindex));
ck_wowrap(Kindex) = filtering(ChannelCoeff,bk_wowrap(1:Kindex));
% ChannelCoeff(1)*bk(1)

ck_woTHP(Kindex) = filtering(ChannelCoeff,ak(1:Kindex));
Q=0;
for Kindex = 2:1:K

    bk(Kindex) = ak(Kindex)/(1+FilterCoeff(1)) - 1/(1+FilterCoeff(1)) * filtering(FilterCoeff(2:end),bk(1:Kindex-1));

    bk_wowrap(Kindex) = ak(Kindex)/(1+FilterCoeff(1)) - 1/(1+FilterCoeff(1)) * filtering(FilterCoeff(2:end),bk_wowrap(1:Kindex-1));
    %     rk(Kindex) = filtering(F

    % Q(Kindex)=(bk(Kindex)+M - mod(bk(Kindex)+M,2*M)) / (2*M);
    Q(Kindex)=(bk(Kindex)-wrapToM(bk(Kindex),M)) / (2*M); % need to check why its correct
    %     [Qsym,~]=quorem(sym(bk(Kindex)+M),sym(2*M));
    %     Q(Kindex) = double(Qsym);
    dk(Kindex) = ak(Kindex) -2*Q(Kindex)*M;
    bkfromD(Kindex) =(ak(Kindex) -2*Q(Kindex)*M )/(1+FilterCoeff(1)) - 1/(1+FilterCoeff(1)) * filtering(FilterCoeff(2:end),bk(1:Kindex-1));

    bk(Kindex) = wrapToM(bk(Kindex),M);
    ck(Kindex) = filtering(ChannelCoeff,bk(1:Kindex));
    ck_wowrap(Kindex) = filtering(ChannelCoeff,bk_wowrap(1:Kindex));
    ck_woTHP(Kindex)= filtering(ChannelCoeff,ak(1:Kindex));
    %     ck(Kindex) = wrapToM(ck(Kindex),2);

end
end