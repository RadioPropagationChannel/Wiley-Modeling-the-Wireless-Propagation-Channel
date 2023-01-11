function RR=BINSequalprobRayleigh(Nbins, sigma)

% Nbins=5;
% sigma=1;

ProbBIN=1/Nbins;

RR=[0];
R1=0;
for ii=1:Nbins-1
    R2=sqrt(-2*sigma^2*log(exp(-R1.^2/(2*sigma^2))-ProbBIN));
    R1=R2;
    RR=[RR; R2];
end
RR=[RR; 99];

