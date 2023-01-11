function DeltaTheta=genLaplacian(S)

stepp=0.05;
MMax=80;
CDF=[];
axisCDF=[];

% Create Laplace CDF

for ii= -MMax : stepp :-stepp
    axisCDF=[axisCDF; ii];
    CDF=[CDF; 0.5*exp(-sqrt(2)*abs(ii)/S)];
end

axisCDF=[axisCDF; 0];
CDF=[CDF;0.5];

for ii=stepp:stepp:MMax
    axisCDF=[axisCDF; ii];
    CDF=[CDF;1-0.5*exp(-sqrt(2)*ii/S)];
end

draw=rand(1,1);
if draw<=0.5
    aa=find(draw>CDF);
    DeltaTheta=axisCDF(aa(end));
end
if draw>0.5
    aa=find(draw<CDF);
    DeltaTheta=axisCDF(aa(1));
end

