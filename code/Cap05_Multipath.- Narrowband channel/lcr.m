Max1=max(stren);
Min1=min(stren);
Max1=ceil(Max1);
Min1=floor(Min1);
%
Lcr=zeros((Max1-Min1).*10+1,1);
Lcrnorm=Lcr;
%
RecTime=length(stren).*1501.5013;
%
for Level=Min1:0.1:Max1,
	for j1=1:length(stren)-1,
		if (stren(j1+1)>stren(j1)) & (stren(j1+1)>=Level) & (stren(j1)<Level),
			Lcr((Level-Min1).*10+1)=Lcr((Level-Min1).*10+1)+1;
		end
	end
end
%
Lcr=Lcr./RecTime;
Xaxis=[Min1:0.1:Max1];
plot(Xaxis,Lcr)
grid
title('Absolute LCR (times/sec)')
ylabel('Crossings per second')
xlabel('Signal level relative to LOS')