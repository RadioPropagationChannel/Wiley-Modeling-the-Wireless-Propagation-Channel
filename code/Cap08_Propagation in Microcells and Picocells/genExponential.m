function T=genExponential(lambda,Narrivals)

sigma=sqrt(1/2*lambda);

I=randn(Narrivals,1);
Q=randn(Narrivals,1);

r=abs((I+j*Q)*sigma);
T=r.^2;