function r=ricetheoretical(r,a,sigma)
A=(r./sigma^2);
B=-(r.^2+a^2);
C=(2*sigma^2);
D=besseli(0,(r.*a)/(sigma^2));
r=A.*exp(B./C).*D;