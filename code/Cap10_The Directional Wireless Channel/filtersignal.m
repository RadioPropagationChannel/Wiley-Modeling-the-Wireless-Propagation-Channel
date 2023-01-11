function [rout,B,A]=filtersignal(rin,Wp,Ws,Rp,Rs)
[N,Wn]=buttord(Wp,Ws,Rp,Rs);
[B,A]=butter(N,Wn);
rout=filter(B,A,rin);