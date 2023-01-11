%Genera datos de la orbita de los satelites de una constelacion, sirve para
%trabajar con las constelaciones incluidas en el directorio "ArchivosTLE", 
%si se quierre otra constelacion distinta basta añadir su fichero TLE a dicho
%directorio.
%SALIDA:se guarda la variable "matriz" en "salida.mat"
%matriz:estructura 
%      -----------satelite1-----------     -----------satelitek-----------     
%      elevacion acimut rango vx vy vz ... elevacion acimut rango vx vy vz 
% t1==>
% t2==> 
% ...
% tn==> 

M_kepler=[];

read_tle;

anho=2008;    % year
mes=1;        % month 
dia=1;        % day #
hora=0;       % hour
min=0;        % minutes
seg=0;        % seconds

lat=42.13;    % latitude (degree)
long=8.44;    % longitude (degree)
alt=0;        % hieght(m)


paso=1/60;  %input('Interval between samples (min): ');
dur=0.02;   %input('Simulation period (days): ');
disp('   ')
disp('Wait ...')

constmatrix=orbits([anho mes dia hora min seg],lat,long,alt,paso,dur,M_kepler);
%whos
%disp('   ')
%disp('Constalletion data is in constallation.mat and constallationaux.mat')

save constellationaux.mat  name anho mes dia hora min seg lat long alt paso dur  
save constellation.mat constmatrix
clear 
disp('   ')
disp('Constalletion data newSatMatrix.mat')
picksats
clear 