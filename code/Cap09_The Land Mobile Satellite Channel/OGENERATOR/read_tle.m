% 
% script read_tle.m
%
% This script reads the NASA 2-line files
% with the parameters of the constellation.
% The name of the file to be read is given
% by the parameter name.
%
% After this function is executed, the Kepler parameteres
% are stored in the matrix M_kepler, with the following format:
%
%  M_kepler = [ Epoch_yr Epoch_day Inc RAAN e ARG_PRG Mean_An Mean_Mot ]
%             [    :         :      :    :  :    :       :       :     ]
%             [    :         :      :    :  :    :       :       :     ]
%             [    :         :      :    :  :    :       :       :     ]
%
%  Each line corresponding to each satellite of the constellation.
%  M_kepler: matrix, n_sats x 8, where n_sats is the number of satellites.
%
%  epoch_yr:     The epoch year (four digits)
%  epoch_day:    Julian day and fractional portion of the day
%  inc:          Inclination [deg]
%  raan:         Right Ascension of the ascending node [deg]
%  e:            Eccentricity (with decimal point)
%  mean_anomaly: [deg]
%  mean_motion:  [revs./day]
%  arg_prg:      Argument of perigee [deg]


% Inicialization
M_kepler = [];
ruta=pwd;
disp('   ')
fichero=input('Intro TLE file (mane.extension): ','s');
name=[ruta '\tleFILES\' fichero]
fid = eval(['fopen(name)'])
if (fid>0)
   end_file=0;
   salir=0;

   while (end_file==0)      
      line0=fgetl(fid);
      if isempty(line0) 
         salir=1;
      end;   
      if isempty(line0)
         salir=1;
      end;   %si line0=vacio devuelve 1 sino devuelve 0
      if (salir==1)
         end_file=1;
      else
         
         line1 = fgets (fid)
         line2 = fgets (fid)
    
         new_sat = [str2num( line1(19:20)) + 2000];          % Epoch Year
      	 new_sat = [new_sat str2num( line1(21:32))];         % Epoch Day
	     new_sat = [new_sat str2num( line2(9:16))];          % Inclination
  	     new_sat = [new_sat str2num( line2(18:25))];         % RAAN
  	     new_sat = [new_sat str2num( ['0.' line2(27:33)])];  % Eccentricity
  	     new_sat = [new_sat str2num( line2(35:42))];	     % Arg. of Perigee
  	     new_sat = [new_sat str2num( line2(44:51))];	     % Mean Anomaly
  	     new_sat = [new_sat str2num( line2(53:63))];	     % Mean Motion
	     M_kepler = [M_kepler; new_sat];
      end;
   end; %while  
   fclose (fid); 
   clear fid line0 line1 line2 new_sat end_file
else  
   error(['File ' name ' cannot be opened.']);
end;
