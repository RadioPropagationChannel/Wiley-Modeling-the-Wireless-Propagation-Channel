function  mat = orbits (start_time, lat, long, ht, delta, tsimul, M_kepler)
%
% For each satellite of aconstellation, this function returns the elevation, azimut
% and range, from the earth station's point of view., plus the velocity coordinates.
%
% Inputs:
% 	yr, mo, day, hour, min, sec: Start time (format: 1997, 30, 1, 14, 26, 34)
% 	lat	: Latitude of the earth station (North lats.: positive numbers (Degrees))
% 	long	: Long. of the earth station (East long.: positive numbers (Degrees))
% 	ht	: height above sea level (meters)
% 	delta	: Time step (Minutes)
% 	tsimul	: Time of simulation (days)
% 	M_kepler: Matrix with the keplerian elements of the satellites.
%
% Outputs:
%	mat	: sampled constellation matrix
%
% By: Manuel Bentosinos-Rico - University of Vigo, Spain
% Release: 6/98

% Constants definition
k1 = 54548481.1;
k2 = 65915.3446;
mu = 398600.5; % [km^3/s^2]
c = pi / 180;  % Conversion from degrees to radians
numb_iterations = 25;   % Number of iterations in the iterative methods

% Number of step for each satellite´s simulation
n_steps = floor ((tsimul* 24 *60) / delta);

if (long < 0)
	long = 360 + long;
end;
lat = lat * c;
long = long * c;
ht = ht / 1000;
delta = delta / 60 / 24; % Calculation of time step in universal time

yr = start_time (1);
mo = start_time (2);
day = start_time (3);
hour = start_time (4);
min = start_time (5);
sec = start_time (6);

mat = [];


% Convertion of angles to radians
M_kepler (:, 3:4) = M_kepler (:,3:4) * c;  % Inc. and raan
M_kepler (:, 6:7) = M_kepler (:,6:7) * c;  % Arg. of Prg. and Mean anom.

% Convertion of rev./day to rad./day
M_kepler (:, 8) = M_kepler (:,8) *  2*pi;

n_sats = size (M_kepler, 1);

for index = 1 : n_sats
 %disp('Satellite Number:');
 %disp(index);

 epoch_yr = M_kepler (index, 1);
 epoch_day = M_kepler (index, 2);
 inc = M_kepler (index, 3);
 raan = M_kepler (index, 4);
 e = M_kepler (index, 5);
 arg_prg = M_kepler (index, 6);
 mean_anomaly = M_kepler (index, 7);
 mean_motion = M_kepler (index, 8);

 % Calculation of a using Newton's method
 k3 = (1 - 1.5 * sin(inc)^2) / (1-e^2)^1.5;
 a = (k1 / mean_motion) ^(2/3);
 i = 1;
 while (i < numb_iterations)
	n = k1 * (1 + k2*k3 / a^2) / a^1.5;
	dn = -(1.5*n/a + k1*k2*k3/a^4.5);
	a = a + (mean_motion - n) / dn;
	if abs(mean_motion - n) > 1e-6
		i = i + 1;
	else
		i = numb_iterations;
	end;
 end

 % Calculation of regresion of nodes and rotation of apsides
 k = mean_motion * k2 / a^2 / (1 - e^2)^2;
 reg_nodes = -k * cos (inc);  % Regresion of nodes
 rot_apsides = k * (2 - 2.5 * sin (inc)^2);

 % Calc of universal time corresponding to start-time as a fraction of a day
 UT = (hour + min/60 + sec/3600) / 24;

 % Calc of Julian date corresponding to start-time
 [jd, d] = jul_date (yr, mo, day, hour, min, sec);


% Calc of Julian date corresponding to epoch
 y = 0+ epoch_yr;%ojo! cambio de siglo antes 1900
 aux = y - 1986;       
 jd_epoch = aux * 365;
 
 [sal]=ch_leap(yr-1);
 if (sal == 1)
    jd_epoch = jd_epoch + floor(aux/4) + 1;   % A day per leap-year is added
 else
    jd_epoch = jd_epoch + floor(aux/4); 
 end

jd_epoch = jd_epoch + epoch_day;
 jd_epoch = jd_epoch + 2446430.5;

% Calculation of time interval (days)
 dt = jd - jd_epoch;

 % Calculation of angles at desired time
 raan = raan + reg_nodes * dt;
 arg_prg = arg_prg + rot_apsides * dt;
 mean_anomaly = mean_anomaly + mean_motion * dt;

 % Starts the simulation for the current satellite
 %step init
 m1 = 0;

 triplets = [];

 for step=1:n_steps
  %d=universal time day
  x = floor (d);
  rm = (d-x) * 24;
  h = floor (rm);
  rm = (rm -h) * 60;
  m = floor(rm);
  s = (rm -m) * 60;
  UT = (h + m/60 + s/3600) / 24;

  % Calculation of angles stepped by delta
  raan = raan + reg_nodes * m1;
  arg_prg = arg_prg + rot_apsides * m1;
  mean_anomaly = mean_anomaly + mean_motion * m1;

  % Calculation of nu and r at desired time using Newton's method

  % Angles mod 2* PI
  mean_anomaly = rem( mean_anomaly, 2*pi);

  ecc_an = pi;  % start value for eccentric anomaly
  i = 1;     % Number of iterations
  while (i <  numb_iterations) 
	dM = 1 - e * cos (ecc_an);
	M = ecc_an - e * sin (ecc_an);
	ecc_an = ecc_an + (mean_anomaly - M) / dM;
	if  abs (mean_anomaly - M) > 1e-8 
		i = i +1;
	else
		i = numb_iterations;
	end;
  end;
  nu = 2 * atan ( sqrt( (1 + e)/(1 - e) ) * tan(ecc_an /2) );
  if (nu < 0)
	nu = 2 * pi + nu;
  end
  % Angles mod 2*PI
  nu = rem (nu, 2*pi);
  r = a * (1 - e * cos (ecc_an));
 

  % PQW position vector
  RP = r * cos (nu);
  RQ = r * sin (nu);

  % PQW velocity vector
  p = a*(1-e^2);
  VP = -1*sqrt(mu/p) * sin( nu);
  VQ = sqrt (mu/p) * (e + cos (nu)); 


  % R-Trans matrix
  cn = cos (raan);
  co = cos (arg_prg);
  ci = cos (inc);
  sn = sin (raan);
  so = sin (arg_prg);
  si = sin (inc);
  R = [ cn*co-sn*so*ci  -cn*so-sn*co*ci ; 
        sn*co+cn*so*ci  -sn*so+cn*co*ci    ;
        so*si           co*si              ];

  % PQW to IJK
  Rijk = R * [RP  RQ]';
  ri = Rijk(1);
  rj = Rijk(2);
  rk = Rijk(3);

  Vijk = R * [VP VQ]';
  vi = Vijk(1);
  vj = Vijk(2);
  vk = Vijk(3);


  % Calculation of GST (in degrees)
  T = (jd - 2415020) / 36525;
%  T = (jd_epoch - 2415020) / 36525;
  gst = 99.691 + 36000.7689*T + .0004*T^2 + UT*360;
  gst = gst * c;
  % Angle mod 2*PI
  gst = rem(gst, 2*pi);
  lst = gst + long;
  lst = rem( lst, 2*pi);


  c_lat = cos (lat);
  s_lat = sin (lat);
  c_lst = cos (lst);
  s_lst = sin (lst);

  % Earth station vector
  r_earth = 6378.14;   % equatorial radius of the earth
  r_earth = r_earth * (1 - s_lat^2 / 298.257); % radius of the earth as a function of latitude
  x = (r_earth + ht) * c_lat;
  z = (r_earth + ht) * s_lat;

  rx = x * c_lst;
  ry = x * s_lst;

  % Range vector
  range_x = ri - rx;
  range_y = rj - ry;
  range_z = rk - z;

  % Calculation of D matrix
   
  D = [ s_lat*c_lst  s_lat*s_lst  -c_lat ;
        -s_lst       c_lst          0    ;
        c_lat*c_lst  c_lat*s_lst  s_lat  ];

  % SEZ components
  rsez = D * [range_x range_y range_z]';
  rs = rsez(1);
  re = rsez(2);
  rz = rsez(3);

  Vsez = D * [vi vj vk]';
  vs = Vsez(1);
  ve = Vsez(2);
  vz = Vsez(3);

  % Range
  range = sqrt (rs^2 + re^2 + rz^2);

  % Look angles
  elevation = asin (rz/range);
  azimut = -atan (re/rs);

  elevation = elevation / c;
  azimut = azimut / c;

  if ((rs<0) & (re>0))
	azimut = azimut;
  elseif ((rs>0) & (re>0))
        azimut = azimut + 180;
  elseif ((rs>0) & (re<0))
        azimut = azimut + 180;
  else
	azimut = azimut + 360;
  end
 
  triplets = [triplets; [elevation azimut range vs ve vz]];

  m1 = delta;
  d = d + m1;
  jd = jd + m1;

 end

 mat = [mat triplets];

end

%save outdata.dat mat -ascii