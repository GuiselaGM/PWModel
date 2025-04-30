% Code with an example to get the significant wave height, dir and period in 
% a buoy location crossed by a Tropical Cyclone using the PWModel 
%_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_
% Reference:
% Grossmann-Matheson et al., 2025, The spatial distribution of ocean wave 
% parameters in tropical cyclones, Ocean Eng, 317, 120091 (2025). 
% https://doi.org/10.1016/j.oceaneng.2024.120091
%_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_
%
% Parameters required to run this code: 
% Vmax = 'Maximum sustained wind speed' 'm/s'
% Vfm = 'Storm translation speed' 'm/s'
% R34 = 'Radius of 34 knot winds (mean of all quadrants)' 'meters'
% Rmax = 'Radius of maximum winds' 'meters'
% Lat = Latitude of Storm center 'degrees'
% Lon = Longitude of Storm center 'degrees' 
% lat_point = latitude of desired point/buoy (degrees)[-90 +90]
% lon_point = longitude of desired point/buoy (degrees)[-180 +180]
% rotationAngle = angle to rotate the storm, negative value is anticlockwise
% from North [-180 +180]. Default is zero (storm moving Northward). It doesn't matter the hemisphere,
% get this value from the storm track positions.
%
% OUTPUT: Hs, dir (including components) and Tp for the required point (ex: buoy position)
% Hs_point = Significant wave Height
% Tp_point = Peak wave Period
% U_point = zonal wave dir component
% V_point = meriodional wave dir component
% Dir_point = Peak wave direction coming from (meteorological convention)
% Dir_point_oc = Peak wave direction going to (oceanographic convention)
%
clear all
close all
%
%addpath('.\wave_diagrams')

% ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
% ENTER PARAMETERS:
% Input TC parameters (Vectors can have more than one TC track point)
Vmax = [50];                                                               % m/s                                                         %m/s
Vfm = [5];                                                                 % m/s
Rmax = [30]*10^3;                                                          % m
R34 = [300]*10^3;                                                          % m
Lat = [29.4]; %Lat = [-29.4];                                              % deg
Lon = [-77.3];                                                             % deg
% Input location of desired point to find Hs and Tp
lon_point = -78.5;                                                         % Buoy longitude
lat_point  = 28.9; %lat_point  = -28.9;                                    % Buoy latitude
% ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>

% Model limits warning
if Vmax > 78 | Vfm > 15 | Rmax > 60*10^3 | R34 > 400*10^3 ...
        | Vmax < 17 | Rmax < 15*10^3 | R34 < 200*10^3 == 1
%disp('Parameter(s) over the limits. Default used instead.');
warning('Parameter(s) over the limits. Default(s) will be used instead.')
end

%Loop to resolve all storm track points available   
for k = 1:length(Vmax)

% Call Fetch model to calculate Hs
[HS,TP,U,V,XX,YY] = PWModel(Vmax(k),Vfm(k),Rmax(k),R34(k),Lat(k));

%Input angle to Rotate the Storm:
prompt = {'Enter rotation relative to North (negative is anticlockwise) [-180^{\circ} +180^{\circ}]'};
dlgtitle = 'Rotation (degrees)';
definput = {'0'};
dims = [1 50];
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
rotationAngle = str2num(answer{1})

% ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
% Call function to rotate and get the Hs, Tp and Dir for buoy location 
% ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>

[Hs_rotated,Tp_rotated,U_rotated,V_rotated] = rotate_get_HsTpdir(HS,TP,U,V,XX,YY,Lat(k),Lon(k),lat_point,lon_point,rotationAngle);

Hs_rotated_all(k) = Hs_rotated;
Tp_rotated_all(k) = Tp_rotated;
U_rotated_all(k) = U_rotated;
V_rotated_all(k) = V_rotated;

end

% Significant Wave Height at the required point (ex: buoy)
Hs_point = Hs_rotated_all';
Tp_point = Tp_rotated_all';
U_point = U_rotated_all';
V_point = V_rotated_all';

%Peak wave dir coming from (meteo convention, starting from N = zero degree)
Dir_point = ((180/pi)*atan2(-U_point,-V_point)); 
Dir_point(Dir_point < 0) =  Dir_point(Dir_point < 0) + 360;                
Dir_point = round(Dir_point);                                              % waves coming from

%Peak wave dir going to (oceanographic convention)
Dir_point_oc = ((180/pi)*atan2(U_point,V_point)); 
Dir_point_oc(Dir_point_oc < 0) =  Dir_point_oc(Dir_point_oc < 0) + 360;            
Dir_point_oc = round(Dir_point_oc);                                        % waves going to

%=======================================================
%Plot Hs and Tp calculated for all given TC track points 
%uncomment if using more than one track point
%=======================================================
% %Hs
% figure;
% plot(Hs_point,'*-k','LineWidth',1);
% grid on
% legend('H_{s} at buoy location');
% ylabel('H_{s} (m)');xlabel('TC track point');
% xticks([1:1:length(Vmax)]);
% title('TC Significant Wave Height');
% xlim([0 length(Hs_point)+1 ])
% ylim([fix(Hs_point(1))-1 round(Hs_point(end))+1])
% %
% %Tp
% figure;
% plot(Tp_point,'*-k','LineWidth',1);
% grid on
% legend('T_{p} at buoy location');
% ylabel('T_{p} (s)');xlabel('TC track point');
% xticks([1:1:length(Vmax)]);
% title('TC Significant Wave Height');
% xlim([0 length(Tp_point)+1 ])
% ylim([fix(Tp_point(1))-1 round(Tp_point(end))+1])

% Save output
fname = ['HsTpDir_TC_' int2str(Vmax(1)) '_' int2str(Vfm(1)) '_' int2str(Rmax(1)/10^3) '_' int2str(R34(1)/10^3)]
save(fname,'Hs_point','Tp_point','U_point','V_point','Dir_point','Dir_point_oc')

% END :)
% Guisela, October/2025
% ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>




