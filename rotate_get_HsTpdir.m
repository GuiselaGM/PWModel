function[Hs_rotated,Tp_rotated,U_rotated,V_rotated] = rotate_get_HsTpdir(HS,TP,U,V,XX,YY,Lat,Lon,lat_point,lon_point,rotationAngle)
%_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_
%Function to rotate the Storm and find wave parameters in a desired position within the
%Tropical Cyclone
%_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_

%==========================================================================
%Input parameters:
%HS = Hs spatial distribution 'meters'[n,m]
%TP = Tp spatial distribution 'seconds'[n,m]
%U and V components spatial distribution [n,m]
%XX = grid position axis x 'km' [n]
%YY = grid position axis y 'km' [m]
%Lat = Latitude of Storm center (degrees)[-90 +90]
%Lon = Longitude of Storm center (degrees)[-180 +180]
%lat_point = latitude of desired point (degrees)[-90 +90]
%lon_point = longitude of desired point (degrees)[-180 +180]
%rotationAngle = angle to rotate the storm, negative value is anticlockwise
%from North [-180 +180], it doesn't matter the hemisphere. (get this value
%prior from two consecutive track positions)

% Output parameters:
% Hs_rotated = Hs for required position in the TC wavefield (ex: buoy
% position)
% Tp_rotated = Tp for required position in the TC wavefield
% U_rotated = U component for required position in the TC wavefield
% V_rotated = V component for required position in the TC wavefield

%==========================================================================
% Reference:
% Grossmann-Matheson et al, 2025
% The spatial distribution of ocean wave parameters in tropical cyclones 
% Ocean Engineering (submitted)
% DOI: (to update)
% ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>

% Distance between desired point and Storm eye
dt_lat = (lat_point - Lat)*111.137;                                        % km
dt_lon = (lon_point - Lon)*111.137;                                        % km
dist_eye = (dt_lat.^2 + dt_lon.^2).^0.5;                                   % km

% y_position=dt_lat;
% x_position=dt_lon.*cosd(lat_point);

% ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
% Rotate the Storm
% ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>

%Rotation
rotation = rotationAngle;                                                  % deg
%hs_rotated = imrotate(Hs,rotation,'nearest','crop');                      % nearest-neighbor interpolation
hs_rotated = imrotate(HS,rotation,'bilinear','crop');                      % bilinear interpolation
tp_rotated = imrotate(TP,rotation,'bilinear','crop');   

%It needs to convert u,v to degrees, imrotate the image and then convert back to u,v
dir_temp = ((180/pi)*atan2(-U,-V));
dir_rotated = imrotate(dir_temp,rotation,'bilinear','crop');
dir_rotated(dir_rotated == 0) = NaN; % to clean the borders of the image 
%converting degrees to u,v (met convention)
u_rotated = -(sin((pi/180)*(dir_rotated+rotation)));
v_rotated = -(cos((pi/180)*(dir_rotated+rotation)));

% %Plot to confirm calculation:
% uu = -cosd(dir_rotated)
% vv = sind(dir_rotated)
% figure, quiversc(XX,YY,uu,vv,'density',20)

% ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
% Evaluate Hs and Tp values on the required position and plot
% ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>

%Define the lat, long grid to plot 
DLONG=(XX/111.137)/cosd(Lat);
DLAT=YY/111.137;
ALONG=DLONG+Lon;
ALAT=DLAT+Lat;

%Hs_rotated = interp2(XX,YY,hs_rotated,x_position,y_position,'nearest');   % cartesian coordinates
Hs_rotated = interp2(ALONG,ALAT,hs_rotated,lon_point,lat_point,'nearest'); % geographical coordinates
%
Tp_rotated = interp2(ALONG,ALAT,tp_rotated,lon_point,lat_point,'nearest'); % geographical coordinates
%
U_rotated = interp2(ALONG,ALAT,u_rotated,lon_point,lat_point,'nearest');   % geographical coordinates
V_rotated = interp2(ALONG,ALAT,v_rotated,lon_point,lat_point,'nearest');   % geographical coordinates
%
Dir_rotated = ((180/pi)*atan2(-U_rotated,-V_rotated));                     % dir coming from
Dir_rotated(Dir_rotated < 0) =  Dir_rotated(Dir_rotated < 0) + 360;        

Hs_rotated(Hs_rotated == 0) = NaN; %taking out zero values; 
Tp_rotated(Tp_rotated == 0) = NaN; %taking out zero values;

% ====================================
% Plot rotated storm and buoy location
% ====================================
%
%--------------------
% Plot Wave field (Hs)
%--------------------
figure;
plot(Lon,Lat,'k+');
axis('equal')
long_min=min(min(ALONG));
long_max=max(max(ALONG));
lat_min=min(min(ALAT));
lat_max=max(max(ALAT));
axis([long_min long_max lat_min lat_max]);
grid on, hold on
%
v=[0,1,2,4,6,8,10,12,14,16,18,20];
[C,h]=contourf(ALONG,ALAT,hs_rotated,v,'k');
clabel(C,h);                                                                       
    q = quiversc(ALONG,ALAT,u_rotated,v_rotated,'density',20);             % plot Peak Dir arrows rotated
    q.AutoScaleFactor=0.6;                                                 % lenght of arrow
    q.Color = 'black';
cb = colorbar;
cb.Label.String = 'H_{s} (m)';
axis([long_min long_max lat_min lat_max]);
axis square %axis equal
colormap(parula);
box on
xlabel('Longitude ({\circ})');ylabel('Latitude ({\circ})');
% draw requested point to find Hs (ex: buoy):
p1 = plot(lon_point,lat_point,'kd','LineWidth',1,'MarkerFaceColor','m');
% draw Storm Eye:
p2 = plot(Lon, Lat, '+k','LineWidth',1);
legend([p1 p2],{'buoy','storm eye'});
if Lat < 0 
    title('Hemisphere South rotated wave field (H_{s} )');
    sname1 = ['Original diagram rotated ', int2str(rotationAngle),char(176), ' (relative to North)'];
    subtitle(sname1,'Interpreter', 'none');   
else
    title('Hemisphere North rotated wave field (H_{s} )');
    sname1 = ['Original diagram rotated ', int2str(rotationAngle),char(176), ' (relative to North)'];
    subtitle(sname1,'Interpreter', 'none'); 
end
caxis([0 round(max(max(hs_rotated)))]);

%TextBox
str = {['H_{s}^{buoy} = ' num2str(Hs_rotated,'%.1f') ' m'],['Dir_{wave}^{buoy} = ' num2str(Dir_rotated,'%.0f') '^{\circ}']};
dim = [.19 .13 .15 .12];
ah = annotation('textbox',dim,'String',str,'FitBoxToText','on');
ah.BackgroundColor = 'w';                                                  % box background color
ah.FontSize = 8;                                                           % box fontsize

%--------------------
% Plot Tp field
%--------------------
figure;
plot(Lon,Lat,'k+');
axis('equal')
long_min=min(min(ALONG));
long_max=max(max(ALONG));
lat_min=min(min(ALAT));
lat_max=max(max(ALAT));
axis([long_min long_max lat_min lat_max]);
grid on,hold on
%
v=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];
[C,h]=contourf(ALONG,ALAT,tp_rotated,v,'k');
clabel(C,h);
cb = colorbar;
cb.Label.String = 'T_{p} (s)';
axis([long_min long_max lat_min lat_max]);
axis square; %axis equal
colormap(parula)
box on
xlabel('Longitude ({\circ})');ylabel('Latitude ({\circ})');
% draw requested point to find Hs (ex: buoy)
p1 = plot(lon_point,lat_point,'kd','LineWidth',1,'MarkerFaceColor','m');
% draw Storm Eye
p2 = plot(Lon, Lat, '+k','LineWidth',1);
legend([p1 p2],{'buoy','storm eye'});
if Lat < 0 
    title('Hemisphere South rotated wave field (T_{p} )');
    sname2 = ['Diagram rotated ', int2str(rotationAngle),char(176), ' (relative to North)'];
    subtitle(sname2,'Interpreter', 'none');   
else
    title('Hemisphere North rotated wave field (T_{p} )');
    sname2 = ['Diagram rotated ', int2str(rotationAngle),char(176), ' (relative to North)'];
    subtitle(sname2,'Interpreter', 'none'); 
end
caxis([0 round(max(max(tp_rotated)))]);

%TextBox
str = {['T_{p}^{buoy} = ' num2str(Tp_rotated,'%.1f') ' s'],['Dir_{wave}^{buoy} = ' num2str(Dir_rotated,'%.0f') '^{\circ}']};
dim = [.19 .13 .15 .12];
ah = annotation('textbox',dim,'String',str,'FitBoxToText','on');
ah.BackgroundColor = 'w';                                                  % box background color
ah.FontSize = 8;                                                           % box fontsize

end

% Guisela, October/2025
% ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>