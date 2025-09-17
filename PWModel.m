function [HS,TP,U,V,XX,YY] = PWModel(Vmax,Vfm,Rmax,R34,Lat)
%_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_
%Function to calculate the fetch using parameterised polynomial equation
%and parameters Hs, peak dir and Tp spatial distribution within a Tropical
%Cyclone
%by Guisela October, 2025
%_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_/(`_
                  
%==========================================================================
%Input parameters:
%Vmax = 'Maximum sustained wind speed' 'm/s'
%Vfm = 'Storm translation speed' 'm/s'
%R34 = 'Radius of 34 knot winds (mean of all quadrants)' 'meters'
%Rmax = 'Radius of maximum winds' 'meters'
%Lat = Latitude of Storm center 'degrees'

%Output parameters:
%Hs = Hs spatial distribution 'meters'
%Tp = Tp spatial distribution 'seconds'
%U = u component peak direction 
%V = v component peak direction 
%XX = grid position axis x 'km'
%YY = grid position axis y 'km'
%==========================================================================

% Reference:
% Grossmann-Matheson et al, 2025
% The spatial distribution of ocean wave parameters in tropical cyclones 
% Ocean Engineering (submitted)
% DOI: (to update)
% ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>

% Aplying limits to the values used in the parameterization
Vmax(Vmax > 78) = 78;Vmax(Vmax < 17) = 17;
Vfm(Vfm > 15) = 15; 
Rmax(Rmax > 60*10^3) = 60*10^3;Rmax(Rmax < 15*10^3) = 15*10^3;
R34(R34 > 400*10^3) = 400*10^3;R34(R34 < 200*10^3) = 200*10^3;

if Vmax > 78 | Vfm > 15 | Rmax > 60*10^3 | R34 > 400*10^3 ...
        | Vmax < 17 | Rmax < 15*10^3 | R34 < 200*10^3 == 1
warning('Parameter(s) over the limits. Default(s) will be used instead.')
end

% %Coefficients:
 a=0.54*10^0; b=-0.169*10^3; c=-1.442*10^3; d=0.3; e=0.143*10^2;  f=-0.043*10^3;  g=0.96*10^4;  h=4.47*10^3; i=1.0*10^5;
 C=0.1;
% 
% %============================
% %Parameterised Fetch Equation
% %============================
 F_P32 = (a.*Vmax.^3 + b.*Vmax.^2 + c.*Vfm.^2 + d.*Vmax.^2.*Vfm + e.*Vmax.*Vfm.^2 + f.*Vmax.*Vfm + g.*Vmax + h.*Vfm + i) .* exp(C.*Vfm);
%  
% %---------------------------------------------------------------------------------
% % Aply correction factors for Rmax (lambda) and for R34  (gamma)
% %---------------------------------------------------------------------------------
 lambda = 0.85*log10((Rmax)/(30*10^3))+1; %Rmax correction factor
% 
 gamma = 0.65*log10((R34)/(300*10^3)) + 1; %R34 correction factor
% 
%==============================================
%Final Model Fetch 
%==============================================
 Fetch = F_P32.*lambda.*gamma; %Final Fetch
% 
%Fetch in km:
 F = Fetch/1000; 
% 
%=============================
% Calculate Hs (maximum value)
%=============================
 gr = 9.81;
% Correction factor defined after validations (overall bias was 11%):
 alpha = 0.89;
 Hs_max = alpha * ((0.0016.*((gr.*F.*1000).^0.5).*Vmax)./gr);

% ===================================================================================
% Find Hs/Hsmax wavefield diagram(s) to get Hs correspondent to the calculated Hs(max)
% ===================================================================================

% Determine run name to load the spatial distributions to interpolate
% OBS: Below are included interpolation between the Hs/Hsmax diagrams when 
% the wind parameters are not exact the same on file names.

% Available runs combinations
 Vmax_runs = [17 30 40 50 65 78];
 Dp_runs = [10 30 50 70 110 150];                                          % runs used delta pressure instead Vmax                                            
 Vfm_runs = [0 2.5 5.0 7.5 10.0 12.5 15.0];                                
 R34_runs = [200 300 400];                                                 % R34 run names in km
 R_runs = [15 30 60];                                                      % Rmax run names in km

 % Determine Vmax run name for spatial distribution
 delVmax=abs(Vmax_runs-Vmax);
 temp = find(delVmax == 0); %  check if there Vmax is exact the same of the run (this case dont need interp)
 if isempty(temp) == 1    
    indice1 = find(delVmax == min(delVmax));  % Index of closest run
    if length(indice1)>1
    indice2 = indice1(2); 
    else    
    out=min(setdiff(delVmax,min(delVmax)));  
    indice2 = find(delVmax==out);          % Index second closest run
        if Vmax_runs(indice1) < (Vmax) & Vmax_runs(indice2) < (Vmax)
        out2 = find(Vmax_runs > (Vmax));
            if isempty(out2) == 1 
            indice2 = indice1;
            else
            indice2 = out2(1);
            end
        elseif Vmax_runs(indice1) > (Vmax) & Vmax_runs(indice2) > (Vmax)
        out3 = find(Vmax_runs < (Vmax));
        indice2 = out3(1);
        end
    end
 else
    indice1 = find(delVmax == 0);
    indice2 = find(delVmax == 0);
 end
 if length(indice2)>1
    indice2 = indice2(1);
 end
 t1b_1=int2str(Vmax_runs(indice1(1)));
 t1_1=int2str(Dp_runs(indice1(1)));
 t1b_2=int2str(Vmax_runs(indice2));
 t1_2=int2str(Dp_runs(indice2));
 
 vmax = [Vmax_runs(indice1(1)) Vmax_runs(indice2)];
 dp = [Dp_runs(indice1(1)) Dp_runs(indice2)];
 clear indice1 indice2 temp

 % Determine Vfm run name for spatial distribution
 delVfm=abs(Vfm_runs-Vfm);
 temp = find(delVfm == 0);%  check if there Vfm is exact the same of the run (this case dont need interp)
 if isempty(temp) == 1 
 indice1 = find(delVfm==min(delVfm));   % Index of closest run
    if length(indice1)>1
    indice2 = indice1(2); 
    else    
    out=min(setdiff(delVfm,min(delVfm)));  
    indice2 = find(delVfm==out);       % Index second closest run
    end
 else
 indice1 = find(delVfm == 0);
 indice2 = find(delVfm == 0);
 end  
 t2_1=int2str(Vfm_runs(indice1(1))*10);
 t2_2=int2str(Vfm_runs(indice2)*10);
 
 vfm = [Vfm_runs(indice1(1))  Vfm_runs(indice2)];
 clear indice1 indice2 temp

 % Determine R34 run name for spatial distribution
 delR34=abs(R34_runs-R34/10^3);
  temp = find(delR34 == 0); %  check if there R34 is exact the same of the run (this case dont need interp)
 if isempty(temp) == 1 
 indice1 = find(delR34==min(delR34)); % Index of closest run
    if length(indice1)>1
    indice2 = indice1(2); 
    else    
    out=min(setdiff(delR34,min(delR34)));  
    indice2 = find(delR34==out);       % Index second closest run
    end
 else
 indice1 = find(delR34 == 0);
 indice2 = find(delR34 == 0);
 end
 t3_1=int2str(R34_runs(indice1(1)));
 t3_2=int2str(R34_runs(indice2));
 
 r34 = [R34_runs(indice1(1))  R34_runs(indice2)];
 clear indice1 indice2 temp
 
%
% Determine Rmax run name for spatial distribution
 delR=abs(R_runs-Rmax/10^3);
 temp = find(delR == 0); %  check if there Rmax is exact the same of the run (this case dont need interp)
 if isempty(temp) == 1 
    indice1 = find(delR==min(delR));  % Index of closest run
     if length(indice1)>1
     indice2 = indice1(2); 
     else  
     out=min(setdiff(delR,min(delR)));  
     indice2 = find(delR==out);        % Index second closest run
        if R_runs(indice1)< (Rmax/10^3) & R_runs(indice2) < (Rmax/10^3)
        indice2 = find(R_runs > (Rmax/10^3));
            if isempty(indice2) == 1
                indice2 = indice1;
            end
        else
        end
     end
 else
 indice1 = find(delR == 0);
 indice2 = find(delR == 0);
 end
 t4_1=int2str(R_runs(indice1(1)));
 t4_2=int2str(R_runs(indice2));
 
 rmax = [R_runs(indice1(1))  R_runs(indice2)];
 clear indice1 indice2 temp
 
%-----------------------------------------------------------------------------------
% Load spatial distribution combinations for posterior wave parameters interpolation 
%-----------------------------------------------------------------------------------

fname1=['wave','_',t1_1,'_',t2_1,'_',t3_1,'_',t4_1];
eval(['load ',fname1]);                                                    % load chosen file
Z1 = Z; clear Z                                                            % Hs
T1 = Tp;clear Tp                                                           % Tp
U1 = upeak; clear upeak                                                    % u component
V1 = vpeak; clear vpeak                                                    % v component
xrm1 = xrm; yrm1 = yrm;                                                    % grid

fname2=['wave','_',t1_2,'_',t2_1,'_',t3_1,'_',t4_1];
eval(['load ',fname2]);
Z2 = Z; clear Z
T2 = Tp;clear Tp
U2 = upeak; clear upeak
V2 = vpeak; clear vpeak

fname3=['wave','_',t1_1,'_',t2_2,'_',t3_1,'_',t4_1];
eval(['load ',fname3]);
Z3 = Z; clear Z
T3 = Tp;clear Tp
U3 = upeak; clear upeak
V3 = vpeak; clear vpeak

fname4=['wave','_',t1_2,'_',t2_2,'_',t3_1,'_',t4_1];
eval(['load ',fname4]);
Z4 = Z; clear Z
T4 = Tp;clear Tp
U4 = upeak; clear upeak
V4 = vpeak; clear vpeak

fname5=['wave','_',t1_1,'_',t2_1,'_',t3_2,'_',t4_1];
eval(['load ',fname5]);
Z5 = Z; clear Z
T5 = Tp;clear Tp
U5 = upeak; clear upeak
V5 = vpeak; clear vpeak

fname6=['wave','_',t1_2,'_',t2_1,'_',t3_2,'_',t4_1];
eval(['load ',fname6]);
Z6 = Z; clear Z
T6 = Tp;clear Tp
U6 = upeak; clear upeak
V6 = vpeak; clear vpeak

fname7=['wave','_',t1_1,'_',t2_2,'_',t3_2,'_',t4_1];
eval(['load ',fname7]);
Z7 = Z; clear Z
T7 = Tp;clear Tp
U7 = upeak; clear upeak
V7 = vpeak; clear vpeak

fname8=['wave','_',t1_2,'_',t2_2,'_',t3_2,'_',t4_1];
eval(['load ',fname8]);
Z8 = Z; clear Z
T8 = Tp;clear Tp
U8 = upeak; clear upeak
V8 = vpeak; clear vpeak

fname9=['wave','_',t1_1,'_',t2_1,'_',t3_1,'_',t4_2];
eval(['load ',fname9]);
Z9 = Z; clear Z
T9 = Tp;clear Tp
U9 = upeak; clear upeak
V9 = vpeak; clear vpeak

fname10=['wave','_',t1_2,'_',t2_1,'_',t3_1,'_',t4_2];
eval(['load ',fname10]);
Z10 = Z; clear Z
T10 = Tp;clear Tp
U10 = upeak; clear upeak
V10 = vpeak; clear vpeak

fname11=['wave','_',t1_1,'_',t2_2,'_',t3_1,'_',t4_2];
eval(['load ',fname11]);
Z11 = Z; clear Z
T11 = Tp;clear Tp
U11 = upeak; clear upeak
V11 = vpeak; clear vpeak

fname12=['wave','_',t1_2,'_',t2_2,'_',t3_1,'_',t4_2];
eval(['load ',fname12]);
Z12 = Z; clear Z
T12 = Tp;clear Tp
U12 = upeak; clear upeak
V12 = vpeak; clear vpeak

fname13=['wave','_',t1_1,'_',t2_1,'_',t3_2,'_',t4_2];
eval(['load ',fname13]);
Z13 = Z; clear Z
T13 = Tp;clear Tp
U13 = upeak; clear upeak
V13 = vpeak; clear vpeak

fname14=['wave','_',t1_2,'_',t2_1,'_',t3_2,'_',t4_2];
eval(['load ',fname14]);
Z14 = Z; clear Z
T14 = Tp;clear Tp
U14 = upeak; clear upeak
V14 = vpeak; clear vpeak

fname15=['wave','_',t1_1,'_',t2_2,'_',t3_2,'_',t4_2];
eval(['load ',fname15]);
Z15 = Z; clear Z
T15 = Tp;clear Tp
U15 = upeak; clear upeak
V15 = vpeak; clear vpeak

fname16=['wave','_',t1_2,'_',t2_2,'_',t3_2,'_',t4_2];
eval(['load ',fname16]);
Z16 = Z; clear Z
T16 = Tp;clear Tp
U16 = upeak; clear upeak
V16 = vpeak; clear vpeak

% Preallocate ZZ
  ZZ=zeros(size(Z1));
  TT=zeros(size(T1));
  UU=zeros(size(U1));
  VV=zeros(size(V1));

%-----------------------------------------------------------------------
% Calculate final Z,T,U,V (Hs,Tp,u and v interpolated from combinations)
%-----------------------------------------------------------------------

if vmax(1)~=vmax(2) & vfm(1)~=vfm(2) & rmax(1)~=rmax(2) & r34(1)~=r34(2)
       
for xx=1:length(xrm1)
    for yy=1:length(yrm1);
      z(1,1,1,1)=Z1(xx,yy);
      z(2,1,1,1)=Z2(xx,yy);
      z(1,2,1,1)=Z3(xx,yy);
      z(2,2,1,1)=Z4(xx,yy);
      z(1,1,2,1)=Z5(xx,yy);
      z(2,1,2,1)=Z6(xx,yy);
      z(1,2,2,1)=Z7(xx,yy);
      z(2,2,2,1)=Z8(xx,yy);
      z(1,1,1,2)=Z9(xx,yy);
      z(2,1,1,2)=Z10(xx,yy);
      z(1,2,1,2)=Z11(xx,yy);
      z(2,2,1,2)=Z12(xx,yy);
      z(1,1,2,2)=Z13(xx,yy);
      z(2,1,2,2)=Z14(xx,yy);
      z(1,2,2,2)=Z15(xx,yy);
      z(2,2,2,2)=Z16(xx,yy);
      ZZ(xx,yy)=interpn(vmax,vfm,r34,rmax,z,Vmax,Vfm,(R34/10^3),(Rmax/10^3));
      t(1,1,1,1)=T1(xx,yy);
      t(2,1,1,1)=T2(xx,yy);
      t(1,2,1,1)=T3(xx,yy);
      t(2,2,1,1)=T4(xx,yy);
      t(1,1,2,1)=T5(xx,yy);
      t(2,1,2,1)=T6(xx,yy);
      t(1,2,2,1)=T7(xx,yy);
      t(2,2,2,1)=T8(xx,yy);
      t(1,1,1,2)=T9(xx,yy);
      t(2,1,1,2)=T10(xx,yy);
      t(1,2,1,2)=T11(xx,yy);
      t(2,2,1,2)=T12(xx,yy);
      t(1,1,2,2)=T13(xx,yy);
      t(2,1,2,2)=T14(xx,yy);
      t(1,2,2,2)=T15(xx,yy);
      t(2,2,2,2)=T16(xx,yy);
      TT(xx,yy)=interpn(vmax,vfm,r34,rmax,t,Vmax,Vfm,(R34/10^3),(Rmax/10^3));
      u(1,1,1,1)=U1(xx,yy);
      u(2,1,1,1)=U2(xx,yy);
      u(1,2,1,1)=U3(xx,yy);
      u(2,2,1,1)=U4(xx,yy);
      u(1,1,2,1)=U5(xx,yy);
      u(2,1,2,1)=U6(xx,yy);
      u(1,2,2,1)=U7(xx,yy);
      u(2,2,2,1)=U8(xx,yy);
      u(1,1,1,2)=U9(xx,yy);
      u(2,1,1,2)=U10(xx,yy);
      u(1,2,1,2)=U11(xx,yy);
      u(2,2,1,2)=U12(xx,yy);
      u(1,1,2,2)=U13(xx,yy);
      u(2,1,2,2)=U14(xx,yy);
      u(1,2,2,2)=U15(xx,yy);
      u(2,2,2,2)=U16(xx,yy);
      UU(xx,yy)=interpn(vmax,vfm,r34,rmax,u,Vmax,Vfm,(R34/10^3),(Rmax/10^3));
      v(1,1,1,1)=V1(xx,yy);
      v(2,1,1,1)=V2(xx,yy);
      v(1,2,1,1)=V3(xx,yy);
      v(2,2,1,1)=V4(xx,yy);
      v(1,1,2,1)=V5(xx,yy);
      v(2,1,2,1)=V6(xx,yy);
      v(1,2,2,1)=V7(xx,yy);
      v(2,2,2,1)=V8(xx,yy);
      v(1,1,1,2)=V9(xx,yy);
      v(2,1,1,2)=V10(xx,yy);
      v(1,2,1,2)=V11(xx,yy);
      v(2,2,1,2)=V12(xx,yy);
      v(1,1,2,2)=V13(xx,yy);
      v(2,1,2,2)=V14(xx,yy);
      v(1,2,2,2)=V15(xx,yy);
      v(2,2,2,2)=V16(xx,yy);
      VV(xx,yy)=interpn(vmax,vfm,r34,rmax,v,Vmax,Vfm,(R34/10^3),(Rmax/10^3));
    end
end

elseif rmax(1)==rmax(2) & vmax(1)~=vmax(2) & vfm(1)~=vfm(2) & r34(1)~=r34(2)
for xx=1:length(xrm1)
    for yy=1:length(yrm1);
      z(1,1,1,1)=Z1(xx,yy);
      z(2,1,1,1)=Z2(xx,yy);
      z(1,2,1,1)=Z3(xx,yy);
      z(2,2,1,1)=Z4(xx,yy);
      z(1,1,2,1)=Z5(xx,yy);
      z(2,1,2,1)=Z6(xx,yy);
      z(1,2,2,1)=Z7(xx,yy);
      z(2,2,2,1)=Z8(xx,yy);
      ZZ(xx,yy)=interpn(vmax,vfm,r34,z,Vmax,Vfm,(R34/10^3));
      t(1,1,1,1)=T1(xx,yy);
      t(2,1,1,1)=T2(xx,yy);
      t(1,2,1,1)=T3(xx,yy);
      t(2,2,1,1)=T4(xx,yy);
      t(1,1,2,1)=T5(xx,yy);
      t(2,1,2,1)=T6(xx,yy);
      t(1,2,2,1)=T7(xx,yy);
      t(2,2,2,1)=T8(xx,yy);
      TT(xx,yy)=interpn(vmax,vfm,r34,t,Vmax,Vfm,(R34/10^3));
      u(1,1,1,1)=U1(xx,yy);
      u(2,1,1,1)=U2(xx,yy);
      u(1,2,1,1)=U3(xx,yy);
      u(2,2,1,1)=U4(xx,yy);
      u(1,1,2,1)=U5(xx,yy);
      u(2,1,2,1)=U6(xx,yy);
      u(1,2,2,1)=U7(xx,yy);
      u(2,2,2,1)=U8(xx,yy);
      UU(xx,yy)=interpn(vmax,vfm,r34,u,Vmax,Vfm,(R34/10^3));
      v(1,1,1,1)=V1(xx,yy);
      v(2,1,1,1)=V2(xx,yy);
      v(1,2,1,1)=V3(xx,yy);
      v(2,2,1,1)=V4(xx,yy);
      v(1,1,2,1)=V5(xx,yy);
      v(2,1,2,1)=V6(xx,yy);
      v(1,2,2,1)=V7(xx,yy);
      v(2,2,2,1)=V8(xx,yy);
      VV(xx,yy)=interpn(vmax,vfm,r34,v,Vmax,Vfm,(R34/10^3));
    end
end

elseif rmax(1)==rmax(2) & r34(1)==r34(2) & vmax(1)~=vmax(2) & vfm(1)~=vfm(2) 
for xx=1:length(xrm1)
    for yy=1:length(yrm1);
      z(1,1,1,1)=Z1(xx,yy);
      z(2,1,1,1)=Z2(xx,yy);
      z(1,2,1,1)=Z3(xx,yy);
      z(2,2,1,1)=Z4(xx,yy);
      ZZ(xx,yy)=interpn(vmax,vfm,z,Vmax,Vfm);
      t(1,1,1,1)=T1(xx,yy);
      t(2,1,1,1)=T2(xx,yy);
      t(1,2,1,1)=T3(xx,yy);
      t(2,2,1,1)=T4(xx,yy);
      TT(xx,yy)=interpn(vmax,vfm,t,Vmax,Vfm);
      u(1,1,1,1)=U1(xx,yy);
      u(2,1,1,1)=U2(xx,yy);
      u(1,2,1,1)=U3(xx,yy);
      u(2,2,1,1)=U4(xx,yy);
      UU(xx,yy)=interpn(vmax,vfm,u,Vmax,Vfm);
      v(1,1,1,1)=V1(xx,yy);
      v(2,1,1,1)=V2(xx,yy);
      v(1,2,1,1)=V3(xx,yy);
      v(2,2,1,1)=V4(xx,yy);
      VV(xx,yy)=interpn(vmax,vfm,v,Vmax,Vfm);
    end
end

elseif rmax(1)==rmax(2) & r34(1)==r34(2) & vfm(1)==vfm(2) & vmax(1)~=vmax(2)  
for xx=1:length(xrm1)
    for yy=1:length(yrm1);
      z(1,1,1,1)=Z1(xx,yy);
      z(2,1,1,1)=Z2(xx,yy);
      ZZ(xx,yy)=interpn(vmax,z,Vmax);
      t(1,1,1,1)=T1(xx,yy);
      t(2,1,1,1)=T2(xx,yy);
      TT(xx,yy)=interpn(vmax,t,Vmax);
      u(1,1,1,1)=U1(xx,yy);
      u(2,1,1,1)=U2(xx,yy);
      UU(xx,yy)=interpn(vmax,u,Vmax);
      v(1,1,1,1)=V1(xx,yy);
      v(2,1,1,1)=V2(xx,yy);
      VV(xx,yy)=interpn(vmax,v,Vmax);
    end
end

elseif rmax(1)==rmax(2) & r34(1)==r34(2) & vfm(1)~=vfm(2) & vmax(1)==vmax(2)  
for xx=1:length(xrm1)
    for yy=1:length(yrm1);
      z(1,1,1,1)=Z1(xx,yy);
      z(2,1,1,1)=Z3(xx,yy);
      ZZ(xx,yy)=interpn(vfm,z,Vfm);
      t(1,1,1,1)=T1(xx,yy);
      t(2,1,1,1)=T3(xx,yy);
      TT(xx,yy)=interpn(vfm,t,Vfm);
      u(1,1,1,1)=U1(xx,yy);
      u(2,1,1,1)=U3(xx,yy);
      UU(xx,yy)=interpn(vfm,u,Vfm);
      v(1,1,1,1)=V1(xx,yy);
      v(2,1,1,1)=V3(xx,yy);
      VV(xx,yy)=interpn(vfm,v,Vfm);
    end
end

elseif rmax(1)==rmax(2) & r34(1)~=r34(2) & vfm(1)~=vfm(2) & vmax(1)==vmax(2)  
for xx=1:length(xrm1)
    for yy=1:length(yrm1);
      z(1,1,1,1)=Z1(xx,yy);
      z(2,1,1,1)=Z3(xx,yy);
      z(1,2,1,1)=Z5(xx,yy);
      z(2,2,1,1)=Z7(xx,yy);
      ZZ(xx,yy)=interpn(vfm,r34,z,Vfm,(R34/10^3));
      t(1,1,1,1)=T1(xx,yy);
      t(2,1,1,1)=T3(xx,yy);
      t(1,2,1,1)=T5(xx,yy);
      t(2,2,1,1)=T7(xx,yy);
      TT(xx,yy)=interpn(vfm,r34,t,Vfm,(R34/10^3));
      u(1,1,1,1)=U1(xx,yy);
      u(2,1,1,1)=U3(xx,yy);
      u(1,2,1,1)=U5(xx,yy);
      u(2,2,1,1)=U7(xx,yy);
      UU(xx,yy)=interpn(vfm,r34,u,Vfm,(R34/10^3));
      v(1,1,1,1)=V1(xx,yy);
      v(2,1,1,1)=V3(xx,yy);
      v(1,2,1,1)=V5(xx,yy);
      v(2,2,1,1)=V7(xx,yy);
      VV(xx,yy)=interpn(vfm,r34,v,Vfm,(R34/10^3));
    end
end

elseif rmax(1)~=rmax(2) & r34(1)==r34(2) & vfm(1)~=vfm(2) & vmax(1)==vmax(2)  
for xx=1:length(xrm1)
    for yy=1:length(yrm1);
      z(1,1,1,1)=Z1(xx,yy);
      z(2,1,1,1)=Z3(xx,yy);
      z(1,2,1,1)=Z9(xx,yy);
      z(2,2,1,1)=Z11(xx,yy);
      ZZ(xx,yy)=interpn(vfm,rmax,z,Vfm,(Rmax/10^3));
      t(1,1,1,1)=T1(xx,yy);
      t(2,1,1,1)=T3(xx,yy);
      t(1,2,1,1)=T9(xx,yy);
      t(2,2,1,1)=T11(xx,yy);
      TT(xx,yy)=interpn(vfm,rmax,t,Vfm,(Rmax/10^3));
      u(1,1,1,1)=U1(xx,yy);
      u(2,1,1,1)=U3(xx,yy);
      u(1,2,1,1)=U9(xx,yy);
      u(2,2,1,1)=U11(xx,yy);
      UU(xx,yy)=interpn(vfm,rmax,u,Vfm,(Rmax/10^3));
      v(1,1,1,1)=V1(xx,yy);
      v(2,1,1,1)=V3(xx,yy);
      v(1,2,1,1)=V9(xx,yy);
      v(2,2,1,1)=V11(xx,yy);
      VV(xx,yy)=interpn(vfm,rmax,v,Vfm,(Rmax/10^3));
    end
end

elseif rmax(1)~=rmax(2) & r34(1)~=r34(2) & vfm(1)==vfm(2) & vmax(1)==vmax(2)  
for xx=1:length(xrm1)
    for yy=1:length(yrm1);
      z(1,1,1,1)=Z1(xx,yy);
      z(2,1,1,1)=Z5(xx,yy);
      z(1,2,1,1)=Z9(xx,yy);
      z(2,2,1,1)=Z13(xx,yy);
      ZZ(xx,yy)=interpn(r34,rmax,z,(R34/10^3),(Rmax/10^3));
      t(1,1,1,1)=T1(xx,yy);
      t(2,1,1,1)=T5(xx,yy);
      t(1,2,1,1)=T9(xx,yy);
      t(2,2,1,1)=T13(xx,yy);
      TT(xx,yy)=interpn(r34,rmax,t,(R34/10^3),(Rmax/10^3));
      u(1,1,1,1)=U1(xx,yy);
      u(2,1,1,1)=U5(xx,yy);
      u(1,2,1,1)=U9(xx,yy);
      u(2,2,1,1)=U13(xx,yy);
      UU(xx,yy)=interpn(r34,rmax,u,(R34/10^3),(Rmax/10^3));
      v(1,1,1,1)=V1(xx,yy);
      v(2,1,1,1)=V5(xx,yy);
      v(1,2,1,1)=V9(xx,yy);
      v(2,2,1,1)=V13(xx,yy);
      VV(xx,yy)=interpn(r34,rmax,v,(R34/10^3),(Rmax/10^3));
    end
end

else
for xx=1:length(xrm1)                                                       % case when none interpolation is needed
    for yy=1:length(yrm1)
      z(1,1,1,1)=Z1(xx,yy);
      ZZ(xx,yy)=Z1(xx,yy); 
      t(1,1,1,1)=T1(xx,yy);
      TT(xx,yy)=T1(xx,yy); 
      u(1,1,1,1)=U1(xx,yy);
      UU(xx,yy)=U1(xx,yy);
      v(1,1,1,1)=V1(xx,yy);
      VV(xx,yy)=V1(xx,yy);
    end
end    
end

Znew = ZZ;                                                                 % Final interpolated values
Tnew = TT;                                                                 
Unew = UU;
Vnew = VV;

%-------------------------------------------------
% Determine wave heights (resultant Hs wave field)
%-------------------------------------------------
 xcontent=xrm;
 ycontent=yrm;
  if Lat<0                                                                 % Southern Hemisphere case
 zcontent = fliplr(Znew');                                                 
 tcontent = fliplr(Tnew');  
 ucontent = fliplr(-Unew');
 vcontent = fliplr(Vnew');
  else                                                                     % Northern Hemisphere case
 zcontent = Znew'; 
 tcontent = Tnew'; 
 ucontent = Unew';
 vcontent = Vnew';
 end
 
HS= zcontent.*Hs_max;                                                      % Significant Wave Height (m)   
TP = tcontent;                                                             % Peak Period (s)
U = ucontent;
V = vcontent;


% Prepare to save result and plot wave field (Storm direction is up to page)
 Req = 30;                                                                 % Spatial distribution Hs/Hsmax was built using R=30km
 XX = xcontent*Req; YY = ycontent*Req;  
 
% Note: The extent of the spatial distribution diagrams is limited to about 300km

% ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
% Plot Interpolated Spatial Distribution 
% ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>

% Prepare the file name
filenamehsdir = ['Hsdir_field_TC_' int2str(Vmax) '_' int2str(Vfm) '_' int2str(Rmax/10^3) '_' int2str(R34/10^3)];
filenametp = ['Tp_field_TC_' int2str(Vmax) '_' int2str(Vfm) '_' int2str(Rmax/10^3) '_' int2str(R34/10^3)];
%
  if Lat < 0 
       sname = ['Hemisphere South'];   
  else
       sname = ['Hemisphere North'];
  end

%==========================================================================
%Plot wave field Tp calculated
%==========================================================================

% fig = figure('name', filename_tp, 'numbertitle', 'off');
figure
    contourf(XX,YY,TP);
    hold on
    cb = colorbar;
    cb.Label.String = 'Tp (s)';
% Plot Contour Lines
    v=[2:2:20];
    [C,h]=contourf(XX,YY,TP,v,'k');
    clabel(C,h,v);
%axis('square')
  xlabel('Distance (km)'); xticks(-300:100:300);
  ylabel('Distance (km)'); yticks(-300:100:300);
  %draw the center of rotation:
  x_center = xrm(141); y_center = yrm(141);
  plot(x_center, y_center, 'k+','LineWidth',2); 
  title('Wave field interpolated diagram (T_{p} )');
tpm = max(max(TP));
subname1 = ['T_{p}^{max} = ', num2str(tpm,'%.1f'),' s'];
subtitle(sname);
%TextBox 1
str1 = {(subname1)};
dim = [.12 .79 .14 .07];
ah = annotation('textbox',dim,'String',str1,'FitBoxToText','on');
ah.BackgroundColor = 'w';                                                  % box background color
ah.FontSize = 8;                                                           % box fontsize
%TextBox 2
str2 = {['V_{max} = ' int2str(Vmax) ' m/s'],['V_{fm} = ' num2str(Vfm,'%.1f') ' m/s'],['R_{max} = ' int2str(Rmax/10^3) ' km'],['R_{34} = ' int2str(R34/10^3) ' km']};
dim = [.12 .12 .14 .18];
ah = annotation('textbox',dim,'String',str2,'FitBoxToText','on');
ah.BackgroundColor = 'w';                                                  
ah.FontSize = 7;   
%
% caxis([0 20]);
caxis([0 round(max(max(TP)))])
colormap(parula);
%
% saveas(gcf,filenametp,'png')

%==========================================================================
%Plot wave field Hs calculated with Peak direction (interpolated)
%==========================================================================
figure;
  v1=[0:2:20];
  [cs,h] = contourf(XX,YY,HS,v1);
  clabel(cs,h);
  xlabel('Distance (km)'); xticks(-300:100:300);
  ylabel('Distance (km)'); yticks(-300:100:300);
  hold on
  %plot Peak Dir arrows
    q = quiversc(XX,YY,U,V,'density',20);
    q.Color = 'black';
  cb = colorbar;
  cb.Label.String = 'H_{s} (m)';
  colormap(parula);
ylim([-280 280]);
xlim([-280 280]);
%draw the center of rotation
  hold on
  x_center = xcontent(141); y_center = ycontent(141); 
  plot(x_center, y_center, 'k+'); 
title('Wave field interpolated diagram (H_{s} )');
hsm = max(max(HS));
subname2 = ['H_{s}^{max} = ', num2str(hsm,'%.1f'),' s'];
%TextBox 1
str3 = {['V_{max} = ' int2str(Vmax) ' m/s'],['V_{fm} = ' num2str(Vfm,'%.1f') ' m/s'],['R_{max} = ' int2str(Rmax/10^3) ' km'],['R_{34} = ' int2str(R34/10^3) ' km']};
dim = [.12 .12 .14 .18];
ah = annotation('textbox',dim,'String',str3,'FitBoxToText','on');
ah.BackgroundColor = 'w';                                                  % box background color
ah.FontSize = 7;                                                           % box fontsize
  %
  if Lat < 0 
       sname = ['Hemisphere South'];   
  else
       sname = ['Hemisphere North'];
  end
subtitle(sname);
%TextBox 2
str4 = {(subname2)};
dim = [.12 .79 .14 .07];
ah = annotation('textbox',dim,'String',str4,'FitBoxToText','on');
ah.BackgroundColor = 'w';                                                  % box background color
ah.FontSize = 8;                                                           % box fontsize
%
% caxis([0 20]);
caxis([0 round(max(max(HS)))])
colormap(parula);
% saveas(gcf,filenamehsdir,'png')
 
end

% Guisela, October/2024
% ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>

