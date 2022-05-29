% Task 1
%----------------------------------------------------------------------
clear all;
clc;

f = 2.4e9;             % 2.4GHz operational frequency
c = 3e8;               % speed of light
lambda = c/f;          % wavelength
beta = 2*pi/lambda;    % Beta
ft = 0.3048;           % 1 feet = 0.3048m
in = 0.0254;           % 1 inch = 0.0254m
erf = 2;               % relative permittivity of facade
sigf = 0.0001;         % conductivity of facade
sigg = 0.0001;         % conductivity of ground
ug = 7;                % relative permittivity of ground
u = 1;                 % permeability
factor = lambda/(4*pi);% parameter related to the antenna effect

% Knife edge dimension, provided in Project outline pdf file Fig.1
L1 = 34*in;
L2 = 6*in;
L3 = 18*in;
D1 = (17+(18/2)+(2.6/2))*in; 
D2 = (25+2.6/2)*in;

%Line of Sight
tx = 0; ty = 0; tz = 3.0525;     %tx(x,y,z) location
rx = 0; ry = 4; rz = 3.0525;     %rx(x,y,z) location
dist1 = ((5*2.6)+(4*25)+(2*17)+18)*8*in;  %Assume 8 sets of periodic wall
dist2 = ((17*5)+(5*2.6)+18/2)*in;
ryv = linspace(ry,dist1,125);
R = sqrt(((rx-tx)^2)+((ryv-ty).^2)+((rz-tz)^2)); %distance from tx to rx
E = (1./R).*exp(-1i*beta*R);    %line of sight
P = 20*log10(abs(E));           

%Ground reflection
tzGImage = -tz;     %transmitter image reflected by ground
Gnd_z = 0;
gradG = (rz-tzGImage)./ryv;      %gradient of reflected ray
pGz = (Gnd_z-tzGImage)./gradG;   %intersection on ground

rG1 = sqrt(((tz).^2)+((pGz).^2));     %Distance of 1st Ray
rG2 = sqrt(((rz).^2)+((ryv-pGz).^2)); %Distance of 2nd Ray
RGtotal = rG1+rG2;                    %Total Distance

anglei = atand((pGz-ty)./tz);
angler = asind((sind(anglei))./(sqrt(ug)));

%Refection Coefficient
ReflectionC = ((cosd(angler))-((sqrt(ug)).*(cosd(anglei))))./(cosd(angler)+((sqrt(ug)).*(cosd(anglei)))); 

%1st order ground reflected ray
Er_G = (1./RGtotal).*exp(-1i*beta*RGtotal).*(ReflectionC); 
%Ground Reflection Electric Field in dB%
PG = 20.*log10(abs(Er_G)); 

%Diffracted rays from wall 
d2r = pi/180;   %Conversion from degree to RAD
r2d = 180/pi;   %Conversion from RAD to degree
n = 2;          %wedge factor consant for for sheet)
btd = 90;       %Right angle

% Longer knife edges
dx1 = (10*ft)-(L3-L2)*in; %distance between the transmitter/receiver and wall
dy1 = (D2*4)+D1 ;
dyv1 = (dy1:dy1:dy1*5); %Range of vector 
% Shorter knife edges
dx2 = 10*ft; %distance between the transmitter/receiver and wall = 10ft
dy2 = 0;
dyv2 = (dy2:D2:D2*36); %Range of vector 
% Left side knife edge
dx3 = (10*ft)-(L1-L2)*in; %distance between the transmitter/receiver and wall
dy3 = -D1;

dz = tz;
phida = [];
Etd = [];
aEtd = [];
Edif = [];

% Calculations of E-field for knife edges #2 to #36
% Assume all the distances between the knife edges are equally distributed
for i = 1:size(dyv2,2)
    
    %edges #7, #13, #19, #25, and #31 are wider than their neighbours.
    if i == '6' or '12' or '18' or '24' or '30' or '36' or '42' or '48'
        sp = sqrt(((dx1-tx)^2)+((dyv2(i)-ty).^2)+((dz-tz)^2)); % incident distance
        s = sqrt(((rx-dx1)^2)+((ryv-dyv2(i)).^2)+((rz-dz)^2)); % observation distance
        phipd = atand(dx1/dyv2(i));  %incident angle
        cc = exp(-sqrt(-1)*btd*sp)/sp;
        pf = exp(-sqrt(-1)*btd*s); %phase factor
        A = sqrt(sp./(s.*(s+sp)));    %spherical spreading factor
        L = ((s.*sp)*(sind(btd))^2)./(s+sp)/lambda; % spherical distance parameter
    else
        sp = sqrt(((dx2-tx)^2)+((dyv2(i)-ty).^2)+((dz-tz)^2)); %incident distance
        s = sqrt(((rx-dx2)^2)+((ryv-dyv2(i)).^2)+((rz-dz)^2)); %observation distance
        phipd = atand(dx2/dyv2(i));  %incident angle
        cc = exp(-sqrt(-1)*btd*sp)/sp;
        pf = exp(-sqrt(-1)*btd*s); %phase factor
        A = sqrt(sp./(s.*(s+sp)));    %spherical spreading factor
        L = ((s.*sp)*(sind(btd))^2)./(s+sp)/lambda; %spherical distance parameter
    end
    
    for j = 1:size(ryv,2)  
        ny2v = ryv(j);
        if i == '6' or '12' or '18' or '24' or '30' or '36' or '42' or '48'
            if ny2v < dyv2
                phi = atand(dx1./(dyv2-ny2v));   %Observation angle
            else
                phi = 90 + atand((dyv2-ny2v)./dx1);
            end
        else
            if ny2v < dyv2
                phi = atand(dx2./(dyv2-ny2v));   %Observation angle
            else
                phi = 90 + atand((dyv2-ny2v)./dx2);
            end
        end
        phida = [phida,phi];
    end
    
    for k = 1:size(ryv,2)
        % Assume at edges #14, #20, #26, #27, #32, #33, #34,and #38, 
        % the incident rays were blocked. 
        % The shadowed edges would be excluded from diffraction field calculation.
        if i == '13' or '19' or '25' or '26' or '31' or '32' or '33' or '37'
            Edrays = 0;
            Etd = [Etd,Edrays];
        else
            [ds,dh,D] =  wdcpub(L(k),phida(k),phipd,btd,n);
            Edrays = cc*A(k)*pf(k)*ds*sqrt(lambda);
            Etd = [Etd,Edrays];
        end
    end
    aEtd(i,:) = Etd;
    phida = [];
    Etd = [];
end

%Calculations of E-field for knife edges #1
sp = sqrt(((dx3-tx)^2)+((dy3-ty).^2)+((dz-tz)^2)); %incident distance
s = sqrt(((rx-dx3)^2)+((ryv-dy3).^2)+((rz-dz)^2)); %observation distance
phipd = atand(dx3/dy3);  %incident angle
cc = exp(-sqrt(-1)*btd*sp)/sp;
pf = exp(-sqrt(-1)*btd*s); %phase factor
A = sqrt(sp./(s.*(s+sp)));    %spherical spreading factor
L = ((s.*sp)*(sind(btd))^2)./(s+sp)/lambda; %spherical distance parameter

for j = 1:size(ryv,2)
    ny2v = ryv(j);
    if ny2v < dy3
        phi = atand(dx3./(dy3-ny2v));   
    else
        phi = 90 + atand((dy3-ny2v)./dx3);
    end
    phida = [phida,phi];
end

for k = 1:size(ryv,2)
    [ds,dh,D] =  wdcpub(L(k),phida(k),phipd,btd,n);
    Edrays = cc*A(k)*pf(k)*ds*sqrt(lambda);
    Etd = [Etd,Edrays];
end

Edif = sum(aEtd);
Edif_total = Edif + Etd;
PKnifeEdge = 20*log10(abs(Edif_total));

figure('Name','TASK 1-Hong Sheng Sing;20018072');
plot(R,P,'c-',ryv,PG,'g-',ryv,PKnifeEdge,'b-');
legend ('Line of Sight','Ground Reflection','Knife Edge Diffraction');
xlabel('Distance (m)'),ylabel('Gain (dB)');
xlim([ry dist1]);
title("Knife Edge Approximation");
