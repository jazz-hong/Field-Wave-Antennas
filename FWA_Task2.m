% Task 2
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

%Line of Sight
tx = 0; ty = 8; tz = 4.5;     %tx(x,y,z) location
rx = 0; ry = 15; rz = 4.5;    %rx(x,y,z) location
dist1 = 116.05;
ryv = (ry:0.5:dist1);        %Took only half distance from the map captured
dist = sqrt(((rx-tx)^2)+((ryv-ty).^2)+((rz-tz)^2)); %distance from tx to rx
E = (1./dist).*exp(-1i*beta*dist);    %line of sight
P = 20*log10(abs(E));               %direct ray alone

%Ground reflection
tzGImage = -tz;     %transmitter image reflected by ground
Gnd_z = 0;
gradG = (rz-tzGImage)./ryv;     %gradient of reflected ray
pGz = (Gnd_z-tzGImage)./gradG;   %intersection on ground

rG1 = sqrt(((tz).^2)+((pGz).^2));     %Distance of 1st Ray
rG2 = sqrt(((rz).^2)+((ryv-pGz).^2)); %Distance of 2nd Ray
RGtotal = rG1+rG2; %Total Distance

anglei = atand((pGz-ty)./tz);
angler = asind((sind(anglei))./(sqrt(ug)));

%Refection Coefficient
rcG = ((cosd(angler))-((sqrt(ug)).*(cosd(anglei))))./(cosd(angler)+((sqrt(ug)).*(cosd(anglei)))); 

%1st order ground reflected ray
Er_G = (1./RGtotal).*exp(-1i*beta*RGtotal).*(rcG); 
%Ground Reflection Electric Field in dB%
EG = 20.*log10(abs(Er_G)); 

%Diffracted rays from wall 
d2r = pi/180;   %Conversion from degree to RAD
r2d = 180/pi;   %Conversion from RAD to degree
n = 2;          %wedge factor consant for for sheet)
btd = 90;       %Right angle

%Knife edge approximation for the pillars of terrace houses
dx1 = 5; %distance between the transmitter/receiver and wall
dy1 = 0;
dyv1 = (0:6:116.05);  %7 pillars, 6 gaps

dz = tz;
phida = [];
Etd = [];
aEtd = [];
Edif = [];

for i = 1:size(dyv1,2)
    sp = sqrt(((dx1-tx)^2)+((dyv1(i)-ty).^2)+((dz-tz)^2)); %incident distance
    s = sqrt(((rx-dx1)^2)+((ryv-dyv1(i)).^2)+((rz-dz)^2)); %observation distance
    phipd = atand(dx1/dyv1(i));  %incident angle
    cc = exp(-sqrt(-1)*btd*sp)/sp;
    pf = exp(-sqrt(-1)*btd*s); %phase factor
    A = sqrt(sp./(s.*(s+sp)));    %spherical spreading factor
    L = ((s.*sp)*(sind(btd))^2)./(s+sp)/lambda; %spherical distance parameter
    
    for j = 1:size(ryv,2)
        ny2v = ryv(j);
        if ny2v < dyv1
            phi = atand(dx1./(dyv1-ny2v));   %observation angle
        else
            phi = 90 + atand((dyv1-ny2v)./dx1);
        end
        phida = [phida,phi];
    end
    
    for k = 1:size(ryv,2)
        [ds,dh,D] =  wdcpub(L(k),phida(k),phipd,btd,n);
        Edrays = cc*A(k)*pf(k)*ds*sqrt(lambda);
        Etd = [Etd,Edrays];
    end
    aEtd(i,:) = Etd;
    phida = [];
    Etd = [];
end

Edif = sum(aEtd);
PKnifeEdge = 20*log10(abs(Edif));

figure('Name','TASK 2-Hong Sheng Sing;20018072');
plot(dist+8,P,'c-',ryv,EG,'g-',ryv,PKnifeEdge,'b-');
legend ('Line of Sight','Ground Reflection','Knife Edge Diffraction');
xlabel('Distance (m)'),ylabel('Gain (dB)');
xlim([ry dist1]);
title("Knife Edge Approximation for Jalan Taman Tasik Semenyih 5 (TTS5) Terrace House");
