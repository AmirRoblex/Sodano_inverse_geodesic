clc
clear all
format long g
%%

% B1 = deg2rad(dms2degrees([36 40 8.09]));
% L1 = deg2rad(dms2degrees([48 30 10.41]));
% B2 = deg2rad(dms2degrees([36 40 45.68]));
% L2 = deg2rad(dms2degrees([48 30 19.8]));
B1 = 0.641048186580116;
L1 = 0.843962686448118;
B2 = deg2rad(dms2degrees([45 00 00]));
L2 = deg2rad(dms2degrees([106 00 00]));

% a0 = 6378388;
% b0 = 6356911.946;
% f = 1-(b0/a0); %Paper's Example Ellipsoid

a0 = 6378137;
b0 = 6356752.31424518;
f = 1-(b0/a0); %wgs84

% a0 = 6378206.4;
% b0 = 6356583.8;
% f = 294.9786982; %Clarke (1866)

L = L2 - L1;

tan_Beta1 = tan(B1)*(1-f);
Beta1 = atan(tan_Beta1);
cos_Beta1 = cos(Beta1);
sin_Beta1 = sin(Beta1);

tan_Beta2 = tan(B2)*(1-f);
Beta2 = atan(tan_Beta2);
cos_Beta2 = cos(Beta2);
sin_Beta2 = sin(Beta2);

a = sin_Beta1 * sin_Beta2;
b = cos_Beta1 * cos_Beta2;
n = (a0 - b0) / (a0 + b0);

sin_L = sin(L);
cos_L = cos(L);
B2_B1 = (B2 - B1) + (2*(sin(Beta2-Beta1))) * (((n+n^2+n^3)*a) - ((n-n^2+n^3)*b));
cos_phi = a + (b * cos_L);
% sin_phi = sqrt(((sin_L * cos_Beta2)^2) + (((sin_Beta2 * cos_Beta1) - (sin_Beta1 * cos_Beta2 * cos_L))^2));     
sin_phi = sqrt(((sin_L * cos_Beta2)^2) + (((sin(B2_B1)) + (2*cos_Beta2*sin_Beta1*(sin(L/2)^2)))^2));
phi = acos(cos_phi);

c = (b * sin_L) / sin_phi;
m = 1 - (c^2);

%
sp1 = (1+f+f^2)*phi;

sp2 = a * (((f+f^2) * sin_phi) - (((f^2/2)*phi^2) * csc(phi)));
sp2 = (((f+f^2)*m*phi) / 2);

sp3 = m * ((-((f+f^2)/2)*phi) - ((f+f^2)/2)*sin_phi*cos_phi + (f^2/2)*phi^2*cos_phi);  %%
sp3 = (((f+f^2)*((2*a)-(m*cos_phi))*sin_phi) / 2);

sp4 = (a^2) * ((-(f^2/2)) * sin_phi * cos_phi);
sp4 = ((f^2*m^2)*(phi+(sin_phi*cos_phi)) / 16);

sp5 = (m^2) * (((f^2/16) * phi) + ((f^2/16) * sin_phi * cos_phi) - ((f^2/2) * phi^2 * cot(phi)) - ((f^2/8) * sin_phi * cos_phi^3));
sp5 = ((f^2 * ((2*a) - (m*cos_phi))^2 * (sin_phi*cos_phi)) / 8);

sp6 = (a * m) * (((f^2/2) * phi^2 * csc(phi)) + ((f^2/2) * sin_phi * cos_phi^2));
sp6 = ((f^2 * (1-m) * (a-(m*cos_phi)) * phi^2 * csc(phi)) / 2);

s0_b0 = sp1 + sp2 + sp3 + sp4 + sp5 + sp6;
s0_b0 = sp1 - sp2 + sp3 + sp4 - sp5 - sp6;
s0_mat = [sp1 ; sp2 ; sp3 ; sp4 ; sp5 ; sp6];
s = s0_b0 * b0;

%
lp1 = ((f+f^2) * phi);
lp2 = a * ((-(f^2/2)*sin_phi) - ((f^2) * phi^2 * csc(phi)));
lp3 = m * ((-((5*f^2)/4) * phi) + ((f^2/4) * sin_phi * cos_phi) + ((f^2) * phi^2 * cot(phi)));
lambda_L_c = lp1 + lp2 + lp3;

lambda = (lambda_L_c * c) + L;
sin_lambda = sin(lambda);
cos_lambda = cos(lambda);
sin2_lambda = sin(lambda/2)^2;


cot_az12 = ((sin_Beta2 * cos_Beta1) - (cos_lambda * sin_Beta1 * cos_Beta2)) / (sin_lambda * cos_Beta2);
az_12 = acot(cot_az12);
cot_az21 = ((sin_Beta2 * cos_Beta1 * cos_lambda) - (sin_Beta1 * cos_Beta2)) / (sin_lambda * cos_Beta1);
az_21 = acot(cot_az21);

if L>=0
    if cot(tan(az_12)) >= 0
        az_12 = az_12;
    elseif cot(tan(az_12)) < 0
        az_12 = az_12+pi;
    end
    if cot(tan(az_21)) >= 0
        az_21 = az_21 + pi;
    elseif cot(tan(az_21)) < 0
        az_21 = az_21+2*pi;
    end
elseif L<0
    if cot(tan(az_12)) >= 0
        az_12 = az_12 + pi;
    elseif cot(tan(az_12)) < 0
        az_12 = az_12 + 2*pi;
    end
    if cot(tan(az_21)) >= 0
        az_21 = az_21;
    elseif cot(tan(az_21)) < 0
        az_21 = az_21 + pi;
    end
end

disp('S = ')
disp(s)

disp('az_12 = ')
disp(degrees2dms(rad2deg(az_12)))

disp('az_21 = ')
disp(degrees2dms(rad2deg(az_21)))


E_wgs84=referenceEllipsoid('wgs84');
[length, az]=distance(dms2degrees([36 40 8.09]), dms2degrees([48 30 10.41]), dms2degrees([36 40 45.68]), dms2degrees([48 30 19.8]), E_wgs84);
% % [length, az]=distance(dms2degrees([20 00 00]), dms2degrees([00 00 00]), dms2degrees([45 00 00]), dms2degrees([106 00 00]), E_wgs84);
% [length, az]=distance(dms2degrees([45 00 00]), dms2degrees([106 00 00]),  dms2degrees([20 00 00]), dms2degrees([00 00 00]), E_wgs84);