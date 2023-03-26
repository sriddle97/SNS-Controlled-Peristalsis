clear
clc
% hmax = 11;              %height max, cm
% hmin = 6.5;             %height min, cm
% l0 = 7.3;               %rhombus side length, cm
% h = linspace(0,2*l0,100);

hmax = 11;              %height max, cm
hmin = 6.5;             %height min, cm
l0 = 0.0888;               %rhombus side length, cm
h = linspace(0,2*l0,100);


l = sqrt(-1*(h.^2)+((2*l0)^2));        %cm

figure
plot(l,h)

l0 = sqrt((.235619-.18326)^2+(.27314-.20141)^2)

h_m = 0.27314-0.11859
l_m = 0.287979-0.18326
l_m_calc = sqrt(-1*(h_m^2)+((2*l0)^2))