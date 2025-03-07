%% IWP Full Scale C.G calculation
clc; close all; clear all

syms xg

%% Inputs

L = 11.10;
R = 1.16;
Mtail = 4.21;

% Mass Items of interest (kg)
m_miscsystems = .95;
m_payload = 5.5; 
m_iwp_batteries = .57;
m_main_batteries = 6.36;
m_ballonets = 3.65; 

m_gondola = m_miscsystems+m_payload+m_iwp_batteries+m_main_batteries
m_solarpanels = 1.51
m_tail = .97
m_iwp = .36
m_motor = .333
m_front_ballonet = m_ballonets/2
m_aft_ballonet = m_ballonets/2
m_point_sum = m_gondola+m_solarpanels+m_tail+m_iwp+m_motor+m_front_ballonet+m_aft_ballonet


% x location 

xCB = L/2;

xsolarpanels = L/2;
xtail = (L/2)+Mtail;
xiwp = (L/2)+Mtail;
xmotor = L;
xfrontballonet = L/4;
xaftballonet = 3*L/4;

% y location (where C.B at 0)

yCB = 0;

ysolarpanels = R;
ytail = 0;
yiwp = 0;
ymotor = 0;
yfrontballonet = -1;
yaftballonet = -1;

% C.G equation

xCG = xCB;
xCG_eqn = ((m_gondola*xg)+(m_solarpanels*xsolarpanels)+(m_tail*xtail)+(m_iwp*xiwp)+(m_motor*xmotor)+(m_front_ballonet*xfrontballonet)+(m_aft_ballonet*xaftballonet))/(m_point_sum) == xCB;


xg = vpasolve(xCG_eqn,xg);

yg = -sqrt((R^2)*(1-(2*xg/L)^2));

yCG = ((m_gondola*yg)+(m_solarpanels*ysolarpanels)+(m_tail*ytail)+(m_iwp*yiwp)+(m_motor*ymotor)+(m_front_ballonet*yfrontballonet)+(m_aft_ballonet*yaftballonet))/(m_point_sum);

Lg = vpa([xg, yg],3)

CG = vpa([xCG, yCG], 3)
CB = vpa([xCB, yCB], 3)