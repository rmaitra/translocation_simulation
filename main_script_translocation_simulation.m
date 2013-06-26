%% 1st Order Simulation of DNA translocation
%Raj Maitra
%Notes:
%Dekker, Fast Translocation through solid-state nanopore
%http://pubs.acs.org/doi/full/10.1021/nl048030d
%drag force ~ 0.3 pN
%voltage force ~ 9.4 pN

close all
clc
clear all
figure(1)

% box and nanopore boundary 
boundary = [800  400;
           -400 -400]; %box boundary (so particles don't leave view)
pore_outer_y = [-15 15]; %pore boundary
pore_outer_x = [-50 0];
pore_inner_y = [-10 10];
pore_inner_x = [-30 -20];
pore_dx = 500;
pore1 = [pore_outer_y; %pore1 boundary
         pore_outer_x;
         pore_inner_y;
         pore_inner_x];
pore2 = [pore_outer_y; %pore2 boundary separated by 500 nm
         pore_outer_x+pore_dx;
         pore_inner_y;
         pore_inner_x+pore_dx];
R = diff(pore_inner_y)*10^-9; %nm
d = diff(pore_outer_x)*10^-9; %nm
         
         
%Universal Parameters
dt = 0.00001; %time step  
kb = 1.3806503*10^(-23)*10^12*10^9; %boltzmann constant pN*nm/K
T = 296.15; %kelvin

%Mass
m = 5; %pg

%DNA Strand Length
a = 0.34*10^-9; %nm
Nbp = 1000; %number of base pairs
R_l = 10*10^-9; %spring (bond) rest length
N=(Nbp*a)/R_l;   %number of Kuhn sections
L=(Nbp*a)*10^9;    %maximum DNA Length (nm)

%Spring Force
k = 1000; %pN/nm spring constant of bonds

%Voltage Force
V = 2000*10^-3; %mV
ld = R_l/(2*a); %linear density
e = 1.6*10^-19; %coulombs
f_v = ld*(e*V/d)*10^12; %voltage force pN: 4*pi*(10^-3 Pa s)*(10*10^-9 m)/(ln((10^-9 m)/(2*10^-9 m))+0.84) to pN*us/nm

%Drag Force
r = 1*10^-9; %nm
eta = 10^-3; %viscosity Pa s
gamma = (4*pi*eta*R_l)/(log(R_l/(2*r))+0.84)*10^9; %pN*us/nm

%Diffusion
D = kb*T/(gamma); %Diffusion constant nm^2/us

%Brownian Force 
omega = sqrt(2*D*dt); %brownian motion

%Store Parameters
parameters = [f_v gamma omega k R_l*10^9 m];

%Initial Settings
position = initialize_dna(N,R_l*10^9)';
velocity = zeros(size(position));
t = 0;
j = 0;
order = 1;
animate_translocate(position,boundary,pore1,pore2)
pause

%Main Loop
while true%position(end,1) < 0
    [position, velocity] = numerical_solve(position, velocity, dt, boundary, pore1, pore2, parameters,order);
    t = t+dt; %total time elapsed
    j = j+1;
    
    if ~rem(j,100)
        animate_translocate(position,boundary,pore1,pore2)%plot position at time t
        pause%pause for a brief moment
    end
end

display(['A ' num2str(Nbp) ' bp strand of DNA translocated through the pore in ' num2str(t) ' us.'])



