clear all
close all
clc

a      =   4.04;
tau    =   16.1/3600;
kappa  =   6.84;
rho_cr =  28.56;
eta    =  32.13;
T      =  10/3600.;
delta  =     0.;
vfree  = 107.11;

L = [ 0.527,  0.494,  0.616,  0.665,  0.635,  0.6  ,  0.688];

Nsegment = 3;
Nhorizon = 2;

%Segment - Time
rho   =  sym('rho',  [Nhorizon,Nsegment]);
v     =  sym('v',    [Nhorizon,Nsegment]);
beta  =  sym('beta', [Nhorizon,Nsegment]);
alpha =  sym('alpha',[Nhorizon,Nsegment]);

C = [];

for k = 1:Nhorizon-1
    for i = 2:Nsegment-1
        theta(k,i) = (1+alpha(k,i)).^a;
        Ve_arg  = -(1/a)*theta(k,i)*(rho(k,i)/rho_cr)^a;
        Ve(k,i)    = vfree*exp(Ve_arg);   
        
        q_im = rho(k,i-1)*v(k,i-1);
        q_i  = rho(k,i  )*v(k,i  );
        
        C = [C;
             rho(k+1,i) - (rho(k,i)   + (T/L(i))*(q_im - q_i));
             v(k+1,i)   -   (v(k,  i) + (T/tau) *(beta(k,i)*Ve(k,i) - v(k,i)) + ...
                                + (T/L(i))*v(k,i)*(v(k,i-1)   - v(k,i)) + ...
                                - beta(k,i)*(1-alpha(k,i))*(eta*T/tau/L(i))*(rho(k,i+1) - rho(k,i))/(rho(k,i) + kappa))];
    end
end

P = [alpha;beta];
dCdP = jacobian(C,P);