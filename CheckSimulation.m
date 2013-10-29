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

Nsegment = 8;
Ntime = 850;

rho =  20*ones(Ntime,Nsegment);
v   = 100*ones(Ntime,Nsegment);
Ve  = 100*ones(Ntime,Nsegment);
beta  =  1*ones(Ntime,Nsegment);
alpha =  0*ones(Ntime,Nsegment);



%Accident
alpha([361:841],5) = 0.1;
beta([361:841],5) = 0.9;

       

for k = 1:Ntime-1
    for i = 2:Nsegment-1
        
        theta(k,i) = (1+alpha(k,i)).^a;
        Ve_arg  = -(1/a)*theta(k,i)*(rho(k,i)/rho_cr)^a;
        Ve(k,i)    = vfree*exp(Ve_arg);   
        
        q_im = rho(k,i-1)*v(k,i-1);
        q_i  = rho(k,i  )*v(k,i  );
        
        rho(k+1,i) = rho(k,i)   + (T/L(i))*(q_im - q_i);
        v(k+1,i)   =   v(k,  i) + (T/tau) *(beta(k,i)*Ve(k,i) - v(k,i)) + ...
                                + (T/L(i))*v(k,i)*(v(k,i-1)   - v(k,i)) + ...
                                - beta(k,i)*(1-alpha(k,i))*(eta*T/tau/L(i))*(rho(k,i+1) - rho(k,i))/(rho(k,i) + kappa);
                            
                                   
                  
    end
end

for i = 1:Nsegment
    figure(i)
    subplot(2,2,1)
    plot(v(:,i))
    subplot(2,2,2)
    plot(rho(:,i))
    subplot(2,2,3)
    plot(Ve(:,i))
    subplot(2,2,4)
    plot(alpha(:,i))
    hold on
    plot(beta(:,i))
    
end

save MatlabSim
% 
% k = 3;i = 2;
% 
% alpha_i = alpha(k,i)
% rho_i   = rho(k,i)
% 
% theta_i = (1 + alpha_i)^a
% Ve_arg  = -(1/a)*theta_i*(rho_i/rho_cr)^a
% Ve      = vfree*exp(Ve_arg)



% for k = 1:Ntime^
%     display('time')
%     k-1
%     for i = 2:Nsegment-1
%         display('rho')
%         rho(k,i)
%         display('v')
%         v(k,i)
%         display('Ve')
%         Ve(k,i)
%         display('theta')
%         theta(k,i)
%     end
%     pause
% end