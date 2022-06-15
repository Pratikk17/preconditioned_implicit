% This script generated fig.5.3 which provides the dependence of 1/phi_0 (see theorem 4.5)
% on gamma. The required data is already stored in "optimization_data_fig_5_3.dat",
% and can be created by running "data_optimization_fig_5_3.m'
%
%

clear all;
close all;

load('optimization_data_fig_5_3.mat')
dt=dT(1)
figure(1); 
plot(gamma(1,:),1./phi0_S(1,:),'r');hold on
xlabel('\gamma'), ylabel('1/\phi_0')
legend('S')

