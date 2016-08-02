clear
close all

rho = dlmread('rho.dat');
o = dlmread('OUT.out');

figure
hold on
for i=1:size(rho,1)    
    plot(rho(i,:))
end
hold off

