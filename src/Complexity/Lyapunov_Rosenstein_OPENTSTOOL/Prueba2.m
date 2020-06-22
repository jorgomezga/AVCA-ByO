% clear all
% clc

function Prueba2
close all
clc
clear all
if isunix
    addpath('/media/Data/Codigos/Toolboxes/x64/ts2/dynsystems/')
    addpath(genpath('/media/Data/Codigos/Toolboxes/NDA/'))
    addpath('/media/Data/Codigos/Lyapunov_Rosenstein/Dynamical/')
    addpath(genpath('/media/Data/Codigos/Toolboxes/OpenTSTOOL/Linux'))
else
    addpath('L:\Codigos\Toolboxes\x64\ts2\dynsystems\')
    addpath(genpath('L:\Codigos\Toolboxes\NDA'))    
%     addpath(genpath('L:\Codigos\Toolboxes\OpenTSTOOL\Windows'))
end

%Henon
% a = 1.4; b = 0.3;
% x(1,1) = 0; x(1,2) = 0;
% for i = 2:501
%     x(i,1) = 1 - a * x(i-1,1)^2 + b * x(i-1,2);
%     x(i,2) = x(i-1,1);
% end
% x(1,:) = [];
% embedded=x;
% embedded=logistic(1:2000,4); embedded=embedded';
% embedded=rosslerpoints(2000,0.1);
embedded=lorenzpoints(5000,0.01);
% [~, embedded]=FOLorenz([16 45.92 4],[0.993 0.993 0.993],500,[0.1 0.1 0.1],0.1);
% embedded=logistic(1:500,4); embedded=embedded';

x=embedded(:,1);
% mean_period = time_series_mean_period(x);
% x=(x-min(x))./(max(x)-min(x));
[iDim,iTao]=embebimiento(x);
% iDim=3;iTao=1;
embedded=embeb(x,iDim,iTao);



Maxdelta_n=3;          %Maximo de iteracion, para el calculo de delta_n
iTau=iTao;        %Exclusion para evitar correlacion temporal
nvecinos=5;
% Lyap=LyapRos1(embedded,iTau,epsilon,Maxdelta_n);
tic
Lyap=LyapRos2(embedded,iTau,Maxdelta_n,nvecinos);
toc

plot(Lyap(:,1),Lyap(:,2))
Lyapa=Lyap(1:3,:);
coef=polyfit(Lyapa(:,1),Lyapa(:,2),1);
LLEmio=coef(1)

%Metodo -1
save temp.dat embedded -ascii -tabs
if ispc
tiseanPath='L:\Codigos\Lyapunov_Rosenstein\';
else
%tiseanPath='/media/Data/Codigos/Lyapunov_Rosenstein/';
tiseanPath='/home/jorge/Documentos/';
end
tic
system([tiseanPath,'lyap_r -V0 ','-s',num2str(Maxdelta_n), ' -t', num2str(iTao),' -m',num2str(iDim),...
    ' -d',num2str(iTao),' -o lyap.dat temp.dat']);
toc
l = load('lyap.dat');
hold on
plot(l(:,1),l(:,2),'r')
l=l(1:3,:);
LLETISEAN=polyfit(l(:,1),l(:,2),1)

%Metodo 3
%Prueba tstool
e=embedded';
%x = largelyap(pointset, query_indices, taumax, maximal_neighbours, k_exclude)
tic
l=largelyap(e', 1:length(e')-Maxdelta_n, Maxdelta_n, nvecinos, iTao); 
toc
plot(0:length(l)-1,l','k');
ltstool=polyfit(1:3,l(1:3)',1)


%Metodo 4
tic
lR=largelyap_R(e', 1:length(e')-Maxdelta_n, Maxdelta_n, nvecinos, iTao); 
toc
plot(0:length(lR)-1,lR','m');
ltstoolR=polyfit(1:3,lR(1:3)',1)
clc