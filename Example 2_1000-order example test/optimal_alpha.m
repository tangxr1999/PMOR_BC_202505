%% alpha optimal choice
%% exa1
clear all ;
clc
k=5;%k>2
ks=0.01;
load('A_bl_09.mat')
tic
lamda = eig(full(A_bl_09));  
save lamda;
alpha_tang = fminimax(@alpha_optimal,0.1, [],[],[],[],10^(-8),[]);
t_ptimal=toc;