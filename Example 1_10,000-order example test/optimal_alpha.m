%% alpha optimal choice
%% exa1
clear all ;
clc
k=5;%k>2
ks=0.01;
load('Matrices_A.mat')
AA = A{1,3};
tic
lamda = eig(full(AA));  
save lamda;
alpha_tang = fminimax(@alpha_optimal,0.1, [],[],[],[],10^(-8),[]);
t_ptimal=toc;