clc;
close all;
clear all;

%load PSSM_DWT_Features_P186;

data=xlsread('PSSM_DPC.csv');

S=dct2(data);%matlab function
    U=S(1:10,1:10);
	feature_PSSM_DCT_P(i,:)=U(:);
    PSSM_DPC_DCT=feature_PSSM_DCT_P;
    