clc;
clear all;
close all;

newDataPSSM=[];  

data=fastaread('Test dataset_148.txt');
[m n]=size(data);

db='D:\NCBI\blast-2.5.0+\db\Swissprot';

    for i=1:148
        i
       arr=[];
       seq=data(i).Sequence;
       newDataPSSM=blastpssm_seq(seq,db);
      
       csvwrite(['pssm' num2str(i) '.xls'] ,newDataPSSM);  
       
       %dlmwrite([num2str(i) '.pssm'], newDataPSSM)
    end
    
%SDNA_PSSM_files_937_107=newDataPSSM;
 
%save PDNA_543_PSSM_files;

