clc;
close all;
clear all;
addpath libsvm-mat-2.88-1;
addpath libsvm-3.20\matlab;

load PSSM_358;
%****************************************
  
Result=0; 
y=[];
Accuracy=0;
TestLabel=[];
Total_Seq_train=358;
DNA_labels=[];
Total_correct=0;
  c1=0;  c2=0; 

DNA_labels(1:179)=1;
DNA_labels(180:358)=2;

  %>>>>>>>>>>
  fold=10;
  single_fold = floor(Total_Seq_train/fold);
test_total = single_fold * fold;
remain = Total_Seq_train - test_total;

aaa=PSSM_358';
for k=1:10
  randperm_data=randperm(size(PSSM_358',2));
 ind=randperm_data;
  feature1=aaa(:,ind);
    permut_labels=DNA_labels(ind);
      aa=find(permut_labels(1:358)==1);
      bb=find(permut_labels(1:358)==2);
  %+++++++++++++++++++++++++++++  train label
  Yout=[];
 Labelstem=[];
 Samplestem=[];
 Samplestem= feature1';
Labelstem= permut_labels';
m = single_fold;
l=1;
A = 1;
        C = 0;

        for T = 1:fold
            C = C + 1;
               T

                if T == 1

                Samples=Samplestem(A + single_fold:end,:)';
                TestSample=Samplestem(1:single_fold,:)';
                Labels=Labelstem(A + single_fold:end,:)';

                TestLabel=Labelstem(1: single_fold,:)';
                A = single_fold;

                else
                    if C == fold
                        s11=Samplestem(1:A,: ); % Jackknifing 
                        l11=Labelstem(1:A,: );

                        Samples=s11';
                        Labels=l11';

                        TestSample=Samplestem(A + 1: end,:)';
                        TestLabel=Labelstem(A + 1: end,:)'; 
                    else
                        s11=Samplestem(1:A,: ); % Jackknifing 
                        l11=Labelstem(1:A,: );
                        A = single_fold;
                        A = T * A;
                        s22=Samplestem((A+1):end,:);
                        l22=Labelstem((A+1):end,:);

                        Samples=[s11;s22]';
                        Labels=[l11;l22]';

                        TestSample=Samplestem((A - single_fold)+ 1: A,:)';

                        TestLabel=Labelstem((A - single_fold)+ 1: A,:)';
                    end  
        end
    
    model = svmtrain(Labels', Samples',' t 0  -c 18.8576 -g 0.21215');
    
   [Predict_label,accuracy, dec_values] = svmpredict(TestLabel', TestSample', model);
 %dec_values=dec_values(A)+A;
    if C == fold
      Yout(l:m + remain) = Predict_label;

    else
        Yout(l:m) = Predict_label;
        l = m + 1;
        m =(T+1)* single_fold;

    end
%PSSM_358_dec_values_10fold_lnr=dec_values;
%save PSSM_358_dec_values_10fold_lnr;
    end
Yout=round(Yout);

yy(1:179)=Yout(aa);
yy(180:358)=Yout(bb);
%yy(193:400)=Yout(cc);
%yy(401:523)=Yout(dd);

Result=find(yy==DNA_labels);
   Total_correct=size(Result,2);
   Accuracy=(Total_correct/Total_Seq_train)*100

if (result< Accuracy)
    result=Accuracy;
    save prediction yy;
end
 end
% per_PSSM_358_dec_10fold_lnr=result;
% save per_PSSM_358_dec_10fold_lnr 

 