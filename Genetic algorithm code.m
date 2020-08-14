function error_percent= Fitness_fun(chromosome_x)

load KNN_RFE_1075;    
load PNN_1075_21;
load SVM_RFE_POL_1075;


  a=find(result_PSSM_DWT_SAAC_DPC_RFE_1075_KNN_11==-1)
  KNN_Pre(a)=1;
  b=find(result_PSSM_DWT_SAAC_DPC_RFE_1075_KNN_11==1)
  KNN_Pre(b)=2;
  
  a=find(result_PSSM_DWT_SAAC_DPC_RFE_1075_PNN_21==1)
  PNN_Pre(a)=1;
  b=find(result_PSSM_DWT_SAAC_DPC_RFE_1075_PNN_21==2)
  PNN_Pre(b)=2;
   
  a=find(result_PSSM_DWT_SAAC_DPC_RFE_1075_SVM_Pol==1)
  SVM_Pre(a)=1;
  b=find(result_PSSM_DWT_SAAC_DPC_RFE_1075_SVM_Pol==2)
  SVM_Pre(b)=2;
  
    
   KNN_Predictions=KNN_Pre';
   PNN_Predictions=PNN_Pre';
   SVM_Predictions=SVM_Pre';
    
   
TotalACC=0;
total=[];
Pred1=[];

c11=[];
c22=[];
c33=[];

c1=0;
c2=0;
c3=0;

T_P=0;
T_N=0;
F_N=0;
F_P=0;

Total_sen=0;
Total_spe=0;
TotalMCC=0;
Total_F_Measure=0;
G_mean=0;
G_mean1=0;
Total_G_mean=0;

    total_sequences=size(SVM_Pre,2)
    total_classifiers=3;
    
    total_predictions(1,:)=KNN_Predictions;
    total_predictions(2,:)=PNN_Predictions;
    total_predictions(3,:)=SVM_Predictions;
    
    
    majority_voting_predictions=[];

    for i=1:total_sequences
        count_class1=0;
        count_class2=0;
            

        for j=1:total_classifiers
            if(total_predictions(j,i)==1)
                count_class1=count_class1+chromosome_x(j);
            elseif(total_predictions(j,i)==2)
                 count_class2=count_class2+chromosome_x(j);                        
            end       
        end
     
     all_class_counts=[count_class1 count_class2];
 
         maximum=all_class_counts(1);
        max_index=1;
        for j=2:3
              if(maximum < all_class_counts(j))
                maximum=all_class_counts(j);
                 max_index=j;
            end
        end
        majority_voting_predictions(i)=max_index;
        chromosome_x
    end

  
     Memb_labels(1:525)=1;
    Memb_labels(526:1075)=2;


    %load D1_features_3prop Final_MultiScale_features GPCRs_labels TrainClassesSize;
    
	Total_Seq_train= total_sequences;
    y=majority_voting_predictions;
      Result=find(y==Memb_labels);
       % Result=find(majority_voting_predictions==labels);
       Total_correct=size(Result,2);
       Accuracy_maj=(Total_correct/Total_Seq_train)*100
       
       c1=find(y(1:525)==1);
   Total_correct=size(c1,2);
   C1=(Total_correct/525)*100
   
   c2=find(y(526:1075)==2);
   Total_correct=size(c2,2);
   C2=(Total_correct/550)*100
   


for i=1:1075
    if(i<=525)
        if(y(i)==1)
            T_P=T_P+1;
        else
            F_N=F_N+1;
        end
    elseif(i>526 && i<=1075)
            if(y(i)==2)
                T_N=T_N+1;
            else
                F_P=F_P+1;
                
            end
    end
end
    TNR=T_N/(T_N+F_P);
    TPR=T_P/(T_P+F_N);
    
    G_mean=sqrt(TNR*TPR); 
    Total_G_mean=Total_G_mean +G_mean;
    
    Sensitivity=T_P/(T_P+F_N);
    Total_sen=Total_sen+Sensitivity;
    
    Specificity=T_N/(T_N+F_P);
    
    a=((T_P+F_P)*(T_P+F_N)*(T_N+F_P)*(T_N+F_N));
    b=sqrt(a);
    d=((T_P*T_N)-(F_P*F_N));
    MCC=d/b;
    
    gh=(T_P/(T_P+F_P));
    rh=(T_P/(T_P+F_N));
    F=(gh*rh)/(gh+rh);
    
    F_Measure=2*F;
    Total_F_Measure=Total_F_Measure+F_Measure;
    
    Total_spe=Total_spe+Specificity;
    
    TotalMCC=TotalMCC+MCC;
    
 for i=1:1075
    if(i<=525)
        if(y(i)==1)
            T_N=T_N+1;
        else
            F_P=F_P+1;
        end
    elseif(i>526 && i<=1075)
            if(y(i)==2)
                T_P=T_P+1;
            else
                F_N=F_N+1;
                
            end
    end
 end
    TNR=T_N/(T_N+F_P);
    TPR=T_P/(T_P+F_N);
    
    G_mean=sqrt(TNR*TPR); 
    Total_G_mean=Total_G_mean +G_mean;
    
    Sensitivity=T_P/(T_P+F_N);
    Total_sen=Total_sen+Sensitivity;
    
    Specificity=T_N/(T_N+F_P);
    
    a=((T_P+F_P)*(T_P+F_N)*(T_N+F_P)*(T_N+F_N));
    b=sqrt(a);
    d=((T_P*T_N)-(F_P*F_N));
    MCC=d/b;
    
    gh=(T_P/(T_P+F_P));
    rh=(T_P/(T_P+F_N));
    
    F=(gh*rh)/(gh+rh);
    F_Measure=2*F;
    
    Total_F_Measure=Total_F_Measure+F_Measure;
    
    Total_spe=Total_spe+Specificity;
    
    TotalMCC=TotalMCC+MCC;   

    
    G_mean1=Total_G_mean/2
    sensi=(Total_sen/2)*100
    speci=(Total_spe/2)*100
    MCCT=(TotalMCC/2)
    T_F_Measure=(Total_F_Measure/2)
    
        error_percent=Accuracy_maj*(-1);
end