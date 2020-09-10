clc;clear;
format long
%% 原始数据
Wrr = load('lncRNAsimilarity.txt');
Wdd = load('diseasesimilarity.txt');
interaction = load('known_lncRNA_disease_interaction.txt');
interaction_ori = interaction;
[nl,~]=size(Wrr);
n=1;
%% 数据处理方法
% for gama=1:1:9
    gama=1;
    Wrr1=reSimilarity(Wdd,interaction_ori);   % reSimilarity：lncRNA 功能相似性
    pkm = cosSim(interaction_ori);
    pkd = cosSim(interaction_ori');
    
   for i=1:nl
        for j=1:nl
           if Wrr(i,j)~=0 
                we(i,j)=1; 
            else
                we(i,j)=0; 
           end
            if Wrr1(i,j)~=0
                wf(i,j)=1;
            else
                wf(i,j)=0;
            end
                km(i,j)=(we(i,j)*Wrr(i,j)+wf(i,j)*Wrr1(i,j)+pkm(i,j))/(we(i,j)+wf(i,j)+1);
        end
   end
    kd = (pkd+Wdd)/2;
    A_ori=WKNKN( interaction_ori, km, kd, 3, 0.9);
 

 
%% 主方法
 [Rt_ori]=BMDA1(km,kd,A_ori,gama);   %Rt: predictive association scores between each miRNA-disease pair
 %%
 Rt_ori_ori=Rt_ori;
 index=find(interaction_ori==1);


%%  5-fold 交叉验证

auc = zeros(1,100);
for k = 1:100
    k
    indices = crossvalind('Kfold', length(index), 5);
    interaction = interaction_ori;
    Rt_ori=Rt_ori_ori;
for cv = 1:5
       cv;
       index_2 = find(cv == indices);
       %%%移除已知关系
       interaction(index(index_2)) = 0;

       %%%计算得分矩阵
%%  数据处理

   Wrr1=reSimilarity(Wdd,interaction);
    pkm = cosSim(interaction);
    pkd = cosSim(interaction');
   
    for i=1:nl
        for j=1:nl
           if Wrr(i,j)~=0 
                we(i,j)=1; 
            else
                we(i,j)=0; 
           end
            if Wrr1(i,j)~=0
                wf(i,j)=1;
            else
                wf(i,j)=0;
            end
                km(i,j)=(we(i,j)*Wrr(i,j)+wf(i,j)*Wrr1(i,j)+pkm(i,j))/(we(i,j)+wf(i,j)+1);
        end
   end
    
    kd = (pkd+Wdd)/2;
    A=WKNKN( interaction, km, kd,3, 0.9);
%% 主方法
    [Rt]=BMDA1(km,kd,A,gama);  
    
%%

      Rt_ori(index(index_2)) = Rt(index(index_2));
       interaction = interaction_ori;
end

%% 画auc曲线
    pre_label_score = Rt_ori(:);
    label_y = interaction_ori(:);
    auc(k) = roc_1(pre_label_score,label_y,'red');
    
end
%% 
auc_ave = mean(auc);
auc_std = std(auc);
% x(n)=gama;    
% y_ave(n)=auc_ave;
% y_std(n)=auc_std;
% n=n+1;
% end
% save auc_ave auc_ave;
% save auc_std auc_std;