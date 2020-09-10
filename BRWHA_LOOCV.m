clc;clear;
% format long
%% ԭʼ����
Wrr = load('lncRNAsimilarity.txt');   %lncRNA ���������
Wdd = load('diseasesimilarity.txt');   % disease ����������
interaction = load('known_lncRNA_disease_interaction.txt');
interaction_ori = interaction;
n=1;
[nl,~]=size(Wrr);
%% ���ݴ�����
% for gama=1:1:9
    gama=0.5;
    Wrr1=reSimilarity(Wdd,interaction_ori);   % reSimilarity��lncRNA ����������
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
%     km = (pkm+Wrr+Wrr1)/3;
    kd = (pkd+Wdd)/2;
    A_ori=WKNKN( interaction_ori, km, kd, 3, 0.9);
 

%% ������
 [Rt_ori]=BMDA(km,kd,A_ori,gama);   %Rt: predictive association scores between each miRNA-disease pair
 index=find(interaction_ori==1);

%%  ��һ������֤��LOOCV��
for u=1:length(index)
    u 
    interaction(index(u))=0;
   

%% ���ݴ�����

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
    
%     km = (pkm+Wrr+Wrr1)/3;
    kd = (pkd+Wdd)/2;
    A=WKNKN( interaction, km, kd, 3, 0.9);

%% ������
    [Rt]=BMDA(km,kd,A,gama);   %Rt: predictive association scores between each miRNA-disease pair
%%

    Rt_ori(index(u))=Rt(index(u));
    interaction = interaction_ori;
    
end
 pre_label_score = Rt_ori(:);
%    save score_loov_ori_L1_cos; 
    label_y = interaction_ori(:);
    auc=roc_1(pre_label_score,label_y,'red');
  
%  aupr=pr_curve(pre_label_score,interaction_ori,colour);
%  x(n) =gama;
%  y(n) = auc;
%  n = n+1;
% end
% plot(x,y,'red')