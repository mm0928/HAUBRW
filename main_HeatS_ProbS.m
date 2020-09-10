function [KK]=main_HeatS_ProbS(trainset9_ratings_matrix, lambda)
trainset9_ratings_matrix=trainset9_ratings_matrix';
m=size(trainset9_ratings_matrix,1); 
n=size(trainset9_ratings_matrix,2);  
% User_uncollected_movies = n-sum(trainset9_ratings_matrix,2); 
Ko = sum(trainset9_ratings_matrix, 1);
Ku = sum(trainset9_ratings_matrix, 2); 
KK=zeros(n,m);

% lambda=0.6;%[0:0.1:1];
beta = 0;%[-2:0.1:1];
num_lambda=size(lambda,2);
num_beta = size(beta,2);
L=20;%[10 50 100];  
num_L = size(L,2);
for h0=1:num_lambda
    W = HeatS_ProbS_recommendations(trainset9_ratings_matrix,lambda(h0)); 
    for h = 1:num_beta
        User_recommender_matrix = zeros(m,n); 
        User_recommender_length = zeros(m,1);
        for k = 1:n
            if Ko(1,k) ~= 0
                Initial_Con(1,k) = Ko(1,k)^beta(h);
            else if Ko(1,k) == 0
                    Initial_Con(1,k) = 0;
                end
                    
            end
        end

       
        for i = 1:m
            KK(:,i) = W*(trainset9_ratings_matrix(i,:).*Initial_Con)'; 
        end
    end
end
end


