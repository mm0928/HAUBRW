function [M]=normFun1(M) 
for i=1:length(M)         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   z=find(M(:,i)~=0);
   M(z,i)=M(z,i)./sum(M(z,i));
end  
% n = size(M);
% L=sum(M);
% for i=1:n
%     for j=1:n
%         M(i,j)=M(i,j)/(((L(i)*L(j))^0.5));
%     end
% end


