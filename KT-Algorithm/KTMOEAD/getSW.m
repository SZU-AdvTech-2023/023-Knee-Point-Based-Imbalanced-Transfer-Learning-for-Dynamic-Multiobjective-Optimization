function [W]=getSW(kneePoints, Hs ,Ht, mu, lambda, dim, kind, p1, p2, p3)
%ITCA����
%Xs��Դ������
%Xt��Ŀ�������ݣ���Xs������ͬ
%mu��ƽ�����ӣ�Խ��Խ����ӳ�������ƶȣ�ԽСԽ����W�ĸ��Ӷ�
%lambda��ƽ�����ӣ�Խ��Խ����ӳ����ͬ�����ݵ������ԣ�ԽСԽ����W�ĸ��Ӷ�
%dim����dimΪ���ڵ���1������ʱ��dimΪ��ά��Ŀ��ά����
%     ��dimΪ����0С��1��С��ʱ����ȡ����������Ӧ������ֵ�ĺ�>=ȫ������ֵ�Ӻ�*dim
%kind���˺���ѡ��:'Gaussian'��'Laplacian'��'Polynomial',����һ�ɷ���-1
%p1,p2,p3���˺�����Ҫ�����Ĳ���
%W���任����n1+n2->dim
%K�����任����
%n1,n2��Դ���ݣ�Ŀ�����ݵ���Ŀ

numOfKnee=size(kneePoints,2);
numOfHs=size(Hs,2);
numOfHt=size(Ht,2);

% 
% s1=size(kneePoints);
% s2=size(Hs);
% n1=s1(2);n2=s2(2);

%%%%%%%%%%% ����K %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X(:,1:numOfKnee)=kneePoints;
X(:,numOfKnee+1:numOfKnee+numOfHs)=Hs;
X(:,numOfKnee+numOfHs+1:numOfKnee+numOfHs+numOfHt)=Ht;
for i=1:numOfKnee+numOfHs+numOfHt
    for j=i:numOfKnee+numOfHs+numOfHt
        K(i,j)=getKernel(X(:,i), X(:,j), kind, p1, p2, p3);
        K(j,i)=K(i,j);
    end
end
    
    
%%%%%%%%%%% ����L %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tempLKH=[];
for i=1:numOfKnee
    tempLKH(1:1, 1:1)=ones(1, 1)/(1*1);
    tempLKH(1+1:1+numOfHs, 1+1:1+numOfHs)=ones(numOfHs, numOfHs)/(numOfHs*numOfHs);
    tempLKH(1:1, 1+1:1+numOfHs)=ones(1, numOfHs)/(-1*numOfHs);
    tempLKH(1+1:1+numOfHs, 1:1)=ones(numOfHs, 1)/(-1*numOfHs);
    LKH{i}=tempLKH;
end
   
    
LHH(1:numOfHs, 1:numOfHs)=ones(numOfHs, numOfHs)/(numOfHs*numOfHs);
LHH(numOfHs+1:numOfHs+numOfHt, numOfHs+1:numOfHs+numOfHt)=ones(numOfHt, numOfHt)/(numOfHt*numOfHt);
LHH(1:numOfHs, numOfHs+1:numOfHs+numOfHt)=ones(numOfHs, numOfHt)/(-numOfHs*numOfHt);
LHH(numOfHs+1:numOfHs+numOfHt, 1:numOfHs)=ones(numOfHt, numOfHs)/(-numOfHs*numOfHt);
    
%%%%%%%%%%% ����H %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H=eye(numOfKnee+numOfHs+numOfHt)-ones(numOfKnee+numOfHs+numOfHt, numOfKnee+numOfHs+numOfHt)/(numOfKnee+numOfHs+numOfHt); 

%%%%%%%%%%% ����W %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% N=(1*eye(n1+n2)+mu*K*L*K+lambda*V1)
TM=zeros(size(mu*K*LKH{1}*K,1),size(mu*K*LKH{1}*K,2));
for i=1:numOfKnee
    TM=TM+mu*K*LKH{i}*K;
end

Temp=(1*eye(numOfKnee+numOfHs+numOfHt)+TM+lambda*K*LHH*K)^(-1)*(K*H*K);
[V,D]=eig(Temp);
D=diag(D);
D=real(D);
[D,I]=sort(D,'descend');

if dim>0 && dim<1
    count=1;
    cur=0;
    s=sum(D);
    while cur/s<dim && D(count)>0
        cur=cur+D(count);
        count=count+1;
    end
else
    count=dim+1;
end

I=I(1:count-1,1);
W=V(:,I');

    
