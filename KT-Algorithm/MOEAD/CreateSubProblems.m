%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA124
% Project Title: Implementation of MOEA/D
% Muti-Objective Evolutionary Algorithm based on Decomposition
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function sp=CreateSubProblems(nObj,nPop,T)

    empty_sp.lambda=[];
    empty_sp.Neighbors=[];

    sp=repmat(empty_sp,nPop,1);%�ظ����鸱����100��lambda���� 100��Neighbors����
    
    %theta=linspace(0,pi/2,nPop);
    
    for i=1:nPop
        lambda=rand(nObj,1);%�������lambda��ֵ
        lambda=lambda/norm(lambda); %norm��lambda���ó�����lambda���������ֵ
        sp(i).lambda=lambda;
        
        %sp(i).lambda=[cos(theta(i))
        %              sin(theta(i))];
    end

    LAMBDA=[sp.lambda]';

    D=pdist2(LAMBDA,LAMBDA);%��LAMBDA֮��ľ���
    %���T���ھ�LAMBDA�ı�ţ�Ȩ��������
    for i=1:nPop
        [~, SO]=sort(D(i,:));%�������������������SO��ס�±�
        sp(i).Neighbors=SO(1:T);%�Ѿ��������T����ݴ洢��sp��
    end

end