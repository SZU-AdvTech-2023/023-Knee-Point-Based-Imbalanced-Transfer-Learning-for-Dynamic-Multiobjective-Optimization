
%parameters  test_suit_name Population_size Iteration_count T_parameter the_group_generation
function res=TrKneeDMOEA(Problem,popSize,MaxIt,T_parameter,group,partNum)
    source=[];
    %partNum=10;%分成十个子区域
    for T = 1:T_parameter(group,3)/T_parameter(group,2) %20
        t= 1/T_parameter(group,1)*(T-1);   
        fprintf(' %d',T);
        
        if T==1 || T==2
            %个体决策变量（或位置），最终 Pareto 前沿，Pareto 前沿数组（包含了每一次）
            [PopX,Pareto,POF_iter]=moead( Problem,popSize,MaxIt, t);    
            %从pareto前沿获取区域拐点
            [kneeF,kneeS]=getKneeGroup(Pareto,partNum);
            LastPopX=PopX;
            LastRank=asignRank(PopX,kneeS);
            kneeArray{T}=kneeS;
        else
   
            kneeS=TPM(kneeArray,T);
 
            [PopX,Pareto,POF_iter]=moead( Problem,popSize,1, t, kneeS);
            Rank=asignRank([PopX kneeS],kneeS);
            
            PopX=[PopX kneeS];
            testPopX=[generateRandomPoints(PopX,Problem) generateRandomPoints(PopX,Problem)];
            [predictPopX,predTime]=predictPopulationKnee(LastPopX',LastRank',PopX',Rank',testPopX',partNum);
            initPopulation=predictPopX;
            if size(initPopulation,2)>popSize
                initPopulation=initPopulation(:,1:floor(popSize/1.2));
            elseif size(predictPopX,2)==0
                initPopulation=PopX;
                if size(initPopulation,2)>=floor(popSize/1.2)
                    initPopulation=initPopulation(:,1:floor(popSize/1.2));
                else
                    initPopulation=initPopulation(:,end);
                end
            end
%              size(initPopulation)
            %initPopulation = zeros(10,50);
            %[PopX,Pareto,POF_iter]=moead( Problem,popSize,1, t, initPopulation);     
            [PopX,Pareto,POF_iter]=moead( Problem,popSize,MaxIt, t, initPopulation); 
            [kneeF,kneeS]=getKneeGroup(Pareto,partNum);
            LastPopX=PopX;
            LastRank=asignRank(PopX,kneeS); 
            kneeArray{T}=kneeS;
        end
    
        res{T}.POF_iter=POF_iter;%
        res{T}.POS=PopX;%
        res{T}.turePOF=getBenchmarkPOF(Problem.Name,group,T,T_parameter );%获取之前采样得到pof
    end

end


function [kneeF,kneeS]=getKneeGroup2(PopX,partNum,t,Problem)
for i=1:size(PopX,2)
    [PopF(:,i),~] = Problem.FObj(PopX(:,i)',t);
end 
    [boundaryS,boundaryF]=getBoundary(PopX,PopF);
    [posArr,pofArr]=partition(PopX,PopF,partNum,boundaryF);
    for partNo=1:partNum
        [kneeS,kneeF]=getKnees(posArr{partNo},pofArr{partNo});
        kneeSArr{partNo}=kneeS;
        kneeFArr{partNo}=kneeF;
    end
    kneeS=cell2mat(kneeSArr);
    kneeF=cell2mat(kneeFArr);
end


function [kneeF,kneeS]=getKneeGroup(Pareto,partNum)
    [boundaryS,boundaryF]=getBoundary(Pareto.X,Pareto.F);
    [posArr,pofArr]=partition(Pareto.X,Pareto.F,partNum,boundaryF);
    for partNo=1:partNum
        [kneeS,kneeF]=getKnees(posArr{partNo},pofArr{partNo});
        kneeSArr{partNo}=kneeS;
        kneeFArr{partNo}=kneeF;
    end
    kneeS=cell2mat(kneeSArr);
    kneeF=cell2mat(kneeFArr);
end


function Rank=asignRank(PopX,KneeX)
for i=1:size(PopX,2)
    for j=1:size(KneeX,2)
        if isequal(PopX(:,i),KneeX(:,j))==1
            Rank(i)=1;
            break;
        else
            Rank(i)=-1;
        end
    end
end
end
