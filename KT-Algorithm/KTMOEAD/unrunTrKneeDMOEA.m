
function res=TrKneeDMOEAv2(Problem,popSize,MaxIt,T_parameter,group,path)
    source=[];
    partNum=10;
    for T = 1:T_parameter(group,3)/T_parameter(group,2)
        t= 1/T_parameter(group,1)*(T-1);   
        fprintf(' %d',T);
        
        if T==1 || T==2
            [PopX,Pareto,POF_iter,Pareto_iter,runTime] = RMMEDA( Problem,popSize,MaxIt, t);       
            [kneeF,kneeS]=getKneeGroup(Pareto,partNum);
            LastPopX=PopX;
            LastRank=asignRank(PopX,kneeS);
            kneeArray{T}=kneeS;
        else
            tic;
            PopX=[addNoise(LastPopX, floor(popSize*4), size(LastPopX,1)) LastPopX];
%             [PopX,Pareto,POF_iter,Pareto_iter,runTime] = RMMEDA( Problem,popSize,MaxIt, t);       
%             [kneeF,kneeS]=getKneeGroup2(PopX,partNum,t,Problem);
            kneeS=TPM(kneeArray,T);
            rt1=toc;
            popsize=popSize;
            if size(kneeS,2) < popsize
                
                NPopX(:,1:size(kneeS,2)) = kneeS;
                %设计一个增加个体的方案
                newpop = addNoise(kneeS, popsize, size(kneeS,1));
                MPopX(:,size(kneeS,2)+1:popsize) = newpop;
            end
            
%             for i=1:NIni
%                 [PopF(:,i),PopV(:,i)] = Problem.FObj(PopX(:,i)',t);
%             end   
            [PopX,Pareto,POF_iter,Pareto_iter,runTime] = RMMEDA( Problem,popSize,2, t, kneeS);
            estiPareto{T}=Pareto;
            Rank=asignRank([PopX kneeS],kneeS);
%             PopX=[newpop kneeS];
%             Rank=[zeros(1,190), ones(1,10)];
            testPopX=[generateRandomPoints(PopX) generateRandomPoints(PopX)];
            [predictPopX,predTime]=predictPopulationV2(LastPopX',LastRank',PopX',Rank',testPopX');
            initPopulation=predictPopX;
            if size(initPopulation,2)>popSize
                initPopulation=initPopulation(:,1:floor(popSize/1.2));
            elseif size(predictPopX,2)==0
                initPopulation=PopX;
                initPopulation=initPopulation(:,1:floor(popSize/1.2));
            end
%              size(initPopulation)


             [PopX,Pareto,POF_iter,Pareto_iter,runTime] = RMMEDA( Problem,popSize,1, t, initPopulation);
             transPareto{T}=Pareto;
             tic
             
            [PopX,Pareto,POF_iter,Pareto_iter,runTime] = RMMEDA( Problem,popSize,MaxIt, t, PopX);
            rt3=toc;
            tic;
            [kneeF,kneeS]=getKneeGroup(Pareto,partNum);
            rt2=toc;
            LastPopX=PopX;
            LastRank=asignRank(PopX,kneeS); 
            kneeArray{T}=kneeS;
            rt{T}.rt1=rt1;
             rt{T}.rt2=rt2;
              rt{T}.predTime=predTime;
              rt{T}.runTime=rt3;
              rt{T}.all=rt1+rt2+predTime+rt3;
              rt{T}.all;
        end
    
        res{T}.POF_iter=POF_iter;
        res{T}.POS=PopX;
        res{T}.turePOF=getBenchmarkPOF(Problem.Name,group,T,T_parameter );
    end
    path
    save([path ,'estimatePOF.mat'],'estiPareto');
    save([path ,'transferPOF.mat'],'transPareto');
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
