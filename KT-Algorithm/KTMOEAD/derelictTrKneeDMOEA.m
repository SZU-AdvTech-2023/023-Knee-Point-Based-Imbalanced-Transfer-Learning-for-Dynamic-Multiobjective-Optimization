
function res=TrKneeDMOEA(Problem,popSize,MaxIt,T_parameter,group)
    source=[];
    partNum=15;
    for T = 1:T_parameter(group,3)/T_parameter(group,2)
        t= 1/T_parameter(group,1)*(T-1);   
        fprintf(' %d',T);
        if T==1
            [PopX,Pareto,POF_iter,Pareto_iter,runTime] = RMMEDA( Problem,popSize,MaxIt, t);       
            [kneeF,kneeS]=getKneeGroup(Pareto,partNum);
            source=[source struct('popX',kneeS','popF',kneeF')];
            [w,individuals]=buildModel(source);
        else
            initPop=predictKnee(w,individuals,popSize);
            [PopX,Pareto,POF_iter,Pareto_iter,runTime] = RMMEDA( Problem,popSize,MaxIt, t, initPop');
            [kneeF,kneeS]=getKneeGroup(Pareto,partNum);
            target=[struct('popX',kneeS','popF',kneeF')];
            source=[source struct('popX',kneeS','popF',kneeF')];
            [w,individuals]=updateModel(source,target);
        end
        
        
        
        res{T}.POF_iter=POF_iter;
        res{T}.POS=PopX;
        res{T}.turePOF=getBenchmarkPOF(Problem.Name,group,T,T_parameter );
    end
end


function res=disnormF(Problem,Fs,Fa,W,y, kind, p1, p2, p3,p,t)
    y=Problem.FObj(y,t);
    res=sum((getNewY(Fs, Fa, y(2), W, kind, p1, p2, p3) - p).^2);
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
