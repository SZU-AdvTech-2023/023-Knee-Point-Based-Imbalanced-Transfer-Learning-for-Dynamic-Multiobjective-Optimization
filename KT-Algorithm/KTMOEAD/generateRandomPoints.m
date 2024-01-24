function resPop=generateRandomPoints(Pop,Problem)
nVar=size(Problem.XLow,1);             % Number of Decision Variables
VarSize=[nVar 1];   % Decision Variables Matrix Size
VarMin = Problem.XLow;         % Decision Variables Lower Bound
VarMax = Problem.XUpp;         % Decision Variables Upper Bound

resPop=Pop;
for j=1:size(Pop,2)
   
    mutIndex=randperm(ceil(size(Pop,1)*0.5));
    %{
    for i=1:mutIndex
        temp=cauchyrnd(Pop(i,j),1);
        if temp <0 || temp>1
           temp =0;
        end
        resPop(i,count)=temp;
        count=count+1;
    end
    %}
    for i=1:mutIndex
        randPop(:,i)=unifrnd(VarMin,VarMax,VarSize);
    end
    resPop = [resPop randPop];
end
