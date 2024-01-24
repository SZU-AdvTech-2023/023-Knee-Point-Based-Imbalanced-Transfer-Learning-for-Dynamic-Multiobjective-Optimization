function pop = TPM(kneeArray,T)
    vec=kneeArray{T-1}-kneeArray{T-2};
    step=vec+normrnd(0,0.1);
    for j=1:size(vec,2)
        for i=1:size(vec,1)-1
            kjl=kneeArray{T-1}(i+1:end,j);
            fenzi=sqrt(sum(kjl.^2));
            pietheta(i)=atan(fenzi/kneeArray{T-1}(i,j));
        end
    end
    num=1;
    while num<11
        theta=ones(1,size(vec,1)-1).*(randi([-1,1],1,1)/10990);
        fi= pietheta+ theta;
        for i=1:size(vec,1)
            if i==1
                u(i)=vec(i)*cos(fi(i));
            elseif i<size(vec,1)
                temp=1;
                for j=1:size(vec,1)-2
                    temp=temp*sin(fi(j));
                end
                u(i)=vec(i)*temp*cos(fi(i));
            else
                temp=1;
                for j=1:size(vec,1)-1
                    temp=temp*sin(fi(j));
                end
                u(i)=vec(i)*temp;
            end
        end
       
        samplePop(:,num)=u;
         num=num+1;
    end
     pop= mod(abs(kneeArray{T-1}+samplePop(:,10)),1);

pop1=mod(abs(kneeArray{T-1}+step.*(kneeArray{T-1}-kneeArray{T-2})),1);
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