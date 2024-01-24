clc
clear
close all
warning('off')
% CreatTruePOF()
con=configure();
repeatMax=con.repeat;
functions=con.TestFunctions;
T_parameter=con.T_parameter;
popSize=con.popSize;

folderPath = 'result'; % 更改为您想要的文件夹名

if ~isfolder(folderPath)
    mkdir(folderPath);
end

for rep=1:20%repeatMax
    for testFuncNo=1:size(functions,2)
        filename1 = fullfile(folderPath, ['MIGD-DF', num2str(testFuncNo), '-Svm', '.txt']);
        fid1 = fopen(filename1,'a');
        partNum=10;
        Problem=TestFunctions(functions{testFuncNo});
        if Problem.NObj==3
            popSize=150;
        end 
        for group=1%:size(T_parameter,1) 

         MaxIt=T_parameter(group,2);
         fprintf('\n KT-DMOEA dec:%d runing on: %s, configure: %d, partnum: %d ,environment:',con.dec,Problem.Name,group, partNum);
         reskt=TrKneeDMOEA(Problem,popSize,MaxIt,T_parameter,group,partNum);        
         [resIGD,resHV]=computeMetrics(reskt,group,rep,testFuncNo,T_parameter);
         fprintf('\n %.3d',resIGD);
         fprintf(fid1,'%f \n',resIGD);
        
        end %configure
    end%testF
end%rep

function [resIGD,resHV]=computeMetrics(resStruct,group,rep,testFuncNo,T_parameter)
     for T=1:size(resStruct,2)
        POFIter=resStruct{T}.POF_iter;%实验POF
        POFbenchmark=resStruct{T}.turePOF;%采样POF
        for it=1:size(POFIter,2)
            pof=POFIter{it};
            pof(imag(pof)~=0) = abs(pof(imag(pof)~=0));
            igd(it)=IGD(pof',POFbenchmark);
            hv(it)=HV(pof',POFbenchmark);
        end
        IGD_T(T)=igd(end); 
        HV_T(T)=hv(end);
     end
     resIGD=mean(IGD_T);
     resHV=mean(HV_T);
     % filename1 = ['KT-DF',num2str(testFuncNo),'-nt',num2str(T_parameter(group,1)),'-taut',num2str(T_parameter(group,2)),'-IGD', '.txt'];
     % filename3 = ['KT-DF',num2str(testFuncNo),'-nt',num2str(T_parameter(group,1)),'-taut',num2str(T_parameter(group,2)),'-HV', '.txt'];
     %{
     if testFuncNo <=6
         filename1 = ['KT-F',num2str(testFuncNo+4),'-nt',num2str(T_parameter(group,1)),'-taut',num2str(T_parameter(group,2)),'-IGD', '.txt'];
         filename3 = ['KT-F',num2str(testFuncNo+4),'-nt',num2str(T_parameter(group,1)),'-taut',num2str(T_parameter(group,2)),'-HV', '.txt'];
     else
         filename1 = ['KT-FDA',num2str(testFuncNo-6),'-nt',num2str(T_parameter(group,1)),'-taut',num2str(T_parameter(group,2)),'-IGD', '.txt'];
         filename3 = ['KT-FDA',num2str(testFuncNo-6),'-nt',num2str(T_parameter(group,1)),'-taut',num2str(T_parameter(group,2)),'-HV', '.txt'];
     end
     %}
     if group == 1 && rep == 5
         fid1 = fopen(filename1,'w');
        
         fid3 = fopen(filename3,'w');
         for i=1:size(IGD_T,2)
              fprintf(fid1,'%f \n',IGD_T(i)); 
         end
         for i=1:size(HV_T,2)
              fprintf(fid3,'%f \n',HV_T(i)); 
         end
         
     end
     
     %{
     if testFuncNo<=6
        filename1 = ['KT-F',num2str(testFuncNo+4),'-nt',num2str(T_parameter(group,1)),'-taut',num2str(T_parameter(group,2)), '-POF', '.txt'];
     else
        filename1 = ['KT-FDA',num2str(testFuncNo-6),'-nt',num2str(T_parameter(group,1)),'-taut',num2str(T_parameter(group,2)), '-POF', '.txt'];
     end
     %}
     % if group == 1 && rep == 3
     %     for T=1:size(resStruct,2)
     %         filename1 = ['KT-DF',num2str(testFuncNo),'-nt',num2str(T_parameter(group,1)),'-taut',num2str(T_parameter(group,2)),'environment',num2str(T),'-POF', '.txt'];
     %         fid1 = fopen(filename1,'w'); 
     %         POFIter=resStruct{T}.POF_iter;
     %          pof = POFIter{size(POFIter,2)};
     %          pof(imag(pof)~=0) = abs(pof(imag(pof)~=0));
     %          for j=1:size(pof,2)
     %              if size(pof,1) == 2
     %                   fprintf(fid1,'%f \t %f \n',pof(1,j),pof(2,j)); 
     %              else
     %                  fprintf(fid1,'%f \t %f \t %f \n',pof(1,j),pof(2,j),pof(3,j)); 
     %              end
     %          end
     %          fprintf(fid1,'\n'); 
     %     end
     %     fclose(fid1);  
     % end
end
