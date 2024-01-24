function [model, beta ] = IKTrAdaBoostTrain(tdX,tdY,tsX,tsY,partNum)
    tX = [tdX ; tsX];
    tY = [tdY ; tsY];
    n = size(tdY,1);
    m = size(tsY,1);
    cnk=1;
    ck=(n-10)/10;
    T = 6;
    w = ones(m+n,1);
    model = cell(2,T); % Cell array to hold both SVM and Decision Tree models
    beta = zeros(1,T);
    for i = 1:n
        if tdY(i)==1
            w(i,:) = w(i)/partNum;
        else
            w(i,:) = w(i)/(n-partNum);
        end
    end
    w(n+1:m+n)=w(n+1:m+n)*1/m;

    for t = 1:T
        % SVM model training
        model{1,t} = svmtrain(w,tY,tX,'-t 0 -q -h 0 -s 2');
        predict_svm = svmpredict(tY,tX,model{1,t},'-q');
        
        % Decision Tree model training
        model{2,t} = fitctree(tX, tY, 'Weights', w);
        predict_tree = predict(model{2,t}, tX);
        
        % Combine predictions from SVM and Decision Tree
        combined_prediction = sign(predict_svm + predict_tree);
        
        if length(combined_prediction) == 0
            break;
        end
        
        Ink=0;
        Ik=0;
        for checkp=1:n
            if checkp > length(combined_prediction) || checkp > length(tdY)
                break;
            end

            if combined_prediction(checkp) ~= tdY(checkp) && tdY(checkp) == 1
                Ik = Ik + 1;
            end
            if combined_prediction(checkp) ~= tdY(checkp) && tdY(checkp) == -1
                Ink = Ink + 1;
            end
        end
        
        sigama=(ck*Ik+0.001)/(cnk*Ink+0.001);
        sW = sum(w(n+1:m+n));
        et = sum(w(n+1:m+n).*(combined_prediction(n+1:m+n)~=tsY)/sW);
        
        if et >= 0.5
            et = 0.499;
        elseif et == 0
            et = 0.001;
        end
        
        bT = et/(1-et);
        beta(t) = bT;
        b = 1/(1+sqrt(2*log(n/T)));
        wUpdate = [(b*ones(n,1)).*sigama.^(combined_prediction(1:n)~=tdY) ; (bT*ones(m,1)).^(-(combined_prediction(n+1:m+n)~=tsY)) ];
        w = w.*wUpdate;
    end
end
