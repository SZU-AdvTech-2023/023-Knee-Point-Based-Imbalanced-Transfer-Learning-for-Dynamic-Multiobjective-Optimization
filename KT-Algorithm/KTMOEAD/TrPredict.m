function Ydash = TrPredict(X, svmmodels, beta)
    % X: features of test data
    N = length(svmmodels);
    start = ceil(N/2);
    l = size(X,1);
    yOne = ones(l,1);
    yTwo = ones(l,1);
    Ydash = ones(l,1);
    for i = start:N
        % SVM model prediction
        predict_svm = svmpredict(yOne, X, svmmodels{1,i}, '-q');

        % Decision Tree model prediction
        predict_tree = predict(svmmodels{2,i}, X);

        % Combine predictions from SVM and Decision Tree
        combined_prediction = sign(predict_svm + predict_tree);

        % Update yOne based on the combined prediction
        yOne = yOne .* ((beta(i) * ones(l, 1)).^(-combined_prediction));

        % Update yTwo (if needed)
        yTwo = yTwo .* ((beta(i) * ones(l, 1)).^(-0.5));
        % predict = svmpredict(yOne,X,svmmodels{i},'-q');
        % %predict = predict == 1;
        % yOne = yOne.*((beta(i)*ones(l,1)).^(-predict));
        % yTwo = yTwo.*((beta(i)*ones(l,1)).^(-0.5));
    end
    Ydash(yOne < yTwo) = -1;
end


