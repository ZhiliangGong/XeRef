function [a,flag] = fitGaussianLikelihood(x,y)

    flag = 0; %fitting good by default
    
    myfun = @(a) (a(1)*exp(-(x-a(2)).^2/2/a(3)^2) - y);
    
    a0 = [0,0,0];
    a0(1) = max(y);
    ind = find(y==max(y),1);
    a0(2) = x(ind);
    if ind == length(y) || ind == 1
        flag = 1;
        a0(3) = range(x)/length(x);
    else
        leftSide = abs(y(1:ind-1) - a0(1)/2);
        ind1 = find(leftSide == min(leftSide),1,'last');
        rightSide = abs(y(ind+1:end) - a0(1)/2);
        ind2 = find(rightSide == min(rightSide),1) + ind;
        a0(3) = (x(ind2) - x(ind1))/2;
    end
    lb = [0,-Inf,0];
    ub = [Inf,Inf,Inf];
    
    options = optimoptions('lsqnonlin','MaxFunEvals',1e25,'MaxIter',1e5,'Display','off');
    a = lsqnonlin(myfun,a0,lb,ub,options);
    a = reshape(a,length(a),1);

end