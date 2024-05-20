% Endsem CH5440 
% Ojas Phadake - CH22B007

clc;
clear all;
close all;

load steamdata.mat

%% A 
nsamples = size(Ftrue, 1);
nvar = size(Ftrue, 2);

Z = Ftrue;
Zs = Z - mean(Z);
Sz = Zs'*Zs/nsamples;

[u,s,v] = svd(Zs/sqrt(nsamples), "econ");
D = s.*s;

fprintf("Using the values which I obtain in D, there is a noticeable difference between the values, as " + ...
    "it goes from 10^-28 to 10^0 which is a very big difference. \nHence, as there are 8 eigenvalues which " + ...
    "are equal, there are 8 constraints.\n")
fprintf("***********************************************************************************************")
%% B

dest = 8;

nfact = nvar - dest;
for i = nfact+1:nvar
    Amat(i-nfact,:) = v(:, i)';% This matrix contains values of constraints which are in total 8
end

% Say the 1st 17 are selected to be independent and the last 8 as dependent
% variables

Aesti = Amat(:, 8:10);
Aesti = [Aesti Amat(:, 12:25)];
Aestd = Amat(:, 1:7);
Aestd = [Aestd Amat(:, 11)];

cond(Aestd)

Restb = -inv(Aestd)*Aesti;

fprintf("\nWe have selected a good set of 17 variables as independent and the other 8 as dependent.\n")

%% C

% Finding out the best set of independent variables so that the condition
% number is minimized

[u1, s1, v1] = svd((Fmeas - mean(Fmeas))/sqrt(nsamples), "econ");
D1= s1.*s1;

for i = nfact+1:nvar
    Amat1(i-nfact,:) = v1(:, i)';% This matrix contains values of constraints which are in total 8
end

index = 1:1:25;
set_independent = nchoosek(index, 8);
cond_num = 10^5; % Intialized as a high value

for i=1:size(set_independent, 1)
    Aestd1 = Amat1(:, set_independent(i, :));
    if cond(Aestd) < cond_num
        min = i;
        cond_num = cond(Aestd1);
    end
    if cond_num < 5 % Set an arbitrary limit for condition number
        break;
    end
end

index(set_independent(min, :)) = [];
Aesti1 = Amat1(:, index);
Restc = -inv(Aestd1)*Aesti1;

fprintf("The true regression matrix is given by Restc: ");
disp(Restc);


% Hypothesis testing begins
lambda = diag(D1);
alpha = 0.05;

Test_statistic = zeros(24,1);
Test_criterion = zeros(24,1);
Degrees_of_freedom = zeros(24,1);
k_value = (2:1:25)';
n = 25;

for d = (n-1):-1:2
   Degrees_of_freedom(d-1) = (d-1)*(d+2)/2;
   lbar = mean(lambda(n-d+1:end));
   Test_statistic(d-1) = (nvar-1)*(d*log(lbar)-sum(log(lambda(n-d+1:end))));
   Test_criterion(d-1) = chi2inv(1-alpha,Degrees_of_freedom(d-1));
end
flag = 1;
d = n-2;
dest = 1;
while flag
    if ( Test_statistic(d) > Test_criterion(d) )
        d = d-1;
        if ( d < 2) 
           flag = 0;
        end
    else
        dest = d + 1;
        flag = 0;
    end
end

disp(table(k_value,Degrees_of_freedom,Test_statistic,Test_criterion));
nfact = nvar - dest;

fprintf("We get dest as 9 which is the true number of constraints by implementing PCA\n");
disp(dest);

L = diag(D1);
err_var = mean(L(end - d + 1:end));

fprintf("The mean error variances are given by: %0.4f\n", err_var)

%% D

maxdiff = max(max(abs(Restc - Restb)));
fprintf("The value of maximum difference comes out as: %0.4f\n", maxdiff)

%% E
theta_pca = subspace(Restb', Restc');
fprintf("The value of theta PCA is: %0.4f", theta_pca)

fprintf("The true regression matrix is given by Restb: ");
disp(Restb);

fprintf("\nAs the value of theta PCA is small, we can say that the solution is correct.")
fprintf("\nAlong with that, even maxdiff is small, which means that in the same way as theta subspace, " + ...
    "we are getting a correct answer. \n ")
