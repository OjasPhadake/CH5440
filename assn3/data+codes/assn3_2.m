% Assignment 3, Multivariate Data Analysis CH5440
% confused zoologist 

clc;
clear;
close all;

%%
xbar = [9 68 129];
S = [7 21 34;
     21 64 102;
     34 102 186];

eigv_max = 250.4009;

[V, D] = eig(S);
normalised_eigvec = V(:, 3);
disp("The normalized eigenvector is: "); 
disp(normalised_eigvec);

remaining_eigvec = V(:, 1:2);

perc_var = D(3,3)/sum(diag(D));
disp("Percentage variance by the first PC is: ");
disp(perc_var);

% One possible set of linear relations must be the eigenvector
% corresponding to the smallest eigenvalue 

lin_reln = V(:, 1)';

data_point = [10.1 73 135.5] - xbar;
score = data_point*normalised_eigvec;
disp("")
disp("The scores matrix is: ")
disp(score)
data_denoised = score*normalised_eigvec' + xbar;

mass_est = (lin_reln(2)*73 + lin_reln(3)*135.5)/-lin_reln(1);
disp("Estimated mass of the lizard is: ");
disp(mass_est)
