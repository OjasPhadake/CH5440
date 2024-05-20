% Assignment 3 Question 3
% CH22B007 Ojas Phadake

clc; clear; close all;
cell_all = struct2table(dir("yalefacespng"));
cell_all = cell_all.name;

matrix = zeros(243, 320, 90); % This will contain all the images
vector_matrix = zeros(77760, 90);% This will contain all the reshaped images in vector form
for i=1:90
    matrix(:, :, i) = imread('yalefacespng/' + string(cell_all(i+2)));
    vector_matrix(:, i) = reshape(matrix(:, :, i), 77760, 1);
end

% In case an image is to be seen, firstly convert it to uint8 and then
% display the reshaped image in size 243 * 320

%% Part 1. Norm Images
vector_norms = vector_matrix(:, 6:6:90); % All normal images
classified_person = zeros(90, 1);
image_num = 1:1:90; image_num = image_num';
true_person_num = floor((image_num-1)/6) + 1;
result = [];
correct_count = 0;

for i=1:90
    % The norm images will be 6*i kind amongst the above dataset
    % After that, use a for loop from 1:15 to check the distances
    dist = zeros(15, 1);
    for j=1:15
        dist(j) = norm(vector_norms(:, j) - vector_matrix(:, i));
    end
    classified_person(i) = find(dist == min(dist, [], "all"));
        if classified_person(i) == true_person_num(i)
            result = [result; "YES"];
            correct_count = correct_count + 1;
        else
            result = [result; "NO"];
        end
end

to_show = table(image_num, true_person_num, classified_person, result);
disp("The following table contains the Image number, Correct label of person, Classified estimate and Result");
disp(to_show)
disp("The number of images obtained correctly by comparing with the norm images is: " + correct_count);
disp("-------------------------------------------------------------------------")

%% Part 2. First Principle Component 
% B - First Principal Component use to make image

PC_comp = zeros(77760, 15);
variance_arr = zeros(15, 1);

for i=1:15
    Z = vector_matrix(:, ((6*i) - 5):(6*i));
    Z_s = Z - mean(Z, 1);
    S_zs = cov(Z_s);
    [V, D] = eig(S_zs);
    D = sort(diag(D), 'descend');
    variance_arr(i) = D(1)/sum(D, "all");
    PC_comp(:, i) = Z*((V(:, end)).^2); % This contains the first PC of every person as each column
end

disp("We multiply by ((V(:, end)).^2) so that condition of orthogonality is " + ...
    "maintained, and also to obtain weights between 0 and 1, so that the " + ...
    "resulting number also lies in between 0 and 255")

% imshow((reshape(uint8(PC_comp(:, 2)), 243, 320)))

classified_person = zeros(90, 1);
image_num = (1:1:90)';
result = [];
correct_count = 0;

for i=1:90

    dist = zeros(15, 1);
    for j = 1:15
        dist(j) = norm(PC_comp(:, j) - vector_matrix(:, i));
    end
    est_person_img = find(dist == min(dist, [], 'all')); % Minimum index is found
        if est_person_img == true_person_num(i)
            result = [result; "YES"];
            correct_count = correct_count + 1;
        else
            result = [result; "NO"];
        end
end

to_show = table(image_num, true_person_num, result);
disp("The following table contains the Image number, Correct label of person, Classified estimate and Result");
disp(to_show)
disp("The number of images obtained correctly by comparing with the first PC is: " + correct_count);
disp("-------------------------------------------------------------------------")

%% Part 3. Two Principle Components

PC_comp_both = zeros(77760, 30);
variance_arr = zeros(15, 1);

for i=1:15
    Z = vector_matrix(:, ((6*i) - 5):(6*i));
    Z_s = Z - mean(Z, 1);
    S_zs = cov(Z_s);
    [V, D] = eig(S_zs);
    D = sort(diag(D), 'descend');
    variance_arr(i) = (D(1) + D(2))/sum(D, "all");
    PC_comp_both(:, 2*i-1) = Z*((V(:, end)).^2);  
    PC_comp_both(:, 2*i) = Z*((V(:, end-1)).^2);
end

classified_person_1 = zeros(90, 1);
classified_person_2 = zeros(90, 1);
image_num = (1:1:90)';
result = [];
correct_count = 0;

for i=1:90

    dist_first = zeros(15, 1);
    dist_second = zeros(15, 1);
    for j = 1:15
        dist_first(j) = norm(PC_comp_both(:, 2*j-1) - vector_matrix(:, i));
        dist_second(j) = norm(PC_comp_both(:, 2*j) - vector_matrix(:, i));
    end
    est_person_img_first = find(dist_first == min(dist_first, [], 'all'));
    est_person_img_second = find(dist_second == min(dist_second, [], 'all'));

    if (est_person_img_first == true_person_num(i)) || (est_person_img_second == true_person_num(i))
        result = [result; "YES"];
        correct_count = correct_count + 1;
    else
        result = [result; "NO"];
    end
    classified_person_1(i) = est_person_img_first;
    classified_person_2(i) = est_person_img_second;
end

to_show = table(image_num, true_person_num, classified_person_1, classified_person_2, result);
disp("The following table contains the Image number, Correct label of person, Classified estimate by PC1, Classified estimate by PC2 and Result");
disp(to_show)
disp("The number of images obtained correctly by comparing with the first 2 PCs is: " + correct_count);
