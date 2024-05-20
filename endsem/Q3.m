% Endsem CH5440 
% Ojas Phadake - CH22B007

clc;
clear all;
close all;

load arx.mat

%% A 
% Scaling appropriately the measurements:
ymeass = ymeas/stdey;
umeass = umeas/stdeu;

L = 10; % lag
nsamples = 1024;

ZL = [];
for i = L+1:-1:1
    ZL = [ZL ymeass(i:nsamples+i-L-1)];
end
for i = L+1:-1:1
    ZL = [ZL umeass(i:nsamples+i-L-1)];
end
[u s v] = svd(ZL/sqrt(nsamples-L),'econ');
lambda = diag(s).^2;

% Hypothesis test
alpha = 0.05;
nvar = size(ZL,2);
tau = zeros(nvar-2,1);
crit = zeros(nvar-2,1);
for d = nvar-1:-1:2
   nu = (d-1)*(d+2)/2;
   nprime = nsamples-nvar - (2*nvar+11)/6;
   lbar = mean(lambda(nvar-d+1:end));
   tau(d-1) = nprime*(d*log(lbar)-sum(log(lambda(nvar-d+1:end))));
   crit(d-1) = chi2inv(1-alpha,nu);
end

disp(table(crit, tau))

% Hypothesis testing
flag = 1;
d = nvar-2;
dest = 1;

while flag
    if ( tau(d) > crit(d) )
        d = d-1;
        if ( d < 2) 
           flag = 0;
        end
    else
        dest = d + 1;
        flag = 0;
    end
end

d1 = dest;
eta1 = L - d1 + 1;

% Last PCA applied using known eta
Zeta1 = [];
for i = eta1+1:-1:1
    Zeta1 = [Zeta1 ymeass(i:nsamples+i-eta1-1)];
end
for i = eta1+1:-1:1
    Zeta1 = [Zeta1 umeass(i:nsamples+i-eta1-1)];
end
[u s v] = svd(Zeta1/sqrt(nsamples-eta1),'econ');
theta1 = v(:,end)';
theta1(1:eta1+1) = theta1(1:eta1+1)/stdey;
theta1(eta1+2:end) = theta1(eta1+2:end)/stdeu;
theta1 = theta1/theta1(1);

%% Now applying for L = 20

fprintf("Now carrying out the same for L = 20. ")
L = 20; % Modified lag

ZL = [];
for i = L+1:-1:1
    ZL = [ZL ymeass(i:nsamples+i-L-1)];
end
for i = L+1:-1:1
    ZL = [ZL umeass(i:nsamples+i-L-1)];
end
[u s v] = svd(ZL/sqrt(nsamples-L),'econ');
lambda = diag(s).^2;

% Hypothesis test
alpha = 0.05;
nvar = size(ZL,2);
tau = zeros(nvar-2,1);
crit = zeros(nvar-2,1);
for d = nvar-1:-1:2
   nu = (d-1)*(d+2)/2;
   nprime = nsamples-nvar - (2*nvar+11)/6;
   lbar = mean(lambda(nvar-d+1:end));
   tau(d-1) = nprime*(d*log(lbar)-sum(log(lambda(nvar-d+1:end))));
   crit(d-1) = chi2inv(1-alpha,nu);
end

disp(table(crit, tau))

% Hypothesis testing
flag = 1;
d = nvar-2;
dest = 1;

while flag
    if ( tau(d) > crit(d) )
        d = d-1;
        if ( d < 2) 
           flag = 0;
        end
    else
        dest = d + 1;
        flag = 0;
    end
end

d = dest;
eta = L - d + 1;

fprintf("We get the value of eta as %0.1d", eta);
fprintf("\nThis means that the system depends on the memory of inputs and outputs of the last 2 time instants. ")

% Last PCA applied using known eta
Zeta = [];
for i = eta+1:-1:1
    Zeta = [Zeta ymeass(i:nsamples+i-eta-1)];
end
for i = eta+1:-1:1
    Zeta = [Zeta umeass(i:nsamples+i-eta-1)];
end
[u s v] = svd(Zeta/sqrt(nsamples-eta),'econ');
theta = v(:,end)';
theta(1:eta+1) = theta(1:eta+1)/stdey;
theta(eta+2:end) = theta(eta+2:end)/stdeu;
theta = theta/theta(1);

fprintf("\nThe obtained coefficient vector is: ")
disp(theta)
fprintf("\nThus we notice that considering the lag as 10 and 20 gave different values of eta.\n")

fprintf("************************************************************************************************************************")

%% B
fprintf("For L = 10:\n")
fprintf("We get the value of eta as %0.1d", eta1);
fprintf("\nThis means that the system depends on the memory of inputs and outputs of the last 2 time instants. ")
fprintf("\nThe obtained coefficient vector is: ")
disp(theta1)

fprintf("For L = 20:\n")
fprintf("We get the value of eta as %0.1d", eta);
fprintf("\nThis means that the system depends on the memory of inputs and outputs of the last 2 time instants. ")
fprintf("\nThe obtained coefficient vector is: ")
disp(theta)

%% Bootstrapping

nboot = 100;
nsub = 700;

[N, nvar] = size(Zeta1); % 1019, 12
theta1 = zeros(nboot,nvar); 

for i = 1:nboot
    ind = randperm(N);
    Zsub = Zeta1(ind(1:nsub),:);
    [u s v] = svd(Zsub/sqrt(nsub),'econ');
    theta1(i,:) = v(:,end)';
    theta1(i,:) = theta1(i,:)/theta1(i,1);  % Normalize theta vector so that first coefficient is unity
end

% Find mean and std of theta elements

thetamean = mean(theta1);
thetastd = std(theta1);
thetal = thetamean - 2*thetastd;
thetau = thetamean + 2*thetastd;
ind = 1:1:6;

fprintf("\nThe 95% confidence intervals using bootstrapping are as follows for L=10:\n")
to_show = table(ind', thetal', thetau');
disp(to_show);

nboot = 100;
nsub = 700;

[N, nvar] = size(Zeta);
theta = zeros(nboot,nvar); 

for i = 1:nboot
    ind = randperm(N);
    Zsub = Zeta(ind(1:nsub),:);
    [u s v] = svd(Zsub/sqrt(nsub),'econ');
    theta(i,:) = v(:,end)';
    theta(i,:) = theta(i,:)/theta(i,1);  % Normalize theta vector so that first coefficient is unity
end

% Find mean and std of theta elements

thetamean = mean(theta);
thetastd = std(theta);
thetal = thetamean - 2*thetastd;
thetau = thetamean + 2*thetastd;
ind = 1:1:16;

fprintf("\nThe 95% confidence intervals using bootstrapping are as follows for L=20:\n")
to_show = table(ind', thetal', thetau');
disp(to_show);

%% 
