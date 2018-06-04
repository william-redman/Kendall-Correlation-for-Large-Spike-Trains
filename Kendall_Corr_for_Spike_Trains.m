%% Kendall Correlation for Spike Trains
%
%   Code written based on method presented in "A Fast Way to Calculate 
%   Kendall Correlations of Large Spike Trains" - William Redman (...). All
%   comments below referring to equations are referring to equations in
%   that paper. 
%   
%   Contact info: wredman@nyu.edu 
%
%   WTR 06/02/2018
%%-----------------------------------------------------------------------%%
%%
function [time, tau] = Kendall_Corr_for_Spike_Trains(X,Y)

%% Initializing 
tic
n = length(X); %it is assumed here that X and Y are of the same size
n_0 = n * (n - 1) /2;

%   Eq. 4
A_X = find(X == 1); 
A_Y = find(Y == 1);

%   Eq. 5
A = intersect(A_X, A_Y);

%   Eq. 6
S_X = 1:n; S_X(A_X) = [];
S_Y = 1:n; S_Y(A_Y) = [];

%   Eq. 7
S = intersect(S_X, S_Y);

%   Eq. 9
delta_X = setdiff(A_X, A_Y); 
delta_Y = setdiff(A_Y, A_X);

%% Calculating K values
K_plus = 0; 
K_minus = 0;

%   First sum in eq. 8
for ii = 1:length(A)
    K_plus  = K_plus + length(find(S > A(ii)));
end

%   Second sum in eq. 8
for ii = 1:length(S) 
    K_plus = K_plus + length(find(A > S(ii)));
end

%   First sum in eq. 10
for ii = 1:length(delta_X)
    K_minus = K_minus + length(find(delta_Y > delta_X(ii)));
end

%   Second sum in eq. 10   
for ii = 1:length(delta_Y)
    K_minus = K_minus + length(find(delta_X > delta_Y(ii)));
end

K = K_plus - K_minus;

%   Eq. 11
ties_X = (sum(X) * (sum(X) - 1) / 2) + ((n - sum(X)) * (n - sum(X) - 1) / 2);
ties_Y = (sum(Y) * (sum(Y) - 1) / 2) + ((n - sum(Y)) * (n - sum(Y) - 1) / 2);

%   Eq. 12
tau = K / (sqrt(n_0 - ties_X) * sqrt(n_0 - ties_Y));

time = toc;

end

