%% Kendall Correlation for Spike Trains 2
%
%   Code written based on method presented in "An O(n) Method of 
%   Calculating Kendall Correlations of Spike Trains" - William Redman. All
%   comments below referring to equations are referring to equations in
%   that paper. 
%   
%   Contact info: wredman@ucsb.edu 
%
%   Written by WTR 11/05/2018 // Last updated by WTR 12/30/2018
%%-----------------------------------------------------------------------%%
%%
function [time, tau] = Kendall_Corr_for_Spike_Trains_2(X,Y)

%% Initializing 
tic
n = length(X); %it is assumed here that X and Y are of the same size
n_0 = n * (n - 1) /2;

A = find((X + Y) == 2);
S = find((X + Y) == 0);

delta_X = find((X - Y) == 1);
delta_Y = find((Y - X) == 1);

%% Calculating K values

K_plus = length(A) * length(S); %Eq. 9

K_minus = length(delta_X) * length(delta_Y); %Eq. 12

K = K_plus - K_minus;

%   Eq. 11
ties_X = (sum(X) * (sum(X) - 1) / 2) + ((n - sum(X)) * (n - sum(X) - 1) / 2); %Eq. 13
ties_Y = (sum(Y) * (sum(Y) - 1) / 2) + ((n - sum(Y)) * (n - sum(Y) - 1) / 2);

%   Eq. 12
tau = K / (sqrt(n_0 - ties_X) * sqrt(n_0 - ties_Y)); %Eq. 14

time = toc;

end

