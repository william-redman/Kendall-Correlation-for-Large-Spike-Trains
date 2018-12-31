%% Knight Kendall Correlation 
%
%   Code written for "An O(n) Method of Calculating Kendall Correlations 
%   of Spike Trains" - William Redman. Discussion of the implementation is
%   in Supplmentary 1 Text. 
%   
%   Contact info: wredman@ucsb.edu 
%
%   Written by WTR 12/20/2018 // Last updated by WTR 12/30/2018
%%-----------------------------------------------------------------------%%
%%
function [time, tau_b] = Knight_Kendall_Corr_2(X, Y)

%% Computing tau
tic 
X_0 = find(X == 0); 
X_1 = find(X == 1); 
Y_0_X_0 = X_0(Y(X_0) == 0); 
Y_1_X_0 = setdiff(X_0, Y_0_X_0); 
Y_0_X_1 = X_1(Y(X_1) == 0); 
Y_1_X_1 = setdiff(X_1, Y_0_X_1);

X = [zeros(1, length(X_0)), ones(1, length(X_1))];
Y = [zeros(1, length(Y_0_X_0)), ones(1, length(Y_1_X_0)), zeros(1, length(Y_0_X_1)), ...
    ones(1, length(Y_1_X_1))]; 

[s_minus, ~] = mergesort_with_inversions2(Y);

T = 0.5 * (sum(X)^2 - sum(X) + sum(abs(X - 1))^2 - sum(abs(X - 1)));
U = 0.5 * (sum(Y)^2 - sum(Y) + sum(abs(Y - 1))^2 - sum(abs(Y - 1))); 

K = length(Y_0_X_0) * length(Y_1_X_1) - s_minus;
n_0 = (length(X)^2 - length(X)) / 2;
tau_b = K / (sqrt(n_0 - T) * sqrt(n_0 - U)); 

time = toc;

end

function [count, sorted] = mergesort_with_inversions2(A)
% Count is the number of inversions; sorted is the sorted array.
% Code taken from https://stackoverflow.com/questions/31557266/counting-inversions-in-matlab-using-mergesort
% and adjusted by WTR. 

n = length(A);
if n == 1
    count = 0;
    sorted = A;
else
    m = ceil(n/2);
    [count1, sorted1] = mergesort_with_inversions2(A(1:m));
    [count2, sorted2] = mergesort_with_inversions2(A(m+1:n));
    [crosscount, sorted] = merge2(sorted1, sorted2);
    count = count1 + count2 + crosscount;
end
end

function [crosscount, z] = merge2(x, y)

n = length(x); 
m = length(y); 
z = zeros(1, n+m);

ix = 1;
iy = 1;
crosscount = 0;
for iz = 1:(n+m)
    if ix > n
        z(iz) = y(iy);
        iy = iy + 1;
    elseif iy > m
        z(iz) = x(ix);
        ix = ix + 1;
    elseif x(ix) <= y(iy)
        z(iz) = x(ix);
        ix = ix + 1;
    elseif x(ix) > y(iy)
        z(iz) = y(iy);
        iy = iy + 1;
        crosscount = crosscount + (n + 1 - ix); 
    end
end
end


