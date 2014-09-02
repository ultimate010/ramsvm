clear
version = 0.1;
function ans = compute_first(k)
    t = (k - 1) ^ (-0.5);
    l = diag(eye(k - 1))';
    ans = t*l;
endfunction

function ans = compute_last(k)
ans = [];
for j = 2:k;
    l = ones(1, k - 1);
    t = -(1 + k ^ 0.5) / ((k - 1) ^ (3 / 2));
    a1 = t * l;
    a2 = zeros(1, k - 1);
    a2(j - 1) = 1;
    t = (k / (k - 1)) ^ 0.5;
    a2 = t * a2;
    t = a1 + a2;
    ans = [ans;t];
endfor
endfunction



function ans = compute_y(k)
    first = compute_first(k);
    last = compute_last(k);
    ans = [first; last];

endfunction


function ans = N_FUNC(i, j)
global Y K;
if j == Y(i)
    ans = K - 1;
else
    ans = 1;
end
endfunction

function ans = H_FUNC(i, j)
global Y;
if j == Y(i)
    ans = 1;
else
    ans = -1;
end
endfunction


%variable define
clear

data =dlmread("data.txt.2", " ", 0, 0);
global X = data(:, 1:6);
global Y = data(:, 7);
%global X = [1, 2, 3, 5, 8, 9; 3, 23, 23, 12, 33, 88; 6, 5, 4, 3, 2, 1; 12, 23, 21, 43, 42, 15; 53, 21, 53, 23, 12, 54; 87, 23, 34, 12, 23, 43;];
%global Y = [2; 2; 2; 2; 1; 1];
%data =dlmread("data.txt.2", " ", 0, 0);
c1 = load("c1.txt");
c2 = load("c2.txt");
c3 = load("c3.txt");

%global X = [c1(1: 20,:); c2(1: 20, :); c3(1: 20, :)]; %[1, 2, 3, 5, 8, 9; 3, 23, 23, 12, 33, 88; 6, 5, 4, 3, 2, 1; 12, 23, 21, 43, 42, 15; 53, 21, 53, 23, 12, 54; 87, 23, 34, 12, 23, 43;]
%global Y = [ones(20, 1); 2 * ones(20, 1); 3 * ones(20, 1)];
%global X = [c1(1: 20,:); c2(1: 20, :); c3(1: 20, :)];
%global Y = [ones(20, 1); 2 * ones(20, 1); 3 * ones(20, 1)];



global N = size(X)(1);
global K = max(Y);
%global ALPHA = zeros(N, K);
global ALPHA = randn(N, K);
global Y_MATRIX = compute_y(K);
global NITER = 1000;
global ITERCOUNT = 0;
global LAMBDA = 0.9;

global XLEN = size(X)(2)
global BETA = zeros(K - 1, XLEN)
global ALPHA_LAST;
global EPSILON = 0.0001;

disp("X:")
X
disp("Y:")
Y

function ans = compute_a_i_j(i, j)
% i to s, j to t
global LAMBDA N Y X K Y_MATRIX ALPHA;
Y_t_q = 0;
for q = 1: K -1
    Y_t_q = Y_t_q + Y_MATRIX(j, q);
end

TEMP1 = zeros(size(X)(2), 1);
for _i = 1:N
    TEMP2 = 0;
    for _j = 1:K
        if _j == j && _i == i
        else
            TEMP2 = TEMP2 + ALPHA(_i, _j) * Y_MATRIX(_j, q) * H_FUNC(_i, _j);
        end
    endfor
TEMP1 = TEMP1 + (X(_i, :)' .* TEMP2);
endfor

fenzi1 = N_FUNC(i, j) * N * LAMBDA;
fenzi2 = (H_FUNC(i, j) * Y_t_q) * (X(i, :) * TEMP1);
fenzi = fenzi1 - fenzi2;
fenmu =  X(i, :) * X(i, :)';
ans = fenzi / fenmu;
endfunction


function iter_onetime()
global ALPHA ALPHA_LAST N K;
ALPHA_LAST = ALPHA;
for i = 1:N
    for j = 1:K
        a_i_j = compute_a_i_j(i, j);
        if a_i_j < 0
            ALPHA(i, j) = 0;
        elseif a_i_j > 0.5
            ALPHA(i, j) = 0.5;
        else
            ALPHA(i, j) = a_i_j;
        end
    endfor
endfor
endfunction

function ans = diff_between()
global ALPHA ALPHA_LAST;
d = ALPHA - ALPHA_LAST;
ans = sum(sum(d.*d));
endfunction

function solve()
global NITER ALPHA K XLEN Y_MATRIX BETA N X Y LAMBDA EPSILON;
printf("Max Iter:%d\n", NITER);
for i = 1:NITER
    iter_onetime();
    diff = diff_between();
    printf("Iter: %d diff: %3f\n", i, diff);
    if diff < EPSILON
        break;
    endif
endfor

% end solve alpha

for q = 1: K - 1;
    TEMP = zeros(1, XLEN);
    for _i = 1 : N
        for _j = 1 : K
            TEMP = TEMP + ALPHA(_i, _j) * Y_MATRIX(_j, q) * H_FUNC(_i, _j) .* X(_i, :);
        endfor
    endfor
    BETA(q,:) = TEMP ./ (N * LAMBDA);
endfor
% end solve beta
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Start compute%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("Begin solving:")
solve();

disp("Alpha:")
ALPHA
disp("Beta:")
BETA
%disp("Y_MATRIX:")
%Y_MATRIX

function ans = avgmin(_x)
global BETA K Y_MATRIX;
fx = zeros(1, K - 1);
for q = 1: K - 1
    fx(q) = BETA(q, :) * _x';
endfor
maxpos = -1;
max_num = -999999999999;

for i = 1: K
    cur = sum(fx .* Y_MATRIX(i, :));
    if cur > max_num
        max_num = cur;
        maxpos = i;
    endif
endfor

ans = maxpos;
endfunction

% begin predict 

predict = zeros(N, 1);
for i = 1:N
 predict(i, 1) = avgmin(X(i, :));
end
disp("Predict result:")
predict

% begin compute accuracy rate
Z = Y - predict;
size(find(Z == 0))(1) / size(Y)(1)
