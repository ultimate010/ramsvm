 function ans = compute_first(k) %j = 1
    t = (k - 1) ^ (-0.5);
    l = diag(eye(k))';
    ans = t*l;
endfunction

function ans = compute_last(k) %j >= 2
ans = [];
for j = 2:k;
    l = diag(eye(k))';
    t = -(1 + k ^ 0.5) / ((k - 1) ^ (3 / 2));
    a1 = t * l;
    a2 = 0 * l;
    a2(j) = 1;
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

% Test
for i = 2:5
    ans = compute_y(i);
    disp("next");
    disp(ans)
endfor
