function x_mod_N = Modified_mod (x, N)

% if x=N then modified mod (x, N)= N

x_mod_N  = mod(x,N);
ind=find(x_mod_N==0);

x_mod_N (ind) = N;

end

