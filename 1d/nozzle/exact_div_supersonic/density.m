function r = density(x, s1, r1)

% define function handle
fun = @(d) densityeq(d, x, s1);

r = fzero(fun, r1);
