syms a b
%mass: UiCAH+ 957.5, Adenine+ 135.13, sugar+ 69.0817, I 126.904

eqns = [a*sqrt(957.5) + b == 105.38, a*sqrt(126.904) + b == 38.533];
%eqns = [a*sqrt(957.5) + b == 105.38, a*sqrt(69.0817) + b == 28.4355];
S = solve(eqns);
sol = [S.a;S.b];
sol = double(sol)
%[X,y] = equationsToMatrix(eqns,a,b);
%z = X\y
%a = 6.4227e+03
%b = 1.2102e+04


