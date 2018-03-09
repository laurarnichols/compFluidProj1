%------------------------------------------------------
% Define parameters of problem
l = 0.5;
start = 0;
n = 10;

temp = 1;
flux = 2;
BCtype = [temp temp];
BCs = [100 500];

% Make k and A a function of space in general
k = 100*ones(n+2,1);
A = 1e-3*ones(n+2,1);

%------------------------------------------------------
% Define grid spacing for each dimension over the
% whole space
dx = l/n*ones(n+1,1);

dx(1) = l/(2*n);
dx(end) = l/(2*n);

%------------------------------------------------------
% Generate grid points
X = zeros(n+2,1);
X(1) = start(1);

for i = 2:n+2
    X(i) = X(i-1)+ dx(i-1);
end

%------------------------------------------------------
% Set initial conditions for T
T0 = zeros(n+2,1);
T0(1) = BCs(1);
T0(end) = BCs(2);

%======================================================
%------------------------------------------------------
% Solve problem
% Inputs: algorithm, initial temp., thermal
% conductivity, area, spacing
% For algorithm: 1 = TDMA, 2 = gaussSeidel
[T1, count1] = solve(1, n, T0, BCtype, BCs, k, A, dx);
[T2, count2] = solve(2, n, T0, BCtype, BCs, k, A, dx);

%------------------------------------------------------
% Set actual solution for comparison
x = linspace(0,l);
TA = BCs(1);
TB = BCs(2);
Tact = TA + (TB - TA)/l *x;

%------------------------------------------------------
% Plot
plot(x,Tact)
hold
plot(X,T1, 'o')
plot(X,T2, '+')