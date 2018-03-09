%------------------------------------------------------
% Define parameters of problem
l = [0.5 0.5];
start = [0 0];
n = [100 100];

temp = 1;
flux = 2;
BCtype = [flux temp; flux temp];
BCs = [-900 100; -900 500];

% Make k and A a function of space in general
k = 100*ones(n(2), n(1));
A = 1e-3*ones(n(2), n(1));

%------------------------------------------------------
% Define grid spacing for each dimension over the
% whole space
dx = l(1)/n(1)*ones(n(1)+1,n(1)+1);

dx(1,:) = l(1)/(2*n(1));
dx(:,1) = l(1)/(2*n(1));
dx(n(1)+1,:) = l(1)/(2*n(1));
dx(:,n(1)+1) = l(1)/(2*n(1));
%---------------------------------
dy = l(2)/n(2)*ones(n(2)+1,n(2)+1);

dy(2,:) = l(2)/(2*n(2));
dy(:,2) = l(2)/(2*n(2));
dy(n(2)+1,:) = l(2)/(2*n(2));
dy(:,n(2)+1) = l(2)/(2*n(2));

%------------------------------------------------------
% Generate grid points
X = zeros(n(1)+2,n(1)+2);
X(:,1) = start(1);

for i = 2:n(1)+2
    for j = 1:n(2)+2
        X(j,i) = X(j,i-1)+ dx(j,i-1);
    end
end

Y = zeros(n(2)+2,n(2)+2);
Y(end,:) = start(2);

for i = 1:n(1)+2
    for j = 2:n(2)+2
        Y(j,i) = Y(j-1,i)+ dy(j-1,i);
    end
end

%------------------------------------------------------
% Set initial conditions for T
T0 = zeros(n(1)+2,n(2)+2);
T0(:,end) = BCs(1,2);
T0(end,:) = BCs(2,2);

%======================================================
%------------------------------------------------------
% Solve problem
% Inputs: algorithm, initial temp., thermal
% conductivity, area, spacing
% For algorithm: 1 = TDMA, 2 = gaussSeidel
[T1, count1] = solve(1, n, T0, BCtype, BCs, k, A, dx, dy);
[T2, count2] = solve(2, n, T0, BCtype, BCs, k, A, dx, dy);

%------------------------------------------------------
% Set actual solution for comparison
x = linspace(0,lX);
Tact = TA + (TB - TA)/lX *x;

%------------------------------------------------------
% Plot
plot(x,Tact)
hold
plot(gridPoints,T1, 'o')
plot(gridPoints,T2, '+')
