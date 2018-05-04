%------------------------------------------------------
% Define parameters of problem
l = 0.5;
start = 0;

N = 100;
error1 = zeros(N-2,1);
count1 = zeros(N-2,1);
error2 = zeros(N-2,1);
count2 = zeros(N-2,1);

for n = 3:N

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
    
    %------------------------------------------------------
    % Set actual solution for comparison
    x = linspace(0,l,n+2);
    TA = BCs(1);
    TB = BCs(2);
    Tact = TA + (TB - TA)/l *x;

    %------------------------------------------------------
    % Solve problem
    % Inputs: algorithm, initial temp., thermal
    % conductivity, area, spacing
    % For algorithm: 1 = TDMA, 2 = gaussSeidel
    
    [T1, count1(n-2)] = solve(1, n, T0, BCtype, BCs, k, A, dx);
    error1(n-2) = getResidual(Tact', T1);
    
    [T2, count2(n-2)] = solve(2, n, T0, BCtype, BCs, k, A, dx);
    error2(n-2) = getResidual(Tact', T2);
end

%------------------------------------------------------
% Plot
plot(x,Tact)
hold
plot(X,T1, 'o')
plot(X,T2, '+')
axis([0 0.5 0 500])
ax = gca;
ax.FontSize = 12;
xlabel('x (m)', 'FontSize', 15, 'Interpreter', 'tex')
ylabel('Temperature (K)', 'FontSize', 15, 'Interpreter', 'tex')

n = 3:N;

figure(2)
plot(n,error1, 'LineWidth', 1.2)
hold
plot(n,error2, 'LineWidth', 1.2)
ax = gca;
ax.FontSize = 12;
xlabel('n', 'FontSize', 20, 'Interpreter', 'tex')
ylabel('Error', 'FontSize', 20, 'Interpreter', 'tex')


figure(3)
plot(n,count1, 'LineWidth', 1.2)
hold
plot(n,count2, 'LineWidth', 1.2)
ax = gca;
ax.FontSize = 12;
xlabel('n', 'FontSize', 20, 'Interpreter', 'tex')
ylabel('Iterations', 'FontSize', 20, 'Interpreter', 'tex')







