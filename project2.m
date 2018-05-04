%u0 = linspace(0, 1,50);
%u0 = 1e-4;
%n0 = 5:5:100;
n0 = 90;
filename = 'epsilon.gif';
u0 = 0.0;
epsilon = [1 5e-1 1e-1 5e-2 1e-2 5e-3 1e-3 5e-4 1e-4 5e-5 1e-5 5e-6 1e-6];
%epsilon = 1e-6;
for imNum = 1:length(epsilon)
    %------------------------------------------------------
    % Define parameters of problem
    l = [0.5 0.5];
    start = [0 0];
    n = [floor(n0) floor(n0)];

    temp = 1;
    flux = 2;
    BCtype = [temp temp; temp temp];
    %BCs = [-100 100; 500 250];
    BCs = [200 200; 200 200];

    %------------------------------------------------------
    % Set initial conditions for T
    T0 = setInitT(n, BCtype, BCs);
    %------------------------------------------------------
    % Make parameters a function of space in general
    

    k = 0.026*ones(n(2)+2, n(1)+2);
    A = 0.5/(n(2)+2)*ones(n(2)+2, n(1)+2);
    u = u0*ones(n(2)+2, n(1)+2);
    rho = 1.2044*ones(n(2)+2, n(1)+2);

    %------------------------------------------------------
    % Define grid spacing for each dimension over the
    % whole space
    dx = l(1)/n(1)*ones(n(2)+2,n(1)+1);

    dx(:,1) = l(1)/(2*n(1));
    dx(:,end) = l(1)/(2*n(1));
    %---------------------------------
    dy = l(2)/n(2)*ones(n(2)+1,n(1)+2);

    dy(2,:) = l(2)/(2*n(2));
    dy(end,:) = l(2)/(2*n(2));

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

    %======================================================
    %------------------------------------------------------
    % Solve problem
    % Inputs: algorithm, initial temp., thermal
    % conductivity, area, spacing
    % For algorithm: 1 = TDMA, 2 = gaussSeidel
    [T, count1] = solve(1, n, T0, BCtype, BCs, k, A, dx, dy, u, rho, epsilon(imNum));

    %------------------------------------------------------
    % Plot
    %plot(x,Tact)
    %hold

    %Tact = getAnalytic(X, Y, T);

    figure(1)
    h = contourf(X,Y,T,20);
    %surf(X,Y,T,'edgecolor', 'none')
    c = colorbar;
    c.Label.String = 'T (K)';
    colormap(jet)
    %az = 45;
    %el = 30;
    %view(az, el);
    %caxis([200, 300]);
    xlabel('X (m)', 'FontSize', 15, 'Interpreter', 'tex')
    ylabel('Y (m)', 'FontSize', 15, 'Interpreter', 'tex')
    %zlabel('T (K)', 'FontSize', 15, 'Interpreter', 'tex')
    annotation('textbox',...
    [0.5 0.75 0.25 0.1],...
    'String',['\epsilon =' num2str(epsilon(imNum), '%.3e')],...
    'FontSize',14, ...
    'Tag' , 'textBox', ...
    'BackgroundColor',[1.0 1.0 1.0])
    axis tight manual % this ensures that getframe() returns a consistent size
    xlim([0 0.5])    % set x limits
    ylim([0 0.5])     % set y limits
    %zlim([200 310])
    pause(0.2)
    frame = getframe(gcf); 
    im = frame2im(frame); 
    [A,map] = rgb2ind(im,256); 

    % Write to the GIF File 
    if imNum == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.3);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.3);
    end
    delete(findall(gcf,'Tag','textBox'));
    
    %input('Waiting...');
    imNum
end

% figure(2)
% h = contourf(X,Y,Tact,10);
% colorbar
% xlabel('X (m)', 'FontSize', 15, 'Interpreter', 'tex')
% ylabel('Y (m)', 'FontSize', 15, 'Interpreter', 'tex')
