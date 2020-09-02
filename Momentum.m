%%
close all
clear all
clc
%% Symbolic function definition and function gradient 
syms X [1 2]
syms Y(X)
% sfunY = sin(X1) + cos(X2) + cos(X1).*sin(X2)
sfunY = 0.5*X1^2 + 0.4*X2^2 + (0.1*X2.^2 + 0.1*X1.^2).*cos(X1).*sin(X2);

gradY = gradient(sfunY);

%% Create anonymous functions for numeric calculations

funY = matlabFunction(sfunY);
funGradY = matlabFunction(gradY);

funY = @(X) funY(X(1),X(2));
funGradY = @(X) funGradY(X(1),X(2));

%%

X = [85;95]; % Initial guess
Y = funY(X);
optimData = table(X',Y,zeros(1,length(X)),zeros(1,length(X)),...
    'VariableNames',{'X','Y','dX','gradY'});

lambda = 0.001; % learning rate
beta = 0.95; % momentum; (0,1)
dX = zeros(length(X),1); 
for i = 1:1000
    gradY = funGradY(X);
    
    dX = gradY + beta.*dX; % dX(k) = gradY + beta*dX(k-1)
        
    newRow = table(X',Y,-lambda.*dX',-gradY',...
        'VariableNames',{'X','Y','dX','gradY'});
    optimData = [optimData; newRow];
    
    X = X - lambda.*dX;
    if dX < 0.005
        X = X + 5*rand(2,1);
    end
    Y = funY(X);
end

%%

% [Xm,Ym] = meshgrid(-100:1:100);
% Zm = arrayfun(@(X,Y)funY([X,Y]),Xm,Ym);
% [DX,DY] = gradient(Zm);
% DX = -DX;
% DY = -DY;

figure

% quiver3(Xm,Ym,Zm,DX,DY,zeros(size(DY)))

sizeLims = 100;

f1 = figure('units','normalized','outerposition',[0 0 1 1],'name','OptimizerFigure');
ax = gca;
hold on;
fs = fsurf(ax,sfunY,[0.5*sizeLims sizeLims],'MeshDensity',50);
fs.FaceAlpha = 0.75;

view(120,60);
pause(2);
for i = 1:height(optimData)
    
    try 
        delete(a1);
        delete(a2);
        delete(p);
    end
    
    p = plot3(ax,optimData.X(i,1),optimData.X(i,2),optimData.Y(i),'r.',...
        'MarkerSize',40);
    plot3(ax,optimData.X(i,1),optimData.X(i,2),optimData.Y(i),'r.',...
        'MarkerSize',10);
    a1 = arrow(ax,[optimData.X(i,:), optimData.Y(i)],...
        optimData.dX(i,:),'k');
    a2 = arrow(ax,[optimData.X(i,:), optimData.Y(i)],...
        optimData.gradY(i,:),'b');
    ax.XLim = [0.5*sizeLims sizeLims];
    ax.YLim = [0.5*sizeLims sizeLims];
    pause(0.2);
end
hold off
%%

function arrowPlot = arrow(ax,P1,P2,clr)
    arrowPlot = plot3(ax,[P1(1);P2(1)+P1(1)],[P1(2);P2(2)+P1(2)],[P1(3);P1(3)],clr,'LineWidth',2);
end