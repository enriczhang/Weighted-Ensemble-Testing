%%%%%%%%%%%%%%%%%%%%%%% plot WE particles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This function plots the evolving WE particles, for visualization purposes

%NOTES:
%This should be commented out in WE_simulation.m for long runs, 
%as it slows down simulations by a large amount

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if t==1

%plot initial particles
x = linspace(0,1,n);
[X,Y] = meshgrid(x,x);
Z = V(X,Y);
close all 
figure
contour(X,Y,Z,40)
hold on
if type == 'u'
z = linspace(0,1,n);
contour(z,z,reshape(RMSD,n,n)',[0 1/4 2/4 3/4 1 5/4 6/4],'LineColor','k');
end
line([0.5 0.5],[0.9 1], 'color','k')
line([0.6 0.6],[0.9 1], 'color','k')
line([0.5 0.6],[0.9 0.9], 'color','k')
line([0.5 0.6],[1 1], 'color','k')
plot(0.1,0.5,'xk','markersize',20)
handle = plot(xs(1,:),xs(2,:),'.','color','b','markersize',16);   %plot particles
axis([0 1 0 1])
pause(1)

else

%delete last WE particle configuration from plot
delete(handle)

%plot colored WE particle configuration: last in A = blue, last in B = red
handle = plot(xs(1,:),xs(2,:),'.','markersize',16,'color','b');
axis([0 1 0 1])
pause(0.01)

end
