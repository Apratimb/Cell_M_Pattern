%% Implementation of CAM in MATLAB
% Based on the Paper: Coherent Motion of Monolayer Sheets under Confinement and Its Pathological Implications
% S S Soumya, Animesh Gupta, Andrea Cugno, Luca Deseri, Kaushik Dayal, Dibyendu Das, Shamik Sen, Mandar M. Inamdar 
% Published: December 21, 2015https://doi.org/10.1371/journal.pcbi.1004670
% Author: Apratim Bajpai, PhD Candidate, Department of Mechanical
% Engineering, NYU Tandon School of Engineering

clearvars;
close all;
filename = 'C1S1.xlsx';
sheet = 1;


R_Pattern = 20;
unnamed = xlsread(filename,sheet);
points = unnamed(:,1:2);
Y0 = (min(unnamed(:,2))+max(unnamed(:,2)))/2;
X0 = (min(unnamed(:,1))+max(unnamed(:,1)))/2;
points = points - [X0*ones(size(points,1),1),Y0*ones(size(points,1),1)];
R_Pattern1 = max([max(points(:,1));max(points(:,2));abs(min(points(:,1)));abs(min(points(:,2)))]);
points = R_Pattern/R_Pattern1*points;
X0 = 0; % X_Coordinate of pattern center
Y0 = 0; %   Y_Coordinate of the pattern center
RI_Pattern= 0;

dij = zeros(size(points,1),size(points,1));
for i = 1:size(points,1)
    for j = 1:size(points,1)
        dij(i,j) = sqrt(sum((points(i,:) - points(j,:)).^2));
    end
end

DT = delaunayTriangulation(points(:,1),points(:,2));

x = DT.Points(:,1);
y = DT.Points(:,2);

% Get all the triangles attached to a particular vertex in the
% triangulation.  
attachedTriangles = vertexAttachments(DT);
for i = 1: size(x,1)
% Use the connectivity list to get the vertex indices of all these
% triangles
    verticesOfTI         = DT.ConnectivityList(attachedTriangles{i},:);
% Find all the unique vertices and remove the current vertex
    neighboursOfInternal{i} = setdiff(unique(verticesOfTI), i);
end
 
% Find the mean a0 of the points.
for vertex = 1:size(DT.Points,1)
    neighboring_vertex = neighboursOfInternal{vertex}';
    temp = 0; 
    ri = [points(vertex,1)-X0,points(vertex,2)-Y0]; %  r of the vertex
    nv = 0;
    for i = 1:length(neighboring_vertex) 
        rj = [points(neighboring_vertex(i),1)-X0,points(neighboring_vertex(i),2)-Y0]; % r of the neighboring vertex
        eij = (rj-ri);  
        mod_rij = abs((eij(1,1)^2+eij(1,2)^2)^0.5);
        temp = temp+mod_rij;
        nv = nv+1;
    end
    a0(vertex) = temp/nv;
end



[theta,r] = cart2pol(points(:,1),points(:,2));
%% Populating the direction preference matrix
ni = unnamed(:,3);
% ni = pi*rand(size(points,1),1);

%% Solving the differential equation
input_solver = [r.*cos(theta);r.*sin(theta);ni];
[sol,vel] = rk4_systems_WORKINGRing(0, 100, 2000, input_solver,R_Pattern,RI_Pattern,a0);
vel = vel';
yv = sol(:,2:end);
ni = yv(:,1+2*size(yv,2)/3:size(yv,2));

ang = 0:2*pi/1000:2*pi;
xcoor = (R_Pattern)*cos(ang);
ycoor = (R_Pattern)*sin(ang);
xicoor = (RI_Pattern)*cos(ang);
yicoor = (RI_Pattern)*sin(ang);

pointsx = yv(:,1:size(yv,2)/3);
pointsy = yv(:,1+size(yv,2)/3:size(yv,2)/3*2);

[theta,r] = cart2pol(pointsx,pointsy);

%% Set up the movie.
F(1:size(pointsx,1)) = struct('cdata',[],'colormap',[]);
v = VideoWriter('CellsMove.avi');
open(v)
for i = 1:size(pointsx,1)
    xlim manual
    ylim manual
    plot(xcoor,ycoor,'Linewidth',5)
    axis equal
    hold on
    plot(xicoor,yicoor,'Linewidth',5)
    axis equal
    hold on
%     DT = delaunayTriangulation(pointsx(:,i),pointsy(:,i));
%     triplot(DT)
%     plot(pointsx(i,:),pointsy(i,:),'o')
%     voronoi(pointsx(:,i),pointsy(:,i))
    plot(pointsx(i,:),pointsy(i,:),'o')
    
    F = getframe(gcf);
    writeVideo(v,F)
    xlim([-R_Pattern-5,R_Pattern+5])
    ylim([-R_Pattern-5,R_Pattern+5])
    hold off
end
close(v)

F(1:size(pointsx,1)) = struct('cdata',[],'colormap',[]);
v = VideoWriter('QuiverVelocity.avi');
open(v)
for i = 1:size(pointsx,1)-1
    xlim manual
    ylim manual
    plot(xcoor,ycoor,'Linewidth',5)
    hold on
    plot(xicoor,yicoor,'Linewidth',5)
    axis equal
    hold on
    quiver(pointsx(i,:),pointsy(i,:),cos(ni(i,:)+theta(i,:)),sin(ni(i,:)+theta(i,:)))
    hold on
    quiver(pointsx(i,:),pointsy(i,:),vel(i,1:size(vel,2)/3),vel(i,size(vel,2)/3+1:2*size(vel,2)/3))
    hold off
    F = getframe(gcf);
    writeVideo(v,F)
end
close(v)
v = VideoWriter('VeronoiQuiver.avi');
open(v)
pointsx(isnan(pointsx)) = 0;
pointsy(isnan(pointsy)) = 0;
for i = 1:size(pointsx,1)-1
    plot(xcoor,ycoor,'Linewidth',5); hold on; xlim manual;
    xlim([-R_Pattern-5,R_Pattern+5]);ylim([-R_Pattern-5,R_Pattern+5]);axis equal;axis off;hold off;
    ylim manual;
    hold on
    plot(xicoor,yicoor,'Linewidth',5)
    axis equal
    voronoi(pointsx(i,:),pointsy(i,:))
    hold on
    pos = [-RI_Pattern -RI_Pattern 2*RI_Pattern 2*RI_Pattern]; 
    rectangle('Position',pos,'FaceColor',[1 1 1],'Curvature',[1 1])
    axis equal
    hold on
    fillout(xcoor,ycoor,[-R_Pattern-5,R_Pattern+5,-R_Pattern-5,R_Pattern+5],'w')
    hold on
    quiver(pointsx(i,:),pointsy(i,:),cos(ni(i,:)+theta(i,:)),sin(ni(i,:)+theta(i,:)))
    hold off
    axis off
    xlim([-R_Pattern-5,R_Pattern+5]);ylim([-R_Pattern-5,R_Pattern+5]);axis equal;axis off;hold off;
        F = getframe(gcf);
    writeVideo(v,F)
end
close(v)