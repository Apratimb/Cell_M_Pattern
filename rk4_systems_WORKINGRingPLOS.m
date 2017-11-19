function [sol,vel] = rk4_systems_WORKINGRingPLOS(a, b, N, alpha,R_Pattern,RI_Pattern,a0)

m = size(alpha,1);
if m == 1
   alpha = alpha';
end

h = (b-a)/N;
t(1) = a;
w(:,1) = alpha;
vel = [];
for i = 1:N
   k1 = h*f(t(i), w(:,i),R_Pattern,R_Pattern/RI_Pattern,a0);
   k2 = h*f(t(i)+h/2, w(:,i)+0.5*k1,R_Pattern,R_Pattern/RI_Pattern,a0);
   k3 = h*f(t(i)+h/2, w(:,i)+0.5*k2,R_Pattern,R_Pattern/RI_Pattern,a0); 
   k4 = h*f(t(i)+h, w(:,i)+k3,R_Pattern,R_Pattern/RI_Pattern,a0);
   w(:,i+1) = w(:,i) + (k1 + 2*k2 + 2*k3 + k4)/6;
   t(i+1) = a + i*h;
   vel = cat(2,vel,f(t(i), w(:,i),R_Pattern,R_Pattern/RI_Pattern,a0));
end


sol = [t' w'];


function dy = f(t,y,R_Pattern,rO2i,a0)
% a0 = ones(size(a0))*mean(a0);
v0 = 1; %   Cell velocity
R0 = 1; %   Radius of sphere of influence
Req = 5/6; %    Radius of cell
Fadh = 1; %  Parameter of Adhesie force
Frep = 3; %    Parameter od Repulsive force
F_Wall = 7;
u = 1;  %   Mobility
T = 1;  %   Relaxation time

X0 = 0; %   X_Coordinate of pattern center
Y0 = 0; %   Y_Coordinate of the pattern center
x = y(1:size(y,1)/3,1);
y1 = y(size(y,1)/3+1:size(y,1)/3*2,1);

points = [x,y1];
x(isinf(x)) = 0;
y1(isnan(y1)) = 0;
x(isnan(x)) = 0;
y1(isinf(y1)) = 0;
ni = y(2*size(y,1)/3+1:size(y,1),1);
dij = zeros(size(points,1),size(points,1));
for i = 1:size(points,1)
    for j = 1:size(points,1)
        dij(i,j) = sqrt((x(i,:) - x(j,:))^2+(y1(i,:)-y1(j,:))^2);
    end
end

DT = delaunayTriangulation(x(:,1),y1(:,1));
xdt = DT.Points(:,1);
ydt = DT.Points(:,2);

% Get all the triangles attached to a particular vertex in the
% triangulation.  
attachedTriangles = vertexAttachments(DT);
for i = 1: size(xdt,1)
% Use the connectivity list to get the vertex indices of all these
% triangles
    verticesOfDT = DT.ConnectivityList(attachedTriangles{i},:);
% Find all the unique vertices and remove the current vertex
    neighboursVERTEXDT{i} = setdiff(unique(verticesOfDT), i);
end

Forces = zeros(size(points));
for vertex = 1:size(DT.Points,1)
    neighboring_vertex = neighboursVERTEXDT{vertex}';
    temp = 0; 
    ri = [points(vertex,1)-X0,points(vertex,2)-Y0]; %  r of the vertex
    n = 0;
    for i = 1:length(neighboring_vertex) 
        rj = [points(neighboring_vertex(i),1)-X0,points(neighboring_vertex(i),2)-Y0]; % r of the neighboring vertex
        eij = (rj-ri);  
        mod_rij = abs((eij(1,1)^2+eij(1,2)^2)^0.5);
        n = n+1;
        if or(mod_rij - 2.3*a0(i) > 0,i == j)
            Forces(vertex,:) = Forces(vertex,:) + [0,0];
        elseif and(mod_rij - a0(i) > 0,mod_rij - a0(i) <1.3*a0(i))
            Forces(vertex,:) = Forces(vertex,:) + Fadh*(mod_rij - a0(i))*eij/mod_rij;
        else
            Forces(vertex,:) = Forces(vertex,:) + Frep*(mod_rij - a0(i))*eij/mod_rij;
        end
    end
end


[theta,r] = cart2pol(real(x),real(y1));

points = [x,y1];
Forceijx = Forces(:,1); 
Forceijy = Forces(:,2); 
Forceix = sum(Forceijx,2);
Forceiy = sum(Forceijy,2);


%%  Getting the wall force
for i = 1:size(x,1)
    if R_Pattern - abs(r(i)) < R0
        Fwall(i,:) = -F_Wall*exp(-2*(R_Pattern - abs(r(i))/R0))*[cos(theta(i)),sin(theta(i))];
    elseif and(-R_Pattern/rO2i + abs(r(i)) < R0,~(isinf(rO2i)))
        Fwall(i,:) = F_Wall*exp(-2*(-R_Pattern/rO2i + abs(r(i))/R0))*[cos(theta(i)),sin(theta(i))];
    else
        Fwall(i,:) = [0,0];
    end
end

%% Total force
Forceix = Forceix+Fwall(:,1);
Forceiy = Forceiy+Fwall(:,2);


%% Solving the differential equations
dy = zeros(3*size(x,1),1);
for i = 1:size(dy,1)/3
    dy(i) = v0*cos(ni(i,1)+theta(i))+Forceix(i);
end
for i = size(dy,1)/3+1:2*size(dy,1)/3
    dy(i) = v0*sin(ni(i-size(y,1)/3,1)+theta(i-size(y,1)/3))+Forceiy(i-size(y,1)/3);
end
for i = 2*size(dy,1)/3+1:size(dy,1)
   dy(i) = asin(0.1*(cos(ni(i-2*size(dy,1)/3)+theta(i-2*size(dy,1)/3))*dy(i-size(dy,1)/3)-sin(ni(i-2*size(dy,1)/3)+theta(i-2*size(dy,1)/3))*dy(i-2*size(dy,1)/3))/(dy(i-2*size(dy,1)/3)^2+(dy(i-size(dy,1)/3))^2)^0.5);
end

