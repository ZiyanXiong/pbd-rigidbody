% Define transformation matrices E1 and E2
E1 = eye(4);  % Example transformation matrix for the first cuboid
E2 = eye(4);  % Example transformation matrix for the second cuboid
%{
E1 = [0.999912777360948	6.44689071252614e-05	0.0132073280444948	0.100000000000000;
6.44689071252614e-05	0.999952349068647	-0.00976193811984401	0;
-0.0132073280444948	0.00976193811984401	0.999865126429595	0.496554687500000;
0	0	0	1];

E2 = [1	1.20242028917953e-36	-3.98492288580733e-18	0.200000000000000;
1.20242028917953e-36	1	6.03484847078002e-19	0;
3.98492288580733e-18	-6.03484847078002e-19	1	1.550000000000;
0	0	0	1];
%}


E1(1,4) = 0.65;
E1(3,4) = 1.6;
E1(1:3,1:3) = [ 0.5403023,  0.0000000, -0.8414710;
            0.0000000,  1.0000000,  0.0000000;
            0.8414710,  0.0000000,  0.5403023 ];

E2(3,4) = 0.5;

E2(3,4) = E2(3,4) + 0.05;

%{
E1 = [    0.9999    0.0000    0.0161    0.8119;
   -0.0000    1.0000    0.0001    0.0014;
   -0.0161   -0.0001    0.9999    1.5013;
         0         0         0    1.0000];

E2 = [    0.9999   -0.0001    0.0165    1.2281;
    0.0001    1.0000    0.0000    0.0014;
   -0.0165   -0.0000    0.9999    2.5263;
         0         0         0    1.0000];
%}
% Define cuboid parameters
cuboid1Size = [1, 1, 1]; % Size of the first cuboid [length, width, height]
cuboid2Size = [1, 1, 1]; % Size of the second cuboid [length, width, height]

collisions = odeBoxBox_mex(E2,cuboid1Size,E1,cuboid2Size);

cla;
% Plotting the cuboids
hold on;

% Plotting the first cuboid
plotCuboid(E1, cuboid1Size, 'r');

% Plotting the second cuboid
plotCuboid(E2, cuboid2Size, 'b');

for i = 1 : collisions.count
	x = [collisions.pos(:,i),collisions.pos(:,i)-collisions.nor*0.5];
	plot3(x(1,:),x(2,:),x(3,:),'g-', 'LineWidth',3);
end

x = [[0 0 0]', collisions.nor];
plot3(x(1,:),x(2,:),x(3,:),'-', 'color', [0.4940 0.1840 0.5560], 'LineWidth',3);

hold off;
axis equal;
% Set custom axis limits
xlim([-2, 2]);  % Example limits for x-axis
ylim([-2, 2]);  % Example limits for y-axis
zlim([0, 3]);  % Example limits for z-axis

grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Cuboids');

% Function to plot a cuboid given transformation matrix and size
function plotCuboid(E, size, color)
    % Vertices of the cuboid in the local coordinate system
    vertices = [-1, -1, -1; % Vertex 1
                 1, -1, -1; % Vertex 2
                 1,  1, -1; % Vertex 3
                -1,  1, -1; % Vertex 4
                -1, -1,  1; % Vertex 5
                 1, -1,  1; % Vertex 6
                 1,  1,  1; % Vertex 7
                -1,  1,  1];% Vertex 8
     vertices = vertices * 0.5;
    % Transform the vertices using the transformation matrix E
    transformedVertices = (E * [vertices, ones(8,1)]')';
    
    % Extracting the transformed vertices
    transformedVertices = transformedVertices(:, 1:3);
    
    % Rescaling the transformed vertices by the given size
    scaledVertices = transformedVertices .* size;
    
    % Faces of the cuboid
    faces = [1, 2, 3, 4; % Bottom face
             5, 6, 7, 8; % Top face
             1, 2, 6, 5; % Side face
             2, 3, 7, 6; % Side face
             3, 4, 8, 7; % Side face
             4, 1, 5, 8];% Side face
     
    % Plotting the cuboid
    patch('Vertices', scaledVertices, 'Faces', faces, 'FaceColor', color, 'FaceAlpha', 0.5);
end
