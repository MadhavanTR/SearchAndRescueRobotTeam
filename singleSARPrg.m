clear all; close all; clc
%% SAR Operation Code adapted from the following Leader-follower Code

% 3-mode Search and Rescue Robot Team
% edited by Madhavan
% 11/2019

%% Experiment Constants

% Run the simulation for a specific number of iterations
iterations = 5000;
radius = 0.33;  % 0.25
% objectResc = randn(1, 2);    %(0.5, 0.5), (0.4, -0.4), (0.4, 0.6)
% if( objectResc(1)<=-2 || objectResc(1)>= 1)
%     objectResc(1)=objectResc(1)/2;
%     if( objectResc(1)<=-1.5 || objectResc(1)>= 1)
%         objectResc(1)=objectResc(1)/2;
%     end
% end
% if( objectResc(2)<=-0.7 || objectResc(2)>= 0.7)
%     objectResc(2)=objectResc(2)/2;
%     if( objectResc(2)<=-1 || objectResc(2)>= 1)
%         objectResc(2)=objectResc(2)/2;
%     end
% end
% objectResc = round(objectResc, 2);
objectResc = [-0.6 0.5];
safeGoal=[-1.2, -0.6];  %(-1.3, -0.7)

%% Set up the Robotarium object

N = 9;
dxiTemp=zeros(2, N-1);
dxiTempNew=zeros(2, N);
initial_positions = generate_initial_conditions(N, 'Width', 1, 'Height', 1, 'Spacing', 0.2);
initial_positions(1, :) = initial_positions(1, :)-1.2;
%initial_positions(2, :) = initial_positions(2, :)-0.5;
[px, py] = lloydsAlgorithm(0.01*rand(N,1),0.01*rand(N,1), [0, 0; 0, 1.8; 2.8, 1.8; 2.8, 0], 50, false);   %[(0, 0), (0, 1), (1, 1), (1, 0)]
goToPoints = [px-1.4, py-0.9];
r = Robotarium('NumberOfRobots', N, 'ShowFigure', true, 'InitialConditions', initial_positions);

%% Initialize for video recording

% vid = VideoWriter('SARExp4_3.mp4', 'MPEG-4');
% vid.Quality = 100;
% vid.FrameRate = 72;
% open(vid);
% writeVideo(vid, getframe(gcf));

%% Create the desired Laplacian

% Inital Graph laplacian, we can change it later depending on the leader
followers = -completeGL(N-1);
L = zeros(N, N);
L(2:N, 2:N) = followers;
L(2, 2) = L(2, 2) + 1;  % 2nd alone is connected to the leader
L(2, 1) = -1;

% Initialize velocity vector
dxi = zeros(2, N);
dxiNew = zeros(2, N);
dxiCopy=dxi;

% Inital assumed state for the system
state = 1;

% These are gains for our formation control algorithm
formation_control_gain = 5; %Originally 10
desired_distance = 0.3;

%% Grab tools we need to convert from single-integrator to unicycle dynamics

% Single-integrator -> unicycle dynamics mapping
si_to_uni_dyn = create_si_to_uni_dynamics('LinearVelocityGain', 0.8);

% Single-integrator barrier certificates
uni_barrier_cert = create_uni_barrier_certificate_with_boundary();

% Single-integrator position controller
leader_controller = create_si_position_controller('XVelocityGain', 0.8, 'YVelocityGain', 0.8, 'VelocityMagnitudeLimit', 0.1);

waypoints = [-1 0.8; -1 -0.8; objectResc; safeGoal]';    %waypoints(3, :) - center of the object to be rescued
                                                         %waypoints(4, :) - safe goal
close_enough = 0.01;   %Originally 0.05

for t = 1:iterations
    
    % Retrieve the most recent poses from the Robotarium.  The time delay is
    % approximately 0.033 seconds
    x = r.get_poses();
    
    %% Algorithm
    
    for i = 1:N

        %Zero velocity and get the topological neighbors of agent i
        if  state ~= 1
            dxi(:, i) = [0 ; 0];
            neighbors = topological_neighbors(L, i);
            for j = neighbors
                dxi(:, i) = dxi(:, i) + ...
                    formation_control_gain*(norm(x(1:2, j) - x(1:2, i))^2 -  desired_distance^2)*(x(1:2, j) - x(1:2, i));
                    %fprintf('2nd bot %d in dx and %d in dy\n', dxi(1, 2), dxi(2, 2));
            end
        end
    end
    
    %% Make the leader travel between waypoints
    
    waypoint = waypoints(:, state);
    if state == 4
        distanceChange=waypoints(:, 4)-waypoints(:, 3);
        waypoint = [waypoints(1, 3)+0.1*distanceChange(1); waypoints(2, 3)+0.1*distanceChange(2)]; % INCREASE AND DECREASE ACCORDING TO GOAL
        waypoints(:, 3) = waypoint;
        if norm(distanceChange)>close_enough
            %fprintf('%d', distanceChange);
            state=3;
        else
            break;
        end
    end
    switch state        
        case 1
            % Write Voronoi here 
            %goToPoints=[linspace(-1.3, 1.3, N); linspace(-0.6, 0.6, N)]'; %(+-1.2, +-0.5)     
            for i=1:N
                dxiNew(:, i) = leader_controller(x(1:2, i), goToPoints(i, :)');
            end
            fprintf('\nMode %d - Exploration', 1);
            if (round(x(1:2, :)-goToPoints', 3)<0.05)
                state = 2;
                distVect=x(1:2, :)-objectResc';
                for i=1:N
                    distVal(i) = norm(distVect(:, i));
                end
                [~, leader_agent] = min(distVal);
            end
            
        case 2
            %  Come to leader
            L = findLaplacian(N, followers, leader_agent);
            waypoint = x(1:2, leader_agent);
            dxi(:, leader_agent) = leader_controller(x(1:2, leader_agent), waypoint); %   Instead of this, leader will have reached the position, and dxi(:, 1) = 0
            %fprintf('In 2 in %d th iteration\n', t);
            fprintf('\nMode 2 - Leader follower');
            dxiCopy=dxi;
            dxiCopy(:, leader_agent)=[];
            if (round(dxiTemp-dxiCopy(:, 1:N-1), 3) == 0)            %abs(dxi(:, 2:N)) < 0.2        %USE THIS, IF ALL BOTS HAVE REACHED LEADER AND NEARLY AT REST, SWITCH TO FORMATION
                state = 3;
            end
            dxiTemp = dxiCopy(:, 1:N-1);

        case 3
             fprintf('\nState 3 - Circle formation around the object and bringing it back');               
             xi = x(1:2,:);
             circularTargets = [ waypoint(1) + radius*cos( 0:2*pi/N:2*pi*(1- 1/N) ) ; waypoint(2) + radius*sin( 0:2*pi/N:2*pi*(1- 1/N) ) ];
             errorToInitialPos = xi - circularTargets;                % Error
             errorNorm = [1,1]*(errorToInitialPos.^2);                % Norm of error
             while max( errorNorm ) > 0.05
                % Update state variables   
                r.step();
                x = r.get_poses();                                    % States of real unicycle robots
                xi = x(1:2,:);                                        % x-y positions
                % Update errors
                errorToInitialPos = xi - circularTargets;
                errorNorm = [1,1]*(errorToInitialPos.^2);
                % Conput control inputs
                dxi = -0.5.*errorToInitialPos;
                dxu = si_to_uni_dyn(dxi, x);
                dxu = uni_barrier_cert(dxu, x);
                % Assing new control inputs to robots
                r.set_velocities(1:N, dxu);                       % Assign dummy zero velocity
                % writeVideo(vid, getframe(gcf));
             end
             state = 4;
             clear errorNorm circularTargets errorToInitialPos
    end

%% Avoid actuator errors

    if(state==1)
        dxi=dxiNew;
    end
    % To avoid errors, we need to threshold dxi
    norms = arrayfun(@(x) norm(dxi(:, x)), 1:N);
    threshold = 3/4*r.max_linear_velocity;
    to_thresh = norms > threshold;
    dxi(:, to_thresh) = threshold*dxi(:, to_thresh)./norms(to_thresh);
    
    %% Use barrier certificate and convert to unicycle dynamics
    
    dxu = si_to_uni_dyn(dxi, x);
    dxu = uni_barrier_cert(dxu, x);
    
    %% Send velocities to agents
    
    %Set velocities
    r.set_velocities(1:N, dxu);
    
    %Iterate experiment
    r.step();
    
    %writeVideo(vid, getframe(gcf)); % Record a video frame every 10 iterations
    
end

% close(vid);
% We can call this function to debug our experiment!  Fix all the errors
% before submitting to maximize the chance that your experiment runs
% successfully.
r.debug();

function L_mat= findLaplacian(N, followers, leader_agent)
    L = zeros(N, N);
    if(leader_agent==1)
        L(2:N, 2:N) = followers;
        L(2, 2) = L(2, 2) + 1;  % 2nd alone is connected to da leader
        L(2, 1) = -1;
    elseif(leader_agent==N)
        L(1:N-1, 1:N-1) = followers;
        L(N-1, N-1) = L(N-1, N-1) + 1;
        L(N-1, N) = -1;
    else
        L(1:leader_agent-1, 1:leader_agent-1) = followers(1:leader_agent-1, 1:leader_agent-1);
        L(leader_agent, :) = zeros(1, N);
        L(leader_agent+1, leader_agent) = -1;
        L(leader_agent+1:N, leader_agent+1:N) = followers(leader_agent:N-1, leader_agent:N-1);
        L(1:leader_agent-1, leader_agent+1:N) = followers(1:leader_agent-1, leader_agent:N-1);
        L(leader_agent+1:N, 1:leader_agent-1) = followers(leader_agent:N-1, 1:leader_agent-1);
        L(leader_agent+1, leader_agent+1)=L(leader_agent+1, leader_agent+1)+1;
    end
    L_mat=L;
end

function [ poses ] = generate_initial_conditions(N, varargin)
    
    poses = zeros(3, N);
    
    parser = inputParser;
    parser.addParameter('Spacing', 0.3);
    parser.addParameter('Width', 3.0);
    parser.addParameter('Height', 1.8);
    parse(parser, varargin{:});
    
    spacing = parser.Results.Spacing;
    width = parser.Results.Width;
    height = parser.Results.Height;

    numX = floor(width / spacing);
    numY = floor(height / spacing);
    values = randperm(numX * numY, N);

    for i = 1:N
       [x, y] = ind2sub([numX numY], values(i));
       x = x*spacing - (width/2); 
       y = y*spacing - (height/2);
       poses(1:2, i) = [x ; y];
    end
    
    poses(3, :) = (rand(1, N)*2*pi - pi);
end

