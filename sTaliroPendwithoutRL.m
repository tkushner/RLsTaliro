function sTaliroPendwithoutRL
%% STaliro without RL
% Things that are weird
% (1) function - I think rewardFunc is the right thought, except we want to
% go "bigger" - like we want the whole pendulum with the pushing aspect too
% (like the function has to be state space & time, and right now its just
% reward value in space) but we want the whole pendulum. (we can do the
% simulink model of the pendulum? bc staliro takes simulink models and it
% seems straightforward enough to make the pendulum in simulink - easy
% google how to)
% (2) cparray - control points that parameterize the input signal. So our
% input is "push left, push right, don't push"  my two guesses are either 
%       (a) array[n,1] where n is the number of times we let an action
%       happen within our time window and 1 because we have 3 different
%       inputs (except I think left is negative, right is positive and no
%       move is zero?) 
%       (b) just array[1] because you can only do one action? (ie push)
% (3) phi - Ok so phi is woreded in terms of falsification which we want to
% find phi = '!(<>[0,0], x)' (except that instead of phi is bad, we want
% phi is good so I think if we just take that phi and then reverse our
% output map, it'll work?)
% (4) preds - ughhh this is trickier... its like, the mapping of positions
% x to their state in space? it has to do with phi, so if we say 'x' in
% phi, i think we need to be like 'this is what x actually is in our model'
% in the preds part (type <help dp_taliro> in command line to read more
% about it)
%       

close all;

%% SETTINGS

%%% What do we call good?
rewardFunc = @(t, y, x)(-(abs(x(1))).^2 + -0.25*(abs(x(2))).^2); % Reward is -(quadratic error) from upright position. Play around with different things!

%%% Confidence in new trials?
learnRate = 0.99; % How is new value estimate weighted against the old (0-1). 1 means all new and is ok for no noise situations.

%%% Exploration vs. exploitation
% Probability of picking random action vs estimated best action
epsilon = 0.5; % Initial value
epsilonDecay = 0.98; % Decay factor per iteration.

%%% Future vs present value
discount = 0.9; % When assessing the value of a state & action, how important is the value of the future states?

%%% Inject some noise?
successRate = 1; % How often do we do what we intend to do?
% E.g. What if I'm trying to turn left, but noise causes
% me to turn right instead. This probability (0-1) lets us
% try to learn policies robust to "accidentally" doing the
% wrong action sometimes.

winBonus = 100;  % Option to give a very large bonus when the system reaches the desired state (pendulum upright).

startPt = [pi, 0]; % Start every episode at vertical down.
input_range = [-pi pi; -pi pi]; % Inputs are bound between these values.

maxEpi = 2000; % Each episode is starting with the pendulum down and doing continuous actions for awhile.
maxit = 1500; % Iterations are the number of actions taken in an episode.
substeps = 2; % Number of physics steps per iteration (could be 1, but more is a little better integration of the dynamics)
dt = 0.05; % Timestep of integration. Each substep lasts this long

% Torque limits -- bang-bang control
tLim = 1;
actions = [0, -tLim, tLim]; % Only 3 options, Full blast one way, the other way, and off.
cp_array = [3; 3]; % Number of 'control points' for each input. Not sure about this.

% Create an S-Taliro options object
opt = staliro_options();
opt.optimization_solver = 'SA_Taliro';
opt.falsification = 0;
opt.parameterEstimation = 1;
opt.runs = 1;

ii = 1;
preds(ii).par = 'r1';
preds(ii).b = 0;
preds(ii).value = 0;
preds(ii).range = [-pi pi];

ii = ii+1;
preds(ii).par = 'r2';
preds(ii).b = 0;
preds(ii).value = 0;
preds(ii).range = [-pi pi];

phi = ['<>(r1 -> r2)'; '<>(r1 /\ r2)'];

% Make the un-updated values on the value map transparent. If not, then
% we see the reward function underneath.
transpMap = true;

% Write to video?
doVid = false;

if doVid
    writerObj = VideoWriter('qlearnVid.mp4','MPEG-4');
    writerObj.FrameRate = 60;
    open(writerObj);
end

%% Discretize the state so we can start to learn a value map
% State1 is angle -- play with these for better results. Faster convergence
% with rough discretization, less jittery with fine.
x1 = -pi:0.05:pi;
%State2 angular rate
x2 = -pi:0.1:pi;

%Generate a state list
states=zeros(length(x1)*length(x2),2); % 2 Column matrix of all possible combinations of the discretized state.
index=1;
for j=1:length(x1)
    for k = 1:length(x2)
        states(index,1)=x1(j);
        states(index,2)=x2(k);
        index=index+1;
    end
end

% Local value R and global value Q -- A very important distinction!
%
% R, the reward, is like a cost function. It is good to be near our goal. It doesn't
% account for actions we can/can't take. We use quadratic difference from the top.
%
% Q is the value of a state + action. We have to learn this at a global
% level. Once Q converges, our policy is to always select the action with the
% best value at every state.
%
% Note that R and Q can be very, very different.
% For example, if the pendulum is near the top, but the motor isn't
% powerful enough to get it to the goal, then we have to "pump" the pendulum to
% get enough energy to swing up. In this case, that point near the top
% might have a very good reward value, but a very bad set of Q values since
% many, many more actions are required to get to the goal.

% Initialize the "cost" of a given state to be quadratic error from the goal state. Note the signs mean that -angle with +velocity is better than -angle with -velocity
R = zeros(length(x1)*length(x2), 1);
index=1;
for j=1:length(x1)
    for k = 1:length(x2)
        R(index,1) = rewardFunc(0, 0, [states(index, 1), states(index, 2)]);
        index=index+1;
    end
end

%R = rewardFunc(0, 0, states(:,1), states(:,2)); % Initialize the "cost" of a given state to be quadratic error from the goal state. Note the signs mean that -angle with +velocity is better than -angle with -velocity
Q = repmat(R,[1,3]); % Q is length(x1) x length(x2) x length(actions) - IE a bin for every action-state combination.

% V will be the best of the Q actions at every state. This is only
% important for my plotting. Not algorithmically significant. I leave
% values we haven't updated at 0, so they appear as blue (unexplored) areas
% in the plot.
V = zeros(size(states,1),1);
Vorig = reshape(max(Q,[],2),[length(x2),length(x1)]);

%% Set up the pendulum plot
panel = figure;
panel.Position = [680 558 1000 400];
panel.Color = [1 1 1];
subplot(1,4,1)

hold on
% Axis for the pendulum animation
f = plot(0,0,'b','LineWidth',10); % Pendulum stick
axPend = f.Parent;
axPend.XTick = []; % No axis stuff to see
axPend.YTick = [];
axPend.Visible = 'off';
axPend.Position = [0.01 0.5 0.3 0.3];
axPend.Clipping = 'off';
axis equal
axis([-1.2679 1.2679 -1 1]);
plot(0.001,0,'.k','MarkerSize',50); % Pendulum axis point

hold off

%% Set up the state-value map plot (displays the value of the best action at every point)
colormap('hot');
subplot(1,4,[2:4]);
hold on
map = imagesc(reshape(R,[length(x2),length(x1)]));
axMap = map.Parent;
axMap.XTickLabels = {'-pi' '0' 'pi'};
axMap.XTick = [1 floor(length(x1)/2) length(x1)];
axMap.YTickLabels = {'-pi' '0' 'pi'};
axMap.YTick = [1 floor(length(x2)/2) length(x2)];
axMap.XLabel.String = 'Angle (rad)';
axMap.YLabel.String = 'Angular rate (rad/s)';
axMap.Visible = 'on';
axMap.Color = [0.3 0.3 0.5];
axMap.XLim = [1 length(x1)];
axMap.YLim = [1 length(x2)];
axMap.Box = 'off';
axMap.FontSize = 14;
caxis([3*min(R),max(R)])
pathmap = plot(NaN,NaN,'.g','MarkerSize',30); % The green marker that travels through the state map to match the pendulum animation
map.CData = V;
hold off

%% Staliro
z1=startPt;
results = staliro(rewardFunc, z1, input_range, cp_array, [], preds, 1, opt);

%% Start learning!
% Number of episodes or "resets"
for episodes = 1:maxEpi
    
    z1 = startPt; % Reset the pendulum on new episode.
    
    % Number of actions we're willing to try before a reset
    for g = 1:maxit
        
        % Stop if the figure window is closed.
        if ~ishandle(panel)
            break;
        end
        results = staliro(rewardFunc, z1, input_range, cp_array, phi, preds, 1, opt);
        disp(results)
        %% PICK AN ACTION
        
        % Interpolate the state within our discretization (ONLY for choosing
        % the action. We do not actually change the state by doing this!)
        [~,sIdx] = min(sum((states - repmat(z1,[size(states,1),1])).^2,2));
        
        % Choose an action:
        % EITHER 1) pick the best action according the Q matrix (EXPLOITATION). OR
        % 2) Pick a random action (EXPLORATION)
        if (rand()>epsilon || episodes == maxEpi) && rand()<=successRate % Pick according to the Q-matrix it's the last episode or we succeed with the rand()>epsilon check. Fail the check if our action doesn't succeed (i.e. simulating noise)
            [~,aIdx] = max(Q(sIdx,:)); % Pick the action the Q matrix thinks is best!
        else
            aIdx = randi(length(actions),1); % Random action!
        end
        
        T = actions(aIdx);
        
        %% STEP DYNAMICS FORWARD
        
        % Step the dynamics forward with our new action choice
        % RK4 Loop - Numerical integration
        for i = 1:substeps
            k1 = Dynamics(z1,T);
            k2 = Dynamics(z1+dt/2*k1,T);
            k3 = Dynamics(z1+dt/2*k2,T);
            k4 = Dynamics(z1+dt*k3,T);
            
            z2 = z1 + dt/6*(k1 + 2*k2 + 2*k3 + k4);
            % All states wrapped to 2pi
            if z2(1)>pi
                z2(1) = -pi + (z2(1)-pi);
            elseif z2(1)<-pi
                z2(1) = pi - (-pi - z2(1));
            end
        end
        
        z1 = z2; % Old state = new state
        
        
        %% UPDATE Q-MATRIX
        
        % End condition for an episode
        if norm(z2)<0.01 % If we've reached upright with no velocity (within some margin), end this episode.
            success = true;
            bonus = winBonus; % Give a bonus for getting there.
        else
            bonus = 0;
            success = false;
        end
        
        [~,snewIdx] = min(sum((states - repmat(z1,[size(states,1),1])).^2,2)); % Interpolate again to find the new state the system is closest to.
        
        if episodes ~= maxEpi % On the last iteration, stop learning and just execute. Otherwise...
            % Update Q
            Q(sIdx,aIdx) = Q(sIdx,aIdx) + learnRate * ( R(snewIdx) + discount*max(Q(snewIdx,:)) - Q(sIdx,aIdx) + bonus );
            
            % Lets break this down:
            %
            % We want to update our estimate of the global value of being
            % at our previous state s and taking action a. We have just
            % tried this action, so we have some information. Here are the terms:
            %   1) Q(sIdx,aIdx) AND later -Q(sIdx,aIdx) -- Means we're
            %      doing a weighting of old and new (see 2). Rewritten:
            %      (1-alpha)*Qold + alpha*newStuff
            %   2) learnRate * ( ... ) -- Scaling factor for our update.
            %      High learnRate means that new information has great weight.
            %      Low learnRate means that old information is more important.
            %   3) R(snewIdx) -- the reward for getting to this new state
            %   4) discount * max(Q(snewIdx,:)) -- The estimated value of
            %      the best action at the new state. The discount means future
            %      value is worth less than present value
            %   5) Bonus - I choose to give a big boost of it's reached the
            %      goal state. Optional and not really conventional.
        end
        
        % Decay the odds of picking a random action vs picking the
        % estimated "best" action. I.e. we're becoming more confident in
        % our learned Q.
        epsilon = epsilon*epsilonDecay;
        
        %% UPDATE PLOTS
        
        if episodes>0
            % Pendulum state:
            set(f,'XData',[0 -sin(z1(1))]);
            set(f,'YData',[0 cos(z1(1))]);
            
            % Green tracer point:
            [newy,newx] = ind2sub([length(x2),length(x1)],snewIdx); % Find the 2d index of the 1d state index we found above
            set(pathmap,'XData',newx);
            set(pathmap,'YData',newy);
            
            % The heat map of best Q values
            V = max(Q,[],2); % Best estimated value for all actions at each state.
            fullV = reshape(V,[length(x2),length(x1)]); % Make into 2D for plotting instead of a vector.
            set(map,'CData',fullV);
            if transpMap
                set(map,'AlphaData',fullV~=Vorig); % Some spots have not changed from original. If not, leave them transparent.
            end
            drawnow;
            
            % Take a video frame if turned on.
            if doVid
                frame = getframe(panel);
                writeVideo(writerObj,frame);
            end
        end
        
        % End this episode if we've hit the goal point (upright pendulum).
        if success
            break;
        end
        
    end
end

if doVid
    close(writerObj);
end

end

function zdot = Dynamics(z,T)
% Pendulum with motor at the joint dynamics. IN - [angle,rate] & torque.
% OUT - [rate,accel]
g = 1;
L = 1;
z = z';
zdot = [z(2) g/L*sin(z(1))+T];
end