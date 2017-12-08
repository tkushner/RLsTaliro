% demo for pendulum

clear

g=1;
l=1;

%model = @(t,x,u) [x(2); -g/l * sin(x(1))];
model = @(t,x,u) [-g/l*sin(x(1))];

init_cond = [-pi, pi];
input_range = [];
cp_array = [];

phi = '[]a';
%phi = '<>_[1,5) a';

ii = 1;
preds(ii).str='a';
preds(ii).A = [1];
preds(ii).b = [-pi/2];

time = 10;

opt = staliro_options();

opt.runs = 1;
opt.spec_space = 'X';
opt.ode_solver = 'ode45';
opt.falsification=0;
opt.optimization_solver = 'UR_Taliro';
opt.optim_params.n_tests = 100;
[results, history] = staliro(model,init_cond,input_range,cp_array,phi,preds,time,opt);

% Get Falsifying trajectory
bestRun = results.optRobIndex;
[T1,XT1] = SimFunctionMdl(model,init_cond,input_range,cp_array,results.run(bestRun).bestSample,time,opt);

figure(1)
clf
rectangle('Position',[-1.6,-1.1,0.2,0.2],'FaceColor','r')
hold on
% if (init_cond(1,1)==init_cond(1,2)) || (init_cond(2,1)==init_cond(2,2))
%     plot(init_cond(1,:),init_cond(2,:),'g')
% else
%     rectangle('Position',[init_cond(1,1),init_cond(2,1),init_cond(1,2)-init_cond(1,1),init_cond(2,2)-init_cond(2,1)],'FaceColor','g')
% end
ntests = results.run(bestRun).nTests;
hist = history(bestRun).samples;
%plot(hist(1:ntests,1),hist(1:ntests,2),'*')
plot(cos(hist(1:ntests,1)),sin(hist(1:ntests,1)),'*')
%plot(XT1(:,1),XT1(:,2))

xlabel('y_1')
ylabel('y_2')

%% plot the pendulum
pendulumtopx = 0;
pendulumtopy = l;

initcond=results.run(bestRun).bestSample;
tend=10;

%deq1=@(t,x) [x(2); -g/l * sin(x(1))]; % Pendulum equations uncoupled
deq1=@(t,x) [-g/l * sin(x(1))];
%[t,sol] = ode45(deq1,[0 tend],[cos(initcond) sin(initcond)]);  % uses a numerical ode solver
[t, sol]=ode45(deq1, [0,tend], initcond);
% solx1 = sol(:,1)'; % takes the transpose for plots
% solx2 = sol(:,2)';

iterations = 1; % Sets initial iteration count to 1
pausetime = 0.1;  % Pauses animation for this time
runtime = tend;  % Runs simulations for this time
tx = 0;  
arraysize = size(t);  % Defines array size of time intervals
timestep = t(runtime) - t(runtime-1); 

% cartesianx = l*sin(solx1);  % Converts angles into cartesian coordinates
% cartesiany = l*cos(solx2); 
cartesianx=l*sin(sol);
cartesiany=l*cos(sol);

for i = 1 : max(arraysize)
    figure(2)
    plotarrayx = [pendulumtopx cartesianx(iterations)];
    plotarrayy = [pendulumtopy cartesiany(iterations)];
    plot(cartesianx(iterations),cartesiany(iterations),'ko',plotarrayx,plotarrayy,'r-')
    %title(['Simple pendulum simulation            \theta = ' num2str(solx1(iterations))],'fontsize',12)
    title(['simple pendulum simulation        \theta =' num2str(sol(iterations))],'fontsize',16)
    xlabel('x [m]','fontsize',12)
    ylabel('y [m]','fontsize',12)
    axis([-1 1 -1 2])
    pause(pausetime)  % Shows results at each time interval
    iterations = iterations + 1;  % increases iteration count by 1
end
