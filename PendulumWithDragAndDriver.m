%pendulum with driver input and drag
%the velocities are what is being controlled in the robustness function


omega = 1;      % resonant frequency = sqrt(k/m)
drag = 0.2;        % drag coeficient per unit mass
Amp = 0.1;        % driving amplitude per unit mass

a = omega^2;    % calculate a coeficient from resonant frequency


g=1;
l=1;

model= @(t,x,u) [x(2); -a * x(1) - drag * x(2) + Amp * sin(1.2*t)];

init_cond = [-pi,pi;-5,5];
input_range = [];
cp_array = [];

phi = '[]!(a\/b)';
%phi = '<>_[1,5) a';

ii = 1;
preds(ii).str='a';
preds(ii).A = [0,-1];
preds(ii).b = [-2];

ii = 2;
preds(ii).str='b';
preds(ii).A = [0,-1];
preds(ii).b = [2];

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
if (init_cond(1,1)==init_cond(1,2)) || (init_cond(2,1)==init_cond(2,2))
    plot(init_cond(1,:),init_cond(2,:),'g')
else
    rectangle('Position',[init_cond(1,1),init_cond(2,1),init_cond(1,2)-init_cond(1,1),init_cond(2,2)-init_cond(2,1)],'FaceColor','g')
end
ntests = results.run(bestRun).nTests;
hist = history(bestRun).samples;
plot(hist(1:ntests,1),hist(1:ntests,2),'*')
%plot(cos(hist(1:ntests,1)),sin(hist(1:ntests,1)),'*')
plot(XT1(:,1),XT1(:,2))

xlabel('y_1')
ylabel('y_2')

%%

tBegin = 0;     % time begin
tEnd = 200;      % time end

[val loc]= max(history.rob)
[initcond]=history.samples;
x0 = history.samples(loc,1);       % initial position
v0 = history.samples(loc,2);       % initial velocitie

a = omega^2;    % calculate a coeficient from resonant frequency

% Use Runge-Kutta 45 integrator to solve the ODE
[t,w] = ode45(model, [tBegin tEnd], [x0 v0]);
x = w(:,1);     % extract positions from first column of w matrix
v = w(:,2);     % extract velocities from second column of w matrix


cartesianx=sin(x);
cartesiany=-cos(x);
for i = 1 : tEnd
    figure(2)
    subplot(2,1,1)
    plotarrayx = [0 cartesianx(i)];
    plotarrayy = [0 cartesiany(i)];
    plot(cartesianx(i),cartesiany(i),'ko',plotarrayx,plotarrayy,'r-')
    %title(['Simple pendulum simulation            \theta = ' num2str(solx1(iterations))],'fontsize',12)
    title('simple pendulum simulation','fontsize',16)
    xlabel('x [m]','fontsize',12)
    ylabel('y [m]','fontsize',12)
    axis([-1 1 -1 1])
    
    subplot(2,1,2)
        plot(i,v(i),'bo')
        title('velocity of penndulum','fontsize',12)
        xlabel('t [seconds]','fontsize',12)
        ylabel('velocity','fontsize',12)
        hold on  % Holds previous values
        axis([0 i+1 min(v)-1 max(v)+1])
    pause(.1)  % Shows results at each time interval
end
