% demo for pendulum

clear

g=1;
l=1;

model = @(t,x,u) [x(2); -g/l * sin(x(1))];

init_cond = [-pi/2, pi/2; -pi/2, pi/2];
input_range = [];
cp_array = [];

phi = '[]!a';
%phi = '<>_[1,5) a';

ii = 1;
preds(ii).str='a';
preds(ii).A = [1,0;0,1];
preds(ii).b = [-pi;-pi];

time = 2;

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
plot(XT1(:,1),XT1(:,2))
xlabel('y_1')
ylabel('y_2')


