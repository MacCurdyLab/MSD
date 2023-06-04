function [t,X,eV] = odeVc(odeFunc,tspan,X0,dt)

%odeV with CONSTANT timestep

%REQUIRED INPUTS
%odefunc    function handle for ode to be integrated
%tspan      time vector to integrate over, with identical format to ode45
%X0         state vector initial condition, [x1 x1dot x2 x2dot ...]'
%dt         constant timestep for simulation

nspan = ceil(tspan(end)/dt);

%initialize u, v, t matrices for solutions
x = zeros(numel(X0)/2,nspan);
v = x;
t = zeros(nspan,1);

%at timestep n=1, the solution in the intial condition
x(:,1) = X0(1:2:end);
v(:,1) = X0(2:2:end);

%we find the solution at timestep n=2 by assuming constant acceleration
Xdot = feval(odeFunc,0,X0);
A = Xdot(2:2:end);
x(:,2) = x(:,1) + v(:,1)*dt + A*dt^2;
v(:,2) = Xdot(1:2:end);
t(2) = dt;

stopFlag = false;   %flag for halting integration
n=2;                %solution increment

%Remaining timesteps are computed using Verlet (leapfrog) integration
while ~stopFlag

    n=n+1;
    
    %form state vector
    X = reshape([x(:,n-1) v(:,n-1)]',[],1);

    %% Else perform constant timestep Verlet integration
    Xdot = feval(odeFunc,t(n-1)+dt,X);    
    A = Xdot(2:2:end);
    x(:,n) = 2*x(:,n-1) - x(:,n-2) + A*dt^2;
    v(:,n) = Xdot(1:2:end);
    t(n) = t(n-1)+dt;

    if max(t)>max(tspan)
        break
    end

end

%% Reshape u and v vectors into state vector X
X = zeros(n,length(X0));
X(:,1:2:end) = x(1:end,1:n)';
X(:,2:2:end) = v(1:end,1:n)';
t=t(1:n);
end
