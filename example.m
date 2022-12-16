%example Running script.
%
%This script computes and plots the controlled solution from a steady state to an other steady state in time T>0.
%Given $\theta^0,b^0$ and $\theta^1,b^1$, the problem is to find controls $u_s$ and $u_l$ such that the solution $\theta(t,x)$ and $b(t)$ of the Stefan problem with initial condition $\theta^0$ and $b^0$ and control $u_s$ and $u_l$ satisfies
%	$\theta(T,x)=\theta^1(x)$ and $b(T)=b^1$.
%The solution is written in each phase as
%	$\theta(t,x)=\sum_{k=0}^\infty \alpha_k(t)\frac{(x-b(t))^{2k+1}}{(2k+1)!}+\sum_{k=1}^\infty \beta_k(t)\frac{(x-b(t))^{2k}}{(2k)!}$
%where $\alpha_k$ and $\beta_k$ for $k>0$ are given by a recurrence formula (in particular they are obtained in term of $\alpha_0$).
%The problem is the to find $\alpha_0$.
%Once this is done, the Dirichlet boundary controls are obtained as:
%	$u_s(t)=\theta(t,0)$ and $u_l(t)=\theta(t,1)$.
%
%See the paper <a href="https://hal.archives-ouvertes.fr/hal-03721544">Controllability of the Stefan problem by the flatness approach</a> for more details.
%
%
%It is organised as follows :
%
% 1-Definition of plotting method (matlab and/or gnuplot):
%  * draw_matlab  - plots using matlab
%  * draw_video   - plotting video
%
% 2-Definition of parameters:
%  * cs, cl - diffusivity coefficients in the solid and liquid phases
%  * ntrap  - number of points for some integral evaluations (using the trapeze method)
%  * N      - truncation order in the sums
%  * nt     - number of points for the time evaluation of the solution
%  * nx     - number of points for the spacial evaluation of the temperature
%  * T      - controllability time (T>0)
%  * sigma  - Gevrey order (1<sigma<2)
%  * b0, v0 - define the initial state (theta^0(x)=v0*(x-b0))
%  * b1, v1 - define the final (target) state (theta^1(x)=v1*(x-b1))
%
% 3-Compute $\alpha_0$ such that the control problem is realised
%  * [alpha0s,alpha0l,db]=SOLVE_FLAT(v0,b0,v1,b1,T,t,sigma);
%
% 4-Compute the $\alpha_k$ and $\beta_k$ for k=1..N
%  * [Ys,Yl]=FLAT_ITER(alpha0s,alpha0l,db,t,N,cs,cl);
%
% 5-Compute times of interest (for the time evaluations of the solution)
%  * t_=COMPUTE_TIMES(alpha0s,alpha0l,t,T,ntrap,nt);
%
% 6-Evaluate the solution and the controls
%  * [bs,Ys_,Yl_]=COMPUTE_Y(t,t_,db,Ys,Yl,b0,b1);
%  * [xs,xl,thetass,thetals]=COMPUTE_THETA(nx,Ys_,Yl_,bs,v0,v1);
%  * [usN,ulN]=COMPUTE_UN(Ys_,Yl_,bs,v0,v1);
%
% 7-Plot the solution
%
% Once the simulation has ended, the results are stored in:
%  * t_               - time evaluations
%  * us, ul           - controls at times t_
%  * bs               - position at of the solid/liquid interface at times t_
%  * xs, xl           - abscissa for the spatial evaluation at times t_, xs(:,i)=linspace(0,bs(i),nx)' and xl(:,i)=linspace(bs(i),1,nx)'
%  * thetass, thetals - temperatures in the two phases, for instance thetass(i,j)=\theta(x(i,j),t_(j))
%
%
%Authors: B. Colle, J. Loheac and T. Takahashi.

p=gcp('nocreate');
if isempty(p), gcp(); end

% plotting parameters
draw_matlab=false;
draw_video=false;
draw_matlab=true;
draw_video=true;
if draw_video, t_video=10; end

%%%%%%%%%%%%%%%%
%% Parameters %%
%%%%%%%%%%%%%%%%
cl=1;
cs=1;
ntrap=1001;
N=2;
nt=101;
T=1;
nx=101;
syms t 'real'; assume(t>=0);
sigma=3/2;
b0=1/2;
v0=1;
v1=4;
b1=3/4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute controlled solution %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tstart0=tic();

fprintf('Computing the alpha_0''s... '); tstart=tic();
[alpha0s,alpha0l,db]=solve_flat(v0,b0,v1,b1,T,t,sigma);
fprintf('Done (%fs)\n',toc(tstart));
fprintf('Computing alpha_k and beta_k (k=1..%d)... \n',N); tstart=tic();
[Ys,Yl]=flat_iter(alpha0s,alpha0l,db,t,N,cs,cl);
fprintf('Done (%fs)\n',toc(tstart));
fprintf('Computing the time discretization (%d points)... ',nt); tstart=tic();
t_=compute_times(alpha0s,alpha0l,t,T,ntrap,nt);
fprintf('Done (%fs)\n',toc(tstart));
fprintf('Ccomputation of theta_s(t,x) and theta_l(t,x)... \n'); tstart=tic();
[bs,Ys_,Yl_]=compute_y(t,t_,db,Ys,Yl,b0,b1);
[xs,xl,thetass,thetals]=compute_theta(nx,Ys_,Yl_,bs,v0,v1);
fprintf('Done (%fs)\n',toc(tstart));
fprintf('Evaluation of u_s and u_l with respect to N... '); tstart=tic();
[usN,ulN]=compute_un(Ys_,Yl_,bs,v0,v1);
fprintf('Done (%fs).\n',toc(tstart));
us=thetass(1,:);
ul=thetals(end,:);

fprintf('Total enlapsed time:\t%fs\n',toc(tstart0));

%%%%%%%%%%%%%%
%% Plotting %%
%%%%%%%%%%%%%%
if draw_video || draw_matlab
	t_i=linspace(0,T,t_video*24);
	b_i=interp1(t_,bs,t_i);
	us_i=interp1(t_,us,t_i);
	ul_i=interp1(t_,ul,t_i);
	xs_i=linspace(0,1,nx)'*b_i;
	xl_i=linspace(0,1,nx)'*(1-b_i)+ones(nx,1)*b_i;
	thetass_i=interp1(t_,thetass',t_i)';
	thetals_i=interp1(t_,thetals',t_i)';

	fprintf('Plotting the solution.\n');
	if draw_video, subplot(2,2,1);
	else, subplot(2,1,1);
	end
	plot(t_i,us_i,t_i,ul_i);
	legend('us','ul')
	xlabel('t');
	if draw_video, subplot(2,2,2);
	else, subplot(2,1,2);
	end
	plot(t_i,b_i);
	legend('b');
	xlabel('t');
	if draw_video
		subplot(2,2,3:4);
		theta_min_max=[min(min(thetass)) max(max(thetals))];
		for i=1:length(t_i)
			plot(xs_i(:,i),thetass_i(:,i),xl_i(:,i),thetals_i(:,i))
			xlabel('x');
			title(strcat('t=',string(t_i(i))));
			ylim(theta_min_max);
			drawnow();
			pause(0.1);
		end
	end
end


%delete(gcp('nocreate'));
