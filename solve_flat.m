function [alpha0s,alpha0l,db]=solve_flat(v0,b0,v1,b1,T,t,sigma,ntrap)
%SOLVE_FLAT Find the expression of the flat outputs.
%See FLAT_ITER for the description of the solution of the Stefan problem.
%This function gives the expression fo alpha_0 in both phases.
% [alpha0s,alpha0l,db]=SOLVE_FLAT(v0,b0,v1,b1,T,t,sigma,ntrap)
%inputs:
% * v0    - initial slope of the temperature
% * b0    - initial position of the solid/liquid interface
% * v1    - final slope of the temperature
% * b1    - final position of the solid/liquid interface
% * T     - control time
% * t     - symbolic time variable
% * sigma - Gevrey order
% * ntrap - number of points for the evaluation of time integrals.
%outputs:
% * alpha0s - function \alpha_0 in the solid phase (symbolic)
% * alpha0l - function \alpha_0 in the liquid phase (symbolic)
% * db      - function \dot{b} (symbolic)
%
%See also STEP_FUNCTION and BUMP_FUNCTION.
%
%Used in EXAMPLE.
%
%Authors: B. Colle, J. Loheac and T. Takahashi.

	alpha0s=v1+(v0-v1)*step_function(t/T,sigma);
	alpha0l=alpha0s;
	if nargin<8, ntrap=1001; end
	t_=linspace(0,T,ntrap);
	Y=bump_function(t_/T,sigma);
	I=trapz(t_,Y);
	db=((b1-b0)/I)*bump_function(t/T,sigma);
	if b1-b0<0,			alpha0l=alpha0l+db;
	elseif b1-b0>0,	alpha0s=alpha0s+db;
	end
end
