function [y]=bump_function(t,sigma)
%BUMP_FUNCTION Gevrey function of compact support in [0,1].
%	[y]=BUMP_FUNCTION(t,sigma)
%inputs:
% * sigma is the gevrey order
% * t is real (symbolic or double array)
%output:
% * y=STEP_FUNCTION(t,sigma)*STEP_FUNCTION(1-t,sigma)			
%
%Used in SOLVE_FLAT.
%
%See also: STEP_FUNCTION.
%
%Authors: B. Colle, J. Loheac and T. Takahashi.

	y=step_function(t,sigma).*step_function(1-t,sigma);
end
