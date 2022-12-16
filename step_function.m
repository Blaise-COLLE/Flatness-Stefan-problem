function [y]=step_function(t,sigma)
%STEP_FUNCTION Gevrey function equal to 1 for t<0 and to 0 for t>1.
%	[y]=STEP_FUNCTION(t,sigma)
%inputs:
% * sigma is the gevrey order
% * t is real (symbolic or double array)
%output:
% * y=0, for t>1,
%     1, for t<0,
%     exp(-(1-t)^(-k))/(exp(-(1-t)^(-k))+exp(-t^(-k))), for 0<t<1, with k=1/(sigma-1).
%
%Used in BUMP_FUNCTION and SOLVE_FLAT.
%
%Authors: B. Colle, J. Loheac and T. Takahashi.

	k=1/(sigma-1);
	switch class(t)
	case 'sym'
		y=piecewise(t<=0,1,t>=1,0,0<t<1,exp(-(1-t)^(-k))/(exp(-(1-t)^(-k))+exp(-t^(-k))));
	case 'double'
		tt=t(t>0 & t<1);
		y=zeros(size(t));
		y(t>0 & t<1)=exp(-(1-tt).^(-k))./(exp(-(1-tt).^(-k))+exp(-tt.^(-k)));
		y(t<=0)=1;
	otherwise
		error('unexpected class: ''%s'' for variable t',class(t));
	end
end

%On calcule en t dans (0,1) la step function.
