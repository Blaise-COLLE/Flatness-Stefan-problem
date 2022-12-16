function t_=compute_times(alpha0s,alpha0l,t,T,ntrap,nt)
%COMPUTE_TIMES Compute the time instants for the evaluation, based on the variation of alpha0s and alpha0l.
%Evaluation is taking lots of time. Hence, we only evaluate the solution at some "well-chosen" time instants.
% [t_]=COMPUTE_TIMES(alpha0s,alpha0l,t,T,ntrap,nt)
%inputs:
% * alpha0s - flat output in the solid phase (symbolic)
% * alpha0l - flat output in the liquid phase (symbolic)
% * t       - symbolic time variable
% * T       - final time.
% * ntrap   - number of points for computing time integrals
% * nt      - number of points for evaluation (length of t_)
%output:
% * t_ - times which will be used for the evaluation of the solution.
%
%Used in EXAMPLE.
%
%Authors: B. Colle, J. Loheac and T. Takahashi.

	t_=linspace(0,T,ntrap);
	dalpha0s=diff(alpha0s,t);
	dalpha0l=diff(alpha0l,t);
	da_=abs(double(subs(dalpha0s,t,vpa(t_))))+abs(double(subs(dalpha0l,t,vpa(t_))));
	da_(1)=0; da_(end)=0;
	da_=da_+T/(nt-1);
	a_=cumtrapz(t_,da_);
	aq_=linspace(0,a_(end),nt); aq_=aq_(2:end-1);
	t_=interp1(a_,t_,aq_);
	t_=[0 t_ T];
	dtmin=5*T/(nt-1);
	ind=find(diff(t_)>dtmin);
	tad=[];
	for i=ind
		dti=t_(i+1)-t_(i);
		nti=ceil(dti/dtmin)+1;
		tadi=linspace(t_(i),t_(i+1),nti);
		tad=[tad tadi(2:end-1)];
	end
	t_=sort([t_ tad]);
end
