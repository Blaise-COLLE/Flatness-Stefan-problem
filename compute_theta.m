function [xs,xl,thetass,thetals]=compute_theta(nx,Ys_,Yl_,bs,v0,v1)
%COMPUTE_THETA Compute the value of the temperatures in each phase.
% [xs,xl,thetass,thetals]=COMPUTE_THETA(nx,Ys_,Yl_,bs,v0,v1)
%inputs:
% * nx  - number of discretisation points in each phase
% * Ys_ - the solution is written as $\theta=\sum_{k=0}^N y_k(t) \frac{(x-b(t))^k}{k!}$, $Y_s(k,l)=y_k(t_l)$.
% * Yl_ - same as Ys_ for the liquid phase
% * bs  - position of the solid/liquid interface
% * v0  - initial slope of the temperature
% * v1  - final solpe of the temperature
% * xs  - abscissa in the solid phase
% * xl  - abscissa in the liquid phase
%outputs:
% * thetass - temperatures in the solid phase
% * thetals - temperatures in the liquid phase
%
%Used in EXAMPLE.
%
%Authors: B. Colle, J. Loheac and T. Takahashi.

	[N,nt]=size(Ys_);
	thetass=zeros(nx,nt); thetals=zeros(nx,nt);
	xs=linspace(0,1,nx)'*bs;
	xl=linspace(0,1,nx)'*(1-bs)+bs;
	thetass(:,1)=v0*(xs(:,1)-bs(1));
	thetals(:,1)=v0*(xl(:,1)-bs(1));
	fprintf('\tEvaluation of theta_s and theta_l... '); tstart=tic();
	if isempty(gcp('nocreate'))
		for i=2:nt-1
			xbs=xs(:,i)-bs(i); xbs_=xbs;
			xbl=xl(:,i)-bs(i); xbl_=xbl;
			thetass(:,i)=Ys_(1,i)*xbs;
			thetals(:,i)=Yl_(1,i)*xbl;
			for k=2:N
				xbs=xbs.*xbs_/k;
				xbl=xbl.*xbl_/k;
				thetass(:,i)=thetass(:,i)+Ys_(k,i)*xbs;
				thetals(:,i)=thetals(:,i)+Yl_(k,i)*xbl;
			end
		end
	else
		parfor i=2:nt-1
			xbs=xs(:,i)-bs(i); xbs_=xbs;
			xbl=xl(:,i)-bs(i); xbl_=xbl;
			Ysi_=Ys_(:,i);
			Yli_=Yl_(:,i);
			ttas=Ysi_(1)*xbs
			ttal=Yli_(1)*xbl;
			for k=2:N
				xbs=xbs.*xbs_/k;
				xbl=xbl.*xbl_/k;
				ttas=ttas+Ysi_(k)*xbs;
				ttal=ttal+Yli_(k)*xbl;
			end
			thetass(:,i)=ttas;
			thetals(:,i)=ttal;
		end
	end
	thetass(:,end)=v1*(xs(:,end)-bs(end));
	thetals(:,end)=v1*(xl(:,end)-bs(end));
	fprintf('Done (%fs).\n',toc(tstart));
end
