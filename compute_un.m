function [us,ul]=compute_un(Ys_,Yl_,bs,v0,v1)
%COMPUTE_UN Evaluate the controls at given times.
% [us,ul]=COMPUTE_UN(Ys_,Yl_,bs,v0,v1)
%inputs:
% * Ys_, Yl_, bs - evaluations of Ys, Yl and b made in COMPUTE_Y
% * v0, v1       - initial and final temperature slope
%outputs:
% * us - evaluation fo the control in the solid phase
% * ul - evaluation fo the control in the liquid phase
%
%See also COMPUTE_Y.
%
%Used in EXAMPLE.
%
%Authors: B. Colle, J. Loheac and T. Takahashi.

	[K,nt]=size(Ys_);
	xbs=-bs(2:end-1);
	xbl=(1-bs(2:end-1));
	Ys__=Ys_(:,2:end-1);
	Yl__=Yl_(:,2:end-1);
	us=zeros(K,nt-2); ul=zeros(K,nt-2);
	us(1,:)=xbs.*Ys__(1,:);
	ul(1,:)=xbl.*Yl__(1,:);
	for k=2:K
		xbs=-bs(2:end-1).*xbs/k;
		xbl=(1-bs(2:end-1)).*xbl/k;
		us(k,:)=xbs.*Ys__(k,:)+us(k-1,:);
		ul(k,:)=xbl.*Yl__(k,:)+ul(k-1,:);
	end
%	us=us(2:2:end,:);
%	ul=ul(2:2:end,:);
%	N=size(us,1);
%	us=[-v0*bs(1)*ones(N,1) us -v1*bs(end)*ones(N,1)];
%	ul=[v0*(1-bs(1))*ones(N,1) ul v1*(1-bs(end))*ones(N,1)];
	us=[-v0*bs(1)*ones(K,1) us -v1*bs(end)*ones(K,1)];
	ul=[v0*(1-bs(1))*ones(K,1) ul v1*(1-bs(end))*ones(K,1)];
end
