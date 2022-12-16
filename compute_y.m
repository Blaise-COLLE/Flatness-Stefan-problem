function [bs,Ys_,Yl_]=compute_y(t,t_,db,Ys,Yl,b0,b1)
%COMPUTE_Y Evaluate Y and b at given time instants.
% [bs,Ys_,Yl_]=COMPUTE_Y(t,t_,db,Ys,Yl,b0,b1)
%inputs:
% * t  - symbolic time variable
% * t_ - evaluation times
% * db - time derivative of b (symbolic)
% * Ys - Coefficients in the solid phase (symbolic), $\theta_s(.,x)=\sum_{k=1}^N Ys(k)\frac{(x-b)^k}{k!}$.
% * Yl - Same as Ys for the liquid phase
% * b0 - initial position of the solid/liquid interface
% * b1 - initial position of the solid/liquid interface
%outputs:
% * bs  - evaluation of $b=b0+\int_0^t db$ at times given in t_
% * Ys_ - evaluation of Ys at times given in t_
% * Yl_ - evaluation of Yl at times given in t_
%
%Used in EXAMPLE.
%
%Authors: B. Colle, J. Loheac and T. Takahashi.

	p=gcp('nocreate');
	nt=length(t_);
	t__=vpa(t_);
	N=length(Ys);
	bs=double(subs(vpa(db),t,t__));
 	bs=cumtrapz(t_,bs)+b0; bs(end)=b1;
	if isempty(p)
		fprintf('\tEvaluation of y_k (k=1..%d), y_{2k+1}=alpha_k and y_{2k}=beta_k... ',N); tstart=tic();
		Ys_=double(subs(Ys,t,t__));
		Yl_=double(subs(Yl,t,t__));
		fprintf('Done (%fs).\n',toc(tstart));
	else
		fprintf('\tEvaluation of y_k (k=1..%d), y_{2k+1}=alpha_k and y_{2k}=beta_k... \n',N); tstart=tic();
		Y_=zeros(2*N,nt);
		Y(2:2:2*N)=Ys;
		Y(1:2:2*N)=Yl;
		numWorkers=p.NumWorkers;
		if 2*N<=numWorkers
			parfor j=1:2*N
				tst=tic();
				Y_(j,:)=double(subs(Y(j),t,t__));
				if mod(j,2)==0, fprintf('\t\tk=%d (solid phase)...  Done (%fs)\n',j/2,toc(tst));
				else, fprintf('\t\tk=%d (liquid phase)... Done (%fs)\n',(j+1)/2,toc(tst));
				end
			end
		else
			done=zeros(1,2*N);
			todo=2*N;
			runs=zeros(1,numWorkers);
			for pd=1:numWorkers
				F(pd)=parfeval(@mysubs,2,Y(todo),t,t__);
				runs(pd)=todo;
				todo=todo-1;
			end
			while ~all(done)
				pd=fetchNext(F);
				j=runs(pd);
				[etime,Y_(j,:)]=fetchOutputs(F(pd));
				if mod(j,2)==0, fprintf('\t\tk=%d (solid phase)...  Done (%fs)\n',j/2,etime);
				else, fprintf('\t\tk=%d (liquid phase)... Done (%fs)\n',(j+1)/2,etime);
				end
				done(j)=true;
				if todo>0
					F(pd)=parfeval(@mysubs,2,Y(todo),t,t__);
					runs(pd)=todo;
					todo=todo-1;
				end
			end
		end
		Ys_=Y_(2:2:end,:);
		Yl_=Y_(1:2:end,:);
		fprintf('\tDone (%fs).\n',toc(tstart));
	end
end

function [etime,Y_]=mysubs(Y,t,t__)
	etime=tic();
	Y_=double(subs(Y,t,t__));
	etime=toc(etime);
end
