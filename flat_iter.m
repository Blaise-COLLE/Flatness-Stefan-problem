function [Ys,Yl]=flat_iter(alpha0s,alpha0l,db,t,N,cs,cl)
%FLAT_ITER Compute the recurrence formula.
%Solution of the Stefan problem is written:
%	$\theta(t,x)=\sum{k=1}^\infty y_k(t)\frac{x^k}{k!}$
%with $y_{2k+1}=\alpha_k$ and $y_{2k}=\beta_k$.
%$\alpha$ and $\beta$ satisfy the recurrences
%	$c\alpha_{k+1}=\dot{\alpha_k}-\dot{b}\beta_{k+1}$ and
%	$c\beta_{k+1}=\dot{\beta_k}-\dot{b}\alpha_k$,
%with $\alpha_0$ given and $\beta_0=0$.
%
% [Ys,Yl]=FLAT_ITER(alpha0s,alpha0l,db,t,N,cs,cl)
%inputs:
% * alpha0s - $\alpha_0$ in the solid phase (symbolic)
% * alpha0l - $\alpha_0$ in the liquid phase (symbolic)
% * db      - $\dot{b}$ (symbolic)
% * t       - symbolic time variable
% * N       - truncation order in the sum (N=size(Ys,1)=size(Yl,1))
% * cs      - diffusivity coefficient in the solid phase
% * cl      - diffusivity coefficient in the liquid phase
%outputs:
% * Ys - $y$ for the solid phase (symbolic array of length N)
% * Yl - $y$ for the liquid phase (symbolic array of length N)
%
%Note: alpha0s, alpha0l and db are computed from SOLVE_FLAT.
%
%Used in EXAMPLE.
%
%Authors: B. Colle, J. Loheac and T. Takahashi.

	p=gcp('nocreate');
	if isempty(p)
		Ys=flat_iter1phase(alpha0s,db,t,N,cs);
		Yl=flat_iter1phase(alpha0l,db,t,N,cl);
		Ys=vpa(Ys);
		Yl=vpa(Yl);
	else
		Asl=sym('Asl',[1 2]);
		Asl(1)=alpha0s; Asl(2)=alpha0l;
		Ysl=sym('Ysl',[2*N+1 2]);
		c=[cs cl];
		parfor i=1:2
			Ysl(:,i)=flat_iter1phase(Asl(i),db,t,N,c(i));
		end
		fprintf('\tConversion vpa... '); tstart=tic();
		numWorkers=p.NumWorkers;
		if numel(Ysl)<=numWorkers
			Ysl=vpa(Ysl);
		else
			done=zeros(1,numel(Ysl));
			todo=numel(Ysl);
			runs=zeros(1,numWorkers);
			for pd=1:numWorkers
				F(pd)=parfeval(@vpa,1,Ysl(todo));
				runs(pd)=todo;
				todo=todo-1;
			end
			while ~all(done)
				pd=fetchNext(F);
				Ysl(runs(pd))=fetchOutputs(F(pd));
				done(runs(pd))=true;
				if todo>0
					F(pd)=parfeval(@vpa,1,Ysl(todo));
					runs(pd)=todo;
					todo=todo-1;
				end
			end
		end
		Ys=Ysl(:,1);
		Yl=Ysl(:,2);
		fprintf('Done (%fs)\n',toc(tstart));
	end
end

function Y=flat_iter1phase(alpha0,db,t,N,c)
	Y=sym('Y',[2*N+1 1]);
	Y(1)=alpha0;
	fprintf('\tk=1... '); tstart=tic();
	Y(2)=simplify(-db*alpha0)/c;
	Y(3)=simplify(diff(Y(1),t)-db*Y(2))/c;
	fprintf('Done (%fs)\n',toc(tstart));
	for k=2:N
		fprintf('\tk=%d... ',k); tstart=tic();
		Y(2*k)=(diff(Y(2*(k-1)),t)-db*Y(2*k-1))/c;
		Y(2*k+1)=(diff(Y(2*k-1),t)-db*Y(2*k))/c;
		fprintf('Done (%fs)\n',toc(tstart));
	end
end
