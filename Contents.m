% Control of the two pahses Stefan problem using flatness.
%
% Files by categories
%  
%Gevrey Functions
%   BUMP_FUNCTION  - Gevrey function of compact support in [0,1].
%   STEP_FUNCTION  - Gevrey function equal to 1 for t<0 and to 0 for t>1.
%
%Solving the control problem using flatness
%   COMPUTE_THETA  - Compute the value of the temperatures in each phase.
%   COMPUTE_TIMES  - Conpute the time instants for the evaluation, based on the variation of alpha0s and alpha0l.
%   COMPUTE_UN     - Evaluate the controls at given times.
%   COMPUTE_Y      - Evaluate Y and b at given time instants.
%   FLAT_ITER      - Compute the recurrence formula.
%   SOLVE_FLAT     - Find the expression of the flat outputs.
%
%example:
%   EXAMPLE        - Running script.
%
%Authors: Blaise Colle, Jerome Loheac and Takeo Takahashi.

%see http://www.fast.u-psud.fr/~moisy/ml/ for the generation of the html documentary files.
makehtmldoc('*.m','upper',... %, 'code'...
	'color', '#ffff00',...
	'title', 'Help for \f',...
	'firstline', '<a href="Contents.html">Back</a>', ...
	'lastline', 'Codes related to the paper: <a href="https://hal.archives-ouvertes.fr/hal-03721544">Controllability of the Stefan problem by the flatness approach</a>.'...
);
if ~exist('Documentation','dir'), mkdir('Documentation'); end
movefile('*.html','./Documentation/');
