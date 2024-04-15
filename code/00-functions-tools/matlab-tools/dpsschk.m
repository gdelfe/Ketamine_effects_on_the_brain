function [e,v]=dpsschk(tapers)
%DPSSCHK Helper function for multitaper routines.
%   E = DPSSCHK(TAPERS) returns the DPSS array based on the input
%     TAPERS.  If TAPERS is a cell containing the DPSS array then E 
%     just returns the array.
%     If TAPERS is in the form [N, NW, K] then E contains the corresponding
%     sequences.
%
%   [E, V] = DPSSCHK(TAPERS) returns the eigenvalues as well as the 
%      sequences.

%       Author: B. Pesaran, 03-12-98
%

nout=max(nargout,1);

flag=0;

if ~iscell(tapers)
	sz=size(tapers);
	if max(sz == 3)	flag=1; end

	if flag 
%   	disp('Calculating tapers ... ');
		N=tapers(1);
		NW=tapers(2);
		K=tapers(3);
        if rem(N,1)
            N = round(N); %deal with rounding errors
        end
		if K < 1 
			error('Error:  K must be greater than or equal to 1');
		end
		if K < 3
%			disp('Warning:  K is less than 3'); 
		end	
		if K > 2*NW-1 
			error('Error:  K must be less than 2*P-1');
		end 
		[e,v]=dpss(N,NW,'calc');
		e=e(:,1:K);
      v=v(1:K);
%       disp('	... Done');
	end
	if ~flag
%   	disp('Tapers already calculated');
		e=tapers;
		v=0;
	end
end

if iscell(tapers)
  e=tapers{1,1};
  v=0;
  if length(tapers) == 2
    v=tapers{1,2};
  end
end

