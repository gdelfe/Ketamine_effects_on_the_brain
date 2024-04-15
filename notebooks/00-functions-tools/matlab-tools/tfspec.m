function [spec, f, ti, err] = tfspec(X, tapers, sampling, dn, fk, pad, pval, flag, contflag, Errorbar)

%TFSPEC  Moving window time-frequency spectrum using multitaper techniques.
%
% [SPEC, F, TI, ERR] = TFSPEC(X, TAPERS, SAMPLING, DN, FK, PAD, PVAL, FLAG, CONTFLAG, ERRORBAR) 
%
%  Inputs:  X		=  Time series array in [Space/Trials,Time] form.
%	    TAPERS 	=  Data tapers in [K,TIME], [N,P,K] or [N,W] form.
%			   	    [N,W] Form:  N = duration of analysis window in s.
%                                W = bandwidth of frequency smoothing in Hz.
%               Defaults to [N,3,5]  where N is NT/10
%				and NT is duration of X. 
%               
%	    SAMPLING 	=  Sampling rate of time series X in Hz. 
%				Defaults to 1.
%	    DN		=  Window step.
%			       	Defaults to N./10;
%	    FK 	 	=  Frequency range to return in Hz in
%                               either [F1,F2] or [F2] form.  
%                               In [F2] form, F1 is set to 0.
%			   	Defaults to [0,SAMPLING/2]
%	    PAD		=  Padding factor for the FFT.  
%			      	i.e. For N = 500, if PAD = 2, we pad the FFT 
%			      	to 1024 points; if PAD = 4, we pad the FFT
%			      	to 2048 points.
%				Defaults to 2.
%	   PVAL		=  P-value to calculate error bars for.
%				Defaults to 0.05 i.e. 95% confidence.
%
%	   FLAG = 0:	calculate SPEC seperately for each channel/trial.
%	   FLAG = 1:	calculate SPEC by pooling across channels/trials. 
%      CONTFLAG = 1; There is only a single continuous signal coming in.
%
%  Outputs: SPEC	=  Spectrum of X in [Space/Trials, Time, Freq] form. 
%	    F		=  Units of Frequency axis for SPEC.
%	    ERR 	=  Error bars in[Hi/Lo, Space, Time, Freq]  
%			   form given by the Jacknife-t or Chi2 interval for PVAL.
%       TI      =
% 
%   See also DPSS, PSD, SPECGRAM.

%   Author: Bijan Pesaran, version date 15/10/98.
%               Optimized when not computing error bars.

sX = size(X);
nt = sX(2);              % calculate the number of points
nch = sX(1);               % calculate the number of channels

if nargin < 3 sampling = 1.; end
n = floor(nt./10)./sampling;
if nargin < 2 tapers = [n,3,5]; end
if length(tapers) == 2
   n = tapers(1);
   w = tapers(2);
   p = n*w;
   k = floor(2*p-1);
   tapers = [n,p,k];
%   disp(['Using ' num2str(k) ' tapers.']);
end
if length(tapers) == 3
   tapers(1) = floor(tapers(1).*sampling);  
   tapers = single(dpsschk(tapers));
end
if nargin < 4 || isempty(dn)
    dn = n./10; 
end
if nargin < 5 || isempty(fk)
    fk = [0,sampling./2.]; 
end
if length(fk) == 1
    fk = [0,fk];
end
if nargin < 6 || isempty(pad)
    pad = 2; 
end
if nargin < 7 || isempty(pval)
    pval = 0.05; 
end
if nargin < 8 || isempty(flag)
    flag = 0; 
end
if nargin < 9 || isempty(contflag) 
    contflag = 0; 
end
if nargin < 10 || isempty(Errorbar)
    Errorbar.Type = 'Chi-squared';
end

nch = single(nch);
K = single(length(tapers(1,:))); 
N = length(tapers(:,1));

if N > nt error('Error: Tapers are longer than time series'); end

% Determine outputs
errorchk = 0;
if nargout > 3 errorchk = 1; end

dn = floor(dn.*sampling);
nf = max(256, pad*2^nextpow2(N+1)); 
nfk = floor(fk./sampling.*nf);
nwin = single(floor((nt-N)./dn));          % calculate the number of windows
f = linspace(fk(1),fk(2),diff(nfk));

if ~flag				% No pooling across trials
    spec = zeros(nch,nwin,diff(nfk),'single');  %Note that the array initialization is single precision
    if ~errorchk				%  Don't estimate error bars
       
        for win = 1:nwin
            % Here the optimized spectral loop starts.
            if contflag
                tmp = detrend(X(:,dn*win:dn*win+N-1))';
                if size(tmp,1)>N   %machine precision work-around? added by alo for weird behavior 181000020
                    tmp = tmp(1:N,:);
                end
            else
                mX = sum(X(:,dn*(win-1)+1:dn*(win-1)+N),1)./nch;
                tmp = (X(:,dn*(win-1)+1:dn*(win-1)+N) - mX(ones(1,nch),:)).';
            end
            for ch = 1:nch
                Xk = fft(tapers.*tmp(:,ch*ones(1,K)),nf);
                Xk = Xk(nfk(1)+1:nfk(2),:);
                spec(ch,win,:) = (sum(Xk.*conj(Xk),2)./(K)).';
            end
            %  The optimized loop ends here
        end
    else 					%  Estimate error bars - this is not optimized
	    % Broken
        [ftmp, dum, err_tmp] = dmtspec(tmp, tapers, sampling, fk, pad, pval);
        spec(ch,win,:) = ftmp;
        err(1,ch,win,:) = err_tmp(1,:);
        err(2,ch,win,:) = err_tmp(2,:);
    end
end


if flag					% Pooling across trials
    spec = zeros(nwin,diff(nfk),'single');

    %disp('Flag = 11');
    %     ind = repmat([1:nch],K,1); ind = ind(:);
    for win = 1:nwin
        %  The optimized loop starts here
        if contflag
            tmp = X(:,dn*(win-1)+1:dn*(win-1)+N);
        else
            mX = sum(X(:,dn*(win-1)+1:dn*(win-1)+N),1)./nch; %finds sum for that window
            tmp = (X(:,dn*(win-1)+1:dn*(win-1)+N) - mX(ones(1,nch),:)).'; %subtract off sum?
        end
        if ~errorchk				%  Don't estimate error bars
            SX = zeros(diff(nfk),1,'single');
            for ch = 1:nch
                Xk = fft(tapers.*tmp(:,ch*ones(1,K)),nf);
                Xk = Xk(nfk(1)+1:nfk(2),:);
                SX = SX + sum(Xk.*conj(Xk),2);
            end
            spec(win,:) = SX'./(K.*nch);
            %  The optimized loop ends here
        else					%  Estimate error bars - This is not optimized
            [ftmp, dum, err_tmp] = dmtspec(tmp', tapers, sampling, fk, pad, pval, flag, Errorbar);
            spec(win,:) = ftmp;
            err(1,win,:) = err_tmp(1,:);
            err(2,win,:) = err_tmp(2,:);
        end
    end
end

ti = linspace(N/2,nt-N/2,nwin);

if size(spec,1) == 1 && length(size(spec)) > 2
    spec = sq(spec);
end
