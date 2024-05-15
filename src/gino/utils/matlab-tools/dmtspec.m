
function [spec, f, err] = dmtspec(X, tapers, sampling, fk, pad, pval, flag, Errorbar)
% DMTSPEC calculates the direct multitaper spectral estimate for time series.
%
% [SPEC, F, ERR] = DMTSPEC(X, TAPERS, SAMPLING, FK, PAD, PVAL, FLAG, ERRORBAR)
%
%  Inputs:  X		=  Time series array in [Space/Trials,Time] form.
%	    TAPERS 	=  Data tapers in [K,TIME], [N,P,K] or [N,W] form.
%			   	Defaults to [N,3,5] where N is duration of X.
%	    SAMPLING 	=  Sampling rate of time series X in Hz.
%				Defaults to 1.
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
%
%       ERRORBAR = Structure. ERRRORBAR.Type - 'Jacknife' or 'Chi-squared'
%
%
%  Outputs: SPEC	=  Spectrum of X in [Space/Trials, Freq] form.
%	    F		=  Units of Frequency axis for SPEC.
%	    ERR 	=  Error bars for SPEC in [Hi/Lo, Space/Trials, Freq]
%			   form given by a Jacknife-t interval for PVAL.
%

% Modification History:
%           Written by:  Bijan Pesaran, 08/97
%           Modified:    Added error bars BP 08/27/98

sX = size(X);
nt = sX(2);
nch = sX(1);

%  Set the defaults

if nargin < 3 sampling = 1; end
nt = nt./sampling;
if nargin < 2 tapers = [nt,3,5]; end
if length(tapers) == 2
    n = tapers(1);
    w = tapers(2);
    p = n*w;
    k = floor(2*p-1);
    tapers = [n,p,k];
    %   disp(['Using ' num2str(k) ' tapers.']);
end
if length(tapers) == 3
    tapers(1) = round(tapers(1).*sampling);
    tapers = dpsschk(tapers);
end
if nargin < 4 fk = [0,sampling./2]; end
if length(fk) == 1
    fk = [0,fk];
end
if nargin < 5 pad = 2; end
if nargin < 6 pval = 0.05;  end
if nargin < 7 flag = 0; end

N = length(tapers(:,1));
nt = round(nt.*sampling);
if N ~= nt error('Error:  Length of time series and tapers must be equal'); end

K = length(tapers(1,:));
nf = max(256,pad*2.^(nextpow2(N+1)));
nfk = floor(fk./sampling.*nf);
dof = 2.*nch.*K;

% Determine outputs
f = linspace(fk(1),fk(2),diff(nfk));
errorchk = 0;
if nargout > 2 errorchk = 1; end
if nargin < 8 && errorchk
    Errorbar.Type = 'Chi-squared';
end

if ~flag		% No pooling across trials
    %disp('No pooling across trials');
    spec = zeros(nch, diff(nfk));
    err = zeros(2, nch, diff(nfk));
    if nch == 1 mX = sum(X)./nt; else mX = sum(X,1)./nch; end
    for ch=1:nch
        tmp = (X(ch,:) - mX)';
        xk = fft(tapers(:,1:K).*tmp(:,ones(1,K)),nf)';
        xk = xk(:,nfk(1)+1:nfk(2));
        Sk = xk.*conj(xk);
        spec(ch,:) = sum(Sk,1)./K;
        
        if errorchk	%  Estimate error bars using Jacknife
            switch Errorbar.Type
                case 'Jacknife'
                    for ik = 1:K
                        indices = setdiff([1:K],[ik]);
                        xj = xk(indices,:);
                        jlsp(ik,:) = log(sum(xj.*conj(xj))./(K-1));
                    end
                    lsig = sqrt(K-1).*std(jlsp,1);
                    crit = tinv(1-pval./2,dof-1);		%   Determine the scaling factor
                    err(1,ch,:) = exp(log(spec(ch,:))+crit.*lsig);
                    err(2,ch,:) = exp(log(spec(ch,:))-crit.*lsig);
                case 'Chi-squared'
                    a = chi2inv(1-pval./2,dof);
                    b = chi2inv(pval./2,dof);
                    err(1,ch,:) = spec(ch,:).*dof./b;
                    err(2,ch,:) = spec(ch,:).*dof./a;
            end
        end
    end
end


if flag			% Pooling across trials
    spec = zeros(1, diff(nfk),'single');
    err = zeros(2, diff(nfk),'single');
    
    Xk = zeros(nch*K, diff(nfk),'single');
    mX = sum(X,1)./nch;
    for ch=1:nch
        tmp = (X(ch,:) - mX)';
        xk = fft(tapers(:,1:K).*tmp(:,ones(1,K)),nf)';
        Xk((ch-1)*K+1:ch*K,:) = xk(:,nfk(1)+1:nfk(2));
    end
    spec = sum(Xk.*conj(Xk),1)./(K.*nch);
    
    if errorchk		%  Estimate error bars
        switch Errorbar.Type
            case 'Jacknife'
                for ik = 1:nch*K
                    if mod(ik,1000)==0
                        sprintf('Jacknife iter: %d of %d', [ik, nch*K])
                    end
                    indices = setdiff([1:K*nch],[ik]);
                    xj = Xk(indices,:);
                    jlsp(ik,:) = log(sum(xj.*conj(xj),1)./(K*nch-1));
                end
                lsig = sqrt(nch*K-1).*std(jlsp,1);
                crit = tinv(1-pval./2,dof-1);		%   Determine the scaling factor
                err(1,:) = exp(log(spec)+crit.*lsig);
                err(2,:) = exp(log(spec)-crit.*lsig);
            case 'Chi-squared'
                a = chi2inv(1-pval./2,dof);
                b = chi2inv(pval./2,dof);
                err(1,:) = spec.*dof./b;
                err(2,:) = spec.*dof./a;
        end
    end
end
