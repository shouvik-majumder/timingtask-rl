function [lp, g, resp, q] = ttPolicy(tau, u, s, family)
%TTPOLICY  Log density of one lick time under a policy and its score with respect to each channel.
%
%   [lp, g, resp, q] = ttPolicy(tau, u, s, family)
%
%   family 'ln'      log-normal                                    u = [z]              s = [sigma]
%          'lnexp'   (1-pi) LN(z, sigma) + pi Exp(r)              u = [z; log r; ell]  s = [sigma]
%          'lnmix2'  (1-pi) LN(z2, s2) + pi LN(z1, s1)            u = [z2; z1; ell]   s = [s2; s1]
%   z is a log median, ell the logit of the early component's weight pi = 1 / (1 + exp(-ell)).
%   The first channel is always the timed component's log median.
%
%   lp    log p(tau)
%   g     d lp / d u, one entry per channel: for a mixture, the responsibility of the component
%         times that component's own score; for the logit, responsibility minus weight
%   resp  responsibility of the early component (0 for 'ln')
%   q     struct with the predictive median, the mean of log tau and a CDF handle
lt = log(tau);
switch family
    case 'ln'
        z = u(1); sg = s(1); d = (lt - z) / sg;
        lp = -lt - log(sg) - 0.9189385332046727 - 0.5 * d^2;
        g = d / sg; resp = 0;
        if nargout > 3, q = struct('median', exp(z), 'meanLog', z, 'cdf', @(x) normcdf((log(x) - z) / sg)); end
    case 'lnexp'
        z = u(1); r = exp(u(2)); ell = u(3); pi1 = 1 / (1 + exp(-ell)); sg = s(1); d = (lt - z) / sg;
        lT = -lt - log(sg) - 0.9189385332046727 - 0.5 * d^2; lE = log(r) - r * tau;
        mx = max(lT, lE); pT = (1 - pi1) * exp(lT - mx); pE = pi1 * exp(lE - mx);
        lp = mx + log(pT + pE); rE = pE / (pT + pE); rT = 1 - rE;
        g = [rT * d / sg; rE * (1 - r * tau); rE - pi1]; resp = rE;
        if nargout > 3
            cdf = @(x) (1 - pi1) * normcdf((log(x) - z) / sg) + pi1 * (1 - exp(-r * x));
            meanLog = (1 - pi1) * z + pi1 * (-0.5772156649015329 - log(r));      % E[log X] of an exponential = -gamma_E - log r
            q = struct('median', medianOf(cdf), 'meanLog', meanLog, 'cdf', cdf);
        end
    case 'lnmix2'
        z2 = u(1); z1 = u(2); ell = u(3); pi1 = 1 / (1 + exp(-ell)); s2 = s(1); s1 = s(2);
        d2 = (lt - z2) / s2; d1 = (lt - z1) / s1;
        l2 = -lt - log(s2) - 0.9189385332046727 - 0.5 * d2^2; l1 = -lt - log(s1) - 0.9189385332046727 - 0.5 * d1^2;
        mx = max(l1, l2); p2 = (1 - pi1) * exp(l2 - mx); p1 = pi1 * exp(l1 - mx);
        lp = mx + log(p1 + p2); r1 = p1 / (p1 + p2); r2 = 1 - r1;
        g = [r2 * d2 / s2; r1 * d1 / s1; r1 - pi1]; resp = r1;
        if nargout > 3
            cdf = @(x) (1 - pi1) * normcdf((log(x) - z2) / s2) + pi1 * normcdf((log(x) - z1) / s1);
            q = struct('median', medianOf(cdf), 'meanLog', (1 - pi1) * z2 + pi1 * z1, 'cdf', cdf);
        end
    otherwise
        error('ttPolicy:family', 'unknown family %s', family);
end
end
function m = medianOf(cdf)
lo = 1e-3; hi = 60;
for i = 1:40, mid = sqrt(lo * hi); if cdf(mid) < 0.5, lo = mid; else, hi = mid; end, end
m = sqrt(lo * hi);
end
