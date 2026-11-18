function S = bbSample(n, P, rs)
%BBSAMPLE  Draw n complete latent units from the population defined by P.
%
%   S = BBSAMPLE(N,P) returns a struct with fields
%       S.Z    N-by-d   latent coordinates
%       S.XB   N-by-b   backbone attributes
%       S.XQ   1-by-M   cell, S.XQ{m} is N-by-pm  (exclusive attributes)
%       S.r    1-by-M   cell, S.r{m}  is N-by-pm  (oracle conditional means)
%       S.sd   1-by-M   cell, S.sd{m} is N-by-pm  (oracle conditional sds)
%
%   S = BBSAMPLE(N,P,RS) uses the RandStream RS, which is how the experiment
%   drivers keep replications reproducible.
%
%   The rows returned here are COMPLETE: every unit carries every attribute.
%   That is the latent population of Section 2, never observed in practice.
%   Use BBSAMPLEBLOCKS to impose Condition D1 and produce what the design
%   actually yields.
%
%   Cross-block dependence is controlled by P.theta.  With theta = 0 the
%   innovations are independent across blocks, so Cond_CI holds and the
%   product coupling P_perp coincides with P.  With theta ~= 0 every
%   block-wise marginal law L(X_B, X_Qm) is unchanged, but Cond_CI fails.
%   That pair is exactly the construction used by experiment E4.
%
%   See also BBPARAMS, BBEMBED, BBSAMPLEBLOCKS.

if nargin < 3 || isempty(rs)
    rs = RandStream.getGlobalStream();
end

Z  = rand(rs, n, P.d);
XB = bbEmbed(Z, P);

% Conditional means and standard deviations, evaluated exactly.
Rmat  = bbMeanMap(Z,P);          % n x q
SDmat = bbSdMap(Z,P);            % n x q

% Correlated standard normal innovations: E = randn * cholE'
E = randn(rs, n, P.q) * P.cholE.';

Xall = Rmat + SDmat .* E;        % n x q

XQ = cell(1,P.M); rC = cell(1,P.M); sdC = cell(1,P.M);
for m = 1:P.M
    idx    = (P.blockId == m);
    XQ{m}  = Xall(:,idx);
    rC{m}  = Rmat(:,idx);
    sdC{m} = SDmat(:,idx);
end

S = struct('Z',Z,'XB',XB,'XQ',{XQ},'r',{rC},'sd',{sdC},'n',n);
end

% -------------------------------------------------------------------------
function R = bbMeanMap(Z,P)
%BBMEANMAP  Oracle conditional mean r(z), n-by-q.
U = Z * P.Wr.' - P.tr.';                         % n x q
switch P.meanType
    case 'weier'
        % Weierstrass series  W(u) = sum_k b^{-k*beta} cos(2 pi b^k u + phi).
        % This is Holder-beta at EVERY point and differentiable at none, so
        % it saturates Cond_A uniformly on the support.  That property is
        % what makes it the right test function for the rate in T2: a mean
        % with an isolated singularity (see 'abspow' below) is Holder-beta
        % only on a null set, is smooth elsewhere, and the smooth part then
        % dominates the average, producing a FASTER decay than 2*beta/d_B.
        % Verifying the exponent with such a function reports a failure that
        % is an artefact of the test function, not of the theorem.
        W = zeros(size(U));
        for k = 0:P.weierTerms
            W = W + P.weierBase^(-k*P.beta) * ...
                cos(2*pi * P.weierBase^k * U + P.phir.');
        end
        R = (P.Lmean * P.sr.') .* W;
    case 'abspow'
        % |.|^beta is Holder-beta with constant 1, but the exponent is
        % attained only at the kink u = 0.  Retained deliberately as the
        % contrast case: it shows that the majorant of T2 is a worst-case
        % statement over the Holder class, not a prediction for every member
        % of it.
        R = (P.Lmean * P.sr.') .* (abs(U).^P.beta);
    case 'smooth'
        R = (P.Lmean * P.sr.') .* sin(2*pi*U) / (2*pi);
    otherwise
        error('bbSample:badMeanType','Unknown meanType "%s".',P.meanType);
end
end

% -------------------------------------------------------------------------
function SD = bbSdMap(Z,P)
%BBSDMAP  Oracle conditional standard deviation s(z), n-by-q, strictly > 0.
if P.hetero == 0
    SD = P.sigma0 * ones(size(Z,1),P.q);
else
    V  = sin(2*pi * (Z * P.Wsig.'));
    SD = P.sigma0 * (1 + P.hetero * V);
end
end
