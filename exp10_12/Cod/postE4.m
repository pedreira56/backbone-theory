function B = postE4(hint, nOracle)
%POSTE4  Bounds of Observation O2 for the population of experiment E4.
%
%   B = POSTE4() finds E4.mat (see BBFINDMAT) and computes, for the pair of
%   attributes followed in E4 (the first exclusive attribute of each block),
%   the interval of cross-block correlations compatible with the observable
%   data, writing E4_bounds.txt next to the .mat.
%   B = POSTE4(FOLDER) searches FOLDER first.
%   B = POSTE4(FOLDER,NORACLE) uses NORACLE draws (default 2e6).
%
%   Needs only this file and bbFindMat.m: the population parameters are read
%   from the saved struct and the oracle maps are evaluated here, so the
%   generator code need not be on the path.
%
%   THE BOUND
%   For attributes a and a' in different blocks,
%       Cov(Xa,Xa') = Cov(ra(XB),ra'(XB)) + E[ Cov(Xa,Xa' | XB) ].
%   The first term is identified.  The second is not, but by Cauchy-Schwarz
%   its integrand is at most sd_a(XB) sd_a'(XB) in absolute value, with
%   equality attainable when the conditional laws are Gaussian.  Hence
%       | Cov(Xa,Xa') - Cov(ra,ra') |  <=  E[ sd_a(XB) sd_a'(XB) ].
%   Every quantity on the right is computable from block-wise data: the
%   conditional mean and dispersion of block 1 from block 1, those of block 2
%   from block 2, and the backbone law from either.  Here they are evaluated
%   from the oracle maps, which is what a consistent block-wise estimator
%   would converge to.

if nargin < 1, hint = ''; end
if nargin < 2 || isempty(nOracle), nOracle = 2e6; end
addpath(fileparts(mfilename('fullpath')));
[matFile, outDir] = bbFindMat('E4', hint);
E4 = load(matFile);
P  = E4.Pprime;                          % the dependent population
pm = P.pm;

rs = RandStream('mt19937ar','Seed',P.seed + 12345);
Z  = rand(rs, nOracle, P.d);             % the backbone law, in latent form
R  = oracleMean(Z,P);                    % conditional means
G  = oracleSd(Z,P);                      % conditional standard deviations

i1 = 1;  i2 = pm + 1;                    % first attribute of each block
r1 = R(:,i1);  r2 = R(:,i2);
g1 = G(:,i1);  g2 = G(:,i2);

c    = cov(r1,r2);  covR = c(1,2);       % identified part
half = mean(g1.*g2);                     % half-width, covariance units
den  = sqrt((var(r1) + mean(g1.^2)) * (var(r2) + mean(g2.^2)));

% The population value, for reference: with correlated innovations the
% conditional covariance is theta * sd1 * sd2.
covTrue = covR + P.theta*half;

B = struct('covIdentified',covR,'halfWidthCov',half, ...
           'corrLower',(covR - half)/den,'corrUpper',(covR + half)/den, ...
           'corrCentre',covR/den,'corrTrue',covTrue/den, ...
           'theta',P.theta,'nOracle',nOracle);

lines = {
 sprintf('E4 - bounds of Observation O2 (oracle sample of %d)',nOracle)
 ''
 sprintf('  identified part  Cov(r1(XB),r2(XB))      : %+.5f',B.covIdentified)
 sprintf('  half-width       E[sd1(XB) sd2(XB)]      : %.5f',B.halfWidthCov)
 ''
 '  in correlation units'
 sprintf('  interval compatible with the data        : [%+.4f, %+.4f]',B.corrLower,B.corrUpper)
 sprintf('  centre (value under Cond_CI)             : %+.4f',B.corrCentre)
 sprintf('  population value (theta = %.2f)          : %+.4f',B.theta,B.corrTrue)
 };
fprintf('\n'); fprintf('%s\n',lines{:});
fid = fopen(fullfile(outDir,'E4_bounds.txt'),'w');
fprintf(fid,'%s\n',lines{:});
fclose(fid);
fprintf('\nWrote E4_bounds.txt to:\n  %s\n',outDir);
end

% =======================================================================
function R = oracleMean(Z,P)
%ORACLEMEAN  Conditional mean r(z), as used by BBSAMPLE.
U = Z * P.Wr.' - P.tr.';
switch P.meanType
    case 'weier'
        R = zeros(size(U));
        for k = 0:P.weierTerms
            R = R + P.weierBase^(-k*P.beta) * ...
                cos(2*pi * P.weierBase^k * U + P.phir.');
        end
        R = (P.Lmean * P.sr.') .* R;
    case 'abspow'
        R = (P.Lmean * P.sr.') .* (abs(U).^P.beta);
    case 'smooth'
        R = (P.Lmean * P.sr.') .* sin(2*pi*U) / (2*pi);
    otherwise
        error('postE4:badMeanType','Unknown meanType "%s".',P.meanType);
end
end

function SD = oracleSd(Z,P)
%ORACLESD  Conditional standard deviation s(z), as used by BBSAMPLE.
if P.hetero == 0
    SD = P.sigma0 * ones(size(Z,1),P.q);
else
    SD = P.sigma0 * (1 + P.hetero * sin(2*pi * (Z * P.Wsig.')));
end
end
