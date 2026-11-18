function out = runE1(quick)
%RUNE1  Driver for experiment E1 (Theorem T1).
%
%   OUT = RUNE1()      full run   (n up to 1e5; use this for the paper)
%   OUT = RUNE1(true)  quick run  (n up to 5e3; smoke test, ~1 min)
%
%   Layout: put every .m file in  ...\Exp1\Cod\  and results appear in
%   ...\Exp1\Result\ , created automatically.  See BBOUTDIR.
%
%   Outputs written to the Result folder:
%     E1.mat            raw output
%     E1_fig.pdf/.png   two-panel figure (rate curves | tail bound)
%     E1_table.txt      the summary table, as printed
%
%   FIGURE DESIGN
%   Panel (a) shows the q = 1 rate curves.  Panel (b) shows the tail bound
%   against its theoretical constant.  The two panels verify different
%   things: (a) checks the exponent across sample sizes, (b) checks the
%   exponential form AND the constant within a single sample size, using the
%   whole distribution rather than a mean.  The q = 2 moment is not plotted:
%   it would repeat panel (a) with every slope doubled.  It goes in the table
%   instead, where its role -- an internal consistency check, since the q = 2
%   slope must be exactly twice the q = 1 slope -- is what it can actually do.
%
%   COLOUR SCHEME
%   Colour encodes the geometric dimension d_B; line style encodes the
%   geometry (solid = cube, dashed = torus, dash-dot = curved embedding).
%   Two curves therefore share a colour precisely when they share a rate.
%   The pair 'cube d=b=2' and 'graph d=2 b=10' comes out same colour,
%   different style, parallel -- which is the claim of T1 rendered directly
%   in the plot: same dimension, same rate, cardinality irrelevant.
%   (MATLAB's default 'lines' palette cycles after 7 entries and would assign
%   those two curves the identical colour by accident, which is worse than
%   useless: it reads as a duplicated plot.)

if nargin < 1 || isempty(quick), quick = false; end

codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);
outDir = bbOutDir(mfilename('fullpath'));
fprintf('[E1] code   : %s\n[E1] results: %s\n\n', codeDir, outDir);

% ---- regimes ----------------------------------------------------------
% Cube and torus are PAIRED at d = 3,5,8.  The cube has a boundary, the flat
% torus does not; the pair isolates the boundary effect and lets the paper
% show, rather than assert, that the deviation seen on the cube in high
% dimension is a boundary and finite-n artefact and not a failure of T1.
%
% The 'graph' regimes use Lg = 6.  With a nearly flat embedding the curve for
% (d=2,b=10) lands on top of the curve for the plain square, because the
% embedding then almost preserves distances and the CONSTANT coincides too.
% A genuinely curved embedding separates the constants while leaving the
% slope untouched.
R = struct('name',{},'params',{});
R(end+1) = struct('name','cube  d=b=1',   'params',{{'d',1,'b',1}});
R(end+1) = struct('name','cube  d=b=2',   'params',{{'d',2,'b',2}});
R(end+1) = struct('name','cube  d=b=3',   'params',{{'d',3,'b',3}});
R(end+1) = struct('name','torus d=b=3',   'params',{{'d',3,'b',3,'geometry','torus'}});
R(end+1) = struct('name','cube  d=b=5',   'params',{{'d',5,'b',5}});
R(end+1) = struct('name','torus d=b=5',   'params',{{'d',5,'b',5,'geometry','torus'}});
R(end+1) = struct('name','cube  d=b=8',   'params',{{'d',8,'b',8}});
R(end+1) = struct('name','torus d=b=8',   'params',{{'d',8,'b',8,'geometry','torus'}});
R(end+1) = struct('name','graph d=2 b=10','params',{{'d',2,'b',10,'embedding','graph','kGraph',2,'Lg',6}});
R(end+1) = struct('name','graph d=4 b=12','params',{{'d',4,'b',12,'embedding','graph','kGraph',2,'Lg',6}});

cfg = struct();
if quick
    cfg.nGrid  = round(logspace(2,3.7,6));
    cfg.repMax = 40;
    cfg.nRecv  = 200;
else
    cfg.nGrid  = round(logspace(2,5,10));
    cfg.repMax = 200;
    cfg.nRecv  = 500;
end

out = E1_nnRate(R,cfg);
out.outDir = outDir;

% ---- consistency check: slope(q=2) must be 2 x slope(q=1) -------------
k1 = find(out.qList == 1, 1);
k2 = find(out.qList == 2, 1);
if ~isempty(k1) && ~isempty(k2)
    out.tbl.ratio_q2_q1 = out.slope(:,k2)./out.slope(:,k1);
end

disp(' ');
disp('=== E1: fitted exponents (theory: slope = -q/d_B, ratio = 2) ===');
disp(out.tbl);

% ---- colour by d_B, line style by geometry ---------------------------
dVals  = [1 2 3 4 5 8];
dCols  = [0.00 0.45 0.74;      % d = 1  blue
          0.85 0.33 0.10;      % d = 2  orange
          0.47 0.67 0.19;      % d = 3  green
          0.49 0.18 0.56;      % d = 4  purple
          0.00 0.60 0.60;      % d = 5  teal
          0.64 0.08 0.18];     % d = 8  dark red

f = figure('Position',[100 100 1000 400]);

% ---- panel (a): rate curves, q = 1 -----------------------------------
subplot(1,2,1); hold on; box on
for g = 1:numel(R)
    P  = bbParams(R(g).params{:});
    c  = dCols(dVals == P.d, :);
    switch P.embedding
        case 'graph',                 ls = '-.';
        otherwise
            if strcmp(P.geometry,'torus'), ls = '--'; else, ls = '-'; end
    end
    plot(out.nGrid, squeeze(out.moment(g,:,k1)), ls, 'Color',c, ...
         'Marker','o', 'MarkerSize',3.5, 'LineWidth',1.1, ...
         'DisplayName',sprintf('%s  (dhat = %.2f)',R(g).name,out.dHat(g,k1)));
end
set(gca,'XScale','log','YScale','log','FontSize',9);
xlabel('n   (donor block size)');
ylabel('E[ R ]');
text(0.03,0.06,'(a)','Units','normalized','FontSize',11,'FontWeight','bold');
legend('Location','southwest','FontSize',7,'Box','off');

% ---- panel (b): tail bound against the theoretical constant ----------
% T1 gives P(R > eps) <= exp(-n c_B eps^{d_B}), so
%     -log P(R > eps)/n  >=  c_B eps^{d_B},
% a straight line through the origin when plotted against eps^{d_B}.  On the
% flat torus with uniform measure the constant is not merely SOME c_B: it is
% the volume of the unit ball, V_d = pi^{d/2}/Gamma(d/2+1).  Overlaying that
% line turns the panel from a check of the exponent into a check of the
% exponent and the constant together.
gSel = find(strcmp({R.name},'torus d=b=3'),1);
if isempty(gSel), gSel = 3; end
jSel = numel(out.nGrid);
Rv   = out.R{gSel,jSel};
n    = out.nGrid(jSel);
dB   = out.dTrue(gSel);
Vd   = pi^(dB/2)/gamma(dB/2 + 1);

eps  = linspace(quantile(Rv,0.02), quantile(Rv,0.98), 40);
surv = arrayfun(@(e) mean(Rv > e), eps);
keep = surv > 0 & surv < 1;
xx   = eps(keep).^dB;
yy   = -log(surv(keep))/n;
cFit = xx(:)\yy(:);

subplot(1,2,2); hold on; box on
plot(xx, Vd*xx, '--', 'Color',[.45 .45 .45], 'LineWidth',1.2, ...
     'DisplayName',sprintf('theoretical  V_d = %.3f',Vd));
plot(xx, yy, 'o', 'Color',dCols(dVals == dB,:), 'MarkerSize',3.5, ...
     'LineWidth',1, 'DisplayName',sprintf('empirical  c_B = %.3f',cFit));
set(gca,'FontSize',9);
xlabel(sprintf('\\epsilon^{d_B}     (torus, d_B = %g, n = %d)',dB,n));
ylabel('-log P(R > \epsilon) / n');
text(0.03,0.94,'(b)','Units','normalized','FontSize',11,'FontWeight','bold');
legend('Location','northwest','FontSize',8,'Box','off');

exportgraphics(f, fullfile(outDir,'E1_fig.pdf'), 'ContentType','vector');
exportgraphics(f, fullfile(outDir,'E1_fig.png'), 'Resolution',200);

out.tail = struct('regime',R(gSel).name,'n',n,'dB',dB, ...
                  'cFit',cFit,'Vtheory',Vd,'eps',eps,'surv',surv);

fprintf('\n=== E1: tail bound (%s, n = %d) ===\n',R(gSel).name,n);
fprintf('  empirical constant   c_B = %.4f\n',cFit);
fprintf('  theoretical constant V_d = %.4f   (volume of the unit ball in R^%g)\n',Vd,dB);
fprintf('  relative discrepancy     = %.2f%%\n',100*(cFit-Vd)/Vd);

% ---- persist ----------------------------------------------------------
save(fullfile(outDir,'E1.mat'),'-struct','out');
fid = fopen(fullfile(outDir,'E1_table.txt'),'w');
fprintf(fid,'E1 - fitted exponents (theory: slope = -q/d_B, ratio = 2)\n\n');
fprintf(fid,'%s\n',evalc('disp(out.tbl)'));
fprintf(fid,'\nTail bound (%s, n = %d)\n',R(gSel).name,n);
fprintf(fid,'  empirical   c_B = %.4f\n  theoretical V_d = %.4f\n  discrepancy     = %.2f%%\n', ...
        cFit,Vd,100*(cFit-Vd)/Vd);
fclose(fid);

fprintf('\nWrote E1.mat, E1_fig.pdf, E1_fig.png, E1_table.txt to:\n  %s\n', outDir);
end
