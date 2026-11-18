function out = runE2(quick)
%RUNE2  Driver for experiment E2 (Theorem T2).
%
%   OUT = RUNE2()      full run   (n up to 1e5; use this for the paper)
%   OUT = RUNE2(true)  quick run  (n up to 5e3; smoke test)
%
%   Layout: put every .m file in  ...\Exp2\Cod\  and results appear in
%   ...\Exp2\Result\ , created automatically.  See BBOUTDIR.
%
%   Outputs: E2.mat, E2_fig.pdf/.png, E2_table.txt
%
%   THE CONTRAST REGIME
%   One regime uses meanType 'abspow', a conditional mean with a single
%   kink.  It is Holder-beta at the kink and infinitely differentiable
%   everywhere else, so it belongs to the Holder-beta class without
%   saturating it.  Its excess decays FASTER than n^{-2*beta/d_B}.  This is
%   not a failure of T2: the majorant is a worst-case statement over the
%   class, and a member that is smooth off a null set does better than the
%   worst case.  The regime is included precisely so the paper can say this
%   with a measurement instead of a caveat.

if nargin < 1 || isempty(quick), quick = false; end

codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);
outDir = bbOutDir(mfilename('fullpath'));
fprintf('[E2] code   : %s\n[E2] results: %s\n\n', codeDir, outDir);

base = {'M',2,'pm',3,'sigma0',0.5,'rhoWithin',0.3,'theta',0};

R = struct('name',{},'params',{});
R(end+1) = struct('name','Weier  beta=0.40 d=2', ...
    'params',{[base {'d',2,'b',2,'meanType','weier','beta',0.40,'hetero',0.4}]});
R(end+1) = struct('name','Weier  beta=0.50 d=2', ...
    'params',{[base {'d',2,'b',2,'meanType','weier','beta',0.50,'hetero',0.4}]});
R(end+1) = struct('name','Weier  beta=0.75 d=2', ...
    'params',{[base {'d',2,'b',2,'meanType','weier','beta',0.75,'hetero',0.4}]});
R(end+1) = struct('name','Lipschitz beta=1 d=2', ...
    'params',{[base {'d',2,'b',2,'meanType','smooth','beta',1.00,'hetero',0.4}]});
R(end+1) = struct('name','Weier  beta=0.50 d=3', ...
    'params',{[base {'d',3,'b',3,'meanType','weier','beta',0.50,'hetero',0.4}]});
R(end+1) = struct('name','kink   beta=0.50 d=2', ...
    'params',{[base {'d',2,'b',2,'meanType','abspow','beta',0.50,'hetero',0.4}]});
R(end+1) = struct('name','Lipschitz beta=1 d=2 homo', ...
    'params',{[base {'d',2,'b',2,'meanType','smooth','beta',1.00,'hetero',0.0}]});

cfg = struct();
if quick
    cfg.nGrid = round(logspace(2,3.7,6)); cfg.repMax = 40; cfg.nRecv = 200;
    cfg.nOracle = 2e5;
else
    cfg.nGrid = round(logspace(2,5,10));  cfg.repMax = 200; cfg.nRecv = 500;
    cfg.nOracle = 2e6;
end

out = E2_riskDecomp(R,cfg);
out.outDir = outDir;

disp(' ');
disp('=== E2: rate of the excess (theory: slope = -2*beta/d_B) ===');
disp('    and relative error of the three identities of T2, in percent');
disp(out.tbl);

% ---- figure -----------------------------------------------------------
f = figure('Position',[100 100 1000 400]);

% (a) the three identities of T2, as relative discrepancy against the exact
% conditional expressions, each with its OWN standard error.
% What must not be plotted here is the empirical risk itself.  The two risks
% converge to floors that differ by a factor of two, and near the floor the
% remaining signal is a fraction of a percent while the sampling error of a
% mean of squared residuals is of the same order or larger; a plot of the
% risks shows noise oscillating around the limit and reads as a failure to
% converge.  The excess is measured precisely in panel (b), in closed form.
% What this panel adds is the part panel (b) cannot show: that the exact
% conditional identities hold at every n.
% The error bars are estimated from the spread across replications, not from
% a closed-form expression.  Receivers within a replication share the donor
% block and are therefore dependent, and the difference identity carries a
% substantially larger variance than the other two -- its discrepancy is
% driven by a cross term in the donor residual that the conditional
% identities do not contain.  A single shared band computed as sqrt(2/N)
% understates all three and understates the third badly.
gA = find(strcmp({R.name},'Lipschitz beta=1 d=2'),1);
if isempty(gA), gA = 4; end
s2 = out.floorS2(gA);
subplot(1,2,1); hold on; box on
yline(0,'-','Color',[.45 .45 .45],'LineWidth',0.8,'HandleVisibility','off');
errorbar(out.nGrid, 100*out.relMed(gA,:), 200*out.seMed(gA,:), ...
     'o-','Color',[0.00 0.45 0.74],'MarkerSize',4,'LineWidth',1.1, ...
     'CapSize',3,'DisplayName','conditional identity, mean');
errorbar(out.nGrid, 100*out.relVal(gA,:), 200*out.seVal(gA,:), ...
     'o-','Color',[0.85 0.33 0.10],'MarkerSize',4,'LineWidth',1.1, ...
     'CapSize',3,'DisplayName','conditional identity, value');
errorbar(out.nGrid, 100*out.relDiff(gA,:), 200*out.seDiff(gA,:), ...
     'o-','Color',[0.47 0.67 0.19],'MarkerSize',4,'LineWidth',1.1, ...
     'CapSize',3,'DisplayName','difference identity');
set(gca,'XScale','log','FontSize',9);
xlabel('n   (donor block size)');
ylabel('discrepancy from the exact identity  (%),  \pm2 s.e.');
text(0.03,0.94,'(a)','Units','normalized','FontSize',11,'FontWeight','bold');
legend('Location','best','FontSize',8,'Box','off');

% (b) the excess and its rate
subplot(1,2,2); hold on; box on
cols = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; ...
        0.49 0.18 0.56; 0.00 0.60 0.60; 0.64 0.08 0.18; 0.40 0.40 0.40];
for g = 1:numel(R)
    if strcmp(R(g).name,'Lipschitz beta=1 d=2 homo'), continue, end
    ls = '-';
    if contains(R(g).name,'kink'), ls = '-.'; end
    plot(out.nGrid, out.excess(g,:), ls, 'Color',cols(g,:), 'Marker','o', ...
         'MarkerSize',3.5,'LineWidth',1.1, ...
         'DisplayName',sprintf('%s  (%.2f / %.2f)', ...
                       R(g).name, out.slope(g), -2*out.beta(g)/out.dB(g)));
end
set(gca,'XScale','log','YScale','log','FontSize',9);
xlabel('n   (donor block size)');
ylabel('E[ ( r_a(Z) - r_a(X_B) )^2 ]');
text(0.03,0.06,'(b)','Units','normalized','FontSize',11,'FontWeight','bold');
legend('Location','southwest','FontSize',7,'Box','off');

exportgraphics(f, fullfile(outDir,'E2_fig.pdf'), 'ContentType','vector');
exportgraphics(f, fullfile(outDir,'E2_fig.png'), 'Resolution',200);

% ---- console summary of the identity checks --------------------------
fprintf('\n=== E2: identity checks (worst case over the sample-size grid) ===\n');
fprintf('  %-28s %10s %10s\n','identity','max |err|','max |z|');
fprintf('  %-28s %8.3f %% %9.2f\n','conditional, mean transfer', ...
        max(out.tbl.errMed_pct), max(out.tbl.zMed));
fprintf('  %-28s %8.3f %% %9.2f\n','conditional, value transfer', ...
        max(out.tbl.errVal_pct), max(out.tbl.zVal));
fprintf('  %-28s %8.3f %% %9.2f\n','difference', ...
        max(out.tbl.errDiff_pct), max(out.tbl.zDiff));
fprintf('  (z = discrepancy in units of its own standard error)\n');
fprintf('\n=== E2: limits at the largest n ===\n');
for g = 1:numel(R)
    f2 = out.floorS2(g);
    fprintf('  %-26s  med/s2 = %.4f   val/(2 s2) = %.4f\n', ...
        R(g).name, out.riskMed(g,end)/f2, out.riskVal(g,end)/(2*f2));
end

save(fullfile(outDir,'E2.mat'),'-struct','out');
fid = fopen(fullfile(outDir,'E2_table.txt'),'w');
fprintf(fid,'E2 - excess rate (theory: slope = -2*beta/d_B)\n');
fprintf(fid,'     identity errors in %%%% and in units of their own standard error (z)\n\n');
fprintf(fid,'%s\n',evalc('disp(out.tbl)'));
fprintf(fid,'\nLimits at the largest n (n = %d)\n',out.nGrid(end));
fprintf(fid,'  %-26s %10s %14s %12s\n','regime','s_a^2(B)','med / s2','val / (2 s2)');
for g = 1:numel(R)
    f2 = out.floorS2(g);
    fprintf(fid,'  %-26s %10.5f %14.4f %12.4f\n', ...
        R(g).name, f2, out.riskMed(g,end)/f2, out.riskVal(g,end)/(2*f2));
end
fclose(fid);

fprintf('\nWrote E2.mat, E2_fig.pdf, E2_fig.png, E2_table.txt to:\n  %s\n', outDir);
end
