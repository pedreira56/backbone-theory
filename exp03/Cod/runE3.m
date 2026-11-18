function out = runE3(quick)
%RUNE3  Driver for experiment E3 (Proposition P4 and Theorem T3).
%
%   OUT = RUNE3()      full run   (n up to 1e5; use this for the paper)
%   OUT = RUNE3(true)  quick run  (n up to 5e3; smoke test)
%
%   Layout: put every .m file in  ...\Exp3\Cod\  and results appear in
%   ...\Exp3\Result\ , created automatically.  See BBOUTDIR.
%
%   Outputs: E3.mat, E3_fig.pdf/.png, E3_table.txt
%
%   WHAT THE FIGURE SHOWS
%   Panel (a) puts Proposition P4 and Theorem T3 side by side: the
%   distributional error of the value rule decays to zero at the rate the
%   theory gives, while the same error for the mean rule stops at a strictly
%   positive floor and never improves, however large the donor block.
%   Panel (b) is the dichotomy itself.  Each point is one value of k, and the
%   two axes are the two criteria; the curve is monotone in both, in
%   opposite directions.  No rule occupies the lower-left corner, and none
%   can: the two criteria are optimised by opposite ends of the family.

if nargin < 1 || isempty(quick), quick = false; end

codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);
outDir = bbOutDir(mfilename('fullpath'));
fprintf('[E3] code   : %s\n[E3] results: %s\n\n', codeDir, outDir);

base = {'M',2,'pm',3,'sigma0',0.5,'rhoWithin',0.3,'theta',0, ...
        'hetero',0.4,'sdShared',true};

R = struct('name',{},'params',{});
R(end+1) = struct('name','Lipschitz beta=1 d=2', ...
    'params',{[base {'d',2,'b',2,'meanType','smooth','beta',1.00}]});
R(end+1) = struct('name','Weier  beta=0.50 d=2', ...
    'params',{[base {'d',2,'b',2,'meanType','weier','beta',0.50}]});
R(end+1) = struct('name','Lipschitz beta=1 d=3', ...
    'params',{[base {'d',3,'b',3,'meanType','smooth','beta',1.00}]});

cfg = struct();
if quick
    cfg.nGrid = round(logspace(2,3.7,6)); cfg.repMax = 40; cfg.nRecv = 200;
    cfg.nOracle = 2e5;
else
    cfg.nGrid = round(logspace(2,5,10));  cfg.repMax = 200; cfg.nRecv = 500;
    cfg.nOracle = 2e6;
end

out = E3_dichotomy(R,cfg);
out.outDir = outDir;

disp(' ');
disp('=== E3: distributional rate for the value rule, and two limits ===');
disp('    riskRatio_val_med should tend to 2 (Theorem T2)');
disp('    w2med_over_floor  should tend to 1 (Proposition P4)');
disp(out.tbl);

gA = 1;                       % regime shown in the figure
jL = numel(out.nGrid);        % largest sample size

f = figure('Position',[100 100 1020 420]);

% ---- (a) P4 floor against T3 decay ------------------------------------
% Drawing order matters.  The claim of P4 is that the mean rule stops exactly
% AT the floor, so the floor must stay visible beneath the curve rather than
% be covered by it: it is drawn first, thick and grey, and the mean-transfer
% points go on top as open markers with no connecting line.  The coincidence
% is then something the reader sees rather than infers.
% No interpretive text is placed inside either panel: what the reader is
% meant to conclude belongs in the caption and in the surrounding prose, not
% painted onto the axes.
% Panel labels sit outside the axes box.  Inside, every corner is either
% occupied by a curve or wanted by the legend, and a label dropped into one
% of them collides with something -- as it did in an earlier version.
subplot(1,2,1); hold on; box on
wV = squeeze(out.w2K(gA,:,1));
plot(out.nGrid([1 end]), out.tau2E(gA)*[1 1], '-', ...
     'Color',[0.74 0.74 0.74],'LineWidth',5, ...
     'DisplayName','floor  E[\tau^2(X_B)]   (Proposition P4)');
plot(out.nGrid, wV, 'o-','Color',[0.85 0.33 0.10], ...
     'MarkerSize',4.5,'LineWidth',1.3,'MarkerFaceColor','w', ...
     'DisplayName',sprintf('value transfer   (slope %.2f / %.2f)', ...
                   out.slopeW(gA), -2*out.beta(gA)/out.dB(gA)));
plot(out.nGrid, out.w2M(gA,:), 'o','Color',[0.00 0.45 0.74], ...
     'MarkerSize',5.5,'LineWidth',1.4,'MarkerFaceColor','w', ...
     'DisplayName','mean transfer');
set(gca,'XScale','log','YScale','log','FontSize',9);
xlabel('n   (donor block size)');
ylabel('W_2^2 to the true conditional law');
ylim([min(wV)/4, out.tau2E(gA)*6]);
xlim([out.nGrid(1)/1.6, out.nGrid(end)*1.6]);
text(-0.02,1.05,'(a)','Units','normalized','FontSize',11, ...
     'FontWeight','bold','HorizontalAlignment','left');
legend('Location','southwest','FontSize',8,'Box','off');

% ---- (b) the dichotomy as a trade-off curve ---------------------------
% The two criteria on the two axes, one point per k.  The curve is monotone
% in both and in opposite directions, so the lower-left corner -- small on
% both criteria at once -- is empty, and no member of the family reaches it.
% That emptiness is the content of the panel and is stated in the caption
% rather than annotated on the axes.
subplot(1,2,2); hold on; box on
rk = squeeze(out.riskK(gA,jL,:));
wk = squeeze(out.w2K(gA,jL,:));
rM = out.riskM(gA,jL);  wM = out.w2M(gA,jL);
spanX = max([rk; rM]) - min([rk; rM]);
spanY = max([wk; wM])  - min([wk; wM]);
plot(rk, wk, '-','Color',[0.49 0.18 0.56],'LineWidth',1.4, ...
     'DisplayName','average of k nearest donors');
plot(rk, wk, 'o','Color',[0.49 0.18 0.56],'MarkerSize',5,'LineWidth',1.2, ...
     'MarkerFaceColor','w','HandleVisibility','off');
% The oracle mean is the limit of the same family, not a separate point: as
% k grows the average of k donors loses its dispersion and approaches a point
% mass at the conditional mean.  Left unconnected it reads as a stray marker,
% which is how an earlier version of this panel looked.
plot([rk(end) rM],[wk(end) wM],':','Color',[0.49 0.18 0.56], ...
     'LineWidth',1.2,'HandleVisibility','off');
plot(rM, wM, 's','Color',[0.00 0.45 0.74],'MarkerSize',9,'LineWidth',1.5, ...
     'MarkerFaceColor','w', ...
     'DisplayName','oracle conditional mean  (k \rightarrow \infty)');
set(gca,'FontSize',9);
xlabel('individual quadratic risk');
ylabel('W_2^2 to the true conditional law');
xlim([min([rk; rM])-0.12*spanX, max([rk; rM])+0.14*spanX]);
ylim([min([wk; wM])-0.10*spanY, max([wk; wM])+0.10*spanY]);
for kk = 1:numel(out.kList)
    if kk == 1, ha = 'right'; dx = -0.015*spanX; else, ha = 'left'; dx = 0.015*spanX; end
    text(rk(kk)+dx, wk(kk), sprintf('k=%d',out.kList(kk)), 'FontSize',7.5, ...
         'HorizontalAlignment',ha,'VerticalAlignment','middle');
end
text(-0.02,1.05,'(b)','Units','normalized','FontSize',11, ...
     'FontWeight','bold','HorizontalAlignment','left');
legend('Location','southwest','FontSize',8,'Box','off');

exportgraphics(f, fullfile(outDir,'E3_fig.pdf'), 'ContentType','vector');
exportgraphics(f, fullfile(outDir,'E3_fig.png'), 'Resolution',200);

% ---- console -----------------------------------------------------------
fprintf('\n=== E3: trade-off at n = %d, regime %s ===\n',out.nGrid(jL),R(gA).name);
fprintf('  %8s %20s %16s\n','k','individual risk','W2^2');
for kk = 1:numel(out.kList)
    fprintf('  %8d %20.6f %16.6f\n',out.kList(kk),rk(kk),wk(kk));
end
fprintf('  %8s %20.6f %16.6f\n','mean',out.riskM(gA,jL),out.w2M(gA,jL));
fprintf('  floor E[tau^2] = %.6f\n',out.tau2E(gA));

save(fullfile(outDir,'E3.mat'),'-struct','out');
fid = fopen(fullfile(outDir,'E3_table.txt'),'w');
fprintf(fid,'E3 - dichotomy between pointwise prediction and distributional reconstruction\n\n');
fprintf(fid,'%s\n',evalc('disp(out.tbl)'));
fprintf(fid,'\nTrade-off at n = %d, regime %s\n',out.nGrid(jL),R(gA).name);
fprintf(fid,'  %8s %20s %16s\n','k','individual risk','W2^2');
for kk = 1:numel(out.kList)
    fprintf(fid,'  %8d %20.6f %16.6f\n',out.kList(kk),rk(kk),wk(kk));
end
fprintf(fid,'  %8s %20.6f %16.6f\n','mean',out.riskM(gA,jL),out.w2M(gA,jL));
fprintf(fid,'  floor E[tau^2] = %.6f\n',out.tau2E(gA));
fclose(fid);

fprintf('\nWrote E3.mat, E3_fig.pdf, E3_fig.png, E3_table.txt to:\n  %s\n', outDir);
end
