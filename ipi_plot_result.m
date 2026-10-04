%IPI_PLOT_RESULT  저장된 최적화 결과(S, R)를 불러와 지정한 노드/방향의 IPI 곡선만 그린다.
%
%   ipi_report 의 그림 규약을 그대로 따른다:
%       IPI = a/F [mm/s²/N],   dB(IPI) = 20·log10(|a/F| / 1e-3)   [dB/N ref 1e-3 mm/s²]
%       K_d 는 |a/F| = ω²/K_d 기준곡선(점선)으로 표시
%
%   사용법 : 아래 [INPUT] 블록만 수정하고 스크립트를 실행한다.
%   결과 파일은 예를 들어  save('result.mat', 'S', 'R')  로 저장된 것이어야 한다.
%   (S = ipi_setup 출력, R = ipi_optimize 출력)
%
%   주의 : 파일명을 plot.m 으로 두면 MATLAB 내장 plot 함수를 가려서
%          이 스크립트와 ipi_report 의 plot(...) 호출이 모두 깨진다.

clear; clc;

%% ===================== [INPUT] =====================
resultFile = 'result.mat';      % 결과 파일 경로 (문자열). load(...) 는 붙이지 않는다
sVar       = 'S';               % 파일 안의 setup 변수 이름
rVar       = 'R';               % 파일 안의 결과 변수 이름

% 출력할 노드 ID (적은 순서대로 그린다).  [] 이면 전체 노드
plotNodes  = [];
% 출력할 방향  1 = X, 2 = Y, 3 = Z.       [] 이면 전체 방향
plotDofs   = [1 2 3];

% subplot 배치 [행 열].  선택한 개수가 칸 수보다 많으면 다음 figure 로 넘어간다
subGrid    = [3 3];

printTable = true;              % K_d 표를 명령창에 출력
figName    = 'IPI optimization';
saveFig    = false;             % true 이면 figure 를 png 로 저장
saveDir    = 'fig';
%% ===================================================

% ---- 결과 불러오기 ----
if isstruct(resultFile)                    % load(...) 결과를 그대로 넣은 경우도 허용
    D = resultFile;
else
    D = load(resultFile);
end
assert(isfield(D, sVar) && isfield(D, rVar), ...
       '%s 안에 변수 ''%s'', ''%s'' 가 모두 있어야 한다 (현재: %s)', ...
       resultFile, sVar, rVar, strjoin(fieldnames(D)', ', '));
S = D.(sVar);  R = D.(rVar);
P = S.P;  nt = size(P.tgt, 1);

% ---- 출력 대상 선택 (노드 순서 -> 방향 순서) ----
if isempty(plotNodes), plotNodes = unique(P.tgt(:,1), 'stable')'; end
if isempty(plotDofs),  plotDofs  = 1:3; end
sel = [];
for n = plotNodes(:)'
    if ~any(P.tgt(:,1) == n)
        warning('노드 %d 는 목표 목록(S.P.tgt)에 없다 -> 건너뜀', n);
        continue;
    end
    for d = plotDofs(:)'
        sel = [sel; find(P.tgt(:,1) == n & P.tgt(:,2) == d)];   %#ok<AGROW>
    end
end
assert(~isempty(sel), '선택된 노드/방향 조합이 없다');

% ---- 공통 계산 (ipi_report 와 동일) ----
DBREF = 1e-3;                         % mm/s²  (사내 dB ref)
G     = 9.80665;                      % N/mm -> kgf/mm
dB    = @(A) 20*log10(max(abs(A), realmin) / DBREF);
dBK   = @(Kd) 20*log10(P.w.^2 ./ Kd / DBREF);

A1  = ipi_apply_bush(P, ipi_bush_impedance(R.Kfull, R.Bfull, R.GEfull, P.w));
Kb  = ipi_dynamic_stiffness(P.Ytt, P.w);          % 부시 없음 (baseline)
Kdb = ipi_log_average(Kb, S.kband, S.minRatio);

% ---- 표 ----
if printTable
    fprintf('\n=== K_d  대역 %g~%g Hz   [N/mm]  (kgf/mm = N/mm / %.5g) ===\n', S.band, G);
    fprintf('%-10s %3s | %9s %9s %9s %9s | %8s %8s %8s | %7s %5s\n', 'node', 'dof', ...
            'baseline', 'init', 'opt', 'target', 'base kgf', 'opt kgf', 'tgt kgf', 'err[%]', 'valid');
    for i = sel'
        fprintf('%-10d %3d | %9.0f %9.0f %9.0f %9.0f | %8.0f %8.0f %8.0f | %+7.2f %5.2f\n', ...
                P.tgt(i,1), P.tgt(i,2), Kdb(i), R.Kd0(i), R.Kd(i), R.T(i), ...
                Kdb(i)/G, R.Kd(i)/G, R.T(i)/G, 100*R.e(i), R.ratio(i));
    end
    fprintf('J : init %.4e -> opt %.4e   (%d / %d 개 출력)\n', R.J0, R.J, numel(sel), nt);
end

% ---- 그림 ----
nr = subGrid(1);  nc = subGrid(2);  npf = nr*nc;
nfig = ceil(numel(sel) / npf);
kb = S.kband;
kp = find(P.w > 0);                   % 0 Hz 는 dB 정의 불가
if saveFig && ~exist(saveDir, 'dir'), mkdir(saveDir); end

for f = 1:nfig
    if nfig > 1, name = sprintf('%s (%d/%d)', figName, f, nfig); else, name = figName; end
    hf = figure('Name', name);
    idx = sel((f-1)*npf+1 : min(f*npf, numel(sel)));
    for k = 1:numel(idx)
        i = idx(k);
        subplot(nr, nc, k);
        yb = dB(P.Ytt(i,:));  yo = dB(A1(i,:));
        rb = dBK(Kdb(i));  ro = dBK(R.Kd(i));  rt = dBK(R.T(i));
        ylo = floor(min([yb(kp) yo(kp) rt(kb)]) / 10) * 10 - 10;
        yhi = ceil (max([yb(kp) yo(kp) rt(kb)]) / 10) * 10 + 10;
        hp = patch([S.band fliplr(S.band)], [ylo ylo yhi yhi], [0.92 0.92 0.92], ...
                   'EdgeColor', 'none', 'FaceAlpha', 0.5);  hold on;
        h1 = plot(P.freq(kp), yb(kp), 'k',  'LineWidth', 1);
        h2 = plot(P.freq(kp), yo(kp), 'r',  'LineWidth', 1);
        h3 = plot(P.freq(kb), rb(kb), 'k:', 'LineWidth', 1);
        h4 = plot(P.freq(kb), ro(kb), 'r:', 'LineWidth', 1);
        h5 = plot(P.freq(kb), rt(kb), 'b--','LineWidth', 1);
        grid on; xlim([P.freq(1) P.freq(end)]); ylim([ylo yhi]);
        xlabel('Freq [Hz]'); ylabel('IPI  [dB/N ref 1e-3 mm/s^2]');
        title(sprintf('N%d ACC %s : KD %.0f -> %.0f  (T %.0f) kgf/mm', P.tgt(i,1), ...
                      char('X'+P.tgt(i,2)-1), Kdb(i)/G, R.Kd(i)/G, R.T(i)/G), 'FontSize', 8);
        if k == 1
            legend([h1 h2 h3 h4 h5 hp], {'baseline (no bush)', 'optimized', 'K_d base', ...
                   'K_d opt', 'K_d target', 'band'}, 'Location', 'southeast', 'FontSize', 7);
        end
    end
    if saveFig
        saveas(hf, fullfile(saveDir, sprintf('%s_%02d.png', strrep(figName, ' ', '_'), f)));
    end
end
