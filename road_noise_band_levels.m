function [OA, Booming, Resonance, Rumble] = road_noise_band_levels(spl, opts)
%ROAD_NOISE_BAND_LEVELS  Road Noise 음압 요약값(O/A, Booming, Resonance, Rumble) 계산
%
%   [OA, Booming, Resonance, Rumble] = road_noise_band_levels(spl)
%   [...] = road_noise_band_levels(spl, Name=Value)
%
%   IDX Road Noise 결과 화면의 2x2 요약값과 같은 방식으로 계산한다.
%     1) Hann 창 보정 : 선형 응답 x 0.815 (= 1.63/2)  ->  -1.777 dB
%     2) A-weighting  : IEC 61672-1 보정값을 각 주파수 bin에 더함
%     3) 대역 에너지 합 : L = 10*log10( sum(10.^(L_k/10)) )   (N으로 나누지 않음)
%
%   입력
%     spl  : 비가중(Z) SPL [dB], 주파수 축 opts.Freq 와 같은 길이.
%            행렬이면 각 열을 별도 응답점(마이크)으로 계산한다.
%
%   옵션 (Name=Value)
%     Freq      : 주파수 축 [Hz]                 (기본 2:2:500)
%     Bands     : 경계 [A B C D] [Hz]            (기본 [20 160 210 500])
%     Weighting : "A" 또는 "Z"                   (기본 "A")
%     Window    : "hann_rss" 또는 "none"         (기본 "hann_rss")
%
%   대역 (내부 경계 B, C 의 샘플은 상위 대역에만 포함)
%     O/A       : A <= f <= D
%     Booming   : A <= f <  B
%     Resonance : B <= f <  C
%     Rumble    : C <= f <= D
%
%   구간에 유효한(finite) 샘플이 없으면 NaN 을 반환한다 (프로그램 표시 '—').
%
%   예)
%     [oa, bm, rs, rb] = road_noise_band_levels(spl);
%     fprintf('%.1f  %.1f  %.1f  %.1f dBA\n', oa, bm, rs, rb);

arguments
    spl {mustBeNumeric, mustBeReal}
    opts.Freq (:,1) double = (2:2:500).'
    opts.Bands (1,4) double {mustBeNonnegative} = [20 160 210 500]
    opts.Weighting (1,1) string {mustBeMember(opts.Weighting, ["A","Z"])} = "A"
    opts.Window (1,1) string {mustBeMember(opts.Window, ["hann_rss","none"])} = "hann_rss"
end

f = opts.Freq;
if isvector(spl)
    spl = spl(:);
end
if size(spl, 1) ~= numel(f)
    error('road_noise_band_levels:size', ...
        'spl 길이(%d)가 주파수 축 길이(%d)와 다릅니다.', size(spl, 1), numel(f));
end
if any(diff(opts.Bands) <= 0)
    error('road_noise_band_levels:bands', 'Bands 는 A < B < C < D 순서여야 합니다.');
end

% 1) Hann 창 보정 (선형 x 0.815 를 dB 로 더함)
L = double(spl);
if opts.Window == "hann_rss"
    L = L + 20*log10(1.63/2);
end

% 2) A-weighting
if opts.Weighting == "A"
    L = L + a_weighting_db(f);
end

% 3) 대역별 에너지 합
A = opts.Bands(1); B = opts.Bands(2); C = opts.Bands(3); D = opts.Bands(4);
OA        = band_energy_sum(L, f >= A & f <= D);
Booming   = band_energy_sum(L, f >= A & f <  B);
Resonance = band_energy_sum(L, f >= B & f <  C);
Rumble    = band_energy_sum(L, f >= C & f <= D);
end


function level = band_energy_sum(L, in_band)
% 열마다 대역 내 finite 샘플만 에너지 합산 (overflow 방지를 위해 최대값 기준 정규화)
level = nan(1, size(L, 2));
for j = 1:size(L, 2)
    v = L(in_band, j);
    v = v(isfinite(v));
    if isempty(v)
        continue
    end
    m = max(v);
    level(j) = m + 10*log10(sum(10.^((v - m)/10)));
end
end


function corr = a_weighting_db(f)
% IEC 61672-1:2013 A-weighting 보정값 [dB] (1 kHz 에서 0 dB)
f  = max(f, 1e-9);
f2 = f.^2;
f1 = 20.598997; fb = 107.65265; fc = 737.86223; f4 = 12194.217;
RA = (f4^2 .* f2.^2) ./ ((f2 + f1^2) .* sqrt(f2 + fb^2) .* sqrt(f2 + fc^2) .* (f2 + f4^2));
corr = 20*log10(RA) + 2.00;
end
