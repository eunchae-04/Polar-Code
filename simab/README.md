# `simab/` — Scenario A/B 민감도 스윕 분석 (core의 확장판)

`core/`에서 다뤘던 "LLR 오추정(Scenario A)", "Design SNR 오추정(Scenario B)" 두 시나리오를,
**단일 값 하나가 아니라 넓은 범위로 스윕**하면서 오차가 커질수록 성능이 어떻게 나빠지는지
민감도(sensitivity) 곡선을 뽑아내는 확장판입니다. 이미지 전송 기능도 함께 유지되어,
스윕하는 각 지점마다 이미지 결과도 같이 저장됩니다.

---

## 파일 구성

| 파일                 | 설명                                                  |
| -------------------- | ----------------------------------------------------- |
| `polar_code_simAB.c` | 이 폴더의 메인 소스                                   |
| `image.png`          | (직접 준비 필요) 전송 테스트용 원본 이미지            |
| `result/`            | 실행 시 자동 생성. 스윕 결과 txt, 시나리오별 이미지들 |

---

## 환경 설정

```bash
gcc -std=c11 -O2 polar_code_simAB.c -lm -o polar_simAB.exe
```

`core/`와 동일하게 이미지 전송을 위해 추가 설정이 필요합니다.

1. Python 3 + Pillow (`pip install Pillow`)
2. 코드 상단 `PYTHON_CMD`를 자신의 Python 경로로 수정
   ```c
   #define PYTHON_CMD "C:/Users/user/anaconda3/python.exe"
   ```
3. 같은 폴더에 `image.png` 준비
4. Gnuplot이 PATH에 있어야 그래프 자동 표시

---

## 실행 방법

```bash
# 단발성 alpha 테스트: baseline 대비 SNR별 배율표 + 해당 alpha의 이미지 1장
./polar_simAB.exe --alpha 2.0

# Scenario A: alpha를 0.2~3.0(0.4 간격, 8개 지점)으로 훑기 — 표 + 이미지 8장 + 그래프
./polar_simAB.exe --sweep-a
./polar_simAB.exe --sweep-a 2.0     # 실제 채널 SNR을 기본값(3.0dB) 대신 다른 값으로 스윕

# Scenario B: design SNR을 -2.0~5.0dB(0.5dB 간격)으로 훑기 — 표 + 이미지 + 그래프
./polar_simAB.exe --sweep-b
./polar_simAB.exe --sweep-b 2.0

# 수동 단발 이미지 테스트 (시나리오 조건 없이 baseline LLR로만 전송)
./polar_simAB.exe --image input.png output.png [snr_db]
```

옵션 없이 실행하면 `core/`와 동일하게 Baseline/Scenario A(alpha=0.5 고정)/Scenario B
(design SNR=0dB 고정) 3종 비교가 수행됩니다.

---

## 핵심 이론

### Scenario A — alpha 스윕 (LLR 오추정 민감도)

```c
#define ALPHA_SWEEP_TRUE_SNR_DB DEFAULT_IMAGE_SNR_DB   // 실제 채널 SNR 고정
#define ALPHA_SWEEP_START 0.2
#define ALPHA_SWEEP_END   3.0
#define ALPHA_SWEEP_STEP  0.4
```

실제 채널 SNR은 고정한 채, 디코더가 채널 잡음 세기를 얼마나 잘못 추정하는지(`alpha =
sigma2_est / sigma2_true`)를 0.2배~3.0배까지 훑습니다. **alpha=1.0(정확한 추정)을 기준점**
으로 삼아, 그보다 작거나 큰 쪽으로 갈수록 BER이 얼마나 나빠지는지를 봅니다.

### Scenario B — Design SNR 스윕 (마스크 설계 오차 민감도)

```c
#define DESIGN_SNR_SWEEP_TRUE_DB 3.0
#define DESIGN_SNR_SWEEP_START (-2.0)
#define DESIGN_SNR_SWEEP_END   5.0
#define DESIGN_SNR_SWEEP_STEP  0.5
```

실제 채널 SNR은 고정한 채, Frozen Mask를 설계할 때 가정한 SNR(`design_snr`)을 -2.0dB~
5.0dB까지 훑습니다. **`design_snr == true_snr`인 지점(이상적으로 정확히 설계된 경우)을
기준점**으로 삼습니다.

### 공유 코어: `run_monte_carlo`

```c
static BerFerResult run_monte_carlo(const int *info_mask, double sigma2_true, double sigma2_est) {
    ...
}
```

`core/`에서는 3가지 시나리오마다 거의 같은 반복문이 중복되어 있었는데, 이 폴더에서는
"인코딩→채널→디코딩→채점" 한 프레임 분량의 로직을 `run_monte_carlo` 함수 하나로 뽑아내어
`run_simulation`, `--sweep-a`, `--sweep-b`가 모두 이 함수를 공유합니다. 스윕 지점이 많아질수록
(alpha 8개, design SNR 15개 등) 코드 중복을 줄이는 게 유지보수에 유리하기 때문입니다.

### 스윕 결과의 기준점(baseline) 찾기

```c
if (fabs(alpha_arr[i] - 1.0) < ALPHA_SWEEP_STEP / 2.0) { baseline_ber = alpha_ber[i]; break; }
```

부동소수점 값은 `==`로 직접 비교하면 반올림 오차 때문에 정확히 일치하지 않을 수 있어서,
"스텝 크기의 절반 이내로 가까우면 같은 값으로 취급"하는 방식으로 기준점(alpha=1.0 또는
design_snr=true_snr)을 찾습니다.

> Polar Code 기본 알고리즘(Frozen Mask, 인코더, SC 디코더)에 대한 설명은 `scd/` 및
> 최상위 `polar_SCD_code_explained.md`와 동일합니다.

---

## 출력 결과물

| 파일                                                                  | 내용                                                                                    |
| --------------------------------------------------------------------- | --------------------------------------------------------------------------------------- |
| `scenario_a_alpha_sweep_true*.txt`                                    | alpha별 BER/FER + baseline(alpha=1) 대비 배율                                           |
| `scenario_b_design_snr_sweep_true*.txt`                               | design SNR별 BER/FER + 기준점 대비 배율                                                 |
| `lena_output_alpha_sweep_*.png`, `lena_output_design_snr_sweep_*.png` | 스윕 지점별 전송 이미지                                                                 |
| Gnuplot 그래프                                                        | alpha(로그 스케일)/design SNR을 x축으로 한 BER/FER 민감도 곡선, 기준점 위치에 점선 표시 |

## 참고

`--alpha`, `--image` 단발 옵션과 기본 실행(Baseline/A/B 3종 비교)은 `core/`와 거의 동일한
동작을 하므로, 시나리오 개념 자체가 궁금하면 `core/README.md`를 먼저 보는 것을 추천합니다.
이 폴더는 그 개념을 "범위 스윕"으로 확장한 버전입니다.
