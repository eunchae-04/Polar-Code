# Polar Code (N=1024, K=512) Simulation Platform
본 플랫폼은 가우시안 근사(Gaussian Approximation, GA) 기반의 부호 구축과 최적화된 Exact SPA 연속 제거(Successive Cancellation, SC) 복호화 알고리즘을 64비트 고속 난수 엔진 기반의 C언어 환경으로 구현한 고정밀 통신 시뮬레이터입니다.

프로젝트는 목적에 따라 실행 파일 두 개로 나뉘어 있습니다.

- **`polar_code.c`**: Baseline / Scenario A / Scenario B 3가지를 고정 조건으로 한 번에 비교하는 기본 실행 파일. BER/FER 그래프, 시나리오별 이미지 전송 테스트, 비교 요약표를 자동으로 만듭니다.
- **`polar_code_simAB.c`**: Scenario A(α)와 Scenario B(Design SNR)를 각각 원하는 값/범위로 훑어보는 심화 분석 전용 실행 파일. 단발 테스트(`--alpha`)와 다중 값 스윕(`--sweep-a`, `--sweep-b`)을 지원합니다.



## 1. 핵심 이론 및 알고리즘 명세

### 1) 가우시안 근사 (Gaussian Approximation) 기반 부호 구축
채널 분극화(Channel Polarization) 현상에 의해 분할되는 서브 채널들의 신뢰도를 정밀하게 추적하기 위해 밀도 진화(Density Evolution)의 가우시안 근사 함수 $\phi(x)$와 이진 탐색 기반의 역함수 $\phi^{-1}(y)$ 알고리즘을 수행합니다.

$$\phi(x) = \begin{cases} 1.0, & \text{if } x \le 0 \\ exp(-0.4527x^{0.86} + 0.0218), & \text{if } 0 < x \le 10 \\ \sqrt{\frac{\pi}{x}}\left(1 - \frac{10}{7x}\right)exp\left(-\frac{x}{4}\right), & \text{if } x > 10 \end{cases}$$

- **Bit-Reversal Index 매핑**: 하향식(Bottom-Up) 나비 구조로 계산된 신뢰도 벡터 mu를 물리 트리 구조와 일치시키기 위해 비트 역전 변환을 거쳐 Natural Order로 재배치합니다.
- **Frozen Mask 확정**: 정렬된 신뢰도 지수를 바탕으로 상위 K=512개의 채널을 정보 비트(Information) 영역으로 지정하고, 나머지 채널은 수렴 능력이 떨어지는 가상 채널로 판단하여 0(Frozen)으로 잠금 처리합니다.

### 2) Top-Down 재귀형 부호화기 (Natural Order Encoder)
생성 행렬 $G_N$의 직접적인 대규모 행렬 곱 연산을 배제하고, 복호화 트리 구조와 완벽한 대칭성을 가지는 상향식 재귀 인코딩을 구현합니다. 상위 블록 입력 u에 대해 좌측 하부 브랜치는 $v = u_{left} \oplus u_{right}$ 커널 연산을 거치며, 최종 부호화된 코드 워드 x는 최하단 레이어에서 단순 연결(Plain Concat)되어 송신 심볼로 출력됩니다.

### 3) Exact SPA 기반 연속 제거(SC) 복호화기
하드웨어 구현 편의성을 위한 Min-Sum 근사 방식의 성능 열화를 극복하기 위해, 수학적 연속성을 완벽히 보장하는 야코비안 로그(Jacobian Logarithm) 정밀 보정 항이 결합된 **Exact SPA(Sum-Product Algorithm)** 연산 레이어를 구현했습니다.

- **f-function (Left Path)**:
  $$LLR_{left} = \text{sign}(l_1) \cdot \text{sign}(l_2) \cdot \min(\vert{}l_1\vert{}, \vert{}l_2\vert{}) + \ln(1 + e^{-\vert{}l_1+l_2\vert{}}) - \ln(1 + e^{-\vert{}l_1-l_2\vert{}})$$
- **g-function (Right Path)**: 이전 좌측 트리에서 복호 완료되어 피드백된 부분 코드 워드($u_{coded\_left}$)를 반영하여 우측 채널 우도를 갱신합니다.
  $$LLR_{right} = l_2 + (1 - 2 \cdot u_{coded\_left}) \cdot l_1$$

### 4) Xorshift64 난수 엔진 및 AWGN 채널
기존 C 표준 라이브러리의 `rand()` 함수가 가진 주기적 한계치(32,767)와 의사 난수 뭉침 현상에 따른 시뮬레이션 데이터 왜곡을 방지하기 위해 **Xorshift64 무작위 수 생성기**를 전면 탑재했습니다.
이를 박스-뮬러(Box-Muller) 변환과 결합하여, 높은 SNR 대역에서도 통계적 무결성과 독립성을 유지하는 정밀한 백색 가우시안 잡음(AWGN) 환경을 제공합니다.



## 2. 연구용 비정상 시나리오 해석 기능

본 시뮬레이터는 완벽한 채널 환경(Baseline) 외에도, 실제 수신기 칩셋 구현 시 마주할 수 있는 하드웨어 제약 및 오차 상황을 정밀 모델링하기 위해 2가지 학술적 분석 시나리오를 지원합니다.

### 1) 시나리오 A: 복호기 입력 LLR의 스케일링 오차 (LLR Scaling Error)
실제 채널에 흐르는 진짜 잡음 분산($\sigma^2_{true}$)과 수신기가 추정 수식에서 믿고 있는 잡음 분산($\sigma^2_{est}$) 사이에 의도적인 오차 계수($\alpha$)를 개입시켜 초기 LLR 스케일을 왜곡하는 상황입니다.

$$\sigma^2_{est} = \alpha \times \sigma^2_{true}$$

- **실측 결과 ($\alpha = 0.5$, `polar_code.c` 고정 조건)**: LLR 전체 크기를 2배로 확대시키는 방향의 오차입니다. 이 경우 Exact SPA의 비선형 보정 항($\ln(1+e^{-|l_1 \pm l_2|})$)은 입력 크기가 커질수록 값이 0에 수렴해, 복호기가 점점 **Min-Sum 근사**와 동일하게 동작하게 됩니다. Min-Sum 결합은 부호(sign) 판정이 양(+)의 균일 스케일링에 불변이므로, 이 방향의 미스매치는 이론적 우려와 달리 **BER을 오히려 개선**시키는 것으로 실측되었습니다. `scenario_comparison_summary.txt` 기준 SNR이 높아질수록 baseline 대비 오차율이 계속 낮아지며(0dB: 0.99배 → 3.0dB: 0.069배), 애초에 우려했던 "2.5~2.7dB 부근 bump(역전) 현상"은 관측되지 않았습니다.
- **α>1(축소) 방향 검증 (`polar_code_simAB.c --alpha 2.0`)**: α=2.0(LLR이 절반으로 축소되는 방향)에서는 정반대로 이미지가 전면적으로 심하게 깨지는 것을 확인했습니다. LLR이 작아질수록 비선형 보정 항의 상대적 비중이 커지고 전체 신뢰도 자체가 낮게 추정되어, SC 복호기의 오류 전파(error propagation)가 훨씬 심해지기 때문입니다. 즉 **α<1(과대평가) 방향과 α>1(과소평가) 방향은 완전히 비대칭적**입니다.
- **α 전 구간 스윕 (`polar_code_simAB.c --sweep-a`)**: α = 0.2~3.0(0.4 간격)을 한 번에 훑으며 BER/FER/이미지를 모두 생성해, 이 비대칭성을 α=1(baseline) 기준 좌우로 한눈에 비교할 수 있습니다.

### 2) 시나리오 B: 부호 구축 시의 Design SNR 미스매치 (Design SNR Mismatch)
채널 변동에 맞춰 프로즌 마스크를 동적으로 변경하지 못하고, 특정 설계 환경(Design SNR)에서 고정 생성된 단 하나의 마스크를 사용하여 채널 환경이 다른 대역까지 강제 복호하는 하드웨어 최적화 제약 상황입니다.

- **실측 결과 (Design SNR = 0.0dB 고정, `polar_code.c` 고정 조건)**: Polar Code가 완전한 채널 의존적 부호(Channel-Dependent Code)임을 뒷받침하는 결과가 확인되었습니다. 0.0 dB 환경에 동결된 마스크는 실제 채널 SNR이 설계 SNR에서 멀어질수록(=신호가 좋아질수록) baseline 대비 성능 열화가 급격히 커집니다. `scenario_comparison_summary.txt` 기준 baseline 대비 오차율이 0dB에서 1.00배 수준이던 것이 3.0dB에서는 38.85배까지 벌어집니다.
- **Design SNR 전 구간 스윕 (`polar_code_simAB.c --sweep-b`, 실제 채널 3.0dB 고정)**: 실제 채널을 3.0dB로 고정하고 design SNR을 -2.0~5.0dB(0.5dB 간격)로 훑은 결과, BER 곡선이 단조 증가가 아니라 **뚜렷한 3단계 패턴**으로 나타났습니다.
  1. **Design SNR -2.0~0.0dB (붕괴 구간)**: BER ≈ 0.5, FER ≈ 1.0으로 사실상 랜덤 판정 수준까지 무너짐.
  2. **Design SNR 0.0~0.5dB (급격한 전환 구간)**: BER이 수직에 가깝게 급락(개선).
  3. **Design SNR 0.5~4.5dB (평탄 구간)**: BER이 0.003~0.006 수준에서 design SNR을 넓게 바꿔도 거의 변하지 않음.
  4. **Design SNR ≈ 5.0dB (재상승)**: 실제 채널(3.0dB)보다 지나치게 낙관적인 마스크를 써서 경계 채널을 잘못 골라, 성능이 다시 소폭 나빠짐.

  이 패턴의 원인은 `construct_frozen_mask()`의 초기값 `mu[0] = 2.0/sigma2`에 있습니다. **design SNR이 너무 낮으면 `mu[0]`이 작아 N=1024, 10단계 재귀만으로는 서브채널들의 극화(polarization)가 충분히 일어나지 못해** 어떤 채널이 진짜 좋은지 거의 구분이 안 되고, 그 결과 정보 비트 위치가 사실상 무작위로 뽑혀 BER이 0.5에 근접합니다. Design SNR이 일정 임계값(이 조건에서는 대략 0dB)을 넘어서면 극화가 충분히 강하게 일어나면서 채널 우선순위가 실제 채널의 우선순위와 크게 다르지 않게 안정화되어, 넓은 범위에서 낮고 평탄한 BER 구간(plateau)이 만들어집니다.

  **즉 "design SNR이 실제 채널과 얼마나 정확히 일치하는가"보다 "design SNR이 극화를 일으키기에 충분히 높은가"가 훨씬 결정적인 임계 조건**입니다. 기존에 보고된 "38.85배 열화"는 정확히는 이 임계값보다 아래(0.0dB 부근)에서 설계된 마스크를 상대적으로 좋은 채널(3.0dB)에 강제로 쓴 극단적 케이스였음을, 스윕 데이터가 구체적으로 뒷받침합니다.



## 3. 프로젝트 구조

```text
Polar-Code/
├─ polar_code.c              # 기본 실행 파일: Baseline/A/B 3개 시나리오 고정 비교 + 그래프
├─ polar_code_simAB.c        # 심화 분석 실행 파일: --alpha, --sweep-a, --sweep-b, --image
├─ README.md
├─ image.png                 # 이미지 전송 테스트용 입력 (그레이스케일 변환되어 사용됨)
└─ result/                   # 모든 실행 결과가 이 폴더 아래에 자동 저장됨
   ├─ simulation_baseline.txt / simulation_scenario_a.txt / simulation_scenario_b.txt
   ├─ scenario_comparison_summary.txt
   ├─ lena_output_baseline.png / _preview.png (Scenario A/B도 동일 패턴)
   ├─ simulation_custom_alpha_*.txt, comparison_alpha_*.txt, lena_output_alpha_*.png       (--alpha)
   ├─ scenario_a_alpha_sweep_true*.txt, lena_output_alpha_sweep_*.png                      (--sweep-a)
   └─ scenario_b_design_snr_sweep_true*.txt, lena_output_design_snr_sweep_*.png            (--sweep-b)
```

## 4. 시뮬레이션 환경 구성

| 파라미터 | 설정 값 | 비고 |
| :--- | :--- | :--- |
| **Block Length ($N$)** | 1024 | Polar Code 블록 크기 |
| **Information Bits ($K$)** | 512 | 유효 정보 데이터 크기 |
| **Code Rate ($R$)** | 0.5 | 부호율 |
| **Modulation** | BPSK | $0 \to +1$, $1 \to -1$ 매핑 |
| **SNR Range ($E_b/N_0$)** | 0.0 dB ~ 3.0 dB | `polar_code.c` 기본 스캔 (0.25dB 간격) |
| **Scenario A Setting (기본)** | $\alpha = 0.5$ | 잡음 과소평가 상황 모델링 |
| **Scenario B Setting (기본)** | 0.0 dB Fixed Mask | 부호 설계 SNR 고정 상황 모델링 |
| **Image Test SNR (기본)** | 3.0 dB (고정) | 시나리오별 이미지 전송 테스트용 |
| **Alpha Sweep 범위** | 0.2 ~ 3.0 (0.4 간격, 8개) | `polar_code_simAB.c --sweep-a` |
| **Design SNR Sweep 범위** | -2.0 ~ 5.0 dB (0.5dB 간격) | `polar_code_simAB.c --sweep-b`, 실제 채널 3.0dB 고정 |



## 5. 실행 방법 및 결과 파일 명세

### 1) 기본 3-시나리오 비교: `polar_code.c`
Baseline / Scenario A(α=0.5 고정) / Scenario B(Design SNR=0.0dB 고정)를 한 번에 비교합니다. 대규모 몬테카를로 BER/FER 연산과, 동일한 시나리오 조건으로 진행되는 이미지(그레이스케일) 전송 테스트를 순서대로 수행한 뒤, Gnuplot으로 세 시나리오를 3분할 가로형 시각화 패널(`multiplot`)로 자동 표시합니다.

```bash
# 컴파일 (실행 파일과 같은 폴더에 image.png가 있어야 함)
gcc -std=c11 -O2 polar_code.c -lm -o polar_sim

# 실행
./polar_sim
```

인자 없이 실행하면 다음이 순서대로 진행됩니다:
1. `Baseline` / `Scenario A` / `Scenario B` BER/FER 몬테카를로 시뮬레이션 (0.0~3.0dB, 0.25dB 간격)
2. Baseline 대비 Scenario A/B의 SNR별 BER 배율 비교표 출력 및 저장
3. 동일한 세 조건(Baseline / A / B)으로 `image.png`를 3.0dB 채널에 전송하는 이미지 테스트 3회 (각각 원본/복원/차이 미리보기 이미지 생성)
4. Gnuplot으로 세 시나리오 BER/FER 그래프를 한 화면에 표시

특정 이미지를 수동으로 한 번만 테스트하고 싶다면:
```bash
./polar_sim --image input.png output.png [snr_db]
```

### 2) A/B 개별 심화 분석: `polar_code_simAB.c`
Scenario A(α)와 Scenario B(Design SNR)를 각각 원하는 단일 값 또는 범위로 자유롭게 테스트합니다.

```bash
# 컴파일
gcc -std=c11 -O2 polar_code_simAB.c -lm -o polar_simAB

# 단발성 alpha 테스트 (baseline 대비 SNR별 배율표 + 해당 alpha의 이미지 1장)
./polar_simAB --alpha 2.0

# Scenario A: alpha를 0.2~3.0(0.4 간격, 8개)으로 훑어 표 + 이미지 8장 + 그래프
./polar_simAB --sweep-a
./polar_simAB --sweep-a 2.0        # 실제 채널 SNR을 3.0dB가 아닌 다른 값으로 바꿔서 스윕

# Scenario B: design SNR을 -2.0~5.0dB(0.5dB 간격)으로 훑어 표 + 이미지 + 그래프
./polar_simAB --sweep-b
./polar_simAB --sweep-b 2.0        # 실제 채널 SNR을 다른 값으로 바꿔서 스윕

# 수동 단발 이미지 테스트 (시나리오 조건 없이 baseline LLR로만 전송)
./polar_simAB --image input.png output.png [snr_db]
```

각 모드는 `result/` 아래에 표(txt) + 이미지(png/_preview.png) + (스윕 모드는) gnuplot 그래프를 자동으로 생성합니다. `--sweep-a`/`--sweep-b`는 실행 중 진행 상황(`[3/8] alpha = 1.00 done...`)을 콘솔에 순차 출력한 뒤, 마지막에 baseline(또는 matched 지점) 대비 배율(`Ratio_to_alpha1`, `Ratio_to_matched`)이 포함된 최종 비교표를 보여줍니다.

### 3) 데이터 출력 파일 명세
- `result/simulation_baseline.txt` / `_scenario_a.txt` / `_scenario_b.txt`: `[Eb/No(dB) BER FER]` 형식의 시나리오별 성능 지표 (`polar_code.c`)
- `result/scenario_comparison_summary.txt`: `[Eb/No(dB) Base_BER ScenA_BER A_over_Base ScenB_BER B_over_Base]` 형식의 baseline 대비 배율 비교표 (`polar_code.c`)
- `result/comparison_alpha_*.txt`: 단발 α 테스트의 SNR별 baseline 대비 배율 (`polar_code_simAB.c --alpha`)
- `result/scenario_a_alpha_sweep_true*.txt`: α별 BER/FER 및 α=1 대비 배율 (`--sweep-a`)
- `result/scenario_b_design_snr_sweep_true*.txt`: design SNR별 BER/FER 및 matched 지점 대비 배율 (`--sweep-b`)
- `result/lena_output_*.png`, `_preview.png`: 각 모드별 이미지 전송 결과 및 원본/복원/차이 비교 미리보기



## 🧪 잡음 추정 오차 및 Design SNR 미스매치 분석 (Robustness Analysis)

실제 무선 수신기(Receiver)는 채널 상황을 100% 완벽하게 추정할 수 없으므로, 잡음 분산($\sigma^2$) 및 채널 환경 추정에 오차가 발생합니다. 본 시뮬레이션에서는 이러한 비이상적(Non-ideal) 수신기 환경이 Polar Code 복호 성능에 미치는 영향을 두 가지 시나리오로 나누어 분석했습니다. 아래는 `polar_code.c`의 고정 조건(α=0.5, Design SNR=0.0dB)으로 얻은 `result/scenario_comparison_summary.txt` 실측치입니다.

| Eb/No (dB) | Baseline BER | Scenario A BER | A/Base | Scenario B BER | B/Base |
| :---: | :---: | :---: | :---: | :---: | :---: |
| 0.00 | 4.879e-01 | 4.847e-01 | 0.99x | 4.880e-01 | 1.00x |
| 0.50 | 4.484e-01 | 4.245e-01 | 0.95x | 4.408e-01 | 0.98x |
| 1.00 | 3.762e-01 | 3.014e-01 | 0.80x | 4.293e-01 | 1.14x |
| 1.50 | 2.232e-01 | 1.363e-01 | 0.61x | 3.222e-01 | 1.44x |
| 2.00 | 9.575e-02 | 3.015e-02 | 0.31x | 2.325e-01 | 2.43x |
| 2.50 | 2.830e-02 | 4.973e-03 | 0.18x | 1.933e-01 | 6.83x |
| 3.00 | 3.000e-03 | 2.077e-04 | **0.07x** | 1.165e-01 | **38.85x** |

전체 데이터는 `result/scenario_comparison_summary.txt` 참고. α 및 Design SNR을 다른 값들로 훑은 결과는 `polar_code_simAB.c --sweep-a` / `--sweep-b` 실행으로 `result/scenario_a_alpha_sweep_true*.txt`, `result/scenario_b_design_snr_sweep_true*.txt`에서 확인할 수 있습니다.

### 1. Scenario A: LLR 미스매치 (LLR Mismatch / Noise Estimation Error)

수신기가 실제 채널의 잡음 분산($\sigma^2_{\text{true}}$)을 오해하여, 오차 계수 $\alpha$가 반영된 추정 분산($\sigma^2_{\text{est}} = \alpha \times \sigma^2_{\text{true}}$)으로 초기 LLR을 구하는 상황입니다.

- **발생 원인**: 수신기 내부의 채널 추정기(Channel Estimator) 오차
- **수학적 영향**:
  $$LLR_{\text{est}} = \frac{2}{\sigma^2_{\text{est}}} \cdot r = \frac{1}{\alpha} \left( \frac{2}{\sigma^2_{\text{true}}} \cdot r \right)$$
  - $\alpha < 1.0$ (과소평가): LLR 스케일이 균일하게 확대됨
  - $\alpha > 1.0$ (과대평가): LLR 스케일이 균일하게 축소됨
- **Exact SPA 복호기에서 실제로 관찰된 메커니즘**:
  - Exact SPA의 비선형 보정 항 $\ln(1+e^{-|l_1+l_2|}) - \ln(1+e^{-|l_1-l_2|})$은 입력 LLR의 절댓값이 커질수록 0에 수렴합니다.
  - $\alpha < 1$로 LLR이 확대되면 복호기는 점점 Min-Sum 근사와 동일하게 동작하게 되고, Min-Sum의 하드 판정은 모든 입력에 동일한 양(+)의 배율이 곱해져도 변하지 않는(scale-invariant) 성질을 가집니다. 그 결과 baseline 대비 BER이 오히려 개선됩니다($\alpha=0.5$ 기준 3.0dB에서 0.07x).
  - 반대로 $\alpha > 1$로 LLR이 축소되면 비선형 보정 항의 상대적 비중과 전체 신뢰도 저평가가 겹쳐 오류가 급격히 늘고 오류 전파가 심해집니다($\alpha=2.0$ 실측: 이미지 전면 붕괴).
  - 즉 Scenario A는 **α=1을 기준으로 완전히 비대칭적**입니다. `--sweep-a`로 α 전 구간을 훑으면 이 비대칭 곡선을 한 번에 확인할 수 있습니다.

---
### 2. Scenario B: GA 부호 구축 미스매치 (Design SNR Mismatch)

실제 운용되는 채널 SNR 대역이 변하더라도, 부호 구축(GA 마스크 생성) 시 특정 고정 SNR($\text{SNR}_{\text{design}}$) 기준의 프로즌 마스크(`info_mask`)를 재사용하는 상황입니다.

- **발생 원인**: SNR 변화에 따라 매번 프로즌 마스크를 동적으로 재계산하지 않고, 단일 계산된 마스크를 고정 사용하는 하드웨어 구조
- **열화 메커니즘 (실측 확인)**:
  - Polar Code의 채널 극화(Polarization) 양상은 $E_b/N_0$ 수준에 따라 변화합니다.
  - `polar_code.c` 고정 조건(Design SNR=0.0dB, 실제 채널 0~3.0dB)에서는 실제 SNR과 설계 SNR의 차이가 커질수록 baseline 대비 BER 배율이 1.00x(0dB) → 38.85x(3.0dB)로 지속적으로 증가했습니다.
  - `polar_code_simAB.c --sweep-b`로 실제 채널을 3.0dB에 고정하고 design SNR을 더 넓게(-2.0~5.0dB) 훑어보면, 이 열화는 단조 증가가 아니라 **극화 임계값을 기준으로 한 계단형 패턴**임이 드러납니다: design SNR이 임계값(이 조건에서 약 0dB) 아래면 극화가 거의 일어나지 않아 BER이 0.5(랜덤 수준)까지 붕괴하고, 임계값을 넘으면 넓은 범위(0.5~4.5dB)에서 낮고 평탄한 BER을 유지하며, 실제 채널보다 지나치게 낙관적인 design SNR(≈5dB)에서는 다시 소폭 열화됩니다.
  - 즉 "design SNR이 실제 채널과 정확히 일치하는가"보다 **"design SNR이 극화를 일으키기에 충분히 높은가"가 성능을 가르는 핵심 임계 조건**입니다.

---

### 📊 시나리오 요약 및 비교

| 항목 | Scenario A (LLR Mismatch) | Scenario B (Design SNR Mismatch) |
| :--- | :--- | :--- |
| **변수 위치** | 수신단 LLR 스케일링 계산부 | 송/수신단 GA 프로즌 마스크 생성부 |
| **핵심 원인** | 잡음 분산 추정 오차 ($\sigma^2_{\text{est}} \neq \sigma^2_{\text{true}}$) | 부호 설계 SNR과 실제 채널 SNR의 불일치 |
| **실측 경향** | α=1 기준 비대칭: α<1은 baseline 대비 **개선**(0.99x→0.07x), α>1은 **급격히 악화** | design SNR이 극화 임계값 이상이면 넓게 **평탄**, 임계값 미만이면 **랜덤 수준까지 붕괴**, 실제보다 지나치게 낙관적이면 소폭 재열화 |
| **실험 도구** | `polar_code.c`(α=0.5 고정) + `polar_code_simAB.c --alpha`/`--sweep-a`(α 자유 테스트) | `polar_code.c`(0.0dB 고정) + `polar_code_simAB.c --sweep-b`(design SNR 전 구간 스윕) |

> 💡 **Key Takeaway**
> Scenario A는 α=1을 기준으로 완전히 비대칭적입니다 — LLR을 과대평가(α<1, 확대)하는 방향은 Exact SPA가 Min-Sum에 가까워지며 오히려 성능이 개선되지만, 과소평가(α>1, 축소)하는 방향은 심각하게 악화됩니다. Scenario B는 "design SNR이 실제 채널과 얼마나 정확히 일치하는가"보다 "design SNR이 채널 극화를 일으키기에 충분히 높은가"라는 임계 조건이 훨씬 결정적이며, 이 임계값 아래로 떨어지면 성능이 완만하게가 아니라 급격한 계단형으로 붕괴합니다. 두 시나리오 모두 "오차의 크기"만이 아니라 "오차의 방향"과 "임계 조건 대비 위치"가 실제 성능에 훨씬 중요하다는 것을 보여줍니다.