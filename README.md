# Polar Code Simulation Platform (N=1024, K=512)

본 플랫폼은 Gaussian Approximation 기반의 부호 구축과 Exact SPA 기반 Successive Cancellation 복호 알고리즘을 64비트 난수 엔진을 활용한 C언어 환경에서 구현한 고정밀 통신 시뮬레이터이다.

이 프로젝트는 목적에 따라 세 가지 실행 파일로 구성하는 것을 목표로 한다.

- **`polar_code.c`**: Baseline, Scenario A, Scenario B를 같은 조건에서 동시에 비교하는 기본 실행 파일이다. BER/FER 그래프, 시나리오별 이미지 전송 실험, 비교 요약표를 자동으로 생성한다.
- **`polar_code_simAB.c`**: Scenario A(α)와 Scenario B(Design SNR)를 원하는 값이나 범위로 탐색하는 심화 분석용 실행 파일이다. 단일 테스트(`--alpha`)와 다중 스윕(`--sweep-a`, `--sweep-b`)을 지원한다.
- **`polar_code_SCD_sim.c`**: 단순 SC/Min-Sum 계열 기준 성능을 검증하는 SCD 시뮬레이터다. SCL 시뮬레이터는 같은 구조의 별도 폴더로 추가할 수 있도록 예약해 두는 것이 좋다.

## 1. 핵심 이론 및 알고리즘 명세

### 1) 가우시안 근사 (Gaussian Approximation) 기반 부호 구축

채널 분극화(Channel Polarization)에 의해 분할되는 서브채널의 신뢰도를 정밀하게 추적하기 위해, 밀도 진화(Density Evolution)의 가우시안 근사 함수 $\phi(x)$와 이진 탐색 기반 역함수 $\phi^{-1}(y)$ 알고리즘을 적용한다.

$$\phi(x) = \begin{cases} 1.0, & \text{if } x \le 0 \\ exp(-0.4527x^{0.86} + 0.0218), & \text{if } 0 < x \le 10 \\ \sqrt{\frac{\pi}{x}}\left(1 - \frac{10}{7x}\right)exp\left(-\frac{x}{4}\right), & \text{if } x > 10 \end{cases}$$

- **Bit-Reversal Index 매핑**: 하향식(Bottom-Up) 나비 구조로 계산된 신뢰도 벡터 mu를 물리 트리 구조와 맞추기 위해 비트 역전 변환을 거쳐 Natural Order로 재배치한다.
- **Frozen Mask 확정**: 정렬된 신뢰도 지수를 바탕으로 상위 K=512개의 채널을 정보 비트(Information) 영역으로 지정하고, 나머지 채널은 수렴 능력이 낮은 가상 채널로 간주하여 0(Frozen)으로 잠근다.

### 2) Top-Down 재귀형 부호화기 (Natural Order Encoder)

생성 행렬 $G_N$의 직접적인 대규모 행렬 곱 연산을 피하고, 복호 트리 구조와 대칭성을 유지하는 상향식 재귀 인코딩을 구현한다. 상위 블록 입력 u에 대해 좌측 하부 브랜치는 $v = u_{left} \oplus u_{right}$ 커널 연산을 거치고, 최종 부호화된 코드워드 x는 최하단 레이어에서 단순 연결(Plain Concat)되어 송신 심볼로 출력된다.

### 3) Exact SPA 기반 연속 제거(SC) 복호화기

하드웨어 구현 편의성을 위해 Min-Sum 근사 방식의 성능 열화를 줄이기 위해, 수학적 연속성을 보장하는 Jacobian Logarithm 기반 보정 항이 결합된 **Exact SPA(Sum-Product Algorithm)** 연산 레이어를 구현한다.

- **f-function (Left Path)**:
  $$LLR_{left} = \text{sign}(l_1) \cdot \text{sign}(l_2) \cdot \min(\vert{}l_1\vert{}, \vert{}l_2\vert{}) + \ln(1 + e^{-\vert{}l_1+l_2\vert{}}) - \ln(1 + e^{-\vert{}l_1-l_2\vert{}})$$
- **g-function (Right Path)**: 이전 좌측 트리에서 복호 완료되어 피드백된 부분 코드 워드($u_{coded\_left}$)를 반영하여 우측 채널 우도를 갱신합니다.
  $$LLR_{right} = l_2 + (1 - 2 \cdot u_{coded\_left}) \cdot l_1$$

### 4) Xorshift64 난수 엔진 및 AWGN 채널

기존 C 표준 라이브러리의 `rand()` 함수는 주기적 한계와 의사난수 뭉침 문제를 가지므로, 이를 보완하기 위해 **Xorshift64 난수 생성기**를 적용한다.
이 난수 생성기는 Box-Muller 변환과 결합되어, 높은 SNR 영역에서도 통계적 무결성과 독립성을 유지하는 정밀한 백색 가우시안 잡음(AWGN) 환경을 제공한다.

## 2. 연구용 비정상 시나리오 해석 기능

본 시뮬레이터는 완벽한 채널 환경(Baseline) 외에도, 실제 수신기 칩셋 구현에서 발생할 수 있는 하드웨어 제약과 오차 상황을 정밀하게 모델링하기 위해 두 가지 학술적 분석 시나리오를 지원한다.

### 1) 시나리오 A: 복호기 입력 LLR의 스케일링 오차 (LLR Scaling Error)

실제 채널의 잡음 분산($\sigma^2_{true}$)과 수신기가 추정하는 잡음 분산($\sigma^2_{est}$) 사이에 의도적인 오차 계수 $\alpha$를 삽입하여 초기 LLR 스케일을 왜곡하는 상황이다.

$$\sigma^2_{est} = \alpha \times \sigma^2_{true}$$

- **실측 결과 ($\alpha = 0.5$, `polar_code.c` 고정 조건)**: LLR 전체 크기를 2배로 확대하는 오차 방향에서는 Exact SPA의 비선형 보정 항이 입력 크기에 따라 점점 0에 수렴하여 복호기가 **Min-Sum 근사**와 유사하게 동작한다. Min-Sum 결합은 부호(sign) 판정이 동일한 양(+)의 비례 스케일링에 불변하므로, 이 방향으로의 오차는 오히려 **BER을 개선**시키는 경향을 보인다. `scenario_comparison_summary.txt`를 기준으로 SNR이 높아질수록 baseline 대비 오차율은 계속 감소하며(0dB: 0.99배 → 3.0dB: 0.069배), 우려되었던 "2.5~2.7dB 근처의 bump"는 관측되지 않았다.
- **α>1(축소) 방향 검증 (`polar_code_simAB.c --alpha 2.0`)**: α=2.0인 경우 LLR이 절반으로 축소되면서 이미지가 전반적으로 심하게 깨지는 현상이 나타난다. LLR이 작아질수록 비선형 보정 항의 상대적 영향이 커지고 전체 신뢰도 자체가 낮게 추정되므로, SC 복호기의 오류 전파가 심해지기 때문이다. 즉 **α<1의 과대평가 방향과 α>1의 과소평가 방향은 비대칭적**이다.
- **α 전 구간 스윕 (`polar_code_simAB.c --sweep-a`)**: α = 0.2~3.0 범위를 0.4 간격으로 훑으면 BER/FER/이미지 결과를 한 번에 생성할 수 있으며, α=1을 기준으로 비대칭성을 직관적으로 비교할 수 있다.

### 2) 시나리오 B: 부호 구축 시의 Design SNR 미스매치 (Design SNR Mismatch)

채널 변동에 대응해 프로즌 마스크를 매번 재계산하지 못하고, 특정 설계 환경(Design SNR)에서 고정된 마스크를 사용하여 다른 채널 조건에서 복호를 수행하는 상황을 의미한다.

- **실측 결과 (Design SNR = 0.0dB 고정, `polar_code.c` 고정 조건)**: Polar Code가 완전한 채널 의존형 부호임을 보여주는 결과로, 설계된 마스크를 실제 채널과 다르게 사용하면 성능 저하가 급격히 발생한다. `scenario_comparison_summary.txt` 기준으로 baseline 대비 오차율은 0dB에서 1.00배 수준이었으나, 3.0dB에서는 38.85배까지 증가한다.
- **Design SNR 전 구간 스윕 (`polar_code_simAB.c --sweep-b`, 실제 채널 3.0dB 고정)**: 실제 채널을 3.0dB로 고정하고 design SNR을 -2.0~5.0dB 범위로 훑으면 BER 곡선이 단조 증가가 아니라 **뚜렷한 3단계 패턴**을 보인다.
  1. **Design SNR -2.0~0.0dB (붕괴 구간)**: BER ≈ 0.5, FER ≈ 1.0으로 사실상 랜덤 판정 수준까지 무너진다.
  2. **Design SNR 0.0~0.5dB (급격한 전환 구간)**: BER이 급격하게 하락하며 성능이 개선된다.
  3. **Design SNR 0.5~4.5dB (평탄 구간)**: BER이 약 0.003~0.006 수준으로 유지되며, design SNR 변화에 크게 민감하지 않다.
  4. **Design SNR ≈ 5.0dB (재상승)**: 실제 채널보다 지나치게 낙관적인 마스크를 사용해 경계 채널을 잘못 선택하고, 성능이 다시 약간 악화된다.

  이 패턴의 원인은 `construct_frozen_mask()`의 초기값 `mu[0] = 2.0/sigma2`에 있다. **Design SNR이 너무 낮으면 `mu[0]`이 작아 N=1024, 10단계 재귀만으로는 서브채널의 극화가 충분히 일어나지 못하고**, 어떤 채널이 우수한지 거의 구분되지 않는다. 그 결과 정보 비트 위치가 사실상 무작위로 선택되어 BER이 0.5에 근접한다. 반대로 Design SNR이 일정 임계값(이 조건에서는 대략 0dB)을 넘으면 극화가 충분히 강해지고, 채널 우선순위가 실제 채널의 우선순위와 크게 다르지 않게 안정화되며, 넓은 범위에서 낮고 평탄한 BER 구간이 형성된다.

  **즉 "design SNR이 실제 채널과 얼마나 정확히 일치하는가"보다 "design SNR이 극화를 일으키기에 충분히 높은가"가 훨씬 더 중요한 임계 조건**이다. 기존에 보고된 "38.85배 열화"는 이 임계값보다 낮은 영역에서 설계된 마스크를 상대적으로 좋은 채널(3.0dB)에 강제로 적용한 극단적인 경우로 볼 수 있으며, 스윕 데이터가 이를 뒷받침한다.

## 3. 프로젝트 구조

```text
Polar-Code/
├─ core/
│  ├─ polar_code.c
│  ├─ polar_sim.exe
│  ├─ image.png
│  └─ result/
├─ scd/
│  ├─ polar_code_SCD_sim.c
│  ├─ polar_code_SCD_sim.exe
│  └─ result/
├─ simab/
│  ├─ polar_code_simAB.c
│  ├─ polar_code_simAB.exe
│  ├─ image.png
│  └─ result/
├─ scl/
│  └─ result/
├─ README.md
└─ image.png
```

현재 코드는 각 실행 파일이 자기 전용 폴더 안의 `result/`를 자동 생성하도록 맞춰 두었다. 즉 `core/polar_code.c`는 `core/result/`, `scd/polar_code_SCD_sim.c`는 `scd/result/`, `simab/polar_code_simAB.c`는 `simab/result/` 아래에 파일을 저장한다.

## 4. 프로젝트 요약표

### 4.1 파일별 역할

| 파일                       | 역할                  | 주로 측정하는 것                            |
| -------------------------- | --------------------- | ------------------------------------------- |
| `core/polar_code.c`        | 기본 비교용 메인 코드 | Baseline, Scenario A, Scenario B 비교       |
| `simab/polar_code_simAB.c` | 심화 분석용 코드      | $\alpha$ 스윕, Design SNR 스윕, 커스텀 실험 |
| `scd/polar_code_SCD_sim.c` | 참조용 기본 코드      | 단순 SC 기반 성능 검증                      |
| `scl/polar_code_SCL_sim.c` | 예정                  | CRC/SCL 기반 확장 시나리오                  |

### 4.2 실행 파일과 소스 매핑

| 실행 파일                | 원본 코드                  | 권장 위치 |
| ------------------------ | -------------------------- | --------- |
| `polar_sim.exe`          | `core/polar_code.c`        | `core/`   |
| `polar_simAB.exe`        | `simab/polar_code_simAB.c` | `simab/`  |
| `polar_code_SCD_sim.exe` | `scd/polar_code_SCD_sim.c` | `scd/`    |
| `polar_code_SCL_sim.exe` | `scl/polar_code_SCL_sim.c` | `scl/`    |

> `.exe` 파일은 소스 파일을 컴파일한 결과일 뿐이며, 실제 성능 차이는 소스 코드의 시나리오 설정과 조건에 따라 결정된다.

### 4.3 시나리오 요약

| 시나리오   | 핵심 변수  | 설명                         | 측정 항목             |
| ---------- | ---------- | ---------------------------- | --------------------- |
| Baseline   | 없음       | 이상적 채널 조건의 기준 성능 | BER, FER              |
| Scenario A | $\alpha$   | LLR 스케일링 오차            | BER, FER, 이미지 복원 |
| Scenario B | Design SNR | frozen mask 설계 불일치      | BER, FER, 민감도      |

### 4.4 결과 파일 요약

| 결과 파일                                                         | 설명                    |
| ----------------------------------------------------------------- | ----------------------- |
| `core/result/simulation_baseline.txt`                             | Baseline BER/FER 결과   |
| `core/result/simulation_scenario_a.txt`                           | Scenario A BER/FER 결과 |
| `core/result/simulation_scenario_b.txt`                           | Scenario B BER/FER 결과 |
| `core/result/scenario_comparison_summary.txt`                     | Baseline 대비 비율 비교 |
| `simab/result/*alpha*.txt`                                        | $\alpha$ 실험 결과      |
| `simab/result/*sweep*.txt`                                        | 스윕 실험 결과          |
| `core/result/lena_output_*.png`, `simab/result/lena_output_*.png` | 이미지 복원 결과        |
| `scd/result/*`                                                    | SCD 기준 검증 결과      |

## 5. 시뮬레이션 환경 구성

| 파라미터                      | 설정 값                    | 비고                                                 |
| :---------------------------- | :------------------------- | :--------------------------------------------------- |
| **Block Length ($N$)**        | 1024                       | Polar Code 블록 크기                                 |
| **Information Bits ($K$)**    | 512                        | 유효 정보 데이터 크기                                |
| **Code Rate ($R$)**           | 0.5                        | 부호율                                               |
| **Modulation**                | BPSK                       | $0 \to +1$, $1 \to -1$ 매핑                          |
| **SNR Range ($E_b/N_0$)**     | 0.0 dB ~ 3.0 dB            | `polar_code.c` 기본 스캔 (0.25dB 간격)               |
| **Scenario A Setting (기본)** | $\alpha = 0.5$             | 잡음 과소평가 상황 모델링                            |
| **Scenario B Setting (기본)** | 0.0 dB Fixed Mask          | 부호 설계 SNR 고정 상황 모델링                       |
| **Image Test SNR (기본)**     | 3.0 dB (고정)              | 시나리오별 이미지 전송 테스트용                      |
| **Alpha Sweep 범위**          | 0.2 ~ 3.0 (0.4 간격, 8개)  | `polar_code_simAB.c --sweep-a`                       |
| **Design SNR Sweep 범위**     | -2.0 ~ 5.0 dB (0.5dB 간격) | `polar_code_simAB.c --sweep-b`, 실제 채널 3.0dB 고정 |

## 6. 실행 방법 및 결과 파일 명세

### 1) 기본 3-시나리오 비교: `polar_code.c`

Baseline / Scenario A(α=0.5 고정) / Scenario B(Design SNR=0.0dB 고정)를 한 번에 비교합니다. 대규모 몬테카를로 BER/FER 연산과, 동일한 시나리오 조건으로 진행되는 이미지(그레이스케일) 전송 테스트를 순서대로 수행한 뒤, Gnuplot으로 세 시나리오를 3분할 가로형 시각화 패널(`multiplot`)로 자동 표시합니다.

Windows에서는 `C:\msys64\ucrt64\bin`과 `C:\msys64\usr\bin`이 PATH에 들어가 있어야 `gcc`가 내부 헬퍼와 DLL을 찾을 수 있습니다. 가장 안전한 방법은 MSYS2 UCRT64 쉘에서 빌드하는 것이고, PowerShell에서 직접 실행할 경우에는 먼저 PATH를 추가한 뒤 컴파일하세요.

```bash
# 컴파일 (core 폴더 안에서 실행한다고 가정)
cd core
gcc -std=c11 -O2 polar_code.c -lm -o polar_sim.exe

# 실행
./polar_sim.exe
```

인자 없이 실행하면 다음이 순서대로 진행됩니다:

1. `Baseline` / `Scenario A` / `Scenario B` BER/FER 몬테카를로 시뮬레이션 (0.0~3.0dB, 0.25dB 간격)
2. Baseline 대비 Scenario A/B의 SNR별 BER 배율 비교표 출력 및 저장
3. 동일한 세 조건(Baseline / A / B)으로 `image.png`를 3.0dB 채널에 전송하는 이미지 테스트 3회 (각각 원본/복원/차이 미리보기 이미지 생성)
4. Gnuplot으로 세 시나리오 BER/FER 그래프를 한 화면에 표시

특정 이미지를 수동으로 한 번만 테스트하고 싶다면:

```bash
./polar_sim.exe --image input.png output.png [snr_db]
```

### 2) A/B 개별 심화 분석: `polar_code_simAB.c`

Scenario A(α)와 Scenario B(Design SNR)를 각각 원하는 단일 값 또는 범위로 자유롭게 테스트합니다.

위와 동일하게 Windows에서는 MSYS2 UCRT64 환경 또는 PATH 설정이 필요합니다.

```bash
# 컴파일 (simab 폴더 안에서 실행한다고 가정)
cd simab
gcc -std=c11 -O2 polar_code_simAB.c -lm -o polar_simAB.exe

# 단발성 alpha 테스트 (baseline 대비 SNR별 배율표 + 해당 alpha의 이미지 1장)
./polar_simAB.exe --alpha 2.0

# Scenario A: alpha를 0.2~3.0(0.4 간격, 8개)으로 훑어 표 + 이미지 8장 + 그래프
./polar_simAB.exe --sweep-a
./polar_simAB.exe --sweep-a 2.0        # 실제 채널 SNR을 3.0dB가 아닌 다른 값으로 바꿔서 스윕

# Scenario B: design SNR을 -2.0~5.0dB(0.5dB 간격)으로 훑어 표 + 이미지 + 그래프
./polar_simAB.exe --sweep-b
./polar_simAB.exe --sweep-b 2.0        # 실제 채널 SNR을 다른 값으로 바꿔서 스윕

# 수동 단발 이미지 테스트 (시나리오 조건 없이 baseline LLR로만 전송)
./polar_simAB.exe --image input.png output.png [snr_db]
```

각 모드는 `result/` 아래에 표(txt) + 이미지(png/\_preview.png) + (스윕 모드는) gnuplot 그래프를 자동으로 생성합니다. `--sweep-a`/`--sweep-b`는 실행 중 진행 상황(`[3/8] alpha = 1.00 done...`)을 콘솔에 순차 출력한 뒤, 마지막에 baseline(또는 matched 지점) 대비 배율(`Ratio_to_alpha1`, `Ratio_to_matched`)이 포함된 최종 비교표를 보여줍니다.

### 3) 데이터 출력 파일 명세

- `core/result/simulation_baseline.txt` / `_scenario_a.txt` / `_scenario_b.txt`: `[Eb/No(dB) BER FER]` 형식의 시나리오별 성능 지표 (`core/polar_code.c`)
- `core/result/scenario_comparison_summary.txt`: `[Eb/No(dB) Base_BER ScenA_BER A_over_Base ScenB_BER B_over_Base]` 형식의 baseline 대비 배율 비교표 (`core/polar_code.c`)
- `simab/result/comparison_alpha_*.txt`: 단발 α 테스트의 SNR별 baseline 대비 배율 (`simab/polar_code_simAB.c --alpha`)
- `simab/result/scenario_a_alpha_sweep_true*.txt`: α별 BER/FER 및 α=1 대비 배율 (`--sweep-a`)
- `simab/result/scenario_b_design_snr_sweep_true*.txt`: design SNR별 BER/FER 및 matched 지점 대비 배율 (`--sweep-b`)
- `core/result/lena_output_*.png`, `simab/result/lena_output_*.png`: 각 모드별 이미지 전송 결과 및 원본/복원/차이 비교 미리보기

## 🧪 잡음 추정 오차 및 Design SNR 미스매치 분석 (Robustness Analysis)

실제 무선 수신기(Receiver)는 채널 상황을 100% 완벽하게 추정할 수 없으므로, 잡음 분산($\sigma^2$) 및 채널 환경 추정에 오차가 발생합니다. 본 시뮬레이션에서는 이러한 비이상적(Non-ideal) 수신기 환경이 Polar Code 복호 성능에 미치는 영향을 두 가지 시나리오로 나누어 분석했습니다. 아래는 `polar_code.c`의 고정 조건(α=0.5, Design SNR=0.0dB)으로 얻은 `result/scenario_comparison_summary.txt` 실측치입니다.

| Eb/No (dB) | Baseline BER | Scenario A BER |  A/Base   | Scenario B BER |   B/Base   |
| :--------: | :----------: | :------------: | :-------: | :------------: | :--------: |
|    0.00    |  4.879e-01   |   4.847e-01    |   0.99x   |   4.880e-01    |   1.00x    |
|    0.50    |  4.484e-01   |   4.245e-01    |   0.95x   |   4.408e-01    |   0.98x    |
|    1.00    |  3.762e-01   |   3.014e-01    |   0.80x   |   4.293e-01    |   1.14x    |
|    1.50    |  2.232e-01   |   1.363e-01    |   0.61x   |   3.222e-01    |   1.44x    |
|    2.00    |  9.575e-02   |   3.015e-02    |   0.31x   |   2.325e-01    |   2.43x    |
|    2.50    |  2.830e-02   |   4.973e-03    |   0.18x   |   1.933e-01    |   6.83x    |
|    3.00    |  3.000e-03   |   2.077e-04    | **0.07x** |   1.165e-01    | **38.85x** |

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

| 항목          | Scenario A (LLR Mismatch)                                                            | Scenario B (Design SNR Mismatch)                                                                                                  |
| :------------ | :----------------------------------------------------------------------------------- | :-------------------------------------------------------------------------------------------------------------------------------- |
| **변수 위치** | 수신단 LLR 스케일링 계산부                                                           | 송/수신단 GA 프로즌 마스크 생성부                                                                                                 |
| **핵심 원인** | 잡음 분산 추정 오차 ($\sigma^2_{\text{est}} \neq \sigma^2_{\text{true}}$)            | 부호 설계 SNR과 실제 채널 SNR의 불일치                                                                                            |
| **실측 경향** | α=1 기준 비대칭: α<1은 baseline 대비 **개선**(0.99x→0.07x), α>1은 **급격히 악화**    | design SNR이 극화 임계값 이상이면 넓게 **평탄**, 임계값 미만이면 **랜덤 수준까지 붕괴**, 실제보다 지나치게 낙관적이면 소폭 재열화 |
| **실험 도구** | `polar_code.c`(α=0.5 고정) + `polar_code_simAB.c --alpha`/`--sweep-a`(α 자유 테스트) | `polar_code.c`(0.0dB 고정) + `polar_code_simAB.c --sweep-b`(design SNR 전 구간 스윕)                                              |

> 💡 **Key Takeaway**
> Scenario A는 α=1을 기준으로 완전히 비대칭적입니다 — LLR을 과대평가(α<1, 확대)하는 방향은 Exact SPA가 Min-Sum에 가까워지며 오히려 성능이 개선되지만, 과소평가(α>1, 축소)하는 방향은 심각하게 악화됩니다. Scenario B는 "design SNR이 실제 채널과 얼마나 정확히 일치하는가"보다 "design SNR이 채널 극화를 일으키기에 충분히 높은가"라는 임계 조건이 훨씬 결정적이며, 이 임계값 아래로 떨어지면 성능이 완만하게가 아니라 급격한 계단형으로 붕괴합니다. 두 시나리오 모두 "오차의 크기"만이 아니라 "오차의 방향"과 "임계 조건 대비 위치"가 실제 성능에 훨씬 중요하다는 것을 보여줍니다.
