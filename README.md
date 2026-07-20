# Polar Code (N=1024, K=512) Simulation Platform
본 플랫폼은 가우시안 근사(Gaussian Approximation, GA) 기반의 부호 구축과 최적화된 Exact SPA 연속 제거(Successive Cancellation, SC) 복호화 알고리즘을 64비트 고속 난수 엔진 기반의 C언어 환경으로 구현한 고정밀 통신 시뮬레이터입니다.



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
이를 뷘커 박스-뮬러(Box-Muller) 변환과 결합하여, 높은 SNR 대역에서도 통계적 무결성과 독립성을 유지하는 정밀한 백색 가우시안 잡음(AWGN) 환경을 제공합니다.



## 2. 연구용 비정상 시나리오 해석 기능

본 시뮬레이터는 완벽한 채널 환경(Baseline) 외에도, 실제 수신기 칩셋 구현 시 마주할 수 있는 하드웨어 제약 및 오차 상황을 정밀 모델링하기 위해 2가지 학술적 분석 시나리오를 지원합니다.

### 1) 시나리오 A: 복호기 입력 LLR의 스케일링 오차 (LLR Scaling Error)
실제 채널에 흐르는 진짜 잡음 분산($\sigma^2_{true}$)과 수신기가 추정 수식에서 믿고 있는 잡음 분산($\sigma^2_{est}$) 사이에 의도적인 오차 계수($\alpha = 0.5$)를 개입시켜 초기 LLR 스케일을 왜곡하는 상황입니다.

$$\sigma^2_{est} = \alpha \times \sigma^2_{true}$$

- **학술적 현상**: 입력 LLR 전체에 상수가 곱해지더라도 부호 판정 결과가 일정한 Min-Sum 방식과 달리, 본 플랫폼의 **Exact SPA 복호기**는 비선형 Jacobian 보정 항을 사용하므로 스케일 왜곡에 매우 민감합니다. 연산 보정 값이 꼬이면서 특정 환경(2.5 dB ~ 2.7 dB)에서 에러율이 일시적으로 역전되는 **요철(Bump) 현상**이 선명하게 관찰됩니다.

### 2) 시나리오 B: 부호 구축 시의 Design SNR 미스매치 (Design SNR Mismatch)
채널 변동에 맞춰 프로즌 마스크를 동적으로 변경하지 못하고, 특정 설계 환경(Design SNR = 0.0 dB)에서 고정 생성된 단 하나의 마스크를 사용하여 채널 환경이 좋은 대역까지 강제 복호하는 하드웨어 최적화 제약 상황입니다.

- **학술적 현상**: Polar Code가 완전한 채널 의존적 부호(Channel-Dependent Code)임을 증명하는 지표입니다. 0.0 dB 환경에 동결된 마스크는 신호 세기가 우수한 3.0 dB 대역에 도달하더라도 채널 분극(Polarization) 효율을 전혀 살리지 못해, FER이 그래프 상단(0.3 ~ 1.0)에 그대로 갇혀 정체되는 치명적인 성능 열화(Floor)가 발생합니다.



## 3. 프로젝트 구조

```text
Polar-Code/
├─ include/
│  ├─ polar_config.h
│  ├─ rng.h
│  ├─ polar_math.h
│  ├─ polar_codec.h
│  ├─ simulation.h
│  └─ plot.h
├─ src/
│  ├─ main.c
│  ├─ rng.c
│  ├─ polar_math.c
│  ├─ polar_codec.c
│  ├─ simulation.c
│  └─ plot.c
├─ README.md
├─ simulation_baseline.txt
├─ simulation_scenario_a.txt
└─ simulation_scenario_b.txt
```

## 4. 시뮬레이션 환경 구성

| 파라미터 | 설정 값 | 비고 |
| :--- | :--- | :--- |
| **Block Length ($N$)** | 1024 | Polar Code 블록 크기 |
| **Information Bits ($K$)** | 512 | 유효 정보 데이터 크기 |
| **Code Rate ($R$)** | 0.5 | 부호율 |
| **Modulation** | BPSK | $0 \to +1$, $1 \to -1$ 매핑 |
| **SNR Range ($E_b/N_0$)** | 0.0 dB ~ 3.0 dB | 0.25 dB 간격 고정 스캔 |
| **Scenario A Setting** | $\alpha = 0.5$ | 잡음 과소평가 상황 모델링 |
| **Scenario B Setting** | 0.0 dB Fixed Mask | 부호 설계 SNR 고정 상황 모델링 |



## 5. 실행 방법 및 결과 파일 명세

### 1) 실행 명령
본 프로그램은 대규모 몬테카를로 연산 완료 후, 연동된 Gnuplot 파이프를 통해 세 시나리오의 Y축 스케일을 일치시킨 3분할 독립 가로형 시각화 패널(`multiplot`)을 화면에 자동으로 팝업합니다.

```
# 1. 컴파일
gcc -std=c11 -O2 polar_code.c -lm -o polar_legacy.exe

# 2. 프로그램 실행
./polar_legacy
```

```bash
# 1. 컴파일
gcc -std=c11 -O2 -Iinclude src\main.c src\rng.c src\polar_math.c src\polar_codec.c src\simulation.c src\plot.c -lm -o polar_code.exe

# 2. 프로그램 실행
./polar_code
```

### 2) 데이터 출력 파일 명세
실행이 완료되면 소스코드와 동일한 디렉터리 내에 시나리오별로 분리된 3개의 텍스트 파일이 자동 생성됩니다. 각 파일은 [Eb/No(dB)  BER  FER] 구조로 이루어져 있습니다.
- `simulation_baseline.txt`: 이상적인 환경에서의 성능 지표 데이터
- `simulation_scenario_a.txt`: LLR 스케일링 오차 반영 데이터 ($\alpha = 0.5$)
- `simulation_scenario_b.txt`: 부호 설계 미스매치 반영 데이터 (0.0 dB 고정)
