# Polar Code (N=1024, K=512) Simulation Platform
본 플랫폼은 가우시안 근사(Gaussian Approximation, GA) 기반의 부호 구축과 최적화된 Exact SPA 연속 제거(Successive Cancellation, SC) 복호화 알고리즘을 64비트 고속 난수 엔진 기반의 C언어 환경으로 구현한 고정밀 통신 시뮬레이터입니다.

---

## 1. 핵심 이론 및 알고리즘 명세

### 1) 가우시안 근사 (Gaussian Approximation) 기반 부호 구축
채널 분극화(Channel Polarization) 현상에 의해 분할되는 서브 채널들의 신뢰도를 정밀하게 추적하기 위해 밀도 진화(Density Evolution)의 가우시안 근사 함수 $\phi(x)$와 이진 탐색 기반의 역함수 $\phi^{-1}(y)$ 알고리즘을 수행합니다.

$$\phi(x) = \begin{cases} 1.0, & \text{if } x \le 0 \\ exp(-0.4527x^{0.86} + 0.0218), & \text{if } 0 < x \le 10 \\ \sqrt{\frac{\pi}{x}}\left(1 - \frac{10}{7x}\right)exp\left(-\frac{x}{4}\right), & \text{if } x > 10 \end{cases}$$

- **Bit-Reversal Index 매핑**: 하향식(Bottom-Up) 나비 구조로 계산된 신뢰도 벡터 $mu$를 물리 트리 구조와 일치시키기 위해 비트 역전 변환을 거쳐 Natural Order로 재배치합니다.
- **Frozen Mask 확정**: 정렬된 신뢰도 지수를 바탕으로 상위 $K=512$개의 채널을 정보 비트(Information) 영역으로 지정하고, 나머지 채널은 수렴 능력이 떨어지는 가상 채널로 판단하여 0(Frozen)으로 잠금 처리합니다.

### 2) Top-Down 재귀형 부호화기 (Natural Order Encoder)
생성 행렬 $G_N$의 직접적인 대규모 행렬 곱 연산을 배제하고, 복호화 트리 구조와 완벽한 대칭성을 가지는 상향식 재귀 인코딩을 구현합니다. 상위 블록 입력 $u$에 대해 좌측 하부 브랜치는 $v = u_{left} \oplus u_{right}$ 커널 연산을 거치며, 최종 부호화된 코드 워드 $x$는 최하단 레이어에서 단순 연결(Plain Concat)되어 송신 심볼로 출력됩니다.

### 3) Exact SPA 기반 연속 제거(SC) 복호화기
하드웨어 구현 편의성을 위한 Min-Sum 근사 방식의 성능 열화를 극복하기 위해, 수학적 연속성을 완벽히 보장하는 야코비안 로그(Jacobian Logarithm) 정밀 보정 항이 결합된 **Exact SPA(Sum-Product Algorithm)** 연산 레이어를 구현했습니다.

- **f-function (Left Path)**:
  $$LLR_{left} = \text{sign}(l_1) \cdot \text{sign}(l_2) \cdot \min(\vert{}l_1\vert{}, \vert{}l_2\vert{}) + \ln(1 + e^{-\vert{}l_1+l_2\vert{}}) - \ln(1 + e^{-\vert{}l_1-l_2\vert{}})$$
- **g-function (Right Path)**: 이전 좌측 트리에서 복호 완료되어 피드백된 부분 코드 워드($u_{coded\_left}$)를 반영하여 우측 채널 우도를 갱신합니다.
  $$LLR_{right} = l_2 + (1 - 2 \cdot u_{coded\_left}) \cdot l_1$$

### 4) Xorshift64 난수 엔진 및 AWGN 채널
기존 C 표준 라이브러리의 `rand()` 함수가 가진 주기적 한계($2^{15}-1$)와 의사 난수 뭉침 현상에 따른 시뮬레이션 왜곡을 방지하기 위해 **Xorshift64 무작위 수 생성기**를 전면 탑재했습니다. 
이를 뷘커 박스-뮬러(Box-Muller) 변환과 결합하여, 높은 SNR 대역에서도 통계적 무결성을 유지하는 정밀한 백색 가우시안 잡음(AWGN) 환경을 제공합니다.

---

## 2. 시뮬레이션 환경 구성

| 파라미터 | 설정 값 | 비고 |
| :--- | :--- | :--- |
| **Block Length ($N$)** | 1024 | Polar Code 블록 크기 |
| **Information Bits ($K$)** | 512 | 유효 정보 데이터 크기 |
| **Code Rate ($R$)** | 0.5 | 부호율 |
| **Modulation** | BPSK | $0 \to +1$, $1 \to -1$ 매핑 |
| **SNR Range (Eb/No)** | 0.0 dB ~ 3.0 dB | 0.25 dB 간격 스캔 |
| **Stop Condition** | 1,000 Bit Errors | 각 SNR별 최소 에러 누적 한계 |

---

## 3.  실행 방법
```
# 1. 컴파일
gcc polar_code.c -o polar -lm

# 2. 프로그램 실행
./polar
```