# Polar Code SC/SCL 디코더 시뮬레이션 프로젝트

Polar Code(N=1024, K=512)의 SC(Successive Cancellation) 및 SCL(Successive Cancellation
List) 디코더를 C언어로 구현하고, AWGN 채널에서의 BER/FER 성능을 몬테카를로 시뮬레이션으로
검증하는 캡스톤 디자인 프로젝트입니다.

---

## 폴더 구조

```
Polar-Code/
├── core/         # 최초 프로토타입 (이미지 전송 + Scenario A/B 비교)
├── scd/          # SC 디코더 핵심 시뮬레이션 (+ RNG 검증용 대조군 포함)
├── scl/          # SCL(List) 디코더 확장
└── simab/        # Scenario A/B 민감도 스윕 분석 (core의 확장판)
```

| 폴더                          | 무엇을 하는가                                                                                                                                                                         | 자세히                      |
| ----------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------- |
| [`core/`](./core/README.md)   | Baseline / LLR 오추정(Scenario A) / Design SNR 오추정(Scenario B) 비교 + 이미지(Lena) 전송 데모                                                                                       | [README](./core/README.md)  |
| [`scd/`](./scd/README.md)     | 이미지·시나리오 기능을 걷어내고 SC 디코더 BER/FER 시뮬레이션에만 집중. 커스텀 RNG(Xorshift64)가 표준 `rand()`와 비교해 결과를 왜곡하지 않는지 검증하는 대조군 코드도 같은 폴더에 포함 | [README](./scd/README.md)   |
| [`scl/`](./scl/README.md)     | 여러 경로를 동시에 유지하는 SCL 디코더로 SC의 한계 보완                                                                                                                               | [README](./scl/README.md)   |
| [`simab/`](./simab/README.md) | LLR 오추정(alpha)·Design SNR 오추정을 넓은 범위로 스윕하며 민감도 분석                                                                                                                | [README](./simab/README.md) |

---

## 진행 히스토리 (한눈에 보기)

```
core (이미지 전송 + 3-시나리오 비교 프로토타입)
   │
   ├─→ scd  (SC 알고리즘만 순수하게 분리 · 정리)
   │     │     └─ 같은 폴더 안에 RNG 검증용 대조군(rand() 버전) 포함
   │     └─→ scl (SC를 List 디코딩으로 확장)
   │
   └─→ simab (core의 시나리오 비교를 넓은 범위 스윕 분석으로 확장)
```

1. **core**: 처음엔 Polar Code 인코딩/디코딩에 이미지 전송(PSNR/MSE)까지 한 파일에 다 넣고,
   "LLR을 잘못 추정하면(Scenario A)", "Frozen Mask를 잘못된 SNR로 설계하면(Scenario B)"
   성능이 얼마나 나빠지는지 보는 프로토타입을 만들었습니다.
2. **scd**: 발표/보고용으로 핵심 알고리즘(SC 디코더)만 깔끔하게 분리했습니다. 이미지 전송,
   커맨드라인 옵션을 모두 제거하고 BER/FER 시뮬레이션 + gnuplot 시각화에만 집중합니다.
   같은 폴더 안에, 커스텀 난수 생성기(Xorshift64)가 표준 라이브러리 `rand()`와 비교했을 때
   결과를 왜곡하지 않는지를 검증하는 대조군 코드도 함께 있습니다.
3. **scl**: scd의 SC 디코더를 Path Metric 기반 List 디코더(SCL)로 확장해서, 단일 경로
   탐색의 구조적 한계를 보완했습니다.
4. **simab**: core에서 다뤘던 "LLR 오추정", "Design SNR 오추정" 시나리오를 단일 값이
   아니라 넓은 범위로 스윕하며, 오차가 커질수록 성능이 얼마나 나빠지는지 민감도 곡선을
   뽑아내는 확장판입니다.

---

## 공통 환경 설정

모든 폴더가 아래 도구를 공통으로 사용합니다.

- **컴파일러**: gcc (Windows는 MSYS2/MinGW-w64 UCRT64 환경 기준)
- **시각화**: [Gnuplot](http://www.gnuplot.org/) — PATH에 등록되어 있어야 `popen("gnuplot -persist", ...)`이 동작합니다.
- **이미지 전송 기능이 있는 폴더(`core`, `simab`)만 추가로 필요**:
  - Python 3 + [Pillow](https://pillow.readthedocs.io/) (`pip install Pillow`)
  - 코드 상단의 `PYTHON_CMD` 매크로를 **자신의 실제 Python 실행 파일 경로**로 수정해야 합니다.
  - 같은 폴더에 전송할 원본 이미지(`image.png`, 기본값)가 있어야 합니다.

공통 빌드 패턴:

```bash
gcc -O2 -Wall -o <실행파일명> <소스파일명>.c -lm
./<실행파일명>
```

> Windows PowerShell에서 파일을 옮기거나 다운로드한 뒤에는, 컴파일 전에
> `(Get-Item 파일명).Length`로 크기를 한 번 확인하는 습관을 들이면 좋습니다.
> 0바이트거나 비정상적으로 작으면 전송 중 파일이 손상된 것입니다.

---

## 공통 핵심 이론

다섯 폴더 모두 아래 요소를 기반으로 합니다 (세부 구현 차이는 각 폴더 README 참고).

- **Polar Code (N=1024, K=512, Rate 1/2)**: 채널 극화(channel polarization) 현상을 이용해
  N개의 서브채널을 "좋은 채널(정보비트)"과 "나쁜 채널(프로즌비트)"로 나누는 부호화 방식
- **Gaussian Approximation 기반 Frozen Mask**: 정확한 밀도진화(Density Evolution) 대신
  빠른 근사식(`phi`/`phi_inv`)으로 각 서브채널의 신뢰도를 계산
- **AWGN 채널 모델**: BPSK 변조 + Box-Muller 변환으로 생성한 가우시안 잡음
- **Xorshift64 난수 생성기**: 표준 `rand()`보다 품질과 속도가 보장되는 커스텀 RNG
- **BER/FER 몬테카를로 시뮬레이션**: 오류가 일정 개수(`TARGET_ERRORS`) 모일 때까지 반복하여
  통계적으로 신뢰할 수 있는 오류율을 산출

각 개념에 대한 자세한 설명(코드 라인별 해설 포함)은 아래 문서를 참고하세요.

- `polar_SCD_code_explained.md` — SC 디코더 전체 코드 설명
- `polar_SCL_code_explained.md` — SCL 디코더 코드 설명 (SC와의 차이점 중심)
- `polar_xorshift64_rng_explained.md` — Xorshift64 난수 생성기 설명

---

## 다음 계획

- `simab`의 스윕 결과와 `scl`의 List 디코딩 결과를 결합해, "오차 상황에서 SCL이 SC보다
  얼마나 더 견고한가"를 비교하는 실험 추가 검토
- `scd`의 RNG 검증 결과를 캡스톤 보고서에 부록으로 첨부
- List Size(L)별 성능 비교, CRC-aided SCL 검토 (기존 주간 보고에서 제시한 방향)
