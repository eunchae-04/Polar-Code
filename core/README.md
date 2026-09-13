# `core/` — Polar Code 최초 프로토타입 (이미지 전송 + Scenario A/B 비교)

이 폴더는 이 프로젝트의 **가장 처음 만들어진 버전**입니다. Polar Code 인코딩/디코딩 알고리즘
자체뿐 아니라, 실제 이미지를 전송해서 눈으로 결과를 확인하는 기능과 "오차 상황"을 흉내 낸
두 가지 시나리오 비교 기능까지 한 파일에 모두 들어 있습니다.

이후 발표/보고용으로는 핵심 알고리즘만 뽑아낸 [`scd/`](../scd/README.md)로 정리했지만,
이 폴더는 "왜 SC 디코더에 이런 결함(LLR 오추정, 마스크 설계 오차)이 위험한지"를 실제
이미지 화질 저하로 시각적으로 보여줄 수 있어 데모/설명용으로 유용합니다.

---

## 파일 구성

| 파일           | 설명                                                                    |
| -------------- | ----------------------------------------------------------------------- |
| `polar_code.c` | 이 폴더의 메인 소스                                                     |
| `image.png`    | (직접 준비 필요) 전송 테스트에 쓸 원본 이미지                           |
| `result/`      | 실행 시 자동 생성. 시뮬레이션 결과 txt, 전송된 이미지, 비교 이미지 저장 |

---

## 환경 설정

```bash
gcc -std=c11 -O2 polar_code.c -lm -o polar_sim.exe
```

**추가로 필요한 것 (이미지 전송 기능 때문에):**

1. Python 3 + Pillow
   ```bash
   pip install Pillow
   ```
2. 코드 상단의 아래 매크로를 **자신의 실제 Python 경로**로 수정해야 합니다.
   ```c
   #define PYTHON_CMD "C:/Users/user/anaconda3/python.exe"
   ```
   (Python 실행 파일 위치는 사람마다 다르므로, 그대로 두면 이미지 변환이 실패합니다.)
3. 같은 폴더에 `image.png`(또는 원하는 이름으로 바꾸고 `DEFAULT_LENA_INPUT` 매크로 수정)가 있어야 합니다.
4. Gnuplot이 설치되어 PATH에 등록되어 있어야 그래프가 자동으로 뜹니다.

---

## 실행 방법

```bash
# 기본 실행: Baseline / Scenario A / Scenario B 3종 시뮬레이션 + 시나리오별 이미지 전송
./polar_sim.exe

# 특정 alpha(LLR 오추정 배율) 하나만 baseline과 비교
./polar_sim.exe --alpha 2.0

# 임의의 이미지 파일을 특정 SNR로 직접 전송해보기
./polar_sim.exe --image input.pgm output.pgm 3.0
```

---

## 핵심 이론

### 1. 세 가지 시나리오

| 시나리오                             | 무엇이 잘못되었나                                                 | 코드에서                                     |
| ------------------------------------ | ----------------------------------------------------------------- | -------------------------------------------- |
| **Baseline**                         | 아무것도 잘못되지 않은 이상적인 경우                              | `llr_scale = 1.0`, 매 SNR에 맞는 마스크 사용 |
| **Scenario A (LLR Mismatch)**        | 디코더가 채널 잡음 세기(sigma²)를 실제보다 다르게(`alpha`배) 추정 | `sigma2_est = sigma2_true * llr_scale`       |
| **Scenario B (Design SNR Mismatch)** | Frozen Mask를 실제 채널 SNR이 아닌 다른 SNR(0dB 고정)로 설계      | `use_fixed_mask = true`, `fixed_mask_snr_db` |

실제 통신 시스템에서는 "채널 상태를 완벽하게 안다"는 가정이 항상 성립하지는 않습니다.
이 두 시나리오는 "그 가정이 깨지면 얼마나 손해를 보는가"를 정량적으로 보여줍니다.

### 2. 이미지 전송과 PSNR/MSE

BER/FER은 "비트가 몇 % 틀렸는가"라는 숫자로만 나오는데, 이걸 실제 이미지에 적용해보면
"그 정도 오류율이 실제로 어떻게 보이는가"를 직관적으로 확인할 수 있습니다.

- **MAE (Mean Absolute Error)**: 원본과 복원 이미지의 픽셀 값 차이 평균
- **MSE (Mean Squared Error)**: 차이를 제곱해서 평균 (큰 오류에 더 민감)
- **PSNR (Peak Signal-to-Noise Ratio)**: `10·log10(255² / MSE)` dB 단위. **값이 클수록
  원본과 더 비슷하다**는 뜻이며, 이미지 품질 평가에서 가장 널리 쓰이는 지표

### 3. 왜 PGM으로 변환하나

C 표준 라이브러리만으로는 PNG/JPG 같은 압축 이미지 포맷을 직접 읽고 쓸 수 없습니다.
그래서 Python + Pillow를 외부 프로세스로 호출(`system()`, `popen()`)해서 PNG를 **PGM(단순
그레이스케일 비트맵)**으로 변환한 뒤, C 코드는 이 단순한 포맷만 다룹니다. 결과를 다시
PNG로 저장할 때도 같은 방식으로 역변환합니다.

### 4. Polar Code 핵심 알고리즘

Frozen Mask 구성(Gaussian Approximation), 재귀형 인코더, SC(Min-Sum) 디코더의 원리는
`scd/` 폴더와 동일합니다. 자세한 설명은 프로젝트 최상위의 `polar_SCD_code_explained.md`를
참고하세요.

---

## 출력 결과물 (`result/` 폴더)

| 파일                                                                                | 내용                                              |
| ----------------------------------------------------------------------------------- | ------------------------------------------------- |
| `simulation_baseline.txt`, `simulation_scenario_a.txt`, `simulation_scenario_b.txt` | 시나리오별 SNR-BER-FER 데이터                     |
| `scenario_comparison_summary.txt`                                                   | Baseline 대비 A/B의 BER 배율(x) 요약 표           |
| `lena_output_*.png`                                                                 | 시나리오별로 전송/디코딩된 이미지                 |
| `*_preview.png`                                                                     | 원본·복원·차이 이미지를 나란히 붙인 비교용 이미지 |
