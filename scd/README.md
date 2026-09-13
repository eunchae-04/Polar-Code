# `scd/` — SC(Successive Cancellation) 디코더 핵심 시뮬레이션

`core/`에서 이미지 전송, 커맨드라인 시나리오 옵션(`--alpha`, `--image`) 등 부가 기능을
모두 제거하고, **"Polar Code SC 디코더의 BER/FER 성능을 시뮬레이션하고 그래프로 본다"**는
목적 하나에만 집중한 버전입니다. 이 폴더에는 본 실험 코드 외에, 난수 생성기 선택이
타당한지 검증하는 대조군 코드도 함께 들어 있습니다.

---

## 파일 구성

| 파일                                      | 설명                                                                                                |
| ----------------------------------------- | --------------------------------------------------------------------------------------------------- |
| `polar_code_SCD_sim.c`                    | **본 실험.** Xorshift64 난수 생성기를 사용하는 SC 디코더 시뮬레이션                                 |
| `polar_code_SCD_sim_libc_rand.c`          | **검증용 대조군.** 난수 생성기만 표준 `rand()`로 교체한 버전 (그 외 알고리즘은 본 실험과 100% 동일) |
| `result/simulation_results.txt`           | `polar_code_SCD_sim` 실행 시 자동 생성                                                              |
| `result/simulation_results_libc_rand.txt` | `polar_code_SCD_sim_libc_rand` 실행 시 자동 생성 (파일명이 달라 서로 덮어쓰지 않고 공존)            |

---

## 환경 설정

두 파일 모두 이미지/Python 의존성이 없어 동일한 방식으로 빌드합니다.

```bash
# 본 실험
gcc -O2 -Wall -o polar_code_SCD_sim polar_code_SCD_sim.c -lm

# 검증용 대조군
gcc -O2 -Wall -o polar_code_SCD_sim_libc_rand polar_code_SCD_sim_libc_rand.c -lm
```

Gnuplot이 PATH에 있어야 실행 후 그래프 창이 자동으로 뜹니다 (없어도 시뮬레이션 자체와
결과 파일 저장은 정상 동작합니다).

컴파일 시 아래와 같은 경고가 뜰 수 있는데, **에러가 아니라 무해한 경고**이니 무시해도
됩니다 (재귀 함수 구조를 gcc 정적 분석기가 완벽히 추적하지 못해서 나는 false positive).

```
warning: 'llr_left' may be used uninitialized [-Wmaybe-uninitialized]
```

---

## 실행 방법

```bash
./polar_code_SCD_sim               # 본 실험
./polar_code_SCD_sim_libc_rand     # 검증용 대조군
```

둘 다 옵션 없이 실행하면 Eb/No 0.0dB~3.0dB(0.25dB 간격)를 스윕하며, 각 SNR에서 오류가
10,000개(`TARGET_ERRORS`) 모일 때까지 몬테카를로 시뮬레이션을 반복합니다.

---

## 핵심 이론 (요약)

이 폴더의 알고리즘은 정확히 4단계로 이루어집니다.

1. **Frozen Mask 구성 (Gaussian Approximation)**: `phi`/`phi_inv`로 각 서브채널의
   신뢰도를 근사 계산하고, 가장 좋은 K=512개를 정보비트 위치로 지정
2. **Polar 인코더**: 재귀적 나비 구조(butterfly structure)로 정보비트+프로즌비트를 부호어로 변환
3. **AWGN 채널**: BPSK 변조 후 Box-Muller로 생성한 가우시안 잡음을 추가
4. **SC 디코더 (Min-Sum)**: f-function(`sign×min(|L1|,|L2|)`)과 g-function(결정 피드백)을
   재귀적으로 적용해, 0번 비트부터 순차적으로 하드 디시전

SC 디코더의 구조적 한계는 **"한 번 잘못 판단하면 그 뒤 비트까지 오류가 전파된다"**는
점이며, 이를 보완한 것이 [`scl/`](../scl/README.md) 폴더의 SCL 디코더입니다.

> **더 자세한 설명**: 이 코드를 한 줄 한 줄 왜 이렇게 작성했는지까지 다룬 문서가
> 프로젝트 최상위에 있습니다 → `polar_SCD_code_explained.md`

---

## 난수 생성기(RNG) 검증 — `polar_code_SCD_sim_libc_rand.c`

### 왜 이 대조군이 필요한가

본 실험은 표준 `rand()` 대신 Xorshift64라는 커스텀 난수 생성기를 씁니다. 이론적으로 더
낫다는 근거는 있지만("주기가 검증됨", "속도가 빠름", "플랫폼 의존성이 없음"), **실제
시뮬레이션 결과에 그 차이가 어떻게 나타나는지 직접 확인**해보기 위해 난수 생성기 부분만
표준 라이브러리 `rand()`/`srand()`로 교체한 버전을 만들었습니다. 그 외 알고리즘(Frozen
Mask, 인코더, SC 디코더, 시뮬레이션 루프)은 본 실험과 완전히 동일합니다.

### 무엇을 바꿨나 — `libc_rand64()`

```c
static uint64_t libc_rand64(void) {
    uint64_t r = 0;
    for (int i = 0; i < 5; i++) {
        r = (r << 15) ^ (uint64_t)(rand() & 0x7FFF);
    }
    return r;
}
```

`rand()`는 최소 보장 범위가 0~32767(15비트)뿐이라, 기존 `xorshift64()`가 주던 64비트
정수와 인터페이스를 맞추기 위해 **5번 뽑아서 이어붙였습니다** (15×4+4=64비트를 정확히
채우는 최소 반복 횟수). 이렇게 해야 `rand_double()`, `generate_info_bits()` 등 나머지
코드를 한 글자도 안 건드리고 재사용할 수 있습니다.

64비트로 확장한 이유는 단순히 "크게 만들자"가 아니라, `rand()`의 낮은 해상도(32768가지
값)가 Box-Muller의 `log(u1)` 계산에서 **잡음의 최댓값에 인위적인 상한선을 만드는 문제**를
피하기 위해서입니다 (자세한 내용은 최상위 `polar_xorshift64_rng_explained.md` 참고).

### 결론 (실험 결과 요약)

실제로 본 실험(Xorshift64)과 대조군(rand()) 결과를 겹쳐보면 BER/FER 곡선이 거의 겹칩니다.
다만 이 결과가 뜻하는 것은 **"이 정도 규모의 몬테카를로 실험에서는 두 RNG의 차이가 최종
통계치에 유의미한 영향을 주지 않는다"**이지, "rand()가 Xorshift64만큼 좋은 난수 생성기다"는
아닙니다. 특히 Windows UCRT의 `rand()`는 `RAND_MAX`가 32767로 작고 하위 비트 패턴 문제가
알려져 있어서, **최종 보고서/결과용 시뮬레이션은 Xorshift64(`polar_code_SCD_sim.c`)를
유지**하는 것을 권장합니다. `polar_code_SCD_sim_libc_rand.c`는 그 선택의 타당성을
뒷받침하는 검증 부록으로 활용합니다.

---

## 출력 결과물

| 파일                                      | 내용                                 |
| ----------------------------------------- | ------------------------------------ |
| `result/simulation_results.txt`           | 본 실험(Xorshift64) SNR-BER-FER 결과 |
| `result/simulation_results_libc_rand.txt` | 대조군(rand()) SNR-BER-FER 결과      |

두 파일 모두 `SNR(dB) BER FER` 형식으로 13개 행(0.0~3.0dB, 0.25dB 간격)이 저장되고,
각 실행 직후 Gnuplot으로 로그 스케일 Waterfall 그래프가 자동으로 표시됩니다. 두 결과를
겹쳐 그리면 RNG 차이에 따른 성능 편차를 시각적으로 비교할 수 있습니다.

## 참고

이 폴더의 본 실험 결과(`simulation_results.txt`)는 [`scl/`](../scl/README.md)과 성능을
비교하는 기준(baseline)으로도 사용됩니다.
