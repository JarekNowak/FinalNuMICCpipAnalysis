# Phase 1: multiclass particle classifier (held-out runs)

Trained on run periods (1, 2, 4) (1329038 tracks, balanced classes, 824 trees), tested on periods (3, 5) (600182 tracks). Track selection: backtracked purity > 0.5, track score >= 0.5, length >= 5 cm.

## One-vs-rest ROC area, test runs

| class | AUC |
|---|---|
| muon | 0.9769 |
| pion | 0.9367 |
| proton | 0.9966 |
| other | 0.9861 |

## Confusion matrix, test runs (rows true, columns predicted: mu, pi, p, other)

- muon: 264202 (85.8%), 39226 (12.7%), 1027 (0.3%), 3477 (1.1%)
- pion: 12755 (15.3%), 64797 (77.9%), 1347 (1.6%), 4238 (5.1%)
- proton: 1189 (0.6%), 6932 (3.6%), 182262 (94.9%), 1672 (0.9%)
- other: 176 (1.0%), 1255 (7.4%), 787 (4.6%), 14840 (87.0%)

## Pion identification on the selection's candidate pool (no muon candidate), test runs

### ROC area, pi vs rest

| group | n (pions) | P(pi), new | mp_pion_bdt split (CC1mu2pi) | mp_pion_bdt single | inclusive pion BDT | MIP BDT | LLR score |
|---|---|---|---|---|---|---|---|
| all | 346987 (54933) | 0.954 | 0.888 | 0.886 | 0.818 | 0.718 | 0.670 |
| FHC | 93630 (13894) | 0.953 | 0.889 | 0.886 | 0.817 | 0.720 | 0.671 |
| RHC | 253357 (41039) | 0.954 | 0.887 | 0.886 | 0.819 | 0.717 | 0.669 |
| pi+ vs rest | 327201 (35147) | 0.953 | 0.884 | 0.882 | 0.825 | 0.725 | 0.680 |
| pi- vs rest | 311840 (19786) | 0.955 | 0.894 | 0.892 | 0.807 | 0.706 | 0.651 |
| contained | 276687 (48666) | 0.960 | 0.900 | 0.898 | 0.879 | 0.835 | 0.783 |
| uncontained | 70300 (6267) | 0.884 | 0.817 | 0.816 | 0.580 | 0.378 | 0.324 |
| length < 20 cm | 118906 (16552) | 0.956 | 0.908 | 0.907 | 0.880 | 0.871 | 0.851 |
| length >= 20 cm | 228081 (38381) | 0.954 | 0.877 | 0.875 | 0.797 | 0.637 | 0.582 |

### ROC area, pi vs p

| group | n (pions) | P(pi), new | mp_pion_bdt split (CC1mu2pi) | mp_pion_bdt single | inclusive pion BDT | MIP BDT | LLR score |
|---|---|---|---|---|---|---|---|
| all | 346987 (54933) | 0.985 | 0.945 | 0.942 | 0.961 | 0.957 | 0.893 |
| FHC | 93630 (13894) | 0.986 | 0.947 | 0.945 | 0.961 | 0.958 | 0.895 |
| RHC | 253357 (41039) | 0.985 | 0.944 | 0.941 | 0.960 | 0.956 | 0.893 |
| pi+ vs rest | 327201 (35147) | 0.985 | 0.942 | 0.940 | 0.963 | 0.960 | 0.900 |
| pi- vs rest | 311840 (19786) | 0.986 | 0.948 | 0.946 | 0.957 | 0.951 | 0.881 |
| contained | 276687 (48666) | 0.988 | 0.951 | 0.949 | 0.962 | 0.958 | 0.907 |
| uncontained | 70300 (6267) | 0.954 | 0.855 | 0.847 | 0.942 | 0.940 | 0.741 |
| length < 20 cm | 118906 (16552) | 0.981 | 0.957 | 0.956 | 0.942 | 0.938 | 0.922 |
| length >= 20 cm | 228081 (38381) | 0.987 | 0.938 | 0.935 | 0.966 | 0.961 | 0.871 |

### ROC area, pi vs mu

| group | n (pions) | P(pi), new | mp_pion_bdt split (CC1mu2pi) | mp_pion_bdt single | inclusive pion BDT | MIP BDT | LLR score |
|---|---|---|---|---|---|---|---|
| all | 346987 (54933) | 0.895 | 0.804 | 0.804 | 0.558 | 0.271 | 0.251 |
| FHC | 93630 (13894) | 0.893 | 0.801 | 0.800 | 0.554 | 0.276 | 0.255 |
| RHC | 253357 (41039) | 0.896 | 0.806 | 0.805 | 0.559 | 0.268 | 0.249 |
| pi+ vs rest | 327201 (35147) | 0.894 | 0.798 | 0.798 | 0.571 | 0.283 | 0.267 |
| pi- vs rest | 311840 (19786) | 0.898 | 0.815 | 0.814 | 0.534 | 0.249 | 0.224 |
| contained | 276687 (48666) | 0.860 | 0.761 | 0.758 | 0.596 | 0.406 | 0.352 |
| uncontained | 70300 (6267) | 0.868 | 0.812 | 0.812 | 0.493 | 0.242 | 0.222 |
| length < 20 cm | 118906 (16552) | 0.747 | 0.594 | 0.592 | 0.479 | 0.452 | 0.510 |
| length >= 20 cm | 228081 (38381) | 0.916 | 0.822 | 0.821 | 0.611 | 0.272 | 0.257 |

## Working point

The split mp_pion_bdt at its CC1mu2pi cuts (soft 0.1, hard -0.2): pion efficiency 91.8%, pion purity 37.3%, protons passing 13.4%, muons passing 52.9%.
The new P(pi) at the same pion efficiency: purity 53.4%, protons passing 4.4%, muons passing 34.9%.

## Largest feature importances (gain share)

- chipr: 0.118
- llr: 0.116
- ts: 0.078
- bragg_mip: 0.057
- contained: 0.055
- n_descendents: 0.047
- e_proton: 0.041
- len: 0.025
- n_shr_daughters: 0.024
- range_mom_mu: 0.022
- n_trk_daughters: 0.019
- chimu: 0.018
- defl_mean: 0.015
- bragg_pion: 0.014
- chika: 0.014
