# Phase 1: multiclass particle classifier (held-out runs)

Trained on run periods (1, 2, 4) (1727696 tracks, balanced classes, 1984 trees), tested on periods (3, 5) (787181 tracks). Track selection: backtracked purity > 0.5, track score >= 0.5, length >= 5 cm; tracks with purity <= 0.5 kept and labelled other.

## One-vs-rest ROC area, test runs

| class | AUC |
|---|---|
| muon | 0.9427 |
| pion | 0.8915 |
| proton | 0.9550 |
| other | 0.8501 |

## Confusion matrix, test runs (rows true, columns predicted: mu, pi, p, other)

- muon: 235095 (76.3%), 48741 (15.8%), 2638 (0.9%), 21458 (7.0%)
- pion: 12997 (15.6%), 59474 (71.5%), 2769 (3.3%), 7897 (9.5%)
- proton: 3442 (1.8%), 14962 (7.8%), 163302 (85.0%), 10349 (5.4%)
- other: 21432 (10.5%), 34736 (17.0%), 37328 (18.3%), 110561 (54.2%)

## Pion identification on the selection's candidate pool (no muon candidate), test runs

### ROC area, pi vs rest

| group | n (pions) | P(pi), new | mp_pion_bdt split (CC1mu2pi) | mp_pion_bdt single | inclusive pion BDT | MIP BDT | LLR score |
|---|---|---|---|---|---|---|---|
| all | 453340 (54933) | 0.903 | 0.868 | 0.866 | 0.813 | 0.721 | 0.667 |
| FHC | 121744 (13894) | 0.904 | 0.870 | 0.867 | 0.813 | 0.721 | 0.667 |
| RHC | 331596 (41039) | 0.903 | 0.868 | 0.865 | 0.813 | 0.720 | 0.667 |
| pi+ vs rest | 433554 (35147) | 0.906 | 0.864 | 0.861 | 0.820 | 0.729 | 0.679 |
| pi- vs rest | 418193 (19786) | 0.897 | 0.876 | 0.873 | 0.801 | 0.707 | 0.645 |
| contained | 364466 (48666) | 0.907 | 0.877 | 0.874 | 0.859 | 0.819 | 0.766 |
| uncontained | 88874 (6267) | 0.842 | 0.811 | 0.808 | 0.616 | 0.433 | 0.344 |
| length < 20 cm | 162190 (16552) | 0.898 | 0.885 | 0.884 | 0.849 | 0.835 | 0.797 |
| length >= 20 cm | 291150 (38381) | 0.904 | 0.861 | 0.859 | 0.798 | 0.651 | 0.584 |

### ROC area, pi vs p

| group | n (pions) | P(pi), new | mp_pion_bdt split (CC1mu2pi) | mp_pion_bdt single | inclusive pion BDT | MIP BDT | LLR score |
|---|---|---|---|---|---|---|---|
| all | 453340 (54933) | 0.956 | 0.945 | 0.942 | 0.961 | 0.957 | 0.893 |
| FHC | 121744 (13894) | 0.958 | 0.947 | 0.945 | 0.961 | 0.958 | 0.895 |
| RHC | 331596 (41039) | 0.955 | 0.944 | 0.941 | 0.960 | 0.956 | 0.893 |
| pi+ vs rest | 433554 (35147) | 0.957 | 0.942 | 0.940 | 0.963 | 0.960 | 0.900 |
| pi- vs rest | 418193 (19786) | 0.953 | 0.948 | 0.946 | 0.957 | 0.951 | 0.881 |
| contained | 364466 (48666) | 0.962 | 0.951 | 0.949 | 0.962 | 0.958 | 0.907 |
| uncontained | 88874 (6267) | 0.884 | 0.855 | 0.847 | 0.942 | 0.940 | 0.741 |
| length < 20 cm | 162190 (16552) | 0.954 | 0.957 | 0.956 | 0.942 | 0.938 | 0.922 |
| length >= 20 cm | 291150 (38381) | 0.954 | 0.938 | 0.935 | 0.966 | 0.961 | 0.871 |

### ROC area, pi vs mu

| group | n (pions) | P(pi), new | mp_pion_bdt split (CC1mu2pi) | mp_pion_bdt single | inclusive pion BDT | MIP BDT | LLR score |
|---|---|---|---|---|---|---|---|
| all | 453340 (54933) | 0.852 | 0.804 | 0.804 | 0.558 | 0.271 | 0.251 |
| FHC | 121744 (13894) | 0.846 | 0.801 | 0.800 | 0.554 | 0.276 | 0.255 |
| RHC | 331596 (41039) | 0.853 | 0.806 | 0.805 | 0.559 | 0.268 | 0.249 |
| pi+ vs rest | 433554 (35147) | 0.857 | 0.798 | 0.798 | 0.571 | 0.283 | 0.267 |
| pi- vs rest | 418193 (19786) | 0.842 | 0.815 | 0.814 | 0.534 | 0.249 | 0.224 |
| contained | 364466 (48666) | 0.790 | 0.761 | 0.758 | 0.596 | 0.406 | 0.352 |
| uncontained | 88874 (6267) | 0.838 | 0.812 | 0.812 | 0.493 | 0.242 | 0.222 |
| length < 20 cm | 162190 (16552) | 0.670 | 0.594 | 0.592 | 0.479 | 0.452 | 0.510 |
| length >= 20 cm | 291150 (38381) | 0.877 | 0.822 | 0.821 | 0.611 | 0.272 | 0.257 |

## Working point

The split mp_pion_bdt at its CC1mu2pi cuts (soft 0.1, hard -0.2): pion efficiency 91.8%, pion purity 27.6%, protons passing 13.4%, muons passing 52.9%.
The new P(pi) at the same pion efficiency: purity 31.5%, protons passing 11.7%, muons passing 44.8%.

## Largest feature importances (gain share)

- llr: 0.248
- contained: 0.124
- n_trk_daughters: 0.096
- n_shr_daughters: 0.081
- llr_y: 0.056
- dist: 0.051
- range_mom_mu: 0.048
- len: 0.042
- start_contained: 0.041
- llr_v: 0.031
- ts: 0.029
- mcs_over_range: 0.027
- e_proton: 0.026
- llr_u: 0.025
- planehits_v: 0.023
