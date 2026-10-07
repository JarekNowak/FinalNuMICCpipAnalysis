# Phase 1: multiclass particle classifier (held-out runs)

Trained on run periods (1, 2, 4) (1727696 tracks, balanced classes, 1996 trees), tested on periods (3, 5) (787181 tracks). Track selection: backtracked purity > 0.5, track score >= 0.5, length >= 5 cm; tracks with purity <= 0.5 kept and labelled other.

## One-vs-rest ROC area, test runs

| class | AUC |
|---|---|
| muon | 0.9583 |
| pion | 0.9177 |
| proton | 0.9680 |
| other | 0.8765 |

## Confusion matrix, test runs (rows true, columns predicted: mu, pi, p, other)

- muon: 246625 (80.1%), 42910 (13.9%), 834 (0.3%), 17563 (5.7%)
- pion: 12393 (14.9%), 62510 (75.2%), 1374 (1.7%), 6860 (8.3%)
- proton: 1511 (0.8%), 8322 (4.3%), 173569 (90.4%), 8653 (4.5%)
- other: 17734 (8.7%), 29967 (14.7%), 37059 (18.2%), 119297 (58.5%)

## Pion identification on the selection's candidate pool (no muon candidate), test runs

### ROC area, pi vs rest

| group | n (pions) | P(pi), new | mp_pion_bdt split (CC1mu2pi) | mp_pion_bdt single | inclusive pion BDT | MIP BDT | LLR score |
|---|---|---|---|---|---|---|---|
| all | 453340 (54933) | 0.930 | 0.868 | 0.866 | 0.813 | 0.721 | 0.667 |
| FHC | 121744 (13894) | 0.931 | 0.870 | 0.867 | 0.813 | 0.721 | 0.667 |
| RHC | 331596 (41039) | 0.930 | 0.868 | 0.865 | 0.813 | 0.720 | 0.667 |
| pi+ vs rest | 433554 (35147) | 0.931 | 0.864 | 0.861 | 0.820 | 0.729 | 0.679 |
| pi- vs rest | 418193 (19786) | 0.928 | 0.876 | 0.873 | 0.801 | 0.707 | 0.645 |
| contained | 364466 (48666) | 0.936 | 0.877 | 0.874 | 0.859 | 0.819 | 0.766 |
| uncontained | 88874 (6267) | 0.864 | 0.811 | 0.808 | 0.616 | 0.433 | 0.344 |
| length < 20 cm | 162190 (16552) | 0.933 | 0.885 | 0.884 | 0.849 | 0.835 | 0.797 |
| length >= 20 cm | 291150 (38381) | 0.929 | 0.861 | 0.859 | 0.798 | 0.651 | 0.584 |

### ROC area, pi vs p

| group | n (pions) | P(pi), new | mp_pion_bdt split (CC1mu2pi) | mp_pion_bdt single | inclusive pion BDT | MIP BDT | LLR score |
|---|---|---|---|---|---|---|---|
| all | 453340 (54933) | 0.981 | 0.945 | 0.942 | 0.961 | 0.957 | 0.893 |
| FHC | 121744 (13894) | 0.982 | 0.947 | 0.945 | 0.961 | 0.958 | 0.895 |
| RHC | 331596 (41039) | 0.980 | 0.944 | 0.941 | 0.960 | 0.956 | 0.893 |
| pi+ vs rest | 433554 (35147) | 0.981 | 0.942 | 0.940 | 0.963 | 0.960 | 0.900 |
| pi- vs rest | 418193 (19786) | 0.980 | 0.948 | 0.946 | 0.957 | 0.951 | 0.881 |
| contained | 364466 (48666) | 0.984 | 0.951 | 0.949 | 0.962 | 0.958 | 0.907 |
| uncontained | 88874 (6267) | 0.942 | 0.855 | 0.847 | 0.942 | 0.940 | 0.741 |
| length < 20 cm | 162190 (16552) | 0.980 | 0.957 | 0.956 | 0.942 | 0.938 | 0.922 |
| length >= 20 cm | 291150 (38381) | 0.981 | 0.938 | 0.935 | 0.966 | 0.961 | 0.871 |

### ROC area, pi vs mu

| group | n (pions) | P(pi), new | mp_pion_bdt split (CC1mu2pi) | mp_pion_bdt single | inclusive pion BDT | MIP BDT | LLR score |
|---|---|---|---|---|---|---|---|
| all | 453340 (54933) | 0.877 | 0.804 | 0.804 | 0.558 | 0.271 | 0.251 |
| FHC | 121744 (13894) | 0.873 | 0.801 | 0.800 | 0.554 | 0.276 | 0.255 |
| RHC | 331596 (41039) | 0.879 | 0.806 | 0.805 | 0.559 | 0.268 | 0.249 |
| pi+ vs rest | 433554 (35147) | 0.879 | 0.798 | 0.798 | 0.571 | 0.283 | 0.267 |
| pi- vs rest | 418193 (19786) | 0.874 | 0.815 | 0.814 | 0.534 | 0.249 | 0.224 |
| contained | 364466 (48666) | 0.830 | 0.761 | 0.758 | 0.596 | 0.406 | 0.352 |
| uncontained | 88874 (6267) | 0.849 | 0.812 | 0.812 | 0.493 | 0.242 | 0.222 |
| length < 20 cm | 162190 (16552) | 0.700 | 0.594 | 0.592 | 0.479 | 0.452 | 0.510 |
| length >= 20 cm | 291150 (38381) | 0.898 | 0.822 | 0.821 | 0.611 | 0.272 | 0.257 |

## Working point

The split mp_pion_bdt at its CC1mu2pi cuts (soft 0.1, hard -0.2): pion efficiency 91.8%, pion purity 27.6%, protons passing 13.4%, muons passing 52.9%.
The new P(pi) at the same pion efficiency: purity 38.1%, protons passing 5.4%, muons passing 38.4%.

## Largest feature importances (gain share)

- llr: 0.210
- contained: 0.087
- bragg_mip: 0.086
- n_trk_daughters: 0.069
- n_shr_daughters: 0.060
- dist: 0.047
- range_mom_mu: 0.043
- llr_y: 0.042
- e_proton: 0.038
- len: 0.033
- bragg_p: 0.029
- start_contained: 0.028
- mcs_over_range: 0.024
- bragg_mu: 0.024
- bragg_p_v: 0.021
