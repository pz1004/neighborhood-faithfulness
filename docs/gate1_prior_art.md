# Gate 1: prior-art audit of the combination

2026-09-28. Plan Section 6, Gate 1. Four searches ran in parallel, one per literature, each
with about 30 queries and each told to open every source it cites. Entries marked
*[not opened]* were seen only in search results or reference lists and are unverified.

**Verdict: Gate 1 passes, with the novelty narrowed.** No search found a work that states
per-point simultaneous confidence statements about a display's neighborhood fidelity
against the noise-free truth (claim C1 together with C2 or C3). Almost every component is
anticipated somewhere. What remains new is the combination, the set-level reach bound, and
the measurements (Section 3).

## 1. Claims audited

| id | claim |
|---|---|
| C1 | Per-point, simultaneous (family-wise error controlled) confidence statements about a display's neighborhood fidelity, stated against the noise-free truth, from replicates. |
| C2 | The reach target (an order statistic of true ranks over displayed neighbors), with intervals uniform over tolerances (ε, k′). |
| C3 | A three-way per-point map (certified faithful / certified distorted / unresolved) with error control. |
| C4 | Finite-sample distance bounds between noisy means (studentized projection bound), union-bounded over pairs and turned into nearest-neighbor rank and reach bounds; the members-at-upper / outsiders-at-lower configuration. |
| C5 | A planning rule linking the certifiable tolerance to intrinsic dimension and noise-to-spacing ratio. |
| C6 | How often naive neighborhood maps computed on noisy data mislabel points against the truth. |
| C7 | Data thinning to create pseudo-replicates for single-cell data (a follow-on). |
| C8 | Validation against a near-truth built from held-out genuine replicates. |

## 2. Verdicts

| claim | verdict | closest work | what it lacks |
|---|---|---|---|
| C1 | **not found** (all four searches) | scDEED (Xia, Lee, Li 2024); EMBEDR (Johnson et al. 2022); DynamicViz; NESS; Liu, Ma, Zhong 2025; Zabel, Hennig, Nieselt 2026; Josse, Wager, Husson 2016 (pointwise confidence ellipses for PCA positions of the noise-free signal) | The reference is the noisy data or run-to-run variability. None controls family-wise error or gives a rank interval. Josse et al. targets positions, not neighborhoods, pointwise and for PCA only. |
| C2 | partly anticipated | Mogstad, Romano, Shaikh, Wilhelm 2024 (simultaneous rank sets from pairwise differences, "for any objects"); Mahamunulu 1967 (select a subset of size s containing at least c of the t best, which is exactly reach ≤ k′); Roy Chaudhuri & Kalyanakrishnan 2019; Ram et al. 2009 (rank-approximate nearest neighbors); Venna et al. 2010 (per-query precision with k retrieved and r relevant: faithful at tolerance k′ is precision ≥ 1−ε with r = k′, so reach inverts it in k′); Lueks et al. 2011 and Mokbel et al. 2013 (tolerated rank error, drawn per point on the map) | The quantity exists in evaluation, computed on observed data. Nobody gives intervals for it, ranks distances from each anchor with inference, or assesses an externally chosen set. Intervals uniform in (ε, k′) were not found, though they follow once rank intervals are simultaneous. |
| C3 | partly anticipated | scDEED (trustworthy / intermediate / dubious from uncorrected permutation-null percentiles, drawn on the display); visual-analytics distortion maps (Aupetit 2007; CheckViz 2011; Martins et al. 2014, 2015; Jeon et al. 2022 reliability maps; Nonato & Aupetit 2019); uncertain databases (certain / possible / impossible answers; true hits / true drops / candidates) | No error control, and not stated against the truth. scDEED's null is random neighbors, so its labels answer "better than random?" rather than "faithful at tolerance k′?". |
| C4 | **anticipated** for one neighbor's rank; set-level reach not found | Soliman & Ilyas 2009 (rank interval = [dominators + 1, n − dominated]); Feng, Glavic, Kennedy 2023; Bernecker et al. 2011 (domination count with every object uncertain); Emrich et al. 2010 (min/max-distance tests are sufficient but not complete); Fassina, Ranzato, Zanella 2024; Mason, Tripathy, Nowak 2019 (distance intervals union-bounded over pairs to recover the nearest-neighbor graph) | Frequentist confidence sets from replicates were not found, nor the order-statistic (reach) extreme over a set, nor the studentized projection bound. All are routine given the rest. |
| C5 | not found (analogues only) | Mahamunulu 1967 (indifference-zone sample sizes); Kannan, Shin, Woodruff 2026 (nearest-neighbor recovery requires σ ≲ k^(−1/4)); Amsaleg et al. 2021 and Furon 2020 (the perturbation that changes a neighbor's rank shrinks as local intrinsic dimension grows); Houle 2013; Snoeck, Bergam, Verma 2025 | No rule maps local dimension and noise-to-spacing ratio to a certifiable tolerance. |
| C6 | partly anticipated (aggregate only); per point not found | Lin & Fukuyama 2024 (trustworthiness against noisy data selects different hyperparameters than against the signal); Ahlmann-Eltze & Huber 2023 (kNN graphs from 10× downsampled counts share about 2.9 of 50 neighbors with deep-sequencing neighbors); DEMaP (Moon et al. 2019, one global correlation with noiseless geodesics); Venna et al. 2010 (mean precision against latent ground truth); Li, Kluger, Landa 2025 (distance deflation under heteroskedastic noise) | No per-point mislabel rates of plug-in fidelity maps against the truth. |
| C7 | not found (ingredients exist) | count splitting, data thinning, generalized thinning (Neufeld, Dharamshi, Gao, Witten and coauthors); Batson et al. 2019; Hivert et al. 2024 | No use of folds to evaluate kNN graphs or display neighborhoods. |
| C8 | practice partly anticipated; claim not found | RSA split-half reliability and noise ceilings; copairs (Kalinin et al. 2025, replicate retrieval with Benjamini–Hochberg FDR) | Nobody builds per-object near-truth neighborhoods from held-out replicates to validate per-point certificates. |

## 3. What the paper can claim

**Not new; must be cited and positioned:**
- *Simultaneous rank inference from pairwise differences*: Mogstad et al. 2024, including a
  stepwise improvement that our Bonferroni engines lack. Klein, Wright, Wieczorek 2020.
  Al Mohamad, Goeman, van Zwet 2022.
- *Rank intervals and three-valued answers from dominance under interval uncertainty*:
  Soliman & Ilyas 2009; Feng et al. 2023; Bernecker et al. 2011. Emrich et al. 2010 showed
  that box (min/max) relaxations are incomplete.
- *The tolerance notion*: Mahamunulu 1967; Roy Chaudhuri & Kalyanakrishnan 2019; Ram et al.
  2009; the co-ranking literature (Lee & Verleysen 2009; Lueks et al. 2011).
- *Distance intervals union-bounded over pairs for nearest-neighbor graphs*: Mason et al.
  2019.
- *Unbiased noise-free distances*: crossnobis (Walther et al. 2016; Diedrichsen et al. 2016).
- *Noisy-data fidelity metrics overfit noise*: Lin & Fukuyama 2024; Ahlmann-Eltze & Huber
  2023.

**New, as far as four searches found:**
1. The statement itself: simultaneous, finite-sample intervals, per displayed point, on
   neighborhood reach against the noise-free truth, valid for any display and uniform over
   tolerances.
2. The set-level reach extreme over a box, and the studentized pairwise distance bound
   turned into reach bounds. Both are routine, so they are presented as tools, not
   contributions.
3. The measurements: how often naive maps mislabel against the truth, which engine
   survives which noise model, and the planning rule.

## 4. Threats to answer in the paper

- **scDEED is the output-form competitor.** C6 should measure scDEED's own labels against
  the truth, not only the plug-in map's. This is added to Gate 3.
- **Mogstad et al. is the engine competitor.** Their stepwise procedure applied to
  distance differences is a baseline to add, after Gate 2.
- **Bernecker / Emrich is the database reading.** A reviewer from databases will see C2–C4
  as domination counting applied to confidence regions. The answer is the statistical
  layer and the set-level target, not the combinatorics.
- **Vishwanath & Arias-Castro 2025** (arXiv:2510.22452; abstract read on 2026-09-28) is
  closer than the search suggested. It builds *uniform* confidence sets for the latent
  configuration of classical MDS on noisy dissimilarities, up to rigid motions. The sets
  come from distributional limits and a multiplier bootstrap that adapts to heteroscedastic
  noise.
  - Rigid motions preserve distances, so neighbor statements about the latent configuration
    could be derived from these sets, as ours are from balls. In that sense it is an
    asymptotic ball engine for dissimilarity data.
  - It is not about a display: no nonlinear display, and no comparison of displayed with
    true neighbors. Its abstract makes no neighbor statement, and its guarantee is
    asymptotic.
  - It must be cited as the closest simultaneous-inference work. Its full text must be read
    before any submission.

## 5. Technical notes from the audit

- **Wording.** The pairwise bound is *valid* in finite samples, not *exact* in the sense of
  tight. Inverting the noncentral χ² with a known scale is never looser than the projection
  bound. At D = 20 and x = 29 (spacing about 3τ) the upper bounds are 5.24τ against 7.03τ
  at α = 0.05, and 9.25τ against 10.58τ at α = 10⁻⁷.
- **The gain in the regime that matters.** For digits at R = 8 the gap is small for the
  upper bound, at about 3%. For the lower bound it is about 1.3τ out of a half-width near
  6.8τ, because the orthogonal term is bounded at its far quantile instead of jointly.
  Making that exact with an estimated per-object scale needs a worst-case argument that
  does not follow from the Jensen step used for the upper bound. It is a follow-on
  refinement (noncentral inversion with a studentized scale), measured against E-pair-t
  once Gate 2 is decided.

## 6. Citations (opened by the searches unless marked)

- **Statistics:**
  - Mogstad, Romano, Shaikh, Wilhelm, *Rev. Econ. Stud.* 91(1):476–518, 2024.
  - Klein, Wright, Wieczorek, *JRSS-C* 69(3):589–606, 2020.
  - Al Mohamad, Goeman, van Zwet, arXiv:1812.05507.
  - Hall & Miller, *Ann. Statist.* 37(6B), 2009.
  - Chandra & Xie, arXiv:2512.00316.
  - Mahamunulu, *Ann. Math. Statist.* 38(4):1079–1091, 1967.
  - Roy Chaudhuri & Kalyanakrishnan, ICML 2019.
  - Mason, Tripathy, Nowak, NeurIPS 2019.
  - Jain, Jamieson, Nowak, NeurIPS 2016.
  - Diedrichsen, Provost, Zareamoghaddam, arXiv:1607.01371.
  - Walther et al., *NeuroImage* 137, 2016.
  - Schütt et al., *eLife* 2023.
  - Nili et al., *PLoS ONE* 2020.
  - Zhang, Lee, Lei, arXiv:2408.02060.
  - Kim & Ramdas, arXiv:2503.21639.
  - Li, Kluger, Landa, arXiv:2507.18520.
  - Kannan, Shin, Woodruff, arXiv:2604.03831.
  - Lin & Fukuyama, *PLoS Comput. Biol.* 2024.
  - *[not opened]*: Xie, Singh, Zhang, *JASA* 2009; Kent & Hainsworth 1995; Reiser 2001;
    Houle 2013.
- **DR quality and visual analytics:**
  - Venna, Peltonen, Nybo, Aidos, Kaski, *JMLR* 11:451–490, 2010.
  - Lueks, Mokbel, Biehl, Hammer, arXiv:1110.3917.
  - Nonato & Aupetit, *IEEE TVCG* 25:2650–2673, 2019.
  - Jeon et al., *IEEE TVCG* 28:551–561, 2022; ZADU, arXiv:2308.00282; CHI 2025 survey, arXiv:2501.10168.
  - Görtler et al., *IEEE TVCG* 26(1):822–831, 2020 (UAPCA).
  - Bachmann, Hennig, Kobak, ECML PKDD 2022 (Wasserstein t-SNE).
  - Josse, Wager, Husson, *JCGS* 25(1):28–48, 2016.
  - Nylund, Correll, Padilla, arXiv:2608.14513.
  - Snoeck, Bergam, Verma, arXiv:2508.07119.
  - Furon, arXiv:2010.00990.
  - *[metadata only or not opened]*: Aupetit, *Neurocomputing* 2007; Lespinats & Aupetit, *CGF* 2011 (CheckViz); Martins et al. 2014, 2015; Lee & Verleysen 2009, 2010; Venna & Kaski 2001, 2006; Mokbel et al. 2013; Espadoto et al. 2021; Hägele et al. 2023 (UAMDS); Moon et al. 2019 (DEMaP; its code was read); Amsaleg et al. 2021; Heulot, Aupetit, Fekete 2013.
- **Applied reliability:**
  - Xia, Lee, Li, *Nat. Commun.* 15:1753, 2024 (scDEED).
  - Johnson et al., *Patterns* 2022 (EMBEDR).
  - Sun, Ma, Zou, *Nat. Comput. Sci.* 3:86–100, 2023 (DynamicViz).
  - Liu, Ma, Zhong, *Nat. Commun.* 16:5037, 2025.
  - Zabel, Hennig, Nieselt, *Front. Bioinform.* 5:1719516, 2026.
  - Ma, Li, Hu, Yu, *Patterns* 7(8):101621, 2026 (NESS).
  - Sankaran et al., *Brief. Bioinform.* 2026 (the distortions package).
  - Chari & Pachter, *PLoS Comput. Biol.* 2023.
  - Heiser & Lau, *Cell Rep.* 2020.
  - Huang, Wang, Rudin, Browne, *Commun. Biol.* 5:719, 2022.
  - Ahlmann-Eltze & Huber, *Nat. Methods* 20:665–672, 2023.
  - Neufeld et al., arXiv:2207.00554 and arXiv:2307.12985.
  - Neufeld, Dharamshi, Gao, Witten, arXiv:2301.07276.
  - Dharamshi et al., *JASA* 120:511–523, 2025.
  - Kalinin et al., *Nat. Commun.* 16:5181, 2025.
  - *[not opened]*: Batson et al. 2019; Hivert et al. 2024.
  - A 2023 blog post by Pullin that the proposal cited could **not be located**, so it is
    dropped from the prior-art list.
- **Uncertain data and robustness:**
  - Agarwal et al., *ACM TALG* 2016.
  - Feng, Glavic, Kennedy, *PVLDB* 16, 2023.
  - Soliman & Ilyas, ICDE 2009.
  - Emrich et al., arXiv:2001.05581 (SIGMOD 2010 version *[not opened]*).
  - Bernecker et al., arXiv:1101.2613.
  - Karlaš et al., *PVLDB* 14(3), 2021.
  - Fan & Koutris, ICDT 2022.
  - Jia et al., AAAI 2022.
  - Ranzato, Shakeel, Zanella, CCS 2025.
  - Fassina, Ranzato, Zanella, *KAIS* 2024.
  - Cheng, Kalashnikov, Prabhakar, SIGMOD 2003.
  - Aronov, Iacono, Sheikhan, CCCG 2017.
  - Ram et al., NIPS 2009.
  - Vishwanath & Arias-Castro, arXiv:2510.22452.
