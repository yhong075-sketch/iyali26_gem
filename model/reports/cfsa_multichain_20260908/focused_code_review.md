# Focused diagnostic code review

Reviewed at 2026-09-08T20:07:35.384789+00:00 against official ArviZ v0.22.0 snapshots in `vendor/`. Scope: `_rhat`, `_ess`, rank/fold/split/quantile wrapper, invalid values and lag-1 shape. No model solves or model sampling. This is a focused implementation check, not a full research-source or sampling audit.

**Result:** no core formula defect found. Fifty finite Rhat/bulk-ESS/tail-ESS comparisons matched reference functions exactly (maximum absolute difference 0): three chains, two variables, normal or integer-tied synthetic arrays, draws 8, 9, 100, 101 and 1000; fixed RNG seed 983. Reference functions were extracted directly from saved official source with optional Numba dispatch replaced by NumPy and valid-input guards bypassed for these valid synthetic inputs. The implementation's deterministic IID, shifted-chain, AR(1), constant and nonfinite checks also passed in this review. Lag-1 output shape was `(3, 2)` for all ten comparison arrays. An initial reference harness import failed before comparisons because `packaging.version` had not been explicitly imported; the corrected harness passed.

Static review confirms rank averaging, Blom transform, pooled split median for folded Rhat, original-array 5%/95% thresholds before indicator splitting, Geyer positive/monotone pairs, final positive even-lag contribution and tau lower bound. Lag-1 is per-original-chain Pearson correlation of adjacent draws, not the biased FFT autocorrelation used inside ESS. The wrapper rejects fewer than two chains/eight draws and preserves per-variable NaN results for nonfinite or globally fixed inputs.

**Documentation finding:** `np.maximum` propagates an undefined folded Rhat, while ArviZ's Python `max` can retain a finite bulk value. A balanced binary alternating array of shape `(3,100,1)` yields NaN here versus ArviZ Rhat 0.9899494936611666 because every folded value is identical. This is a conservative fail-closed departure; document it with the already stated constant-ESS departure. `np.minimum` likewise propagates any undefined tail-indicator ESS. The tied comparisons intentionally compared only finite values on both sides; this does not claim all tied or degenerate inputs equal ArviZ.

Coverage: 8 focused properties checked; 7 supported without qualification, 1 supported with the conservative-departure documentation requirement; broad source audit and actual retained-chain diagnostics unchecked. Core functions are private and reviewed through `diagnose`; arbitrary alternate `_autocov(axis=...)` use was not reviewed.

Code SHA256: `b6c3569f2e165f08ead700740e80faff17a77b8680ed62fa888fe393471eb616`.
Reference diagnostics SHA256: `08e4f2d9ce513808ab0802bf50da3b9109e080def599a0513f0cb9e1975a6ca1`.
Reference stats utilities SHA256: `5c6f082cf56cf096e1fefee3ccb7f03249b61eb693da3dbe1f17bf9046252e34`.

Official sources: [diagnostics.py](https://raw.githubusercontent.com/arviz-devs/arviz/v0.22.0/arviz/stats/diagnostics.py), [stats_utils.py](https://raw.githubusercontent.com/arviz-devs/arviz/v0.22.0/arviz/stats/stats_utils.py). Both independently opened during this review. No implementation edits made.
