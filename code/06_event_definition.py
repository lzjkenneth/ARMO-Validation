"""Step 6 - the wetting-event definition, checked three ways (answers "how do you decide what a wetting event is?").

1. Example events: hourly 5 cm soil moisture and rain for 6 events picked at fixed quantiles of window length, with the
   pre-storm baseline, the peak, the 50 % line (the report's rule), the 1/e line and the 0.03 m3/m3 line, and the fitted
   exponential dry-down. Shows what "50 % of the rise" looks like on real data.
2. Window length under the three rules (half / e-fold / absolute 0.03), with the share of windows that are cut short by the
   next storm (censored) - a censored window is ended by new rain, not by drying, so it is a lower bound.
3. Fitted dry-down time constant tau, compared with the half-amplitude window: for a pure exponential the anomaly halves
   after tau * ln 2 = 0.69 tau.
Outputs: figures/fig_event_examples.png, figures/fig_event_windows.png, data/processed/event_definition_summary.csv"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROC = os.path.join(ROOT, "data", "processed")
FIG = os.path.join(ROOT, "figures"); os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})

ev = pd.read_parquet(os.path.join(PROC, "events.parquet"))
hourly = pd.read_parquet(os.path.join(PROC, "uscrn_hourly.parquet"))
st = pd.read_csv(os.path.join(PROC, "station_seasons.csv")).drop_duplicates("wban").set_index("wban")

# ---- 2. summary table
rows = []
for name, lab in (("half", "50 % of the rise (report)"), ("half_p3", "50 %, must stay below 3 h"), ("efold", "1/e of the rise (37 %)"),
                  ("abs", "anomaly >= 0.03 m3/m3")):
    W = ev[f"W_{name}_h"].dropna() / 24.0; c = ev[f"cens_{name}"].dropna().astype(bool)
    rows.append(dict(rule=lab, events=len(W), median_d=W.median(), q25_d=W.quantile(0.25), q75_d=W.quantile(0.75),
                     mean_d=W.mean(), share_lt_1d=(W < 1).mean(), share_lt_2d=(W < 2).mean(), share_censored=c.mean()))
summ = pd.DataFrame(rows)
tau = ev.tau_pre_h.dropna() / 24.0; tau_free = ev.tau_free_h.dropna() / 24.0
fit_ok = ev.tau_pre_h.notna() & ~ev.cens_half.astype(bool)
half_after_peak = (ev.we_half - ev.t_peak) / pd.Timedelta(hours=1)
ratio = (half_after_peak[fit_ok] / (ev.loc[fit_ok, "tau_pre_h"] * np.log(2)))
extra = pd.DataFrame([dict(rule="dry-down tau, floor = pre-storm level", events=len(tau), median_d=tau.median(), q25_d=tau.quantile(0.25),
                           q75_d=tau.quantile(0.75), mean_d=tau.mean(), share_lt_1d=(tau < 1).mean(), share_lt_2d=(tau < 2).mean(),
                           share_censored=np.nan),
                      dict(rule="dry-down tau, free floor (literature style)", events=len(tau_free), median_d=tau_free.median(),
                           q25_d=tau_free.quantile(0.25), q75_d=tau_free.quantile(0.75), mean_d=tau_free.mean(),
                           share_lt_1d=(tau_free < 1).mean(), share_lt_2d=(tau_free < 2).mean(), share_censored=np.nan),
                      dict(rule="peak-to-50 % time / (tau_pre ln2), uncensored", events=len(ratio), median_d=ratio.median(),
                           q25_d=ratio.quantile(0.25), q75_d=ratio.quantile(0.75), mean_d=ratio.mean(), share_lt_1d=np.nan,
                           share_lt_2d=np.nan, share_censored=np.nan)])
summ = pd.concat([summ, extra], ignore_index=True)
summ.to_csv(os.path.join(PROC, "event_definition_summary.csv"), index=False)
print(summ.round(3).to_string(index=False))
print("rise time (onset to peak) median h:", ev.rise_h.median(), "| amplitude median m3/m3:", round(ev.amp.median(), 3),
      "| rain median mm:", round(ev.rain_mm.median(), 1))

# ---- 1. example events at fixed quantiles of the half window
ev_s = ev[~ev.cens_half.astype(bool) & ev.tau_pre_h.notna()].sort_values("W_half_h").reset_index(drop=True)
picks = [ev_s.iloc[int(q * (len(ev_s) - 1))] for q in (0.1, 0.3, 0.5, 0.6, 0.8, 0.95)]
fig, axs = plt.subplots(2, 3, figsize=(11, 6.2))
for ax, e in zip(axs.flat, picks):
    t0 = e.t_start - pd.Timedelta(hours=12); t1 = max(e.we_half, e.we_efold if pd.notna(e.we_efold) else e.we_half) + pd.Timedelta(hours=24)
    g = hourly[(hourly.wban == e.wban) & (hourly.time >= t0) & (hourly.time <= t1)]
    th = (g.time - e.t_start) / pd.Timedelta(hours=1)
    ax.plot(th, g.sm5, color="#1f4e79", lw=1.4, label="5 cm soil moisture (hourly)")
    ax2 = ax.twinx(); ax2.bar(th, g.rain.fillna(0), width=1.0, color="#7fb3d5", alpha=0.6); ax2.set_ylim(0, max(4 * g.rain.max(), 1))
    ax2.invert_yaxis(); ax2.set_yticks([]); ax2.spines["right"].set_visible(False); ax2.spines["top"].set_visible(False)
    hp = (e.t_peak - e.t_start) / pd.Timedelta(hours=1)
    for lev, col, lab in ((e.pre, "0.4", "pre-storm level"), (e.pre + 0.5 * e.amp, "#c0392b", "50 % of rise (report rule)"),
                          (e.pre + e.amp / np.e, "#e67e22", "1/e of rise"), (e.pre + 0.03, "#27ae60", "+0.03 m3/m3")):
        ax.axhline(lev, color=col, lw=0.9, ls="--" if col != "0.4" else ":", label=lab)
    ws, we = [(x - e.t_start) / pd.Timedelta(hours=1) for x in (e.ws_half, e.we_half)]
    ax.axvspan(ws, we, color="#c0392b", alpha=0.08, label="wetting window (50 % rule)")
    tt = np.linspace(hp, hp + min(7 * 24, th.max() - hp), 200)
    ax.plot(tt, e.pre + e.amp * np.exp(-(tt - hp) / e.tau_pre_h), color="#8e44ad", lw=0.9, alpha=0.8,
            label="fitted dry-down (floor = pre-storm level)")
    ax.set_title(f"{st.loc[e.wban, 'name'].replace('_', ' ')} {e.t_start:%d %b %Y}\nrain {e.rain_mm:.0f} mm, rise {e.amp:.3f}, "
                 f"window {e.W_half_h / 24:.1f} d, tau {e.tau_pre_h / 24:.1f} d", fontsize=8)
    ax.set_xlabel("hours after storm start"); ax.set_ylabel("m3/m3")
    ax.set_ylim(e.pre - 0.02, e.pre + e.amp + 0.03)
handles, labels = axs.flat[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, fontsize=8)
fig.suptitle("Six real wetting events (USCRN, picked at the 10th-95th percentile of window length): where the 50 % rule cuts", fontsize=10)
fig.tight_layout(rect=(0, 0.07, 1, 0.97)); fig.savefig(os.path.join(FIG, "fig_event_examples.png")); plt.close(fig)

# ---- 2b. window distributions
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
bins = np.linspace(0, 10, 41)
for name, col, lab in (("half", "#c0392b", "50 % of rise (report)"), ("half_p3", "#922b21", "50 %, 3 h below to end"),
                       ("efold", "#e67e22", "1/e of rise"), ("abs", "#27ae60", "anomaly >= 0.03")):
    W = np.sort(ev[f"W_{name}_h"].dropna() / 24.0)
    ax[0].plot(W, np.arange(1, len(W) + 1) / len(W), color=col, label=f"{lab}: median {np.median(W):.1f} d")
ax[0].axvline(1.4, color="k", lw=0.8, ls=":"); ax[0].text(1.45, 0.05, "1.4 d", fontsize=8)
ax[0].set_xlim(0, 10); ax[0].set_xlabel("wetting window (days)"); ax[0].set_ylabel("share of events with a shorter window")
ax[0].legend(frameon=False, fontsize=8); ax[0].set_title(f"Window length under three rules ({len(ev):,} events, 120 stations)", fontsize=9)
ax[1].hist(tau.clip(upper=15), bins=np.linspace(0, 15, 31), color="#8e44ad", alpha=0.75)
ax[1].set_xlabel("dry-down time constant tau, floor = pre-storm level (days, clipped at 15)"); ax[1].set_ylabel("events")
ax[1].set_title(f"Dry-down tau: median {tau.median():.1f} d (IQR {tau.quantile(.25):.1f}-{tau.quantile(.75):.1f})", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_event_windows.png")); plt.close(fig)
print("figures written")
