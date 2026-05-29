## One-session benchmark: Mooncake (reverse-mode) vs ForwardDiff on the
## per-case latent-joint bdbv_model. Each backend is warmed up (to pay
## the one-off AD compilation) and then timed on a representative fit so
## the reported sampling time reflects steady-state, not compilation.

using BdbvLinelist
using ADTypes: AutoForwardDiff, AutoMooncake

const WARM_SAMPLES = 25
const WARM_CHAINS  = 1
const BENCH_SAMPLES = 1000
const BENCH_CHAINS  = 4

ll = load_linelist(BdbvLinelist.LINELIST_PATH)
d  = build_data(ll)
mk() = bdbv_model(d; family = :gamma)

function run_backend(name, adtype)
    @info "[$name] warm-up ($WARM_SAMPLES samples, $WARM_CHAINS chain) — pays AD compilation"
    tc = @elapsed sample_fit(mk(); samples = WARM_SAMPLES, chains = WARM_CHAINS,
                             adtype = adtype, progress = false)
    @info "[$name] timed fit ($BENCH_SAMPLES samples, $BENCH_CHAINS chains)"
    ts = @elapsed chn = sample_fit(mk(); samples = BENCH_SAMPLES, chains = BENCH_CHAINS,
                                   adtype = adtype, progress = false)
    diag = diagnostics(chn)
    @info "[$name] result" compile_plus_warmup_s=round(tc, digits=1) sampling_s=round(ts, digits=1) max_rhat=round(diag.rhat, digits=3) min_ess=round(Int, diag.ess) ndiv=diag.ndiv
    return (; name, tc, ts, diag)
end

fd = run_backend("ForwardDiff", AutoForwardDiff())
mc = run_backend("Mooncake",    AutoMooncake(; config = nothing))

println("\n=== AD backend comparison (gamma, $BENCH_SAMPLES×$BENCH_CHAINS) ===")
println(rpad("backend", 14), rpad("warmup+compile(s)", 20), rpad("sampling(s)", 14),
        rpad("max_rhat", 10), rpad("min_ess", 10), "ndiv")
for r in (fd, mc)
    println(rpad(r.name, 14), rpad(round(r.tc, digits=1), 20),
            rpad(round(r.ts, digits=1), 14), rpad(round(r.diag.rhat, digits=3), 10),
            rpad(round(Int, r.diag.ess), 10), r.diag.ndiv)
end
println("\nspeedup (sampling, ForwardDiff/Mooncake): ", round(fd.ts / mc.ts, digits=2), "×")
