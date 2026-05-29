# LAMP-Forge Snakemake DAG.
#
# Why both a Snakefile AND a single-process pipeline.py?
#   - pipeline.py is the friendly path: one command, runs locally, easy to
#     read top-to-bottom in the notebook walkthrough.
#   - The Snakefile is for cluster/HPC use, partial re-runs ("just rerun
#     specificity with a new off-target panel"), and for users who want
#     standard Snakemake reporting (DAG, benchmarks, conda envs per rule).
#
# Invocation:
#     snakemake --cores 4 --configfile config/example_config.yaml
#
# Or pass overrides on the command line:
#     snakemake --cores 4 --configfile config/example_config.yaml \
#         --config output_dir=results/myrun

from pathlib import Path

# ----------------------------------------------------------------------------
# Config plumbing
# ----------------------------------------------------------------------------

# Snakemake loads --configfile into the global `config` dict. We resolve the
# few paths we need locally; everything else is read directly from the YAML
# by the Python module functions.
OUTPUT_DIR = Path(config.get("output", {}).get("dir", "results"))
WORK_DIR = OUTPUT_DIR / "work"
CONFIG_PATH = config.get("_config_path", "config/example_config.yaml")

# Allow CLI override: --config _config_path=path/to/my.yaml
if "_config_path" not in config:
    # Snakemake passes the --configfile through implicitly; this file path
    # only matters to the lamp_forge.config loader. If the user wants to
    # change it, they should set _config_path in --config.
    pass


# ----------------------------------------------------------------------------
# Rules
# ----------------------------------------------------------------------------

rule all:
    input:
        OUTPUT_DIR / "lamp_forge_report.html",
        OUTPUT_DIR / "primer_sets.json",
        OUTPUT_DIR / "primer_sets.csv",
        OUTPUT_DIR / "run_manifest.json",


rule fetch:
    output:
        fasta=WORK_DIR / "sequences.fasta",
    params:
        config_path=CONFIG_PATH,
    log:
        OUTPUT_DIR / "logs" / "fetch.log",
    run:
        from lamp_forge import fetch as fetch_mod
        from lamp_forge.config import load_config
        cfg = load_config(params.config_path)
        records = fetch_mod.fetch_for_config(cfg)
        fetch_mod.write_fasta(records, Path(output.fasta))


rule align:
    input:
        fasta=WORK_DIR / "sequences.fasta",
    output:
        alignment=WORK_DIR / "alignment.fasta",
    log:
        OUTPUT_DIR / "logs" / "align.log",
    threads: 4
    shell:
        "mafft --auto --thread {threads} --quiet {input.fasta} > {output.alignment} 2> {log}"


rule conserve:
    input:
        alignment=WORK_DIR / "alignment.fasta",
    output:
        track=OUTPUT_DIR / "conservation.tsv",
        regions=WORK_DIR / "regions.json",
    params:
        config_path=CONFIG_PATH,
    log:
        OUTPUT_DIR / "logs" / "conserve.log",
    run:
        import json
        from dataclasses import asdict
        from lamp_forge import align as align_mod, conserve as conserve_mod
        from lamp_forge.config import load_config
        from lamp_forge.report import write_conservation_tsv

        cfg = load_config(params.config_path)
        msa = align_mod.load_alignment(Path(input.alignment))
        track = conserve_mod.compute_track(msa, window_size=cfg.window_size)
        regions = conserve_mod.find_conserved_regions(
            track,
            entropy_threshold=cfg.entropy_threshold,
            min_region_length=cfg.min_region_length,
        )
        write_conservation_tsv(track, Path(output.track))
        Path(output.regions).parent.mkdir(parents=True, exist_ok=True)
        Path(output.regions).write_text(json.dumps([asdict(r) for r in regions], indent=2))


rule design:
    input:
        regions=WORK_DIR / "regions.json",
    output:
        sets=WORK_DIR / "primer_sets_unscreened.json",
    params:
        config_path=CONFIG_PATH,
    log:
        OUTPUT_DIR / "logs" / "design.log",
    run:
        import json
        from dataclasses import asdict
        from lamp_forge import primer_design
        from lamp_forge.config import load_config
        from lamp_forge.types import ConservedRegion

        cfg = load_config(params.config_path)
        raw = json.loads(Path(input.regions).read_text())
        regions = [ConservedRegion(**r) for r in raw]
        sets = primer_design.design_all(regions, cfg)
        Path(output.sets).write_text(
            json.dumps([_set_to_dict(s) for s in sets], indent=2, default=str)
        )


def _set_to_dict(s):
    """Local helper: serialise a LampPrimerSet for the intermediate JSON."""
    from dataclasses import asdict
    d = asdict(s)
    return d


rule specificity:
    input:
        sets=WORK_DIR / "primer_sets_unscreened.json",
    output:
        sets=WORK_DIR / "primer_sets_screened.json",
        tsv=OUTPUT_DIR / "specificity.tsv",
    params:
        config_path=CONFIG_PATH,
    log:
        OUTPUT_DIR / "logs" / "specificity.log",
    run:
        import json
        from lamp_forge import specificity as spec_mod
        from lamp_forge.config import load_config
        from lamp_forge.types import LampPrimerSet, Primer, SpecificityHit

        cfg = load_config(params.config_path)
        raw = json.loads(Path(input.sets).read_text())
        sets = [_dict_to_set(d) for d in raw]
        off_targets = spec_mod.discover_off_targets(cfg.off_target_dir)
        if off_targets:
            db = Path(WORK_DIR) / "off_targets_db"
            spec_mod.build_blast_db(off_targets, db)
            spec_mod.screen_all(sets, db, cfg)
            spec_mod.write_specificity_tsv(sets, Path(output.tsv))
        else:
            Path(output.tsv).write_text("# No off-target genomes supplied.\n")
        Path(output.sets).write_text(
            json.dumps([_set_to_dict(s) for s in sets], indent=2, default=str)
        )


def _dict_to_set(d):
    """Reverse of _set_to_dict — rebuild a LampPrimerSet from nested dicts."""
    from lamp_forge.types import LampPrimerSet, Primer, SpecificityHit

    def mk_primer(pd):
        if pd is None:
            return None
        return Primer(**pd)

    d = dict(d)
    d["f3"] = mk_primer(d["f3"])
    d["b3"] = mk_primer(d["b3"])
    d["fip"] = mk_primer(d["fip"])
    d["bip"] = mk_primer(d["bip"])
    d["lf"] = mk_primer(d.get("lf"))
    d["lb"] = mk_primer(d.get("lb"))
    d["cross_reactivity_hits"] = [SpecificityHit(**h) for h in d.get("cross_reactivity_hits", [])]
    return LampPrimerSet(**d)


rule report:
    input:
        sets=WORK_DIR / "primer_sets_screened.json",
        alignment=WORK_DIR / "alignment.fasta",
        regions=WORK_DIR / "regions.json",
        fasta=WORK_DIR / "sequences.fasta",
    output:
        html=OUTPUT_DIR / "lamp_forge_report.html",
        json=OUTPUT_DIR / "primer_sets.json",
        csv=OUTPUT_DIR / "primer_sets.csv",
        manifest=OUTPUT_DIR / "run_manifest.json",
    params:
        config_path=CONFIG_PATH,
    log:
        OUTPUT_DIR / "logs" / "report.log",
    run:
        import json
        from lamp_forge import align as align_mod, conserve as conserve_mod, fetch as fetch_mod, report as report_mod, specificity as spec_mod
        from lamp_forge.config import load_config
        from lamp_forge.types import ConservedRegion

        cfg = load_config(params.config_path)
        raw_sets = json.loads(Path(input.sets).read_text())
        sets = [_dict_to_set(d) for d in raw_sets]
        raw_regions = json.loads(Path(input.regions).read_text())
        regions = [ConservedRegion(**r) for r in raw_regions]
        msa = align_mod.load_alignment(Path(input.alignment))
        track = conserve_mod.compute_track(msa, window_size=cfg.window_size)
        records = fetch_mod.load_fasta(Path(input.fasta))
        n_off_targets = len(spec_mod.discover_off_targets(cfg.off_target_dir))
        report_mod.write_json(sets, Path(output.json))
        report_mod.write_csv(sets, Path(output.csv))
        report_mod.render_html(
            sets, regions, track, cfg,
            n_input_seqs=len(records),
            n_off_targets=n_off_targets,
            path=Path(output.html),
        )
        report_mod.write_manifest(
            cfg, fetch_mod.manifest(records), len(regions), len(sets), Path(output.manifest)
        )
