"""Summarise and visualise BAMITA (``MultiwayImputation``) evaluation results.

The BAMITA counterpart of ``report_deepmicrogen``. All the reporting machinery
(loading, summary table, box/violin/degradation, mechanism/condition overlays,
heatmap, facet grid, diversity scatter) is shared -- this script only points the
shared ``run`` at the ``MultiwayImputation`` results, which it reads from
``results/MultiwayImputation/`` and plots under ``results/MultiwayImputation/plot/``.

BAMITA is untuned and has no NN hyperparameter columns, so the per-parameter-set
grouping collapses to one set per (dataset, condition) automatically.
"""
from report_deepmicrogen import run

if __name__ == "__main__":
    run("MultiwayImputation")
