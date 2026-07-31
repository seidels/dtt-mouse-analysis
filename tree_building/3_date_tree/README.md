# Dating the rooted NJ tree

Two steps, run in order, per side (B1/B2) then merged:

1. `run_dating.sh` — BAT-correct negative NJ branches, LSD2-date, merge B1+B2
2. `run_constrained_dating.sh` — re-date under a literature min-age ladder, so
   early splits can't imply more coexisting lineages than the embryo had cells

Both call `run_lsd2.sh` for the actual LSD2 invocation (and the day-unit
newick it emits via `time_tree.R`); constrained re-dating just adds a
minimum-age datefile (`build_min_age_datefile.R`'s output) on top.

## 1. Unconstrained dating

```
ROOTED_TREE_B1=results/2-rooted-nj/B1/nj_rooted_ingroup.nwk \
ROOTED_TREE_B2=results/2-rooted-nj/B2/nj_rooted_ingroup.nwk \
STAGEDIR=results/3-dated-tree MERGE_TOKEN=ge7 \
  bash 3_date_tree/run_dating.sh
```

Takes each side's rooted tree from `2_root_tree/`, corrects negative branch
lengths (an NJ artefact -- LSD2 needs non-negative substitution counts), runs
LSD2 with root/tip calibration dates, then merges B1+B2 on a shared day-0
zygote root.

## 2. Constrained re-dating

Plain LSD2 reads "zero edits" as "~zero elapsed time" and can pack more early
splits into the first few days than the embryo actually had cells. This step
builds an LSD2 minimum-age constraint per node from a literature cell-count
ceiling (`processed_data/sample_matched_ceiling_sourced.csv`) and re-dates,
iterated 3x (each pass re-ranks nodes on the previous pass's output):

```
STAGEDIR=results/3-dated-tree MERGE_TOKEN=ge7 \
RANK_TREE_B1=results/3-dated-tree/B1/lsd2/nj99478/minB2h/time_tree.nwk \
RANK_TREE_B2=results/3-dated-tree/B2/lsd2/nj99478/minB2h/time_tree.nwk \
  bash 3_date_tree/run_constrained_dating.sh
```

This is the tree used downstream: `results/3-dated-tree/merged/merged_time_tree_ge7_attempt1_minage_sourced_minB2h_l0.01.nwk`.

## External prerequisites

- R + `ape`, `BAT`
- LSD2 source, built from `src/`:
  `git clone https://github.com/tothuhien/lsd2.git lsd2 && (cd lsd2/src && make)`
