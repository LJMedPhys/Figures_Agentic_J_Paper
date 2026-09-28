# Figure captions

Ready to paste. Each figure also exists as an `.svg` beside the `.png`.

---

## figure_task1_overlays

**Task 1: the same six nuclei, segmented three ways.**

Top row, all on the identical DAPI crop of position 40 (608 × 608 px at
x = 565, y = 952; 0.322 µm per pixel; scale bar 20 µm): the laboratory's six
hand-drawn ROIs, the six objects the Advanced configuration segmented for
itself, and the seven the Fast configuration segmented. Outline colour is that
source's own TUNEL call — red for TUNEL-positive, blue for TUNEL-negative — and
numbers are each source's own object IDs, which do not correspond between
panels. Fast mode split the laboratory's largest nucleus (Cell 3, a 37 × 30 px
ROI) into two objects; its own merge audit tested all 21 object pairs and merged
none.

Bottom, the background-corrected TUNEL signal of each nucleus, matched between
sources by centroid (all matches within 4 px). A red ring marks every point the
source that produced it called TUNEL-positive. The three traces track one
another, so the configurations agree with the laboratory on the *measurement*;
Cell 6 carries a ring for the laboratory alone, so they disagree only on where
to place the cutoff. No threshold derived from these six values recovers the
laboratory's second positive — applied to the laboratory's own numbers, a
largest-gap split also returns Cell 3 alone — so the laboratory's call for Cell 6
rests on evidence outside the corrected signal.

---

## figure_tasks23_dotplots

**Tasks 2 and 3: the published result and what each configuration reproduced.**

Proliferation index (10 − corrected red / corrected green) of single parasites,
grouped by whether the host macrophage was TUNEL-positive, plotted as in the
source publication. Each column is one source; within a column the three
positions are fixed (TUNEL+, TUNEL−, Unassigned) so that a configuration missing
a group leaves a visible gap rather than shifting its remaining groups, and `n`
is given for every position including the empty ones. Horizontal bars are
medians; jitter is horizontal only.

The published columns are JCI Insight 169020 figure 4C (right) and figure 4D,
read from the article's source data. They match this project's independent parse
of the laboratory's own workbooks to the last digit — 15 values for task 2 and
all 828 for task 3 — so the published column is simultaneously the figure in the
paper and the reference every other column is scored against.

**Task 2** (y axis 0–10): the agent-plus-human repeats reproduce the published
distribution, median for median. The Advanced configuration found 3 of the 15
parasites and produced no TUNEL-positive parasite at all; the Fast configuration
found none of the 15, reported 22 objects elsewhere in the field, and could
assign none of them to a host, so both groups are empty and the comparison
cannot be made. One point of each falls outside the axis.

**Task 3** (y axis 4–10, stars from a two-sided Mann–Whitney U test on the
parasite-level indices, computed identically for every column): both
agent-plus-human repeats recover the direction and the medians of the published
result. Advanced and Fast produced no per-parasite table at all, having been
stopped by the run deadline. The fixed axis clips a comparable minority from
every column that has data — 32 of 828 published, 29 of 583 and 32 of 729 for
the repeats.

---

## figure_workflow_task1

**Task 1: how each configuration reached its answer.**

One column per configuration, read top to bottom; each block is a step and the
grey line beneath it names the tool or method that performed it. The badge under
each configuration is time to result.

Steps are taken from each run's own record — its state ledger for the plan and
the decisions, the quality-control report written beside each output for the
rule actually applied, and the scripts the run saved — not reconstructed after
the fact. The manual column is the laboratory's description of how the published
figure was produced.

The two agent-plus-human repeats ran the same channel-2 pipeline on this task
and produced the same six corrected values to thirteen significant digits, so
one is shown as an example repeat and carries its 32 minutes; the other
additionally measured a second channel first and published that instead, which
differs in what it adopted rather than in the steps drawn here. Advanced and Fast
are shown against the same 85-minute budget used for every automatic run.

Note the divergence at the segmentation step: Advanced tried StarDist, rejected
the result as ambiguous, and fell back to a classical difference-of-Gaussians,
Otsu and watershed pipeline that recovered exactly the laboratory's six nuclei;
Fast kept StarDist and returned seven.
