# agent-based-flocking-model

# Agent-Based Flocking Model (Vicsek with a Leader)

An agent-based simulation of collective motion (flocking): particles move
at constant speed and align their heading with nearby neighbors, plus
noise (the classic Vicsek model). This version adds a **leader** particle
whose heading counts more heavily to neighbors (`weight`) but who is less
influenced by them, with its own noise level (`beta`), to study whether a
leader helps or hurts overall flock alignment.

## Method

- 2D periodic Vicsek dynamics with a spatial hash grid for fast neighbor
  lookup.
- Followers: standard Vicsek update (average neighbor heading + noise).
- Leader: persists more in its own heading, influence from neighbors
  shrinks as `weight` grows, noise scaled by `beta`.
- Compared against a no-leader baseline using the global polar order
  parameter (fraction of particles moving in a common direction), averaged
  over multiple trials per (`weight`, `beta`).

## Results

Small `beta` (leader quieter than followers) raises alignment above
baseline; large `beta` (leader noisier) drags it below baseline — an
unreliable strong signal hurts more than having no leader. See
`presentation.pdf` for full plots.

## Contents

- `vicsek_leader.c` — simulation and parameter sweep
- `presentation.pdf` — slides with background and results

## Running

```bash
gcc -O2 -o vicsek_leader vicsek_leader.c -lm
./vicsek_leader
```

Produces `baseline_alignment.dat` and one `hope_<beta>.dat` file per beta
value.
