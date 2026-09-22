#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "time.h"

#define PI 3.14159265359
#define t_max 1E4              // total simulation time (in units of dt)
#define v_magnitude 0.03       // fixed speed of every particle (Vicsek particles move at constant speed)
#define MAX_PARTICLES_PER_CELL 100

// ---------------------------------------------------------------------------
// Spatial hash grid
//
// Instead of checking every pair of particles (O(N^2)) to find neighbors
// within the interaction radius, we bin particles into a grid of square
// cells. To find a particle's neighbors we only need to look in its own
// cell plus the 8 cells around it, which is much faster for large N.
// ---------------------------------------------------------------------------

// One grid cell: a small fixed-size list of particle indices that currently
// live inside it.
typedef struct {
    int particle_ids[MAX_PARTICLES_PER_CELL];
    int count;
} Cell;

// Convert a particle's (x, y) position into (cell_x, cell_y) grid coordinates,
// wrapping around at the box edges since the simulation uses periodic
// boundary conditions.
void get_cell_indices(float x, float y, float cell_size, float L, int *cell_x, int *cell_y, int grid_size) {
    *cell_x = (int)(x / cell_size);
    *cell_y = (int)(y / cell_size);

    // Handle periodic boundaries
    if (*cell_x < 0) *cell_x += grid_size;
    else if (*cell_x >= grid_size) *cell_x -= grid_size;
    if (*cell_y < 0) *cell_y += grid_size;
    else if (*cell_y >= grid_size) *cell_y -= grid_size;
}

// Flatten 2D cell coordinates into a single 1D index so cells can be stored
// in a plain array.
int get_hash_index(int cell_x, int cell_y, int grid_size) {
    return cell_y * grid_size + cell_x;
}

// Re-bin every particle into the grid. Called once per timestep since
// particles move every step.
void build_spatial_hash(float x[], float y[], int N, float L, float cell_size, Cell *grid, int grid_size) {
    // Clear grid
    int total_cells = grid_size * grid_size;
    for (int i = 0; i < total_cells; i++) {
        grid[i].count = 0;
    }

    // Insert particles into grid
    for (int i = 0; i < N; i++) {
        int cell_x, cell_y;
        get_cell_indices(x[i], y[i], cell_size, L, &cell_x, &cell_y, grid_size);
        int hash_idx = get_hash_index(cell_x, cell_y, grid_size);

        if (grid[hash_idx].count < MAX_PARTICLES_PER_CELL) {
            grid[hash_idx].particle_ids[grid[hash_idx].count] = i;
            grid[hash_idx].count++;
        }
    }
}

// ---------------------------------------------------------------------------
// Core Vicsek update rule (with leader/follower extension)
//
// For a given particle, look at every neighbor within r_cut (using the
// spatial hash to only scan the 3x3 block of cells around it), average
// their heading directions, then add noise. Leaders and followers use two
// different update rules (see below).
// ---------------------------------------------------------------------------
float neighbor_list_spatial(float x[], float y[], float x_part, float y_part, float angles[],
                            float eta, int N, float L, Cell *grid, int grid_size, float cell_size,
                            int leader_list[], int self_index, float vx[], float vy[], int weight, float beta) {

    float r_cut = 1.0;
    float r_cut_sq = r_cut * r_cut;
    float S_x = 0.0, S_y = 0.0;   // running sum of neighbor direction vectors (x and y components)
    int neighbor = 0;              // (weighted) neighbor count, used to normalize S_x, S_y
    float mean_angle = 0.0;

    // Which grid cell does the current particle sit in?
    int center_cell_x, center_cell_y;
    get_cell_indices(x_part, y_part, cell_size, L, &center_cell_x, &center_cell_y, grid_size);

    // Scan the current cell and its 8 neighboring cells (3x3 block)
    for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
            int cell_x = center_cell_x + dx;
            int cell_y = center_cell_y + dy;

            // Handle periodic boundaries
            if (cell_x < 0) cell_x += grid_size;
            else if (cell_x >= grid_size) cell_x -= grid_size;
            if (cell_y < 0) cell_y += grid_size;
            else if (cell_y >= grid_size) cell_y -= grid_size;

            int hash_idx = get_hash_index(cell_x, cell_y, grid_size);

            // Check all particles in this cell
            for (int p = 0; p < grid[hash_idx].count; p++) {
                int j = grid[hash_idx].particle_ids[p];

                // If we (self_index) are a leader, exclude ourself from our
                // own neighbor sum -- the leader's update rule (below) adds
                // its own previous heading back in separately, weighted
                // differently from the neighbor contribution.
                if (leader_list[self_index] == 1 && j == self_index) {
                    continue;
                }

                float delta_x = x[j] - x_part;
                float delta_y = y[j] - y_part;

                // Minimum-image convention for periodic boundary conditions
                if (delta_x >  L/2) delta_x -= L;
                else if (delta_x < -L/2) delta_x += L;
                if (delta_y >  L/2) delta_y -= L;
                else if (delta_y < -L/2) delta_y += L;

                float r_sq = delta_x * delta_x + delta_y * delta_y;

                if (r_sq < r_cut_sq) {
                    // A follower neighbor contributes with weight 1.
                    // A leader neighbor contributes with weight `weight`,
                    // i.e. leaders count as `weight` ordinary particles when
                    // influencing others -- this is what makes them "loud".
                    if (leader_list[j] == 0) {
                        neighbor++;
                        S_x += cosf(angles[j]);
                        S_y += sinf(angles[j]);
                    }
                    else {
                        neighbor += weight;
                        S_x += (float)weight * cosf(angles[j]);
                        S_y += (float)weight * sinf(angles[j]);
                    }
                }
            }
        }
    }

    if (neighbor == 0) {
        // No neighbors found: keep current heading (no direction to average
        // toward).
        S_x = cosf(angles[self_index]);
        S_y = sinf(angles[self_index]);
    } else {
        S_x /= (float)neighbor;
        S_y /= (float)neighbor;
    }

    float u = (float)rand() / RAND_MAX;  // uniform random number in [0, 1)

    if (leader_list[self_index] == 0) {
        // --- FOLLOWER UPDATE (standard Vicsek rule) ---
        // New heading = average neighbor direction + uniform noise in
        // [-eta/2, eta/2].
        mean_angle = atan2f(S_y, S_x) +  eta * (u - 0.5f);
    }
    else {
        // --- LEADER UPDATE ---
        // The leader partially ignores its neighbors and partially keeps
        // moving in its own previous direction, then adds its own
        // (possibly different) noise strength beta*eta.

        // Social direction = average heading of neighbors only (self was
        // excluded above)
        float theta_soc = atan2f(S_y, S_x);
        float nx_soc = cosf(theta_soc);
        float ny_soc = sinf(theta_soc);

        // Leader's own current heading (self term)
        float nx_self = vx[self_index] / v_magnitude;
        float ny_self = vy[self_index] / v_magnitude;

        float w_nei  = 1.0 / (float) weight;   // neighbors' influence on the leader shrinks as weight grows
        float w_self = 1.0f;                   // leader always weights its own previous heading fully

        // Combine neighbor influence and self-persistence
        float v_x = w_nei * nx_soc + w_self * nx_self;
        float v_y = w_nei * ny_soc + w_self * ny_self;

        // Normalize back to a unit direction vector
        float norm = sqrtf(v_x*v_x + v_y*v_y);
        if (norm > 0.0f) {
            v_x /= norm;
            v_y /= norm;
        }

        // Leader's noise strength is scaled by beta relative to followers'
        // eta -- beta < 1 makes the leader "quieter"/more decisive, beta > 1
        // makes the leader noisier/less reliable as a signal.
        mean_angle = atan2f(v_y, v_x) + beta * eta * (u - 0.5f);
    }

    return mean_angle;
}

// ---------------------------------------------------------------------------
// Full simulation run: Vicsek model with exactly one leader particle.
// Returns the time-averaged global polar order parameter Va (see Alignment
// section of the presentation), i.e. how aligned the whole flock is,
// averaged over the second half of the run (after a transient warm-up
// period of t > 2000).
// ---------------------------------------------------------------------------
float viscek(float eta, float rho, float weight, float beta) {
    int N = 40;                    // number of particles
    int leader = 1;                // number of leaders to designate
    int is_leader[N];
    int count_leader = 0;
    float L = sqrtf(N / rho);      // box size set by fixed density rho = N / L^2
    float va_total = 0.0;
    float total = 0.0;
    float vx_leader=0.0,vy_leader=0.0;

    // Spatial hash parameters
    float r_cut = 1.0;
    float cell_size = r_cut;
    int grid_size = (int)(L / cell_size);
    Cell *grid = (Cell*)malloc(grid_size * grid_size * sizeof(Cell));

    float x_pos[N], y_pos[N];
    float vx[N], vy[N], direction[N];
    float new_direction[N];
    float dt = 1.0;
    float t = 0.0;

    float nLx=0.0, nLy=0.0;   // (currently unused leftover from the commented-out
                              // leader-tracking alignment metric below)

    // Random initial positions and headings, uniformly distributed in the box
    for (int i = 0; i < N; i++) {
        x_pos[i] = ((float)rand() / RAND_MAX) * L;
        y_pos[i] = ((float)rand() / RAND_MAX) * L;
        direction[i] = ((float)rand() / RAND_MAX) * 2.0 * PI;
        is_leader[i] = 0;
        vx[i] = v_magnitude * cosf(direction[i]);
        vy[i] = v_magnitude * sinf(direction[i]);
    }

    // Randomly designate `leader` particles as leaders
    while (count_leader < leader){
        float p = (float)rand() / RAND_MAX;
        int i = (int) (p * N);
        if (is_leader[i] == 0) {
            is_leader[i] = 1;
            count_leader ++;
        }
    }

    // ---- Time evolution loop ----
    while (t < t_max) {
        float A_sum = 0.0;
        float sum_x = 0.0, sum_y = 0.0;
        int n_f = 0;

        // Rebuild the spatial hash since particles have moved since last step
        build_spatial_hash(x_pos, y_pos, N, L, cell_size, grid, grid_size);

        // Step 1: compute everyone's new heading based on current positions
        // (synchronous update -- all particles look at the *same* old state)
        for (int i = 0; i < N; i++) {
            new_direction[i] = neighbor_list_spatial(x_pos, y_pos,x_pos[i], y_pos[i],direction,eta,
                                            N, L,grid, grid_size, cell_size,is_leader,i,vx,vy, weight,beta);
        }

        // Step 2: apply new headings, advance positions, wrap around box (PBC)
        for (int i = 0; i < N; i++) {
            direction[i] = new_direction[i];

            vx[i] = v_magnitude * cosf(direction[i]);
            vy[i] = v_magnitude * sinf(direction[i]);

            x_pos[i] += vx[i] * dt;
            y_pos[i] += vy[i] * dt;

            if (x_pos[i] < 0) x_pos[i] += L;
            else if (x_pos[i] > L) x_pos[i] -= L;
            if (y_pos[i] < 0) y_pos[i] += L;
            else if (y_pos[i] > L) y_pos[i] -= L;

            // After the transient warm-up period, accumulate velocity
            // components for the global order parameter Va, and keep track
            // of the leader's current heading (nLx, nLy) -- computed but not
            // currently used, see commented-out block below.
            if (t > 2000){
                if (is_leader[i] == 1){
                    vx_leader = vx[i];
                    vy_leader = vy[i];
                    nLx = vx_leader / v_magnitude;
                    nLy = vy_leader / v_magnitude;
                    sum_x += vx[i];
                    sum_y += vy[i];
                }
                else{
                    sum_x += vx[i];
                    sum_y += vy[i];
                }
            }
        }

        // NOTE: this commented-out block is the "Alignment with Leader"
        // metric A(t) = <cos(theta_k - theta_L)> shown in the presentation
        // slides (dot product of each follower's heading with the leader's
        // heading, averaged over followers). It is currently disabled --
        // the active code below computes the standard *global* Vicsek order
        // parameter Va instead (average heading of the whole flock,
        // leader included), which is a different quantity from the
        // slide plots labeled "Alignment with Leader".
        // if (t > 2000){
        //     for(int k=0; k<N; k++){
        //         if(is_leader[k] == 0){
        //             float dot = vx[k] * nLx + vy[k] * nLy;
        //             float alignment = dot / v_magnitude;
        //             A_sum += alignment;
        //             n_f++;
        //         }
        //     }
        //     float A_global = A_sum / (float)n_f;
        //     total += A_global;
        // }

        // Global polar order parameter: Va = |sum of unit velocities| / (N * v_magnitude).
        // Va = 1 means perfect alignment (everyone moving the same direction),
        // Va = 0 means random/disordered motion.
        if (t > 2000) {
            float Va = sqrtf(sum_x*sum_x + sum_y*sum_y) / (v_magnitude * N);
            total += Va;
        }
        t += dt;
    }

    free(grid);

    // Time-average over the post-warm-up window (t = 2000 .. t_max)
    return total / (t_max - 2000.0);
}

// ---------------------------------------------------------------------------
// Same simulation but with NO leaders at all -- pure classic Vicsek model.
// Used to compute the "Baseline (No Leader)" reference line in the plots.
// ---------------------------------------------------------------------------
float viscek_no_leader(float eta, float rho) {
    int N = 40;
    int is_leader[N];
    float L = sqrtf(N / rho);
    float total = 0.0;

    // Spatial hash parameters
    float r_cut = 1.0;
    float cell_size = r_cut;
    int grid_size = (int)(L / cell_size);
    Cell *grid = (Cell*)malloc(grid_size * grid_size * sizeof(Cell));

    float x_pos[N], y_pos[N];
    float vx[N], vy[N], direction[N];
    float new_direction[N];
    float dt = 1.0;
    float t = 0.0;

    // Initial condition - NO LEADERS
    for (int i = 0; i < N; i++) {
        x_pos[i] = ((float)rand() / RAND_MAX) * L;
        y_pos[i] = ((float)rand() / RAND_MAX) * L;
        direction[i] = ((float)rand() / RAND_MAX) * 2.0 * PI;
        is_leader[i] = 0;  // everyone is a follower
        vx[i] = v_magnitude * cosf(direction[i]);
        vy[i] = v_magnitude * sinf(direction[i]);
    }

    // Time evolution (identical structure to viscek(), but weight=0 and
    // beta=1.0 are passed in as placeholders since there are no leaders to
    // apply them to)
    while (t < t_max) {
        build_spatial_hash(x_pos, y_pos, N, L, cell_size, grid, grid_size);

        for (int i = 0; i < N; i++) {
            new_direction[i] = neighbor_list_spatial(x_pos, y_pos, x_pos[i], y_pos[i],
                                direction, eta, N, L, grid, grid_size, cell_size,
                                is_leader, i, vx, vy, 0, 1.0);  // weight=0, beta=1.0 (unused, no leaders)
        }

        for (int i = 0; i < N; i++) {
            direction[i] = new_direction[i];
            vx[i] = v_magnitude * cosf(direction[i]);
            vy[i] = v_magnitude * sinf(direction[i]);
            x_pos[i] += vx[i] * dt;
            y_pos[i] += vy[i] * dt;

            if (x_pos[i] < 0) x_pos[i] += L;
            else if (x_pos[i] > L) x_pos[i] -= L;
            if (y_pos[i] < 0) y_pos[i] += L;
            else if (y_pos[i] > L) y_pos[i] -= L;
        }

        if (t > 2000){
            float sum_vx = 0.0, sum_vy = 0.0;
            for(int k = 0; k < N; k++){
                sum_vx += vx[k];
                sum_vy += vy[k];
            }
            float Va = sqrtf(sum_vx*sum_vx + sum_vy*sum_vy) / (v_magnitude * N);
            total += Va;
        }

        t += dt;
    }

    free(grid);
    return total / (t_max - 2000.0);
}


int main() {
    float etaa = 2.0;   // fixed follower noise strength for all sweeps below
    float rho = 4.0;    // fixed density
    int n_trials = 5;   // independent repeats per (weight, beta) pair, for averaging + error bars

    // beta = leader's own noise multiplier (relative to eta). Swept over a
    // range from "quieter than followers" (0.7) to "much noisier than
    // followers" (3.5) -- corresponds to the "Below Baseline" / combined
    // sweep plots in the presentation.
    float beta_values[] = {0.7,0.8,0.9,1.0,1.5,2,2.5,3,3.5};
    int n_beta = sizeof(beta_values) / sizeof(beta_values[0]);

    // weight = how strongly a leader's heading counts relative to a single
    // follower when neighbors average directions. Swept from 1 (leader
    // counts as an ordinary particle) up to 50 (leader dominates any local
    // average it's part of).
    int weights[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,
    32,34,36,38,40,42,44,46,48,50};
    int n_weights = sizeof(weights) / sizeof(weights[0]);

    printf("Starting beta parameter sweep...\n");
    printf("eta = %.2f, rho = %.2f, trials per point = %d\n\n", etaa, rho, n_trials);

    // --- Baseline: no leader at all ---
    printf("Computing baseline (no leader)...\n");
    float baseline_sum = 0.0;
    int baseline_trials = 20;
    for (int trial = 0; trial < baseline_trials; trial++) {
        srand(5678 + trial);  // different seed range from the leader runs below
        baseline_sum += viscek_no_leader(etaa, rho);
    }
    float baseline = baseline_sum / baseline_trials;
    printf("Baseline alignment (no leader) = %.6f\n\n", baseline);

    // Save baseline to file for later comparison in plots
    FILE *baseline_file = fopen("baseline_alignment.dat", "w");
    fprintf(baseline_file, "%.6f\n", baseline);
    fclose(baseline_file);

    // --- Main sweep: for each beta, scan over all weight values ---
    for (int b = 0; b < n_beta; b++) {
        float beta = beta_values[b];

        // One output file per beta value, e.g. hope_0.700.dat, hope_1.500.dat, ...
        char filename[100];
        sprintf(filename, "hope_%.3f.dat", beta);
        FILE *outData = fopen(filename, "w");

        if (outData == NULL) {
            printf("Error: Could not open file %s\n", filename);
            continue;
        }

        fprintf(outData, "# Beta = %.3f, Eta = %.2f, Rho = %.2f\n", beta, etaa, rho);
        fprintf(outData, "# Weight  Alignment  StdError\n");

        printf("Processing beta = %.3f...\n", beta);

        for (int w = 0; w < n_weights; w++) {
            int weight = weights[w];
            float sum_align = 0.0;
            float sum_align_sq = 0.0;

            // Run n_trials independent simulations at this (weight, beta)
            // and compute mean + standard error across them
            for (int trial = 0; trial < n_trials; trial++) {
                srand(1234 + trial);
                float align = viscek(etaa, rho, weight, beta);
                sum_align += align;
                sum_align_sq += align * align;
            }

            float mean_align = sum_align / n_trials;
            float variance = (sum_align_sq / n_trials) - (mean_align * mean_align);
            float std_error = sqrtf(variance / n_trials);

            fprintf(outData, "%d  %.6f  %.6f\n", weight, mean_align, std_error);

            if (weight % 5 == 0) {
                printf("  Weight %2d: Alignment = %.4f +/- %.4f\n", weight, mean_align, std_error);
            }
        }

        fclose(outData);
        printf("Completed beta = %.3f -> %s\n\n", beta, filename);
    }

    printf("All simulations complete!\n");
    printf("Output files: hope_*.dat (one per beta value), baseline_alignment.dat\n");
    return 0;
}
