/* Inside Rev.7's output branch, after the physical update. Read-only fields.
   Frame records are little-endian {int32 cell; float32 depth, v1, v3}. */
static int stable_intervals = 0, previous_area = -1;
static double first_arrival = -1, window_start = -1, window_distance = 0;
static int window_area = 0;
static unsigned char *valid = NULL;
if (!valid) {
    size_t count = (size_t)sico.M * sico.N;
    const char *path = getenv("NEPAL_VALID");
    FILE *mask = path ? fopen(path, "rb") : NULL;
    valid = malloc(count);
    if (!mask || !valid || fread(valid, 1, count, mask) != count || fclose(mask)) exit(71);
}
static double closest_sq = 1.e30;
const char *output = getenv("NEPAL_FRAME_DIR");
if (!output || sico.MODEL != 7 || sizeof(float) != 4 || sizeof(int) != 4) exit(71);
char frame_path[PATH_MAX];
if (snprintf(frame_path, sizeof(frame_path), "%s/%05d.bin", output, nout) >= sizeof(frame_path)) exit(71);
FILE *frame = fopen(frame_path, "wb");
if (!frame) exit(72);
int cumulative_area = 0, outlet_cells = 0;
for (int cell = 0; cell < sico.IMAX; cell++) {
    if (cdomain2[cell] != 0) continue;
    int index = px[cell] * sico.N + py[cell];
    double dx = sico.BDWEST + (py[cell]+0.5)*sico.CSZ - outlet_x;
    double dy = sico.BDNORTH - (px[cell]+0.5)*sico.CSZ - outlet_y;
    double distance_sq = dx*dx + dy*dy;
    if (valid[index] && awt[cell][32] >= 0.1f) {
        cumulative_area++;
        if (distance_sq < closest_sq) closest_sq = distance_sq;
    }
    float depth = awt[cell][16];
    if (depth < 0.1f) continue;
    float fields[] = {depth, awt[cell][13], awt[cell][15]};
    if (fwrite(&index, sizeof(index), 1, frame) != 1 || fwrite(fields, sizeof(float), 3, frame) != 3) exit(72);
    if (valid[index] && distance_sq <= 300.0*300.0) outlet_cells++;
}
if (fclose(frame)) exit(72);
if (outlet_cells >= 3 && first_arrival < 0) first_arrival = tsum;
if (first_arrival >= 0 && previous_area > 0 && cumulative_area >= previous_area
    && cumulative_area - previous_area <= 0.005 * previous_area) stable_intervals++;
else stable_intervals = 0;
previous_area = cumulative_area;
const char *reason = "running";
if (!csuccess) reason = "numerical_failure";
else if (!ccontinue) reason = "flow_stopped";
else if (first_arrival >= 0 && tsum-first_arrival >= 600 && stable_intervals >= 5) {
    ccontinue = 0;
    reason = "outlet_stable";
} else if (tsum >= tmax) {
    ccontinue = 0;
    reason = "safety_cap";
} else if (first_arrival < 0 && tsum >= 7200 && cumulative_area > 0) {
    double distance = sqrt(closest_sq);
    if (window_start < 0 || cumulative_area > 1.005 * window_area ||
        window_distance - distance >= sico.CSZ) {
        window_start = tsum;
        window_distance = distance;
        window_area = cumulative_area;
    } else if (tsum - window_start >= 3600) {
        ccontinue = 0;
        reason = "stalled";
    }
}
if (snprintf(frame_path,sizeof(frame_path),"%s/timeline.csv",output) >= sizeof(frame_path)) exit(71);
FILE *timeline = fopen(frame_path,"a");
if (!timeline || fprintf(timeline,"%d,%.6f,%d,%d,%.6f,%s\n",nout,tsum,cumulative_area,outlet_cells,first_arrival,reason)<0 || fclose(timeline)) exit(72);
