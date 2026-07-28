#include "../include/all.hxx"

#include <cerrno>
#include <climits>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <string.h>

using namespace std;

/*
  Restart-state readers for standalone LASAM output.

  The standalone driver writes one restart row per model output time.  These
  readers intentionally use the last non-header data row from each file, which
  lets a later standalone run restart from the final state saved by an earlier
  run.  Values are read in the model's natural units, which are cm for depths,
  pressure heads, storage, and flux-like accumulated depths.

  The public functions in this file exit with an error message if a requested
  restart file is missing or malformed.  That is deliberate: a partial or
  silently misread restart state would be much harder to diagnose than a hard
  initialization failure.
*/

/* One parsed wetting front tuple from data_layers.csv. */
typedef struct {
  double depth_cm;        // Wetting front depth in model-native cm.
  double theta;           // Volumetric water content for the front.
  int layer_num;          // 1-based soil layer containing this front.
  int front_num;          // 1-based wetting front ordering index.
  bool to_bottom;         // Whether this front is in contact with a layer bottom.
  double psi_cm;          // Capillary pressure head in model-native cm.
  double dzdt_cm_per_h;   // Saved front velocity in cm/h, restored for continuity.
} WFRecord;

/* qsort callback used to rebuild the linked list in front_num order. */
static int cmp_record_by_front_num(const void *a, const void *b)
{
  const WFRecord *ra = (const WFRecord*)a;
  const WFRecord *rb = (const WFRecord*)b;
  return (ra->front_num - rb->front_num);
}

/*
  Return true for rows that can contain restart data.

  The standalone files include human-readable headers and, for wetting fronts,
  timestep comment lines.  The restart reader skips those rows and keeps only
  actual saved-state records.
*/
static bool is_restart_data_line(const char *buf)
{
  const char *p = buf;
  while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;

  if (*p == '\0' || *p == '#') return false;
  if (strncmp(p, "Time,", 5) == 0) return false;

  return true;
}

/*
  Read through the whole file and keep the final restart data row.

  This is intentionally not "first data row" behavior: data_layers.csv,
  data_non_vadose_state.csv, and data_giuh_state.csv are time series of saved
  states, so the final valid row is the state at the end of the previous run.
*/
static bool read_last_data_line(FILE *fp, char *buf, size_t buflen)
{
  char line[65536];
  bool found = false;

  while (fgets(line, sizeof(line), fp)) {
    if (!is_restart_data_line(line)) continue;

    snprintf(buf, buflen, "%s", line);
    found = true;
  }

  return found;
}

/*
  Parses state lines like:
    [(18,0.197519,1,1,1,2000,1e-05)|(...)]

  Each tuple is:
    (depth_cm,theta,layer_num,front_num,to_bottom,psi_cm,dzdt_cm_per_h)

  The parser only translates text into WFRecord values.  Physical consistency
  checks, sorting, linked-list insertion, and hydraulic-property reconstruction
  happen later in InitializeWettingFrontsFromCSV().
*/
static int parse_state_line_to_records(const char *line_in, WFRecord **recs_out, int *n_out)
{
  /* Give callers predictable outputs even if parsing fails early. */
  *recs_out = NULL;
  *n_out = 0;

  /*
    strdup lets strtok_r and bracket replacement modify a private copy of the
    row.  The input line may point into a caller-owned buffer.
  */
  char *line = strdup(line_in);
  if (!line) return 1;

  /* Wetting front records are enclosed in square brackets. */
  char *lbr = strchr(line, '[');
  char *rbr = strrchr(line, ']');
  if (!lbr || !rbr || rbr <= lbr) {
    free(line);
    return 2;
  }
  *rbr = '\0';
  char *p = lbr + 1;

  /*
    The number of wetting fronts changes during a simulation, so grow the
    record array dynamically rather than assuming one front per soil layer.
  */
  int cap = 16;
  int n = 0;
  WFRecord *recs = (WFRecord*)malloc(sizeof(WFRecord) * cap);
  if (!recs) {
    free(line);
    return 3;
  }

  char *saveptr = NULL;
  for (char *tok = strtok_r(p, "|", &saveptr); tok != NULL; tok = strtok_r(NULL, "|", &saveptr)) {
    /* Be tolerant of spaces and a leading comma before each tuple. */
    while (*tok == ' ' || *tok == '\t' || *tok == ',') tok++;

    double depth_cm, theta, psi_cm, dzdt_cm_per_h;
    int layer_num, front_num, to_bottom_int;
    int consumed = 0;

    int matched = sscanf(tok, " ( %lf , %lf , %d , %d , %d , %lf , %lf ) %n",
                         &depth_cm, &theta, &layer_num, &front_num,
                         &to_bottom_int, &psi_cm, &dzdt_cm_per_h, &consumed);

    /* A restart written by the current writer must contain all seven fields. */
    if (matched != 7 || consumed == 0 || tok[consumed] != '\0') {
      free(recs);
      free(line);
      return 4;
    }

    if (!std::isfinite(depth_cm) || !std::isfinite(theta) ||
        !std::isfinite(psi_cm) || !std::isfinite(dzdt_cm_per_h)) {
      free(recs);
      free(line);
      return 5;
    }

    if (to_bottom_int != 0 && to_bottom_int != 1) {
      free(recs);
      free(line);
      return 6;
    }

    if (n == cap) {
      cap *= 2;
      WFRecord *tmp = (WFRecord*)realloc(recs, sizeof(WFRecord) * cap);
      if (!tmp) {
        free(recs);
        free(line);
        return 7;
      }
      recs = tmp;
    }

    recs[n].depth_cm = depth_cm;
    recs[n].theta = theta;
    recs[n].layer_num = layer_num;
    recs[n].front_num = front_num;
    recs[n].to_bottom = (to_bottom_int == 1);
    recs[n].psi_cm = psi_cm;
    recs[n].dzdt_cm_per_h = dzdt_cm_per_h;
    n++;
  }

  free(line);
  *recs_out = recs;
  *n_out = n;
  return 0;
}

extern void InitializeWettingFrontsFromCSV(
    int num_layers,
    const char *data_layers_csv_path,
    int *layer_soil_type,
    double *cum_layer_thickness_cm,
    double *frozen_factor,
    struct wetting_front **head,
    struct soil_properties_ *soil_properties)
{
  /*
    Restart initialization owns the wetting front linked list from scratch.
    If anything fails after this point, the function exits rather than leaving
    a partly initialized model.
  */
  *head = NULL;

  FILE *fp = fopen(data_layers_csv_path, "r");
  if (!fp) {
    fprintf(stderr, "ERROR: could not open data_layers file: %s\n", data_layers_csv_path);
    exit(1);
  }

  /* Use the final saved state, not the first state in the file. */
  char state_line[65536];
  if (!read_last_data_line(fp, state_line, sizeof(state_line))) {
    fclose(fp);
    fprintf(stderr, "ERROR: no wetting front restart state lines found in %s\n",
            data_layers_csv_path);
    exit(1);
  }
  fclose(fp);

  WFRecord *recs = NULL;
  int nrecs = 0;
  int perr = parse_state_line_to_records(state_line, &recs, &nrecs);
  if (perr != 0 || !recs || nrecs <= 0) {
    fprintf(stderr,
            "ERROR: failed to parse state line in %s (err=%d). Expected "
            "(depth_cm,theta,layer_num,front_num,to_bottom,psi_cm,dzdt_cm_per_h) tuples.\n",
            data_layers_csv_path, perr);
    exit(1);
  }

  /*
    The linked list order is hydrologically meaningful.  Sort by the saved
    front_num before rebuilding it, then check that no front numbers are
    missing or duplicated.
  */
  qsort(recs, nrecs, sizeof(WFRecord), cmp_record_by_front_num);

  if (recs[0].front_num != 1) {
    fprintf(stderr, "ERROR: first front_num in file is %d, but restart requires numbering from 1\n",
            recs[0].front_num);
    free(recs);
    exit(1);
  }

  for (int i = 1; i < nrecs; i++) {
    if (recs[i].front_num != recs[i - 1].front_num + 1) {
      fprintf(stderr, "ERROR: front_num values are not contiguous at %d -> %d\n",
              recs[i - 1].front_num, recs[i].front_num);
      free(recs);
      exit(1);
    }
  }

  for (int i = 0; i < nrecs; i++) {
    WFRecord r = recs[i];

    /* layer_num is 1-based everywhere in this model. */
    if (r.layer_num < 1 || r.layer_num > num_layers) {
      fprintf(stderr, "ERROR: record layer_num=%d out of range 1..%d\n", r.layer_num, num_layers);
      free(recs);
      exit(1);
    }

    /*
      A front must fall inside the layer it says it belongs to.  This also
      catches old restart files written in mm/x10 units before this reader and
      writer were changed to use cm directly.
    */
    double layer_top_cm = cum_layer_thickness_cm[r.layer_num - 1];
    double layer_bottom_cm = cum_layer_thickness_cm[r.layer_num];
    if (r.depth_cm <= 0.0 ||
        r.depth_cm < layer_top_cm - 1.0e-8 ||
        r.depth_cm > layer_bottom_cm + 1.0e-8) {
      fprintf(stderr,
              "ERROR: restart front_num=%d has depth_cm=%.17g outside layer %d "
              "bounds [%.17g, %.17g] cm.\n",
              r.front_num, r.depth_cm, r.layer_num,
              layer_top_cm, layer_bottom_cm);
      free(recs);
      exit(1);
    }

    /*
      listInsertFront restores depth, theta, front number, layer number, and
      to_bottom.  to_bottom is restart state and is intentionally loaded from
      the file rather than inferred from layer ordering.
    */
    struct wetting_front *current =
        listInsertFront(r.depth_cm, r.theta, r.front_num, r.layer_num, r.to_bottom, head);

    if (current == NULL) {
      fprintf(stderr, "ERROR: listInsertFront returned NULL inserting front_num=%d\n", r.front_num);
      free(recs);
      exit(1);
    }

    /* These state variables are not set by listInsertFront, so restore them. */
    current->psi_cm = r.psi_cm;
    current->dzdt_cm_per_h = r.dzdt_cm_per_h;

    /*
      Hydraulic conductivity is derived state.  Recompute it from the restored
      theta and the soil properties for the front's current layer instead of
      trusting an additional saved value.
    */
    int soil = layer_soil_type[r.layer_num];
    double Se = calc_Se_from_theta(current->theta,
                                   soil_properties[soil].theta_e,
                                   soil_properties[soil].theta_r);
    double Ksat_cm_per_h = frozen_factor[r.layer_num] * soil_properties[soil].Ksat_cm_per_h;
    current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, soil_properties[soil].vg_m);
  }

  free(recs);
}

/*
  Restart variables that are not part of the wetting front linked list.

  These values carry conceptual-reservoir storage, the surface-runoff memory
  needed for wetting front creation, flux caching state used during dry
  timesteps, and the model clock.  The writer stores them as key=value pairs
  so this parser can be insensitive to column position after the leading time
  column.
*/
typedef struct {
  double CR_fast_storage_cm;
  double CR_slow_storage_cm;
  double volon_timestep_cm;
  bool runoff_in_prev_step;
  double precip_previous_timestep_cm;
  bool cache_fluxes;
  int cache_count;
  double previous_AET;
  double previous_PET;
  double previous_recharge;
  double accumulated_PET;
  double accumulated_free_drainage;
  double time_s;
  int timesteps;
} nonvadoseRestartState;

/* Booleans are written as 0/1 in the restart CSV files. */
static bool parse_bool_01(const char *s, bool *out)
{
  if (strcmp(s, "1") == 0) {
    *out = true;
    return true;
  }
  if (strcmp(s, "0") == 0) {
    *out = false;
    return true;
  }
  return false;
}

/* Parse a complete finite floating-point token instead of silently accepting
   malformed text as zero through strtod(). */
static bool parse_finite_double(const char *s, double *out)
{
  char *end = NULL;
  errno = 0;
  double value = strtod(s, &end);

  if (end == s || errno == ERANGE || !std::isfinite(value)) return false;
  while (*end == ' ' || *end == '\t') end++;
  if (*end != '\0') return false;

  *out = value;
  return true;
}

/* Parse a complete integer token and reject overflow or trailing characters. */
static bool parse_int(const char *s, int *out)
{
  char *end = NULL;
  errno = 0;
  long value = strtol(s, &end, 10);

  if (end == s || errno == ERANGE || value < INT_MIN || value > INT_MAX) return false;
  while (*end == ' ' || *end == '\t') end++;
  if (*end != '\0') return false;

  *out = (int)value;
  return true;
}

static int parse_non_vadose_state_kv_line(const char *line_in, nonvadoseRestartState *rst)
{
  /*
    Older restart rows may not have all cache history or model clock fields.
    Defaults here keep those optional values cold-started unless explicit
    values are found below.
  */
  rst->cache_fluxes = false;
  rst->cache_count = 1;
  rst->previous_AET = 0.0;
  rst->previous_PET = 0.0;
  rst->previous_recharge = 0.0;
  rst->accumulated_PET = 0.0;
  rst->accumulated_free_drainage = 0.0;
  rst->time_s = 0.0;
  rst->timesteps = 0;

  char *line = strdup(line_in);
  if (!line) return 1;

  /*
    These fields are required because they affect immediate restart behavior:
    reservoir storage and the previous-step surface-water state.
  */
  bool got_fast = false;
  bool got_slow = false;
  bool got_volon = false;
  bool got_runoff = false;
  bool got_precip_prev = false;

  char *saveptr = NULL;
  for (char *tok = strtok_r(line, ",\r\n", &saveptr);
       tok != NULL;
       tok = strtok_r(NULL, ",\r\n", &saveptr)) {
    while (*tok == ' ' || *tok == '\t') tok++;

    /* Each non-vadose state token is key=value. */
    char *eq = strchr(tok, '=');
    if (!eq) {
      free(line);
      return 2;
    }

    *eq = '\0';
    const char *key = tok;
    const char *val = eq + 1;

    /*
      Unknown keys are ignored.  That keeps the reader forward-tolerant if a
      future writer adds more restart metadata that older code can safely skip.
    */
    if (strcmp(key, "CR_fast_storage_cm") == 0) {
      if (!parse_finite_double(val, &rst->CR_fast_storage_cm)) {
        free(line);
        return 6;
      }
      got_fast = true;
    }
    else if (strcmp(key, "CR_slow_storage_cm") == 0) {
      if (!parse_finite_double(val, &rst->CR_slow_storage_cm)) {
        free(line);
        return 6;
      }
      got_slow = true;
    }
    else if (strcmp(key, "volon_timestep_cm") == 0) {
      if (!parse_finite_double(val, &rst->volon_timestep_cm)) {
        free(line);
        return 6;
      }
      got_volon = true;
    }
    else if (strcmp(key, "runoff_in_prev_step") == 0) {
      if (!parse_bool_01(val, &rst->runoff_in_prev_step)) {
        free(line);
        return 3;
      }
      got_runoff = true;
    }
    else if (strcmp(key, "precip_previous_timestep_cm") == 0) {
      if (!parse_finite_double(val, &rst->precip_previous_timestep_cm)) {
        free(line);
        return 6;
      }
      got_precip_prev = true;
    }
    else if (strcmp(key, "cache_fluxes") == 0) {
      if (!parse_bool_01(val, &rst->cache_fluxes)) {
        free(line);
        return 5;
      }
    }
    else if (strcmp(key, "cache_count") == 0) {
      if (!parse_int(val, &rst->cache_count)) {
        free(line);
        return 7;
      }
      if (rst->cache_count < 1) rst->cache_count = 1;
    }
    else if (strcmp(key, "previous_AET") == 0) {
      if (!parse_finite_double(val, &rst->previous_AET)) {
        free(line);
        return 6;
      }
    }
    else if (strcmp(key, "previous_PET") == 0) {
      if (!parse_finite_double(val, &rst->previous_PET)) {
        free(line);
        return 6;
      }
    }
    else if (strcmp(key, "previous_recharge") == 0) {
      if (!parse_finite_double(val, &rst->previous_recharge)) {
        free(line);
        return 6;
      }
    }
    else if (strcmp(key, "accumulated_PET") == 0) {
      if (!parse_finite_double(val, &rst->accumulated_PET)) {
        free(line);
        return 6;
      }
    }
    else if (strcmp(key, "accumulated_free_drainage") == 0) {
      if (!parse_finite_double(val, &rst->accumulated_free_drainage)) {
        free(line);
        return 6;
      }
    }
    else if (strcmp(key, "time_s") == 0) {
      if (!parse_finite_double(val, &rst->time_s)) {
        free(line);
        return 6;
      }
    }
    else if (strcmp(key, "timesteps") == 0) {
      if (!parse_int(val, &rst->timesteps)) {
        free(line);
        return 7;
      }
    }
  }

  free(line);

  if (rst->time_s < 0.0 || rst->timesteps < 0) return 8;

  /* Refuse restart if any required non-vadose state is absent. */
  if (!(got_fast && got_slow && got_volon && got_runoff && got_precip_prev)) {
    return 4;
  }

  return 0;
}

extern void InitializenonvadoseStateFromCSV(
    const char *non_vadose_state_csv_path,
    struct model_state *state)
{
  FILE *fp = fopen(non_vadose_state_csv_path, "r");
  if (!fp) {
    fprintf(stderr, "ERROR: could not open non-vadose state file: %s\n",
            non_vadose_state_csv_path);
    exit(1);
  }

  /* Use the final saved non-vadose state from the previous run. */
  char line[65536];
  if (!read_last_data_line(fp, line, sizeof(line))) {
    fclose(fp);
    fprintf(stderr, "ERROR: no data lines found in non-vadose state file: %s\n",
            non_vadose_state_csv_path);
    exit(1);
  }

  nonvadoseRestartState rst;
  char *parse_start = line;

  /*
    data_non_vadose_state.csv starts with a Time column.  If a comma appears
    before the first key=value pair, skip that leading time field and parse only
    the restart-state tokens.
  */
  char *comma = strchr(line, ',');
  char *first_equal = strchr(line, '=');
  if (comma && first_equal && comma < first_equal) parse_start = comma + 1;

  int perr = parse_non_vadose_state_kv_line(parse_start, &rst);

  fclose(fp);

  if (perr != 0) {
    fprintf(stderr, "ERROR: failed to parse non-vadose restart state in %s (err=%d)\n",
            non_vadose_state_csv_path, perr);
    exit(1);
  }

  /*
    Restore reservoir storage and the surface-water memory used on the next
    timestep.  lgar_initialize() will use these reservoir values to initialize
    volCRstart/volCRend consistently for the restart mass balance.
  */
  state->lgar_mass_balance.CR_fast_storage_cm = rst.CR_fast_storage_cm;
  state->lgar_mass_balance.CR_slow_storage_cm = rst.CR_slow_storage_cm;
  state->lgar_mass_balance.volon_timestep_cm = rst.volon_timestep_cm;
  state->lgar_bmi_params.runoff_in_prev_step = rst.runoff_in_prev_step;
  state->lgar_bmi_params.precip_previous_timestep_cm = rst.precip_previous_timestep_cm;
  state->lgar_bmi_params.time_s = rst.time_s;
  state->lgar_bmi_params.timesteps = rst.timesteps;

  /*
    Flux-cache values are only meaningful when the current config allows flux
    caching.  If the restart file contains cached fluxes but the current config
    disables caching, intentionally reset those fields.
  */
  if (state->lgar_bmi_params.allow_flux_caching) {
    state->lgar_mass_balance.cache_fluxes = rst.cache_fluxes;
    state->lgar_bmi_params.cache_count = rst.cache_count;
    state->lgar_mass_balance.previous_AET = rst.previous_AET;
    state->lgar_mass_balance.previous_PET = rst.previous_PET;
    state->lgar_mass_balance.previous_recharge = rst.previous_recharge;
    state->lgar_mass_balance.accumulated_PET = rst.accumulated_PET;
    state->lgar_mass_balance.accumulated_free_drainage = rst.accumulated_free_drainage;
  }
  else {
    state->lgar_mass_balance.cache_fluxes = false;
    state->lgar_bmi_params.cache_count = 1;
    state->lgar_mass_balance.previous_AET = 0.0;
    state->lgar_mass_balance.previous_PET = 0.0;
    state->lgar_mass_balance.previous_recharge = 0.0;
    state->lgar_mass_balance.accumulated_PET = 0.0;
    state->lgar_mass_balance.accumulated_free_drainage = 0.0;
  }

  /*
    These are set here for completeness; lgar_initialize() also recomputes the
    reservoir start/end bookkeeping immediately after config initialization.
  */
  state->lgar_mass_balance.volCRend_cm = rst.CR_fast_storage_cm + rst.CR_slow_storage_cm;
  state->lgar_mass_balance.volCRend_timestep_cm = rst.CR_fast_storage_cm + rst.CR_slow_storage_cm;
}

/* Parsed representation of one data_giuh_state.csv row. */
typedef struct {
  int num_giuh_ordinates;
  double *queue_vals;
} GIUHRestartState;

/*
  Parse one GIUH restart row, for example:
    num_giuh_ordinates=5,queue=[0,0.1,0.2,0,0,0]

  The queue has num_giuh_ordinates + 1 entries because the existing GIUH code
  uses indices 0..num_giuh_ordinates.
*/
static int parse_giuh_state_line(const char *line_in, GIUHRestartState *rst)
{
  /* Initialize outputs so the caller can safely inspect/free after failure. */
  rst->num_giuh_ordinates = -1;
  rst->queue_vals = NULL;

  char *line = strdup(line_in);
  if (!line) return 1;

  char *num_ptr = strstr(line, "num_giuh_ordinates=");
  char *queue_ptr = strstr(line, "queue=[");
  if (!num_ptr || !queue_ptr) {
    free(line);
    return 2;
  }

  char *num_value = num_ptr + strlen("num_giuh_ordinates=");
  char *num_end = strchr(num_value, ',');
  if (!num_end) {
    free(line);
    return 3;
  }
  *num_end = '\0';
  if (!parse_int(num_value, &rst->num_giuh_ordinates) ||
      rst->num_giuh_ordinates < 0) {
    free(line);
    return 3;
  }

  /* Queue values are enclosed in square brackets and separated by commas. */
  char *lbr = strchr(queue_ptr, '[');
  char *rbr = strrchr(queue_ptr, ']');
  if (!lbr || !rbr || rbr <= lbr) {
    free(line);
    return 4;
  }

  *rbr = '\0';
  char *vals = lbr + 1;

  int expected_n = rst->num_giuh_ordinates + 1;
  rst->queue_vals = (double*)malloc(sizeof(double) * expected_n);
  if (!rst->queue_vals) {
    free(line);
    return 5;
  }

  int nread = 0;
  char *saveptr = NULL;
  for (char *tok = strtok_r(vals, ",", &saveptr);
       tok != NULL;
       tok = strtok_r(NULL, ",", &saveptr)) {
    while (*tok == ' ' || *tok == '\t') tok++;

    /* More values than expected means this file does not match the config. */
    if (nread >= expected_n) {
      free(rst->queue_vals);
      free(line);
      return 6;
    }

    if (!parse_finite_double(tok, &rst->queue_vals[nread])) {
      free(rst->queue_vals);
      free(line);
      return 8;
    }
    nread++;
  }

  free(line);

  /* Fewer values than expected is also a malformed restart row. */
  if (nread != expected_n) {
    free(rst->queue_vals);
    rst->queue_vals = NULL;
    return 7;
  }

  return 0;
}

extern void InitializeGIUHRunoffQueueFromCSV(
    const char *giuh_state_csv_path,
    double *giuh_runoff_queue,
    int num_giuh_ordinates)
{
  FILE *fp = fopen(giuh_state_csv_path, "r");
  if (!fp) {
    fprintf(stderr, "ERROR: could not open GIUH state file: %s\n", giuh_state_csv_path);
    exit(1);
  }

  /* Use the final saved GIUH queue from the previous run. */
  char line[65536];
  if (!read_last_data_line(fp, line, sizeof(line))) {
    fclose(fp);
    fprintf(stderr, "ERROR: no data lines found in GIUH state file: %s\n", giuh_state_csv_path);
    exit(1);
  }

  GIUHRestartState rst;
  int perr = parse_giuh_state_line(line, &rst);

  fclose(fp);

  if (perr != 0) {
    fprintf(stderr, "ERROR: failed to parse GIUH restart state in %s (err=%d)\n",
            giuh_state_csv_path, perr);
    exit(1);
  }

  /*
    The saved queue length must match the current GIUH configuration.  Loading
    a queue from a different GIUH shape would shift runoff timing.
  */
  if (rst.num_giuh_ordinates != num_giuh_ordinates) {
    fprintf(stderr,
            "ERROR: GIUH restart file has num_giuh_ordinates=%d but model expects %d\n",
            rst.num_giuh_ordinates, num_giuh_ordinates);
    free(rst.queue_vals);
    exit(1);
  }

  /* Copy the parsed queue into the BMI-owned queue array. */
  for (int i = 0; i <= num_giuh_ordinates; i++) {
    giuh_runoff_queue[i] = rst.queue_vals[i];
  }

  free(rst.queue_vals);
}
