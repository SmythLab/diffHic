#include "diffhic.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

extern "C" { 

static const R_CallMethodDef all_call_entries[] = {
	CALLDEF(check_input, 2),
	CALLDEF(cap_input, 3),
	CALLDEF(cluster_2d, 6),
	CALLDEF(split_clusters, 6),
	CALLDEF(get_bounding_box, 3),
	CALLDEF(quadrant_bg, 9),
	CALLDEF(count_background, 13),
	CALLDEF(count_connect, 6),
	CALLDEF(count_patch, 5),
	CALLDEF(iterative_correction, 9),
	CALLDEF(report_hic_pairs, 11),
	CALLDEF(report_hic_binned_pairs, 12),
	CALLDEF(test_parse_cigar, 2),
	CALLDEF(test_fragment_assign, 6),
    CALLDEF(pair_stats, 9),
    CALLDEF(get_missing_dist, 4),
    CALLDEF(directionality, 9),
  	{NULL, NULL, 0}
};

void attribute_visible R_init_diffHic(DllInfo *info)
{
	R_registerRoutines(info, NULL, all_call_entries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
	R_forceSymbols(info, TRUE);
}

}
