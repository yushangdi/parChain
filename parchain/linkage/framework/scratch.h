#ifdef PERF_RANGE
  double range_time = 0;
  double find_nns_time = 0;
  // intT cand_num = 0;
  double max_round_time = 0;
  double max_range_time = 0;
  double max_nns_time = 0;

inline void report_perf_range(intT round){
    if(round > 1 && (round < 5 || round % PRINT_FREQ == 0)){
      UTIL::PrintFunctionItem("PERF_RANGE", "nns", find_nns_time);
      UTIL::PrintFunctionItem("PERF_RANGE", "init_tree_time", init_tree_time);
      // UTIL::PrintFunctionItem("PERF_RANGE", "tree_bb_time", tree_bb_time);
      UTIL::PrintFunctionItem("PERF_RANGE", "range+nnc", range_time);
      UTIL::PrintFunctionItem("PERF_RANGE", "cand_num", cand_num);
  }
}

inline void reset_perf_range(){
    range_time = 0;
    find_nns_time = 0;
    init_tree_time =  0;
    // tree_bb_time = 0;
    cand_num = 0;
  }
#endif
