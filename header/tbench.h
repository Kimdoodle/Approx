#pragma once
#include "seal/seal.h"
#include "CKKS_params.h"
#include "parse.h"
#include "evaluate.h"
#include "measure_time.h"


int find_best_pmd(int s, int depth);
void test_helut_time(int p_num, int e_num, int s_num, int N_num, bool print);
void test_remez_time(int p_num, int e_num, int s_num, int N_num, bool print);