#ifndef BENCH_DECORATOR_H
#define BENCH_DECORATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "nvcl_structures.h"
struct BenchWrapper
{
	size_t nb_fct_call;
	size_t iter_per_fct;
  	// Define function arguments here
  	int (*func)(const OVSEI *sei, OVFrame **frame);
};

void bench_decorator(struct BenchWrapper bw, char* name, const OVSEI *sei, OVFrame **frame)
{
	volatile int dummy;
	double total_avg_time = 0;

	for (size_t i = 0; i < bw.nb_fct_call; ++i) {
		double avg_time = 0;

		for (size_t j = 0; j < bw.iter_per_fct; ++j) {
			double t1 = clock();
			volatile const int r = bw.func(sei, frame);
			double t2 = clock();

			dummy = r;
			double dt = (t2 - t1) * 1000 / CLOCKS_PER_SEC; // millisecond

			avg_time += dt;
		}

		avg_time /= (double)bw.iter_per_fct;
		total_avg_time += avg_time;
	}

	total_avg_time /= (double)bw.nb_fct_call;
	printf("%-20s: %f ms (average time / %lu calls of %lu iteration each)]\n",
			name,
			total_avg_time,
			bw.nb_fct_call,
			bw.iter_per_fct);
}

#endif