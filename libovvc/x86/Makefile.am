if HAVE_SSE4_1
noinst_LTLIBRARIES = libx86optim.la
libx86optim_la_CPPFLAGS = -I${srcdir}/../
libx86optim_la_SOURCES =	rcn_transform_sse.c         \
							vvc_utils_sse.c             \
							rcn_intra_dc_planar_sse.c   \
							rcn_mc_sse.c                \
							rcn_transform_add_sse.c     \
							rcn_alf_sse.c               \
							rcn_intra_cclm_sse.c        \
							rcn_lfnst_sse.c             \
							rcn_intra_mip_sse.c         \
							rcn_sao_sse.c               \
							rcn_dmvr_sse.c              \
							rcn_prof_bdof_sse.c         \
							rcn_df_sse.c
endif
