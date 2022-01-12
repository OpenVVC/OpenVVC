if HAVE_X86_OPTIM
noinst_LTLIBRARIES = libx86optim.la
libx86optim_la_CPPFLAGS = -I${srcdir}/../ -DBITDEPTH=10
noinst_HEADERS =
libx86optim_la_SOURCES =

if HAVE_SSE4_1
libx86optim_la_SOURCES +=	rcn_transform_sse.c         \
							vvc_utils_sse.c             \
							rcn_intra_dc_planar_sse.c   \
							rcn_mc_sse.c                \
							rcn_transform_add_sse.c     \
							rcn_alf_sse.c               \
							rcn_intra_cclm_sse.c        \
							rcn_lfnst_sse.c             \
							rcn_intra_mip_sse.c         \
							rcn_intra_angular_sse.c     \
							rcn_sao_sse.c               \
							rcn_dmvr_sse.c              \
							rcn_prof_bdof_sse.c         \
							rcn_df_sse.c

noinst_HEADERS += rcn_sse.h
endif

if HAVE_AVX2
libx86optim_la_SOURCES +=	rcn_alf_avx2.c              \
							rcn_sao_avx2.c              \
							rcn_intra_mip_avx2.c        \
							rcn_dmvr_avx2.c             \
							rcn_mc_avx2.c               \
							rcn_transform_add_avx2.c


noinst_HEADERS += rcn_avx2.h
endif
endif
