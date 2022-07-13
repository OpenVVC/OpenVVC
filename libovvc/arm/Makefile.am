if HAVE_NEON
noinst_LTLIBRARIES = libarmoptim.la
libarmoptim_la_CPPFLAGS = -I${srcdir}/../ -DBITDEPTH=10

libarmoptim_la_SOURCES =


if HAVE_NEON_SIMDE
libarmoptim_la_SOURCES +=	rcn_transform_sse.c                   \
							vvc_utils_sse.c                       \
							rcn_mc_sse.c                          \
							rcn_transform_add_sse.c               \
							rcn_alf_sse.c                         \
							rcn_intra_cclm_sse.c                  \
							rcn_lfnst_sse.c                       \
							rcn_intra_mip_sse.c                   \
							rcn_sao_sse.c                         \
							rcn_dmvr_sse.c                        \
							rcn_prof_bdof_sse.c                   \
							rcn_intra_angular_sse.c
if !HAVE_EXCLUDE_CLANG
libarmoptim_la_SOURCES +=	rcn_dequant_sse.c
endif
endif

if HAVE_NEON_ASM
libarmoptim_la_SOURCES +=	rcn_mc_neon.S                         \
							rcn_mc_neon_wrapper.c                 \
							rcn_intra_dc_planar_neon.S            \
							rcn_intra_dc_planar_neon_wrapper.c    \
							rcn_dequant_neon.S                    \
							rcn_dequant_neon_wrapper.c			  \
							rcn_transform_add_neon.S              \
							rcn_transform_add_neon_wrapper.c      \
							rcn_transform_neon.S                  \
                            rcn_transform_neon_wrapper.c
endif
endif

noinst_HEADERS = rcn_neon.h