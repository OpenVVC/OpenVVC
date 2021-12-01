if HAVE_NEON
noinst_LTLIBRARIES = libarmoptim.la
libarmoptim_la_CPPFLAGS = -I${srcdir}/../

libarmoptim_la_SOURCES =


if HAVE_NEON_SIMDE
libarmoptim_la_SOURCES +=	rcn_transform_sse.c         \
							vvc_utils_sse.c             \
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
if !HAVE_EXCLUDE_ANDROID
libarmoptim_la_SOURCES +=	rcn_intra_dc_planar_sse.c
endif
endif

if HAVE_NEON_ASM
libarmoptim_la_SOURCES +=	rcn_mc_neon.S               \
							rcn_mc_neon_wrapper.c       
endif
endif
