lib_LTLIBRARIES = libovvc.la
libovvc_la_LDFLAGS = -version-info 0:3:0
libovvc_la_LIBADD=

pkgconfigdir = $(libdir)/pkgconfig
nodist_pkgconfig_DATA = libopenvvc.pc

if HAVE_SSE4_1
SUBDIRS = x86/
libovvc_la_LIBADD += ./x86/libx86optim.la
endif

if HAVE_NEON
SUBDIRS = arm/
libovvc_la_LIBADD += ./arm/libarmoptim.la
endif

dsp_template_SOURCES = \
							rcn_alf.c                \
							rcn_ctu.c                \
							rcn_df.c                 \
							rcn_intra.c              \
							rcn_inter.c              \
							rcn_intra_cclm.c         \
							rcn_prof_bdof.c          \
							rcn_transform.c          \
							rcn_transform_scale.c    \
							rcn_transform_tree.c     \
							rcn_fill_ref.c           \
							rcn_intra_angular.c      \
							rcn_intra_dc_planar.c    \
							rcn_lmcs.c               \
							rcn_residuals.c          \
							rcn_dequant.c            \
							rcn_intra_mip.c          \
							rcn_mc.c                 \
							rcn_sao.c

noinst_LTLIBRARIES = libovvc8bit.la
libovvc_la_LIBADD += libovvc8bit.la
libovvc8bit_la_CPPFLAGS = -DBITDEPTH=8
libovvc8bit_la_SOURCES = $(dsp_template_SOURCES)


noinst_LTLIBRARIES += libovvc10bit.la
libovvc_la_LIBADD += libovvc10bit.la
libovvc10bit_la_CPPFLAGS = -DBITDEPTH=10
libovvc10bit_la_SOURCES = $(dsp_template_SOURCES)


libovvc_la_SOURCES = 	ctudec.c                     \
						data_rcn_angular.c           \
						data_rcn_mip.c               \
						data_rcn_transform.c         \
						data_scan_lut.c              \
						dec_init.c                   \
						dpb.c                        \
						dpb_internal.c               \
						mempool.c                    \
						mvpool.c                     \
						nvcl.c                       \
						nvcl_dpb.c                   \
						nvcl_nal_ph.c                \
						nvcl_nal_aps.c               \
						nvcl_nal_pps.c               \
						nvcl_nal_sps.c               \
						nvcl_nal_sei.c               \
						nvcl_ptl.c                   \
						nvcl_rpl.c                   \
						ovdec.c                      \
						ovutils.c                    \
						ovmem.c                      \
						overror.c                    \
						ovdmx.c                      \
						ovio.c                       \
						ovlog.c                      \
						ovannexb.c                   \
						ovunits.c                    \
						ovframe.c                    \
						ovframepool.c                \
						pp_film_grain.c              \
						post_proc.c                  \
						ovthreads.c                  \
						drv_affine_mvp.c             \
						drv_lines.c                  \
						drv_lfnst.c                  \
						drv_intra.c                  \
						drv_mvp.c                    \
						rcn.c                        \
						rcn_lfnst.c                  \
						rcn_gpm.c                    \
						compat_old.c                 \
						slicedec.c                   \
						vcl_alf.c                    \
						vcl_sao.c                    \
						vcl_sh.c                     \
						vcl_coding_unit.c            \
						vcl_coding_tree.c            \
						vcl_transform_unit.c         \
						vcl_residual_coding.c        \
						vcl_cabac.c


include_HEADERS = 	ovdefs.h               \
					ovunits.h              \
					ovlog.h                \
					ovdec.h                \
					ovframe.h              \
					ovio.h                 \
					ovdmx.h                \
					ovversion.h

ovversion.h:
	$(AM_V_GEN)$(srcdir)/../version.sh $(srcdir)/../VERSION ovversion.h ./
