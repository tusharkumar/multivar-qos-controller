#!/bin/sh

# summarizes the SR improvements over all benchmarks, datasets and XY-config: srt over fixed
export SCORE=score_sr

export LOG_M=43.1

export LOG_F=43.4N2W1bud10

export LOG_R_FR=43.8Accufr
export LOG_R_AC=43.8Accuac
export LOG_R_ACFR=43.8Accuacfr

export LOG_X=46.4

export LOG_B=46.1

perl 1pm/scripts/exp/graph_scores.perl summaryfr \
       1pm/dolbycanyon_encode/${SCORE}_${LOG_M} "mpeg2enc 1X" mpeg2enc_1p canyon 1pm/dolbycity640x480_encode/${SCORE}_${LOG_M} "mpeg2enc 1X" mpeg2enc_1p city640x480 1pm/dolbycity_encode/${SCORE}_${LOG_M} "mpeg2enc 1X" mpeg2enc_1p city320x240 \
       2pm/dolbycanyon_encode/${SCORE}_${LOG_M} "mpeg2enc 2X" mpeg2enc_2p canyon 2pm/dolbycity640x480_encode/${SCORE}_${LOG_M} "mpeg2enc 2X" mpeg2enc_2p city640x480 2pm/dolbycity_encode/${SCORE}_${LOG_M} "mpeg2enc 2X" mpeg2enc_2p city320x240 \
       f/0obj/exp320x240/${SCORE}_${LOG_F} "ferns" ferns_0obj exp320x240 f/0obj/exp640x480/${SCORE}_${LOG_F} "ferns" ferns_0obj exp640x480 \
       r/exp480x360/minrect/${SCORE}_${LOG_R_FR} "rtftr mr" rtftr_minrect exp480x360 r/exp640x480/minrect/${SCORE}_${LOG_R_FR} "rtftr mr" rtftr_minrect exp640x480 \
       r/exp480x360/scalefactor/${SCORE}_${LOG_R_FR} "rtftr sf" rtftr_scalefactor exp480x360 r/exp640x480/scalefactor/${SCORE}_${LOG_R_FR} "rtftr sf" rtftr_scalefactor exp640x480 \
       r/exp480x360/both/${SCORE}_${LOG_R_FR} "rtftr mr,sf" rtftr_both exp480x360 r/exp640x480/both/${SCORE}_${LOG_R_FR} "rtftr mr,sf" rtftr_both exp640x480 \
       b/numcores/${SCORE}_${LOG_B} "bodytrack nc" bt_nc v b/numparticles/${SCORE}_${LOG_B} "bodytrack np" bt_np v b/numcores_numparticles/${SCORE}_${LOG_B} "bodytrack ncnp" bt_ncnp v \
       x/numcores/${SCORE}_${LOG_X} "x264 nc" x264_nc d x/subpel/${SCORE}_${LOG_X} "x264 sp" x264_sp d x/merange/${SCORE}_${LOG_X} "x264 me" x264_me d \
       x/numcores_subpel/${SCORE}_${LOG_X} "x264 nc,sp" x264_ncsp d x/numcores_merange/${SCORE}_${LOG_X} "x264 nc,me" x264_ncme d x/subpel_merange/${SCORE}_${LOG_X} "x264 sp,me" x264_spme d \
       x/numcores_subpel_merange/${SCORE}_${LOG_X} "x264 nc,sp,me" x264_ncspme d

perl 1pm/scripts/exp/graph_scores.perl summaryac \
       f/1obj/exp320x240/${SCORE}_${LOG_F} "ferns" ferns_1obj exp320x240 f/1obj/exp640x480/${SCORE}_${LOG_F} "ferns" ferns_1obj exp640x480 \
       r/exp480x360/minrect/${SCORE}_${LOG_R_AC} "rtftr mr" rtftr_minrect exp480x360 r/exp640x480/minrect/${SCORE}_${LOG_R_AC} "rtftr mr" rtftr_minrect exp640x480 \
       r/exp480x360/scalefactor/${SCORE}_${LOG_R_AC} "rtftr sf" rtftr_scalefactor exp480x360 r/exp640x480/scalefactor/${SCORE}_${LOG_R_AC} "rtftr sf" rtftr_scalefactor exp640x480 \
       r/exp480x360/both/${SCORE}_${LOG_R_AC} "rtftr mr,sf" rtftr_both exp480x360 r/exp640x480/both/${SCORE}_${LOG_R_AC} "rtftr mr,sf" rtftr_both exp640x480

perl 1pm/scripts/exp/graph_scores.perl summaryacfr \
       f/2obj/exp320x240/${SCORE}_${LOG_F} "ferns" ferns_2obj exp320x240 f/2obj/exp640x480/${SCORE}_${LOG_F} "ferns" ferns_2obj exp640x480 \
       r/exp480x360/minrect/${SCORE}_${LOG_R_ACFR} "rtftr mr" rtftr_minrect exp480x360 r/exp640x480/minrect/${SCORE}_${LOG_R_ACFR} "rtftr mr" rtftr_minrect exp640x480 \
       r/exp480x360/scalefactor/${SCORE}_${LOG_R_ACFR} "rtftr sf" rtftr_scalefactor exp480x360 r/exp640x480/scalefactor/${SCORE}_${LOG_R_ACFR} "rtftr sf" rtftr_scalefactor exp640x480 \
       r/exp480x360/both/${SCORE}_${LOG_R_ACFR} "rtftr mr,sf" rtftr_both exp480x360 r/exp640x480/both/${SCORE}_${LOG_R_ACFR} "rtftr mr,sf" rtftr_both exp640x480
