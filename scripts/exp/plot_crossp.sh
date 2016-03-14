#!/bin/sh

# compares the improvements in SR by changing number of X
export LOG_M=43.1
export LOG_F=43.4N2W1bud10
export LOG_R_FR=43.8Accufr
export LOG_R_AC=43.8Accuac
export LOG_R_ACFR=43.8Accuacfr

perl 1pm/scripts/exp/graph_scores.perl crosspmpeg \
       2pm/score_crossp_dolbycanyon_encode_${LOG_M} "2X > 1X" mpeg2enc_2X canyon 2pm/score_crossp_dolbycity640x480_encode_${LOG_M} "2X > 1X" mpeg2enc_2X city640x480 2pm/score_crossp_dolbycity_encode_${LOG_M} "2X > 1X" mpeg2enc_2X city320x240 \
       4pm/score_crossp_dolbycanyon_encode_${LOG_M} "4X > 1X" mpeg2enc_4X canyon 4pm/score_crossp_dolbycity640x480_encode_${LOG_M} "4X > 1X" mpeg2enc_4X city640x480 4pm/score_crossp_dolbycity_encode_${LOG_M} "4X > 1X" mpeg2enc_4X city320x240

perl 1pm/scripts/exp/graph_scores.perl crossprtftr \
       r/score_crossp_exp480x360_minrect_${LOG_R_FR} "fr: mr,sf > mr" rtftr_fr_2X_minrect exp480x360 r/score_crossp_exp640x480_minrect_${LOG_R_FR} "fr: mr,sf > mr" rtftr_fr_2X_minrect exp640x480 \
       r/score_crossp_exp480x360_scalefactor_${LOG_R_FR} "fr: mr,sf > sf" rtftr_fr_2X_scalefactor exp480x360 r/score_crossp_exp640x480_scalefactor_${LOG_R_FR} "fr: mr,sf > sf" rtftr_fr_2X_scalefactor exp640x480 \
       r/score_crossp_exp480x360_minrect_${LOG_R_AC} "ac: mr,sf > mr" rtftr_ac_2X_minrect exp480x360 r/score_crossp_exp640x480_minrect_${LOG_R_AC} "ac: mr,sf > mr" rtftr_ac_2X_minrect exp640x480 \
       r/score_crossp_exp480x360_scalefactor_${LOG_R_AC} "ac: mr,sf > sf" rtftr_ac_2X_scalefactor exp480x360 r/score_crossp_exp640x480_scalefactor_${LOG_R_AC} "ac: mr,sf > sf" rtftr_ac_2X_scalefactor exp640x480 \
       r/score_crossp_exp480x360_minrect_${LOG_R_ACFR} "ac,fr: mr,sf > mr" rtftr_acfr_2X_minrect exp480x360 r/score_crossp_exp640x480_minrect_${LOG_R_ACFR} "ac,fr: mr,sf > mr" rtftr_acfr_2X_minrect exp640x480 \
       r/score_crossp_exp480x360_scalefactor_${LOG_R_ACFR} "ac,fr: mr,sf > sf" rtftr_acfr_2X_scalefactor exp480x360 r/score_crossp_exp640x480_scalefactor_${LOG_R_ACFR} "ac,fr: mr,sf > sf" rtftr_acfr_2X_scalefactor exp640x480
