file="minrectangle_controller_"
filelog="minrectangle_controller_log_"
dir1="minrect_ANN/minrect_ANN_"
slash="/"
mkdir minrect_ANN
for k in -2 -1 0 1 2 3 4 5 6 7 8 9 10 11
          	do
                    dir2=$dir1$k
                    file_name=$dir2$slash$file$k
                    log_name=$dir2$slash$filelog$k
                   # echo $file_name
                    rm -rf $dir2
                    mkdir $dir2

                    if [ $k -eq -2 ]; then
                             FILE1=$file_name USE_SCALE_FACTOR=false DEFAULT_VALUE=-1 mean_exec_time=0.1 lower_frac=0.1 upper_frac=0.1 ADVANCED_CONTROLLER=true SOURCE_FILE=test.avi DEST_DIR=$dir2 ./rtftr > $log_name
                    else
                             FILE1=$file_name USE_SCALE_FACTOR=false DEFAULT_VALUE=$k mean_exec_time=0.1 lower_frac=0.1 upper_frac=0.1 ADVANCED_CONTROLLER=false SOURCE_FILE=test.avi DEST_DIR=$dir2 ./rtftr > $log_name
                    fi
                done
      



