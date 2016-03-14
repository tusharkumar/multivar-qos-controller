dir1="maxnumiters/maxnumiters_"
file="maxnumiters_controller"
slash="/"
filelog="maxnumiters_controller_log_"

mkdir maxnumiters

for k in -2 -1 0 1 2 3 4 5 6 7 8 9 10 11
      do
          dir2=$dir1$k
          file_name=$dir2$slash$file$k
          log_name=$dir2$slash$filelog$k

          rm -rf $dir2
          mkdir $dir2

          if [ $k -eq -2 ]; then
	            MEAN_EXEC_TIME=0.04 DEST_DIR=$dir2 ADVANCED_CONTROLLER=true LOWER_FRAC=0.1 UPPER_FRAC=0.1 DEFAULT_VALUE=-1 FILENAME=$file_name ./ferns_demo -v mousepad.mp4 -k 100 > $log_name
          else
                    MEAN_EXEC_TIME=0.04 DEST_DIR=$dir2 ADVANCED_CONTROLLER=false LOWER_FRAC=0.1 UPPER_FRAC=0.1 DEFAULT_VALUE=$k FILENAME=$file_name ./ferns_demo -v mousepad.mp4 -k 100 > $log_name
          fi
       done
 

