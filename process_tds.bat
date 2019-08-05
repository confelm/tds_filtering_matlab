ECHO OFF
ECHO  "process tds data"
matlab -r "process_tds; quit" -logfile process_tds.log
ECHO "Done"