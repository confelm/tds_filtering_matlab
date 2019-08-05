#!/bin/bash
set mypath=%cd%
eval matlab -r "process_tds; quit" -logfile process_tds.log