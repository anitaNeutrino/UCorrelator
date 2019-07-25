### Grab ANITA-4 spectrum averages
#!/bin/bash
UCORRELATOR_SPECAVG_DIR=/unix/anita4/SpectrumAverages

echo "Grabbing ANITA-4 spectrum averages..."

FIRSTRUN=41
LASTRUN=367
CHECK_MARK="\033[0;32m\xE2\x9C\x94\033[0m"
BASE_URL=https://users.rcc.uchicago.edu/~abl/timeAverages/

for (( RUN=$FIRSTRUN; RUN<=$LASTRUN; RUN++ ))
do
    ### Unavailable
    if ([ $RUN -eq 46 ] || [ $RUN -eq 60 ] || [ $RUN -eq 61 ] )
    then
	continue
    fi
    echo -n "Getting averages for run..." $RUN
    wget -q $BASE_URL/$RUN\_10.root -O $UCORRELATOR_SPECAVG_DIR/$RUN\_10.root
    if [ -f $UCORRELATOR_SPECAVG_DIR/$RUN\_10.root ]; then
	echo -e "\\r${CHECK_MARK} Grabbed the spectrum averages for" $RUN
    else
	echo -e "!!!!! Could not grab the spectrum averages for" $RUN
    fi
done
