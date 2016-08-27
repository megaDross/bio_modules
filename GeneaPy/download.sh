#!/bin/bash

question() {
    while true; do

        echo -e "$1"
        read ANSWER

        if [ "$ANSWER" == "y" ]; then
            $2
            break
        elif [ "$ANSWER" == "n" ]; then
            break
        fi
    done
}


ttuner() {
    URL="http://downloads.sourceforge.net/project/tracetuner/tracetuner/tracetuner_3.0.6beta/tracetuner_3.0.6beta.tar.bz2"
    wget -P ~/bin/ $URL 
    tar -xjf ~/bin/tracetuner_3.0.6beta.tar.bz2 -C ~/bin
    cd ~/bin/tracetuner_3.0.6beta/src/; make
}


genome() {
    
    while [ "$GENOME" != "hg19" ] ||  ["$GENOME" != "hg38" ]; do

        echo -e "\nWhich version of the human genome to you wish to download (hg19 or hg38)?"
        read GENOME

        if [ "$GENOME" == "hg19" ]; then
            URL="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz"
        elif [ "$GENOME" == "hg38" ]; then
            URL="http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz"
        else
            echo "$GENOME is invalid. Type hg38 or hg19"
        fi
    done
    
    mkdir ~/.config/genome
    wget -P ~/.config/genome/ $URL
    cat ~/.config/genome/chr*.fa > ~/.config/genome/$GENOME.fa
    rm ~/.config/genome/chr*.fa
}


question "Do you want to download ttuner (y or n)?" ttuner

question "\nDo you want to download the human genome (y or n)?" genome




