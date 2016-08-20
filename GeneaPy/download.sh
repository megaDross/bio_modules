#!/bin/bash

question() {
    while [[ "$ANSWER" != "y" && "$ANSWER" != "n" ]]; do
        if [ "$ANSWER" == "y" ]; then
            $2
        fi
        echo "$1"
        read ANSWER
    done
}


ttuner() {
    URL="http://downloads.sourceforge.net/project/tracetuner/tracetuner/tracetuner_3.0.6beta/tracetuner_3.0.6beta.tar.bz2"

    wget -P ~/bin/ $URL 
    tar -xjf ~/bin/tracetuner_3.0.6beta.tar.bz2 -C ~/bin
    cd ~/bin/tracetuner_3.0.6beta/src/; make
}


genome() {
    echo "Which version of the human genome to you wish to download (hg19 or hg38)?"
    read GENOME

    if [ "$GENOME" == "hg19" ]; then
        URL="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz"
    elif [ "$GENOME" == "hg38" ]; then
        URL="http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz"
    else
        echo "$GENOME is invalid. Type hg38 or hg19"
        exit 1
    fi

    mkdir ~/.configure/genome
    wget -P ~/.configure/genome/ $URL
    cat ~/.configure/genome/chr*.fa > ~/.configure/genome/$GENOME.fa
    rm ~/.configure/genome/chr*.fa
}


question "Do you want to download ttuner (y or n)?" ttuner

question "Do you want to download the human genome (y or n)?" genome




