# run velocyto 

-----
### Example run at NIG supercomputer (05/18/24)
```sh
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l medium
#$ -l s_vmem=96G,mem_req=96G
#$ -o velocyto.txt
#$ -e velocyto_Err.txt


source ~/.bashrc

# move to working directoryへ移動
wd=$HOME/../path_to_DIR
cd $wd

conda activate velo3.6

## Step1. gtfファイルの準備：  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# wdにgene.gtfファイルを用意する（wdにないとダメっぽい）
LOCK_FILE="genes.gtf"
if [ -f ${LOCK_FILE} ]; then
  
    echo "genes.gtf file used in Cellranger Count command exists in current directory ! "

  else

    echo "genes.gtf is copied to current directoy, copy genes.gtf files in current directory"

    # ファイルがない場合，working directoryへコピーする
    cp $HOME/../genes.gtf .
    #gunzip genes.gtf.gz
fi


## Step2. sample listの準備: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## cellranger countのフォルダを読み出してリスト作成．
## すでにリストファイルがあった場合，削除して新たに作成．
LOCK_FILE2="SampleList.txt"

if [ -f ${LOCK_FILE2} ]; then

    rm ${LOCK_FILE2}

fi

ls $HOME/../fastq -F | grep / | sed s/\\/$//g >> ./${LOCK_FILE2}


## Step3. run velocyto>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cat ./${LOCK_FILE2} | while read line
do
  velocyto run10x \
    --mask $HOME/Genome/maskedGTF/mm39_UCSC_rmsk.gtf \
    $wd/${line} \
    ./genes.gtf
done

```
