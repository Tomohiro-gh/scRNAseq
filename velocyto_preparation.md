# 環境構築後のvelocytoの準備

### Step1. mappingに使用したgtfファイルを用意する
gtfファイルをworking directoryへコピーしておく

### Step2. masked gtfの準備
繰り返し配列をはじめとした特定の領域をマスクしたgtf ,  genome dataを得る
-> velocytoではこの情報を使用することを推奨している

これを作るためにはrepeat maskerというものが有名らしい
http://www.repeatmasker.org/


取得の仕方は下記URLを参照
https://labs.wsu.edu/winuthayanon/scrna/how-to-analyze-single%e2%80%90cell-rna%e2%80%90seq/explaining-velocyto-command-line/

実際には
UCSC genome browser へアクセス　
https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=Mouse&db=mm10&hgta_group=allTracks&hgta_track=rmsk&hgta_table=0&hgta_regionType=genome&position=chr12%3A56694976-56714605&hgta_outputType=primaryTable&hgta_outputType=gff&hgta_outFileName=mm10_rmsk.gtf


clade : Vertebrate ，genome : Zebrafish を指定
assmebly: GRCz11/DanRer11

Track: RepeatMasker
tableで　rmsk 

region : genome
output format : GTF

output filename : 名前を指定

get outputで作成
