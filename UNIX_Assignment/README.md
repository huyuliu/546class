# UNIX Assignment

## Data Inspection

### Attributes of `fang_et_al_genotypes`

```linux
$ less fang_et_al_genotypes.txt
$ tail fang_et_al_genotypes.txt
$ wc fang_et_al_genotypes.txt
  2783  2744038 11051939 fang_et_al_genotypes.txt
$ du -h fang_et_al_genotypes.txt
  6.1M	fang_et_al_genotypes.txt
$ tail -n +6 fang_et_al_genotypes.txt | awk -F "\t" '{print NF; exit}'
  986
```

By inspecting this file I learned that:

* It is a file with header on first line, and the first column is the sample ID,while the snp data is start from the forth column.
* It has 2783 lines, 2744038 words and 11051939 byte counts.
* The size of the file is 6.1M.
* It has 986 columns.

### Attributes of `snp_position.txt`

```linux
$ less snp_position.txt
$ tail snp_position.txt
$ wc snp_position.txt 
  984 13198 82763 snp_position.txt
$ du -h snp_position.txt 
  38K	snp_position.txt
$ tail -n +6 snp_position.txt | awk -F "\t" '{print NF; exit}'
  15
```

By inspecting this file I learned that:

* It is a file with header on first line with the snp ID on the first column.
* It has 984 lines, 13198 words and 82763 byte counts.
* The size of the file is 38K.
* It has 15 columns.

## Data Processing

```linux
$ sort -k3 fang_et_al_genotypes.txt > genotypes_sort.txt
$ less -S genotypes_sort.txt
$ awk -f transpose.awk genotypes_sort.txt > transposed_genotypes.txt
$ sort -k1 transposed_genotypes.txt > genotypes_sort_transposed.txt 
$ sort -k1 snp_position.txt > snp_sort.txt
$ join -1 1 -2 1 -a 1 genotypes_sort_transposed.txt snp_sort.txt > all_genotype.txt
  join: genotypes_sort_transposed.txt:479: is not sorted: PZA00710.1 ......
  join: snp_sort.txt:651: is not sorted: PZA03081.1 ......
$ join -1 1 -2 1 genotypes_sort_transposed.txt snp_sort.txt | wc -l
  join: genotypes_sort_transposed.txt:479: is not sorted: PZA00710.1 ......
  join: snp_sort.txt:651: is not sorted: PZA03081.1 ......
  982
$ less -S all_genotype.txt
$ less -S snp_sort.txt
$ less -S genotypes_sort_transposed.txt
$ vi snp_sort.txt
$ vi genotypes_sort_transposed.txt
$ rm all_genotype.txt
$ join -1 1 -2 1 -a 1 -t $'\t' genotypes_sort_transposed.txt snp_sort.txt > all_genotype.txt
  join: genotypes_sort_transposed.txt:479: is not sorted: PZA00710.1 ......
  join: snp_sort.txt:477: is not sorted: PZA00710.1 ......
  join: snp_sort.txt:651: is not sorted: PZA03081.1 ......
$ less -S all_genotype.txt
$ vi all_genotype.txt
$ cut -f3 genotypes_sort.txt | uniq -c 
      1 Group
     22 TRIPS
     15 ZDIPL
     17 ZLUXR
     10 ZMHUE
    290 ZMMIL
   1256 ZMMLR
     27 ZMMMR
    900 ZMPBA
     41 ZMPIL
     34 ZMPJA
     75 ZMXCH
     69 ZMXCP
      6 ZMXIL
      7 ZMXNO
      4 ZMXNT
      9 ZPERR
$ tail -n +6 all_genotype.txt | awk -F " " '{print NF; exit}'
  2795
$ head -n 3 all_genotype.txt | awk  -F "\t" '{for(i=1;i<=NF;i++){if($i~/ZMMIL|ZMMLR|ZMMMR/)print i}}' > maize_column.txt
$ wc -l maize_column.txt
  1573 maize_column.txt
$ head -n 3 all_genotype.txt | awk  -F "\t" '{for(i=1;i<=NF;i++){if($i~/ZMPBA|ZMPIL|ZMPJA/)print i}}' > teosinte_column.txt
$ wc -l teosinte_column.txt
  975 teosinte_column.txt
$ cat maize_column.txt
$ cut -f 1,66-1638,2784- all_genotype.txt > maize_all.txt
$ cat teosinte_column.txt
$ cut -f 1,1639-2613,2784- all_genotype.txt > teosint_all.txt
```

I sorted genotype data by group, and transposed the genotype data, so the columns become rows.
Then I sorted both snp file and transposed genotype file by first column and tried to merged them. But it gave me warning, so I inspected the trouble lines and found out `PZA00710.1` and `PZA00710.16` cannot be sorted correctly in SNP file.
I tested different code and cannot solve it, so I just manually adjusted them.
After that, I changed `Sample_ID` in genotype file into `SNP_ID` to match with SNP file, and merged them again. Even it still gave me warning, it actually worked.
In the end, I moved headers to the first three lines and counted how many maize samples and how many teosinte samples inside the file, and separated them.

### Maize Data

```linux
$ head -n 1 maize_all.txt | awk  -F "\t" '{for(i=1;i<=NF;i++){if($i~/Chromosome/)print i}}'
1576
$ head -n 1 maize_all.txt | awk  -F "\t" '{for(i=1;i<=NF;i++){if($i~/mult_positions/)print i}}'
1579
$ head -n 3 maize_all.txt | tee maize_chr1.txt maize_chr2.txt maize_chr3.txt maize_chr4.txt maize_chr5.txt maize_chr6.txt maize_chr7.txt maize_chr8.txt maize_chr9.txt maize_chr10.txt maize_multi.txt maize_unknown.txt
$ cat maize_all.txt | awk -F "\t" '{if($1579~/^C/)print $0}' >> maize_multi.txt
$ cat maize_all.txt | awk -F "\t" '{if($1576~/^u/)print $0}' >> maize_unknown.txt
$ cat maize_all.txt | awk -F "\t" '{if($1576~/^1$/){if($1579~/^C/);else print $0}}' >> maize_chr1.txt
$ cat maize_all.txt | awk -F "\t" '{if($1576~/2/){if($1579~/^C/);else print $0}}' >> maize_chr2.txt
$ cat maize_all.txt | awk -F "\t" '{if($1576~/3/){if($1579~/^C/);else print $0}}' >> maize_chr3.txt
$ cat maize_all.txt | awk -F "\t" '{if($1576~/4/){if($1579~/^C/);else print $0}}' >> maize_chr4.txt
$ cat maize_all.txt | awk -F "\t" '{if($1576~/5/){if($1579~/^C/);else print $0}}' >> maize_chr5.txt
$ cat maize_all.txt | awk -F "\t" '{if($1576~/6/){if($1579~/^C/);else print $0}}' >> maize_chr6.txt
$ cat maize_all.txt | awk -F "\t" '{if($1576~/7/){if($1579~/^C/);else print $0}}' >> maize_chr7.txt
$ cat maize_all.txt | awk -F "\t" '{if($1576~/8/){if($1579~/^C/);else print $0}}' >> maize_chr8.txt
$ cat maize_all.txt | awk -F "\t" '{if($1576~/9/){if($1579~/^C/);else print $0}}' >> maize_chr9.txt
$ cat maize_all.txt | awk -F "\t" '{if($1576~/10/){if($1579~/^C/);else print $0}}' >> maize_chr10.txt
$ wc -l maize_chr1.txt maize_chr2.txt maize_chr3.txt maize_chr4.txt maize_chr5.txt maize_chr6.txt maize_chr7.txt maize_chr8.txt maize_chr9.txt maize_chr10.txt maize_multi.txt maize_unknown.txt
$ awk 'NR == 1; NR > 1 {print $0 | "sort -k1577V"}' maize_chr1.txt > maize_chr1_sort.txt
$ awk 'NR == 1; NR > 1 {print $0 | "sort -k1577V"}' maize_chr2.txt > maize_chr2_sort.txt
$ awk 'NR == 1; NR > 1 {print $0 | "sort -k1577V"}' maize_chr3.txt > maize_chr3_sort.txt
$ awk 'NR == 1; NR > 1 {print $0 | "sort -k1577V"}' maize_chr4.txt > maize_chr4_sort.txt
$ awk 'NR == 1; NR > 1 {print $0 | "sort -k1577V"}' maize_chr5.txt > maize_chr5_sort.txt
$ awk 'NR == 1; NR > 1 {print $0 | "sort -k1577V"}' maize_chr6.txt > maize_chr6_sort.txt
$ awk 'NR == 1; NR > 1 {print $0 | "sort -k1577V"}' maize_chr7.txt > maize_chr7_sort.txt
$ awk 'NR == 1; NR > 1 {print $0 | "sort -k1577V"}' maize_chr8.txt > maize_chr8_sort.txt
$ awk 'NR == 1; NR > 1 {print $0 | "sort -k1577V"}' maize_chr9.txt > maize_chr9_sort.txt
$ awk 'NR == 1; NR > 1 {print $0 | "sort -k1577V"}' maize_chr10.txt > maize_chr10_sort.txt
$ awk 'NR <= 3; NR > 3 {print $0 | "sort -rnk1577"}' maize_chr1.txt >  maize_chr1_sort_r.txt
$ awk 'NR <= 3; NR > 3 {print $0 | "sort -rnk1577"}' maize_chr2.txt >  maize_chr2_sort_r.txt
$ awk 'NR <= 3; NR > 3 {print $0 | "sort -rnk1577"}' maize_chr3.txt >  maize_chr3_sort_r.txt
$ awk 'NR <= 3; NR > 3 {print $0 | "sort -rnk1577"}' maize_chr4.txt >  maize_chr4_sort_r.txt
$ awk 'NR <= 3; NR > 3 {print $0 | "sort -rnk1577"}' maize_chr5.txt >  maize_chr5_sort_r.txt
$ awk 'NR <= 3; NR > 3 {print $0 | "sort -rnk1577"}' maize_chr6.txt >  maize_chr6_sort_r.txt
$ awk 'NR <= 3; NR > 3 {print $0 | "sort -rnk1577"}' maize_chr7.txt >  maize_chr7_sort_r.txt
$ awk 'NR <= 3; NR > 3 {print $0 | "sort -rnk1577"}' maize_chr8.txt >  maize_chr8_sort_r.txt
$ awk 'NR <= 3; NR > 3 {print $0 | "sort -rnk1577"}' maize_chr9.txt >  maize_chr9_sort_r.txt
$ awk 'NR <= 3; NR > 3 {print $0 | "sort -rnk1577"}' maize_chr10.txt >  maize_chr10_sort_r.txt
$ sed -i "s/\?/-/g" maize_chr*_sort_r.txt
```

Find out the columns that record chromosome information and multiposition information.
Then creat and add headers to separate chromosome files.
Separate SNPs into each file based on chromosome and multiposition.
Sort each file based on position.
Sort file in reverse order based on position and replace `?` by `-`.

### Teosinte Data

```linux
$ head -n 1 teosinte_all.txt | awk  -F "\t" '{for(i=1;i<=NF;i++){if($i~/Chromosome/)print i}}'
978
$ head -n 1 teosinte_all.txt | awk  -F "\t" '{for(i=1;i<=NF;i++){if($i~/mult_positions/)print i}}'
981
$ head -n 3 teosinte_all.txt | tee teosinte_chr1.txt teosinte_chr2.txt teosinte_chr3.txt teosinte_chr4.txt teosinte_chr5.txt teosinte_chr6.txt teosinte_chr7.txt teosinte_chr8.txt teosinte_chr9.txt teosinte_chr10.txt teosinte_multi.txt teosinte_unknown.txt
$ cat teosinte_all.txt | awk -F "\t" '{if($981~/^C/)print $0}' >> teosinte_multi.txt
$ cat teosinte_all.txt | awk -F "\t" '{if($978~/^u/)print $0}' >> teosinte_unknown.txt
$ cat teosinte_all.txt | awk -F "\t" '{if($978~/^1$/){if($981~/^C/);else print $0}}' >> teosinte_chr1.txt
$ cat teosinte_all.txt | awk -F "\t" '{if($978~/2/){if($981~/^C/);else print $0}}' >> teosinte_chr2.txt
$ cat teosinte_all.txt | awk -F "\t" '{if($978~/3/){if($981~/^C/);else print $0}}' >> teosinte_chr3.txt
$ cat teosinte_all.txt | awk -F "\t" '{if($978~/4/){if($981~/^C/);else print $0}}' >> teosinte_chr4.txt
$ cat teosinte_all.txt | awk -F "\t" '{if($978~/5/){if($981~/^C/);else print $0}}' >> teosinte_chr5.txt
$ cat teosinte_all.txt | awk -F "\t" '{if($978~/6/){if($981~/^C/);else print $0}}' >> teosinte_chr6.txt
$ cat teosinte_all.txt | awk -F "\t" '{if($978~/7/){if($981~/^C/);else print $0}}' >> teosinte_chr7.txt
$ cat teosinte_all.txt | awk -F "\t" '{if($978~/8/){if($981~/^C/);else print $0}}' >> teosinte_chr8.txt
$ cat teosinte_all.txt | awk -F "\t" '{if($978~/9/){if($981~/^C/);else print $0}}' >> teosinte_chr9.txt
$ cat teosinte_all.txt | awk -F "\t" '{if($978~/10/){if($981~/^C/);else print $0}}' >> teosinte_chr10.txt
$ wc -l teosinte_chr1.txt teosinte_chr2.txt teosinte_chr3.txt teosinte_chr4.txt teosinte_chr5.txt teosinte_chr6.txt teosinte_chr7.txt teosinte_chr8.txt teosinte_chr9.txt teosinte_chr10.txt teosinte_multi.txt teosinte_unknown.txt
$ awk 'NR == 1; NR > 1 {print $0 | "sort -k979V"}' teosinte_chr1.txt > teosinte_chr1_sort.txt
$ awk 'NR == 1; NR > 1 {print $0 | "sort -k979V"}' teosinte_chr2.txt > teosinte_chr2_sort.txt
$ awk 'NR == 1; NR > 1 {print $0 | "sort -k979V"}' teosinte_chr3.txt > teosinte_chr3_sort.txt
$ awk 'NR == 1; NR > 1 {print $0 | "sort -k979V"}' teosinte_chr4.txt > teosinte_chr4_sort.txt
$ awk 'NR == 1; NR > 1 {print $0 | "sort -k979V"}' teosinte_chr5.txt > teosinte_chr5_sort.txt
$ awk 'NR == 1; NR > 1 {print $0 | "sort -k979V"}' teosinte_chr6.txt > teosinte_chr6_sort.txt
$ awk 'NR == 1; NR > 1 {print $0 | "sort -k979V"}' teosinte_chr7.txt > teosinte_chr7_sort.txt
$ awk 'NR == 1; NR > 1 {print $0 | "sort -k979V"}' teosinte_chr8.txt > teosinte_chr8_sort.txt
$ awk 'NR == 1; NR > 1 {print $0 | "sort -k979V"}' teosinte_chr9.txt > teosinte_chr9_sort.txt
$ awk 'NR == 1; NR > 1 {print $0 | "sort -k979V"}' teosinte_chr10.txt > teosinte_chr10_sort.txt
$ awk 'NR <= 3; NR > 3 {print $0 | "sort -rnk979"}' teosinte_chr1.txt >  teosinte_chr1_sort_r.txt
$ awk 'NR <= 3; NR > 3 {print $0 | "sort -rnk979"}' teosinte_chr2.txt >  teosinte_chr2_sort_r.txt
$ awk 'NR <= 3; NR > 3 {print $0 | "sort -rnk979"}' teosinte_chr3.txt >  teosinte_chr3_sort_r.txt
$ awk 'NR <= 3; NR > 3 {print $0 | "sort -rnk979"}' teosinte_chr4.txt >  teosinte_chr4_sort_r.txt
$ awk 'NR <= 3; NR > 3 {print $0 | "sort -rnk979"}' teosinte_chr5.txt >  teosinte_chr5_sort_r.txt
$ awk 'NR <= 3; NR > 3 {print $0 | "sort -rnk979"}' teosinte_chr6.txt >  teosinte_chr6_sort_r.txt
$ awk 'NR <= 3; NR > 3 {print $0 | "sort -rnk979"}' teosinte_chr7.txt >  teosinte_chr7_sort_r.txt
$ awk 'NR <= 3; NR > 3 {print $0 | "sort -rnk979"}' teosinte_chr8.txt >  teosinte_chr8_sort_r.txt
$ awk 'NR <= 3; NR > 3 {print $0 | "sort -rnk979"}' teosinte_chr9.txt >  teosinte_chr9_sort_r.txt
$ awk 'NR <= 3; NR > 3 {print $0 | "sort -rnk979"}' teosinte_chr10.txt >  teosinte_chr10_sort_r.txt
$ sed -i "s/\?/-/g" teosinte_chr*_sort_r.txt
```

Same with maize data, but work on teosinte data.

### Organize Files

I moved maize files and teosinte files into separate folders, and the final files are in the folder called `result` inside both `maize` and `teosinte` folders.


