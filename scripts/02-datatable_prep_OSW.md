Datatable preparation eDNA metabarcoding data OSW project
================

This script takes your Blast output from the GMGI database, Mitofish
database, and NCBI database to create one datatable with read counts and
taxonomic assignment.

**Workflow summary:**  
1. Load libraries  
2. Load metadata  
3. Load BLAST output from GMGI, Mitofish, and NCBI  
4. Load DADA2 ASV Table  
5. Taxonomic Assignment 6. ASV read count filtering 7. Collapsing read
counts by species 8. Exporting results df

# Load libraries

``` r
library(ggplot2) ## for plotting
library(dplyr) ## for data table manipulation
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyr) ## for data table manipulation
library(readr) ## for reading in tsv files
library(readxl) ## for reading in excel files
library(stringr) ## for data transformation
library(strex) ## for data transformation
library(writexl) ## for excel output
library(purrr) ## for data transformation
library(funrar) ## for make_relative()
library(tidyverse) ## for data transformation
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ forcats   1.0.0     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

# Metadata input

## Identify paths for metadata and project data

Each user needs to write in their specific directory outputs prior to
the file name.

Change wkdir?

``` r
### User edits:
### 1. change paths of input and output as desired 

## GMGI Fish database
path_GMGIdb = "C:/BoxDrive/Box/Science/Fisheries/Projects/eDNA/Metabarcoding Lab Resources/Reference Databases/GMGI_Vert_Ref.xlsx"
path_fishbase_tax = "C:/BoxDrive/Box/Science/Fisheries/Projects/eDNA/Metabarcoding Lab Resources/Reference Databases/taxonomic_classification_fishbase.csv"
path_mitofish_tax = "C:/BoxDrive/Box/Science/Fisheries/Projects/eDNA/Metabarcoding Lab Resources/Reference Databases/taxonomic_classification_mitofish.csv"

## BLAST results
path_blast_gmgi = "data/BLASToutput/BLASTResults_GMGI.txt"
path_blast_mito = "data/BLASToutput/BLASTResults_Mito.txt"
path_blast_ncbi_taxassigned = "data/BLASToutput/NCBI_taxassigned.txt"
path_blast_ncbi = "data/BLASToutput/BLASTResults_NCBI.txt"

## ASV table results 
## confirm that the ASV_table.len.tsv name is correct for user's project
path_asv_table = "data/ASV_table.len.tsv"

# output paths 
path_choice_required = "data/BLASToutput/Taxonomic_assignment/Choice_required_GMGI_multiplehits.xlsx"
path_disagree_list = "data/BLASToutput/Taxonomic_assignment/SampleReport_taxonomic_ID.xlsx"
```

## Load project metadata

Metadata specific to each project. This contains information about each
sample (e.g., month, site, time, sample type, etc.). Confirm that sample
IDs match those used in the ASV_table.len.tsv file.

``` r
### User edits:
### 1. change path of metadata file

### Project specific information
meta <- read.csv("data/metadata/samplesheet.csv") %>% 
  ### selecting only sample ID column (we don't need the paths of raw files here)
  dplyr::select(sampleID) %>%
  
  ## adding a sample Type column based on contents of SampleID column 
  mutate(SampleType = case_when(grepl("BK", sampleID) ~ "Blank", 
                                grepl("C", sampleID) ~ "Control",
                                TRUE ~ "Inside WEA")) %>%
  
  ## creating a new column for site
      ## this is slightly round about way to do this, but works 
  separate(sampleID, c("Month", "Site1", "Site2", "Depth"), sep = "_", remove=F) %>%
  unite(Site, Site1, Site2, sep = "-", remove=T) %>%
  
  ## editing notation of B/S to Bottom and Surface 
  mutate(Depth = case_when(grepl("B", Depth) ~ "Bottom", grepl("S", Depth) ~ "Surface", grepl("NA", Depth) ~ NA),
  ## adding lease area information based on site information       
         Lease_area = case_when(grepl("REV", Site) ~ "Revolution Wind", grepl("SF", Site) ~ "South Fork", 
                                grepl("501", Site) ~ "Vineyard Wind 1", grepl("VW", Site) ~ "Vineyard Wind 1", 
                                grepl("DI", Site) ~ NA),

  ## adding full month name 
         Month = case_when(Month == "Aug" ~ "August", Month == "Sep" ~ "September", Month == "Oct" ~ "October",
                           Month == "Nov" ~ "November", Month == "July" ~ "July"))

sample_date <- read_xlsx("data/metadata/OSW_master_samplelist.xlsx", sheet=1) %>% dplyr::select(-SampleType) %>%
  mutate(Site = gsub("_", "-", Site))
latlon <- read_xlsx("data/metadata/OSW_master_samplelist.xlsx", sheet=2) %>% 
  mutate(Site = gsub("_", "-", Site))
## add code to manipulate sample ID if needed
## change samplesheet back to metadata.csv

meta <- meta %>% left_join(., latlon, by = "Site") %>% left_join(., sample_date, by = c("Month", "Site", "Depth"))

meta %>% write_xlsx("data/metadata/full_metadata.xlsx")
```

## Load database metadata

No user edits in this section because paths have already been set above.

``` r
# Load GMGI database information (common name, species name, etc.)
gmgi_db <- read_xlsx(path_GMGIdb, sheet = 1) %>% dplyr::rename(sseqid = Ref) %>%
  ## removing > from beginning of entires within Ref column
  mutate(sseqid = gsub(">", "", sseqid))

# Create a variable that defines taxonomic assignment levels
phylo_classifications = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "species")

# Load taxonomic assignment information from fishbase and mitofish
tax_fishbase <- read.csv(path_fishbase_tax, header=T, col.names = c("Phylo", "Species_name")) %>%
   ## creating taxonomic assignment columns
  separate(Phylo, phylo_classifications, sep = ";", remove=T)

tax_mito <- read.csv(path_mitofish_tax, header=T, col.names = c("Species_name", "Phylo")) %>%
   ## creating taxonomic assignment columns
  separate(Phylo, phylo_classifications, sep = ";", remove=T)
```

    ## Warning: Expected 7 pieces. Missing pieces filled with `NA` in 6 rows [4190, 4191, 5390,
    ## 5391, 7206, 7207].

# BLAST data input

No user edits unless user changed blastn parameters from fisheries team
default.

``` r
## Setting column header names and classes
blast_col_headers = c("ASV_ID", "sseqid", "pident", "length", "mismatch", "gapopen",
                                        "qstart", "qend", "sstart", "send", "evalue", "bitscore")
blast_col_classes = c(rep("character", 2), rep("numeric", 10))
```

## GMGI database

No user edits.

``` r
Blast_GMGI <- read.table(path_blast_gmgi, header=F, col.names = blast_col_headers, colClasses = blast_col_classes) %>%
  ## blast changes spaces to hyphons so we need to change that back to match our metadata
  mutate(sseqid = gsub("-", " ", sseqid)) %>%
  ## join with GMGI database information
  left_join(., gmgi_db, by = "sseqid")

## Check how many ASVs were identified with the GMGI Database
length(unique(Blast_GMGI$ASV_ID)) 
```

    ## [1] 667

``` r
### OSW 667
```

## Mitofish database

No user edits.

``` r
Blast_Mito <- read.table(path_blast_mito, header=F, col.names = blast_col_headers, colClasses = blast_col_classes) %>%
  # renaming sseqid to species name
  dplyr::rename(Species_name = sseqid) %>%
  
  # replacing _ with spaces
  mutate(Species_name = gsub("_", " ", Species_name))
```

## NCBI database

No user edits.

``` r
NCBI_taxassigned <- read.delim2(path_blast_ncbi_taxassigned, header=F, col.names = c("staxid", "Phylo")) %>%
  ## creating taxonomic assignment columns
  separate(Phylo, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species_name"), sep = ";") %>%
  ## creating species column based on Species_name
  mutate(., species = str_after_nth(Species_name, " ", 1))

Blast_NCBI <- read.table(path_blast_ncbi, header=F,
                           col.names = c("ASV_ID", "sseqid", "sscinames", "staxid", "pident", "length", "mismatch",
                                         "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"),
                           colClasses = c(rep("character", 3), "integer", rep("numeric", 9))) %>%
  left_join(., NCBI_taxassigned, by = "staxid")
```

# Load DADA2 ASV Table

The column headers will be the Sample IDs and the first column is the
ASV ID. ASVs are given a “rank” based on sum of reads from that ASV
(pre-filtering). ‘Random’ indicates that if ASVs are tied, then the code
will randomly assign a rank for those tied. Because we don’t need an
exact rank here, ‘random’ will do for a tie-breaker.

No user edits.

OSW: 727 ASVs with 667 of those identified with GMGI database.

``` r
ASV_table <- read_tsv(path_asv_table, show_col_types = FALSE) %>%
  ## calculate the sum of all reads for each ASV
  mutate(., ASV_sum = rowSums(across(where(is.numeric)))) %>% 
  
  ## calculate a ranking based on those sum calculated above
  mutate(ASV_rank = rank(-ASV_sum, ties.method='random')) %>%
  
  ## move the sum and rank columns to after ASV_ID and arrange by rank
  relocate(c(ASV_sum,ASV_rank), .after = ASV_ID) %>% arrange((ASV_rank))

## creating list of rankings
ASV_rank_list <- ASV_table %>% dplyr::select(ASV_ID, ASV_sum, ASV_rank)
```

# Taxonomic Assignment

Identifying where NCBI, Mito, and GMGI disagree on tax assignment. With
the hierarchial approach, ASVs that match to GMGI and several other
databases will only result in GMGI assignment. By reviewing this df, we
can be sure we aren’t missing an assignment in our GMGI curated
database.

**Sub-workflow:**  
1. Identify any ASVs that contain multiple hits within the GMGI
database. 2. Identify entries that mismatch between GMGI, Mitofish, and
NCBI databases. 3. Assign taxonomy based on hierarchical approach. 4.
Edit taxonomy annotations based on mismatch table.  
5. Adjusting common name for those entries that don’t have one (from
Mito or GMGI).

## Identify any ASVs that contain multiple hits within the GMGI database

At this point, a fisheries team member needs to make choices about which
taxonomic assignment to accept.

#### Create list of those ASVs with multiple hits

No user edits.

``` r
multiple_hit_choice <- Blast_GMGI %>% group_by(ASV_ID) %>%
  ## take top percent identity hit, count the number of top hits, and filter to those with more than 1 top hit 
  slice_max(pident, n=1) %>% count() %>% filter(n>1) %>%
  
  ## adding BLAST_GMGI information with these ASVs and ASV rank and sum
  left_join(., Blast_GMGI, by = "ASV_ID") %>%
  left_join(., ASV_rank_list, by = "ASV_ID") %>%
  
  ## moving database percent ID to be next to Blast percent ID
  relocate(c(db_percent_ID, ASV_sum, ASV_rank), .after = pident); multiple_hit_choice
```

    ## # A tibble: 77 × 20
    ## # Groups:   ASV_ID [21]
    ##    ASV_ID         n sseqid pident db_percent_ID ASV_sum ASV_rank length mismatch
    ##    <chr>      <int> <chr>   <dbl> <chr>           <dbl>    <int>  <dbl>    <dbl>
    ##  1 1535a44c6…     2 3.00_…   99.1 100             22221       66    107        1
    ##  2 1535a44c6…     2 2.50_…   99.1 100             22221       66    107        1
    ##  3 1535a44c6…     2 3.00_…   98.1 100             22221       66    107        2
    ##  4 1535a44c6…     2 3.00_…   98.1 100             22221       66    107        2
    ##  5 1535a44c6…     2 3.00_…   98.1 100             22221       66    107        2
    ##  6 1535a44c6…     2 3.00_…   98.1 100             22221       66    107        2
    ##  7 1535a44c6…     2 3.00_…   98.1 100             22221       66    107        2
    ##  8 15e97dd14…     2 2.00_…   98.1 99               1847      255    108        2
    ##  9 15e97dd14…     2 2.00_…   98.1 100              1847      255    108        2
    ## 10 1b7bef208…     2 3.50_…  100   100              3409      180    106        0
    ## # ℹ 67 more rows
    ## # ℹ 11 more variables: gapopen <dbl>, qstart <dbl>, qend <dbl>, sstart <dbl>,
    ## #   send <dbl>, evalue <dbl>, bitscore <dbl>, Common_name <chr>,
    ## #   Species_name <chr>, Category <chr>, `Kingdom, Phyl, …` <lgl>

``` r
## export this data frame as excel sheet 
multiple_hit_choice %>% write_xlsx(path_choice_required)
```

Based on the output above, user needs to make some choices.

#### Choosing one of several hits.

Insert choice sheet for columns with x? instead of manually

``` r
### User edits:
### 1. Add filtering cases using the format ASV_ID == "" ~ sseqid == ""
### example: ASV_ID == "1535a44c66f8850e6d30284f8ddeb38d" ~ sseqid == "3.00_Human_mito1"
### Keep the TRUE ~ TRUE at the end

Blast_GMGI_edited <- Blast_GMGI %>% 
  ### picking one of several hits
  filter(case_when(
    ASV_ID == "15e97dd1453cc2ecb8b5082f9b514004" ~ 
      sseqid == "2.00_common eider_Somateria mollissima_OR_Bufflehead_bucephala albeola_OR_other eiders sea ducks",
    ASV_ID == "edb60a770bbe58d0ee147058b2a1b13f" ~ 
      sseqid == "2.00_common eider_Somateria mollissima_OR_Bufflehead_bucephala albeola_OR_other eiders sea ducks",
    ASV_ID == "b29a188a6e4589bbbcc9ae9c5e5570bd" ~ 
      sseqid == "2.00_common eider_Somateria mollissima_OR_Bufflehead_bucephala albeola_OR_other eiders sea ducks",
    ASV_ID == "80c040f8fa7ca320845742c52669d5e9" ~ 
      sseqid == "2.00_common eider_Somateria mollissima_OR_Bufflehead_bucephala albeola_OR_other eiders sea ducks",
    ASV_ID == "46db97c5565c94eb6c74f3e91483addd" ~ sseqid == "1.00_Northern_sand_lance_Ammodytes_dubius",
    ASV_ID == "911894fe810a42fd66c3c43018d29753" ~ sseqid == "1.00_Atlantic_menhaden_LS16_or_river_herrings_Clupeidae_sp",
    ASV_ID == "1535a44c66f8850e6d30284f8ddeb38d" ~ sseqid == "3.00_Human_mito1",
    ASV_ID == "1b7bef208071964dff913b3e58cc6deb" ~ sseqid == "3.00_Human_mito1",
    ASV_ID == "267d361028b5d98f0610470de52f3135" ~ sseqid == "3.50_Human_chromo4_AABC161701",
    ASV_ID == "282aac0ec13fc675ac82d31f2b3e70e9" ~ sseqid == "3.50_Human_chromosome1_659M14",
    ASV_ID == "7950b1078efc076defba9c936b970ef7" ~ sseqid == "3.50_Human_chromo17_RP1113L8",
    ASV_ID == "20a9b66022e64e1ed2b8c4527d0ff2ac" ~ sseqid == "2.00_Great_black_backed_gull_and_other_Larus_gulls",
    ASV_ID == "4450a6fa10b56881617cff33c5585aa8" ~ sseqid == "1.10_Mummichog_Fundulus_heteroclitus_cluster3",
    ASV_ID == "97692e77675cee377312fb71ecc51b4e" ~ sseqid == "2.00_Great Shearwater_Ardenna gravis",
    ASV_ID == "9f0c60d28f7c1bc1c462a9092bf3cda7" ~ sseqid == "2.00_Great_black_backed_gull_and_other_Larus_gulls",
    ASV_ID == "ba9d084edcf15b5ebdecedeccda0c8c7" ~ sseqid == "2.00_Great_black_backed_gull_and_other_Larus_gulls",
    ASV_ID == "c0a3f3ed23f04247d92740a9502f8b57" ~ sseqid == "2.00_Great_black_backed_gull_and_other_Larus_gulls",
    ASV_ID == "e0848e25382836806b28673b01f892a2" ~ sseqid == "2.00_Great_black_backed_gull_and_other_Larus_gulls",
    ASV_ID == "d66441776b9986f2d4a95cf5895dc693" ~ sseqid == "3.00_Human_mito1",
    ASV_ID == "d66441776b9986f2d4a95cf5895dc693" ~ sseqid == "3.50_Human_chromo17_RP1113L8b",
    ASV_ID == "d9aff8e5fde346fe4ae2aa2c21952512" ~ sseqid == "3.50_Human_chromo17_RP1113L8b",
    TRUE ~ TRUE)) %>%
  ## filter out to be unknown
  filter(!ASV_ID == "2744f56d73d9591aa2ddedd05637bb27")
```

#### Confirming that all entries have been dealth with

No user edits.

``` r
### Check the below output to confirm the filtering steps above worked (if it worked, it won't be in output)
Blast_GMGI_edited %>% group_by(ASV_ID) %>% slice_max(pident, n=1) %>% count() %>% filter(n>1)
```

    ## # A tibble: 0 × 2
    ## # Groups:   ASV_ID [0]
    ## # ℹ 2 variables: ASV_ID <chr>, n <int>

## Identify entries that mismatch between GMGI, Mitofish, and NCBI databases

Creating a df called “Disagree”. Review the output before moving onto
the next section.

No user edits. Come back to this section, it’s not including entries
that disagree because of an NA entry.

686 ASV that are identified through three databases.

``` r
# Disagree2 <- Blast_GMGI_edited %>% group_by(ASV_ID) %>% 
#   dplyr::rename(., GMGI_db_ID = db_percent_ID, GMGI_pident = pident) %>%
#   ## Creating new columns with species name based on pident information
#   mutate(
#     GMGI_100 = if_else(GMGI_pident == 100, Species_name, NA),
#     GMGI_lessthan100 = if_else(GMGI_pident < 100, Species_name, NA)) %>%
#   
#   ## taking only the top hit per ASV ID
#   slice_max(GMGI_pident, n = 1, with_ties = FALSE) %>% ungroup() %>%
# 
#   ## filtering to distinct rows with selected columns
#   distinct(ASV_ID, GMGI_db_ID, GMGI_pident, GMGI_100, GMGI_lessthan100) %>%
#   
#   ## adding Mitofish and editing the Blast_Mito df in the process
#   full_join(Blast_Mito %>% dplyr::select(ASV_ID, Species_name) %>%
#               dplyr::rename(Mitofish = Species_name) %>%
#               distinct() %>% group_by(ASV_ID) %>%
#               mutate(Mitofish = paste0(Mitofish, collapse = ";")),
#             by = "ASV_ID") %>%
#   
#   ## adding NCBI and editing the Blast_NCBI df in the process
#   full_join(Blast_NCBI %>% dplyr::select(ASV_ID, Species_name) %>%
#               dplyr::rename(NCBI = Species_name) %>%
#               distinct() %>% group_by(ASV_ID) %>%
#               mutate(NCBI = paste0(NCBI, collapse = ";")),
#             by = "ASV_ID") %>%
#   
#   ## adding ASV rank and sum information
#   left_join(., ASV_rank_list, by = "ASV_ID") %>%
# 
#   ## filtering out duplicate rows
#   distinct() %>%
#   
#   ## filtering to those entries that mismatch between GMGI, Mitofish, and NCBI
#   filter(!(GMGI_100 %in% GMGI_lessthan100) | !(GMGI_100 %in% Mitofish) | !(GMGI_100 %in% NCBI))
# 
# ## export this data frame as excel sheet 
# Disagree %>% write_xlsx(path_disagree_list)
```

## Assign taxonomy based on hierarchical approach

Taxonomic identification is taken from GMGI 100%, then GMGI \<100%, then
Mitofish 100%, and finally NCBI 100%.

No user edits.

``` r
ASV_table_taxID <- ASV_table %>% 
  
  ## 1. Top hit from GMGI's database
  left_join(Blast_GMGI_edited %>% group_by(ASV_ID) %>% slice_max(pident, n = 1) %>%
                            dplyr::select(ASV_ID, Species_name) %>% ungroup(),
            by = join_by(ASV_ID)) %>%
  
  ## 2. Mitofish database
  ### join df, select ASV_ID and Species_name columns, rename Species_name to Mito, call only distinct rows
  left_join(., Blast_Mito %>% dplyr::select(ASV_ID, Species_name) %>% dplyr::rename(Mito = Species_name) %>% distinct() %>%
              
              ### group by ASV_ID, and collapse all species names separated by ;, then take only distinct rows
              group_by(ASV_ID) %>% mutate(Mito = paste0(Mito, collapse = ";")) %>% distinct(), by = "ASV_ID") %>%
  
  ### if GMGI annotation is NA, then replace with Mitofish 
  mutate(., Species_name = ifelse(is.na(Species_name), Mito, Species_name)) %>%

  ## 3. NCBI database; same functions as above
  left_join(., Blast_NCBI %>% dplyr::select(ASV_ID, Species_name) %>% dplyr::rename(NCBI = Species_name) %>% distinct() %>%
              group_by(ASV_ID) %>% mutate(NCBI = paste0(NCBI, collapse = ";")) %>% distinct(), by = "ASV_ID") %>%
  mutate(., Species_name = ifelse(is.na(Species_name), NCBI, Species_name)) %>%
  
  ## 4. if Species name is STILL not filled, call it "unassigned"
  mutate(., Species_name = ifelse(is.na(Species_name), "unassigned", Species_name)) %>%  

  ## removing Mito spp and NCBI spp
  dplyr::select(-Mito, -NCBI) %>%
  
  ## move species name to be after ASV_ID
  relocate(., c(Species_name), .after = ASV_ID)

## confirming number of ASVs in ASV_table and ASV_table_taxID
## the output should be TRUE
nrow(ASV_table) == nrow(ASV_table_taxID)
```

    ## [1] TRUE

## Edit taxonomy annotations based on mismatch table

Override any annotations:

``` r
### User edits:
### 1. Add mutate cases using the format ASV_ID == "" ~ ""
### example: ASV_ID == "abb58e582fcc5bd9d2526b4bf98ed7a3" ~ "Ardenna griseus or Ardenna gravis",
### Keep the TRUE ~ Species_name at the end

### 2. Add ifelse() cases using the format ifelse(grepl('', Species_name), "", Species_name
### example: ifelse(grepl('Homo sapiens', Species_name), "Homo sapiens", Species_name
### example: ifelse(grepl('Phyllozelus siccus', Species_name), "unassigned", Species_name)

ASV_table_taxID <- ASV_table_taxID %>%
  ## Use mutate case_when() for specific ASV IDs
  mutate(Species_name = case_when(
    ASV_ID == "d33d1d19c918f501785c1b1c4c550b91" ~ "Myrophis punctatus",
    ASV_ID == "01f9f22b222897bd5c0372d8640c15e5" ~ "Rhomboplites aurorubens",
    TRUE ~ Species_name)) %>%

  ## OR Use ifelse() function to find entries with specific content that might apply to mulitple ASVs
  mutate(Species_name = ifelse(grepl('Homo sapiens', Species_name), "Homo sapiens", Species_name)) %>%

  ## changing all bacterium or Staphylococcus groups to bacterium and homo sapiens to human
  mutate(Species_name = ifelse(grepl('bacterium|Staphylococcus|Bordetella|Ralstonia|Nitzschia', 
                                     Species_name), "unassigned", Species_name),
         Species_name = ifelse(grepl('Homo sapiens', Species_name), "Homo sapiens", Species_name),
         Species_name = ifelse(grepl('Sinacroneuria|Peltoperla', Species_name), "unassigned", Species_name),
         Species_name = ifelse(grepl('Delphinus delphis;Lagenorhynchus albirostris;Orcinus orca', 
                                     Species_name), "unassigned", Species_name),
         Species_name = ifelse(grepl('Limanda aspera;Cleisthenes herzensteini;Nesiarchus nasutus;Hippoglossoides elassodon', Species_name), "unassigned", Species_name),
         Species_name = ifelse(grepl('Opsanus tau;Opsanus phobetron', Species_name), "unassigned", Species_name),
         Species_name = ifelse(grepl('Stenella coeruleoalba;Delphinus delphis', Species_name), "unassigned", Species_name),
         
         ## below are entries that would be taken out with future common offender list
         Species_name = ifelse(grepl('Protoginella maestratii', Species_name), "unassigned", Species_name),
         Species_name = ifelse(grepl('Podilymbus podiceps', Species_name), "unassigned", Species_name),
         Species_name = ifelse(grepl('Cololabis saira', Species_name), "unassigned", Species_name))
```

### Confirm all entries are dealt with

No user edits.

``` r
## Output will be blank
ASV_table_taxID %>% dplyr::select(Species_name) %>% distinct() %>% 
  filter(., grepl(";", Species_name)) %>% arrange(Species_name) 
```

    ## # A tibble: 0 × 1
    ## # ℹ 1 variable: Species_name <chr>

## Adjusting common name for those entries that don’t have one (from Mito or NCBI)

No user edits.

``` r
### add common name column to df
ASV_table_taxID <- ASV_table_taxID %>%
  left_join(., gmgi_db %>% dplyr::select(Species_name, Common_name, Category) %>% distinct(), by = "Species_name") %>%
  relocate(., c(Common_name, Category), .after = Species_name)

### print entries with no common name
ASV_table_taxID %>% dplyr::select(Species_name, Common_name) %>% filter(is.na(Common_name)) %>% distinct()
```

    ## # A tibble: 4 × 2
    ##   Species_name            Common_name
    ##   <chr>                   <chr>      
    ## 1 unassigned              <NA>       
    ## 2 Rhomboplites aurorubens <NA>       
    ## 3 Myrophis punctatus      <NA>       
    ## 4 Sphyrna lewini          <NA>

Editing common names and category when needed.

``` r
### User edits:
### 1. Add mutate cases using the format ifelse(grepl('', Species_name), "", Common_name
### example: ifelse(grepl('unassigned', Species_name), "unassigned", Common_name)
### 2. Add mutate cases using the format ifelse(grepl('', Species_name), "", Category

ASV_table_taxID <- ASV_table_taxID %>% 
  # changing specific entries for Common name
  mutate(Common_name = ifelse(grepl('unassigned', Species_name), "unassigned", Common_name),
         Common_name = ifelse(grepl('Rhomboplites aurorubens', Species_name), "Vermilion snapper", Common_name),
         Common_name = ifelse(grepl('Sphyrna lewini', Species_name), "Scalloped hammerhead", Common_name),
         Common_name = ifelse(grepl('Myrophis punctatus', Species_name), "Speckled worm-eel", Common_name)
         ) %>%
  
  # changing specific entries for category
  mutate(Category = ifelse(grepl('unassigned', Species_name), "unassigned", Category),
         Category = ifelse(grepl('Rhomboplites aurorubens', Species_name), "Teleost Fish", Category),
         Category = ifelse(grepl('Sphyrna lewini', Species_name), "Elasmobranch", Category),
         Category = ifelse(grepl('Myrophis punctatus', Species_name), "Teleost Fish", Category)
         )

## printing list of species name without common names 
## after additions to mutate function above, this output should be zero 
ASV_table_taxID %>% dplyr::select(Species_name, Common_name) %>% filter(is.na(Common_name)) %>% distinct()
```

    ## # A tibble: 0 × 2
    ## # ℹ 2 variables: Species_name <chr>, Common_name <chr>

# Filtering: Filter ASV by less than 0.1% reads and then collapse by group

## Filter out reads that are less than 0.1% of ASV (row) total per sample.

Create an output of what you’re losing with filtering Can sumVar be
ASV_sum

``` r
### User edits:
### 1. In mutate(sumVar = ), change FIRST_SAMPLE:LAST_SAMPLE to sample 1 column and the last sample column 
### 2. In mutate(across(.cols = )), change FIRST_SAMPLE:LAST_SAMPLE to sample 1 column and the last sample column
### example: C1_bottom:Sed_Neg2

ASV_table_taxID_filtered <- ASV_table_taxID %>%
  ## telling the df we are doing the following function by rows (ASVs)
  rowwise() %>%

  ## filtering out any values that are less than 0.001 of the total ASV read # in each sample
  mutate(across(.cols = (7:ncol(.)),
                .fns = ~ ifelse((.x/ASV_sum)<0.001, NA, .x))) %>% ungroup()
```

# Collapsing read counts by species name

Col 7: last col?

``` r
### User edits:
### 1. In mutate(sumVar = ), change FIRST_SAMPLE:LAST_SAMPLE to sample 1 column and the last sample column 
### example: C1_bottom:Sed_Neg2

## Expected number of species from all our ASVs
length(unique(ASV_table_taxID_filtered$Species_name))
```

    ## [1] 88

``` r
ASV_table_taxID_collapsed <- ASV_table_taxID_filtered %>% 
  # removing original ASV_ID to collapse
  dplyr::select(-ASV_ID) %>%  
  
  ## group by Species_name and sample
  group_by(Species_name, Common_name, Category) %>%
  
  ## sum down column by species name and sample to collapse
  summarise(across(c(Aug_501_1_B:Sep_VW_S7_S), ~ sum(., na.rm = TRUE))) %>% ungroup()
```

    ## `summarise()` has grouped output by 'Species_name', 'Common_name'. You can
    ## override using the `.groups` argument.

``` r
## conforming that value is as expected
nrow(ASV_table_taxID_collapsed)
```

    ## [1] 88

# Creating results output with metadata

Raw reads

``` r
## Raw reads matrix (wide format)
ASV_table_taxID_collapsed %>% write_xlsx("data/results/Rawreads_matrix.xlsx")

### User edits:
### 1. In gather("sampleID", "reads", ..), change FIRST_SAMPLE:LAST_SAMPLE to sample 1 column and the last sample column 
### example: C1_bottom:Sed_Neg2

## Raw reads long format and filtering out entries with zero reads
ASV_table_taxID_collapsed %>% 
  gather("sampleID", "reads", c(Aug_501_1_B:Sep_VW_S7_S)) %>%
  left_join(., meta, by = "sampleID") %>%
  filter(reads > 0) %>% write_xlsx("data/results/Rawreads_longformat.xlsx")
```

    ## Warning in left_join(., meta, by = "sampleID"): Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 4753 of `x` matches multiple rows in `y`.
    ## ℹ Row 1 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.

Relative Abundance of reads (per sample)

``` r
relative_ab <- ASV_table_taxID_collapsed %>%
  ## calculating relative abundance
  mutate(across(c(Aug_501_1_B:Sep_VW_S7_S), ~ .x / sum(.x, na.rm = TRUE))) %>%
  ## removing samples with no left left (sum = 0)
  select(where(~ !any(is.nan(.x))))

relative_ab %>% write_xlsx("data/results/Relative_abundance_matrix.xlsx")

relative_ab %>%
  gather("sampleID", "rel_ab", c(Aug_501_1_B:Sep_VW_S7_S)) %>%
  left_join(., meta, by = "sampleID") %>% write_xlsx("data/results/Relative_abundance_longformat.xlsx")
```

    ## Warning in left_join(., meta, by = "sampleID"): Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 4753 of `x` matches multiple rows in `y`.
    ## ℹ Row 1 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.