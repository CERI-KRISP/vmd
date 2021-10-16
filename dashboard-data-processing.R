#Process and update dashboard db

#Load libraries
library(readr)
library(dplyr)
library(DBI)
library(sqldf)
library(data.table)
library(tidyr)
library(writexl)
library(readxl)

## Import nextclade
sa_nextclade_210731 <- read_delim("~/temp/Collaborations/covid/wave3_data/SA_gisaid/sa_gisaid_auspice_input_hcov-19_2021_03_01-2021-07-31/sa.sequences-nextclade.200330-210731.tsv", 
                                                   "\t", escape_double = FALSE, trim_ws = TRUE)
sa_nextclade_2108 <- read_delim("~/temp/Collaborations/covid/wave3_data/SA_gisaid/sa_gisaid_auspice_input_hcov-19_2021_03_01-2021-07-31/sa.sequences-nextclade.210801-210831.tsv", 
                                                   "\t", escape_double = FALSE, trim_ws = TRUE)

## clades
df_nc_clades_210731 <- sa_nextclade_210731 %>% select(seqName, clade)
df_nextclade_2108 <- sa_nextclade_2108 %>% select(seqName, clade)
df_nc_clades <- rbind(df_nc_clades_210731, df_nextclade_2108)

## Mutations


## Pangolin
sa_seqs_210731_pangolin <- read_csv("~/temp/Collaborations/covid/wave3_data/SA_gisaid/sa_gisaid_auspice_input_hcov-19_2021_03_01-2021-07-31/sa.sequences.200330-210731.pangolin.csv")
sa_seqs_210831_pangolin <- read_csv("~/temp/Collaborations/covid/wave3_data/SA_gisaid/sa_gisaid_auspice_input_hcov-19_2021_03_01-2021-07-31/sa.sequences.210801-210831.pangolin.csv")
sa_seqs_pangolin <- rbind(sa_seqs_210731_pangolin, sa_seqs_210831_pangolin)

dim(sa_seqs_pangolin)#

## import PHE/PHA variants file
variants <- read_delim("variants.tsv", "\t",escape_double = FALSE, trim_ws = TRUE)

## import augur metadata
sa_metadata_200301_210228 <- read_delim("~/temp/Collaborations/covid/wave3_data/SA_gisaid/sa_gisaid_auspice_input_hcov-19_2021_03_01-2021-07-31/sa.metadata.200301-210228.tsv", 
                                        "\t", escape_double = FALSE, trim_ws = TRUE)
sa_metadata_210301_210615 <- read_delim("~/temp/Collaborations/covid/wave3_data/SA_gisaid/sa_gisaid_auspice_input_hcov-19_2021_03_01-2021-07-31/sa.metadata.210301-210615.tsv", 
                                        "\t", escape_double = FALSE, trim_ws = TRUE)
sa_metadata_210616_210731 <- read_delim("~/temp/Collaborations/covid/wave3_data/SA_gisaid/sa_gisaid_auspice_input_hcov-19_2021_03_01-2021-07-31/sa.metadata.210616-210731.tsv", 
                                        "\t", escape_double = FALSE, trim_ws = TRUE)
sa_metadata_210801_210831 <- read_delim("~/temp/Collaborations/covid/wave3_data/SA_gisaid/sa_gisaid_auspice_input_hcov-19_2021_03_01-2021-07-31/sa.metadata.210801-210831.tsv", 
                                        "\t", escape_double = FALSE, trim_ws = TRUE)
sa_metadata <- rbind(sa_metadata_200301_210228, sa_metadata_210301_210615, sa_metadata_210616_210731, sa_metadata_210801_210831)

## SA sanger sequences,
sanger_sequenced <- read_delim("~/temp/Collaborations/covid/wave3_data/SA_gisaid/sa_gisaid_auspice_input_hcov-19_2021_03_01-2021-07-31/sanger_sequenced.tsv", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
sanger_sequenced <- sanger_sequenced %>% mutate(
  coverage = round(((29903-as.numeric(totalMissing))/29903)*100,0)
)
write_xlsx(sanger_sequenced, "../assembly_auto/merge_sanger/manuscript/submission/Table S3.xlsx")

sa_metadata_sanger <- sa_metadata %>% filter(strain %in% sanger_sequenced$seqName)
write_xlsx(sa_metadata_sanger, "../assembly_auto/merge_sanger/manuscript/submission/Table S2.xlsx")

### Master table
sa_metadata <- sa_metadata %>% left_join(df_nc_clades, by = c("strain"="seqName"))
sa_metadata <- sa_metadata %>% left_join(sa_seqs_pangolin[, c("taxon","lineage")], by = c("strain"="taxon"))
sa_metadata <- sa_metadata %>% mutate(variant = clade,
                                      variant = ifelse(lineage == "A.23.1", "VUI-202102/01", variant),
                                      variant = ifelse(lineage == "B.1.1.318", "VUI-21FEB-04", variant),
                                      variant = ifelse(lineage == "C.1", "VUI-21MAY-02", variant),
                                      variant = ifelse(lineage == "B.1.1.7", "Alpha", variant),
                                      variant = ifelse(clade == "20H (Beta, V2)", "Beta", variant),
                                      variant = ifelse(clade == "21A (Delta)", "Delta+", variant),
                                      variant = ifelse(lineage == "C.1.2", "VUI-26JUN-01", variant),
                                      variant = ifelse(clade == "21D (Eta)", "Eta", variant),
                                      week = (Sys.Date() - 30))



sa_metadata_variants <- sa_metadata %>% filter(variant %in% variants$name) %>% select("strain", "date", week, date_submitted, "clade", "lineage", variant)
               
sa_metadata_variants$date <- as.character(sa_metadata_variants$date)
sa_metadata_variants$date_submitted <- as.character(sa_metadata_variants$date_submitted)
sa_metadata_variants$week <- as.character(sa_metadata_variants$week)

#sa_variant_counts <- sqldf('SELECT   variant, count(*), min(date) as mrca 
#                           FROM sa_metadata_variants GROUP BY variant')

#sa_variant_counts_new <- sqldf('SELECT variant, date(date) date, date(week) as week, COUNT(variant) 
#                           FROM sa_metadata_variants WHERE date >= week GROUP BY variant')

sa_variant_info <- sqldf('SELECT 
                            * 
                          FROM 
                    (SELECT variant as name, count(variant) as "number of cases", min(date) as "first case" FROM sa_metadata_variants GROUP BY name) as a  LEFT JOIN 
                    (SELECT variant as name, COUNT(variant) as "new cases in last 1 month" FROM sa_metadata_variants  WHERE date_submitted >= week GROUP BY name) as b ON a.name = b.name')
sa_variant_info <- sa_variant_info[, -4] 
sa_variant_info <- sa_variant_info %>% left_join(variants, by = "name" )

#, count(*) 

### Create database
# Create an RSQLite database
con <- dbConnect(RSQLite::SQLite(), "dashboardDB")

dbWriteTable(con, "tbl_lineages", sa_seqs_pangolin)
dbWriteTable(con, "tbl_nc_clades", df_nc_clades)
dbWriteTable(con, "tbl_variants", variants) 
dbWriteTable(con, "tbl_metadata", sa_metadata)
dbWriteTable(con, "tbl_metadata_variants", sa_metadata_variants)
#dbListTables(con)
dbDisconnect(con)


## Africa variants
variants_Africa <- read_delim("variants_Africa.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

#Africa_df <- read_excel("Africa_all_data_15September_gooddates.xlsx")
Africa_df <- read_excel("Africa_all_data_14October2021.xlsx")

sadc_countries <- c("Angola", "Botswana", "Democratic Republic of the Congo", "Eswatini", "Lesotho", "Madagascar", "Malawi", "Mauritius", "Mozambique", "Namibia", "South Africa", "Union of the Comoros", "Zambia", "Zimbabwe")

Africa_df <- Africa_df %>% mutate(variant = Nextstrain_clade,
                                   variant = ifelse(pango_lineage == "A.23.1", "A.23.1", variant),
                                   variant = ifelse(pango_lineage == "B.1.1.318" | pango_lineage %like% "AZ", "B.1.1.318", variant),
                                   variant = ifelse(pango_lineage == "C.1", "C.1", variant),
                                   variant = ifelse(Nextstrain_clade == "20I (Alpha, V1)", "Alpha", variant),
                                   variant = ifelse(Nextstrain_clade == "20H (Beta, V2)", "Beta", variant),
                                   variant = ifelse(Nextstrain_clade == "21A (Delta)", "Delta", variant),
                                   variant = ifelse(pango_lineage == "C.1.2", "C.1.2", variant),
                                   variant = ifelse(Nextstrain_clade == "21D (Eta)", "Eta", variant),
                                   variant = ifelse(pango_lineage == "A.23.1", "A.23.1", variant),
                                   variant = ifelse(pango_lineage == "C.36.3", "C.36.3", variant),
                                   #week = (Sys.Date() - 7),
                                   week = (as.Date("2021-10-14") - 7),
                                   #month = (Sys.Date() - 30),
                                   month = (as.Date("2021-10-14") - 30),
                                   sadac = ifelse(country %in% sadc_countries, 1, 0))


Africa_df <- Africa_df %>% filter(variant %in% variants_Africa$name) %>% 
  select("strain", "date", week, month, date_submitted, Nextstrain_clade, pango_lineage, variant, country, sadac)

Africa_df$date <- as.character(Africa_df$date)
Africa_df$date_submitted <- as.character(Africa_df$date_submitted)
Africa_df$week <- as.character(Africa_df$week)
Africa_df$month <- as.character(Africa_df$month)

con <- dbConnect(RSQLite::SQLite(), "dashboardDB")
dbWriteTable(con, "tbl_metadata_Africa", Africa_df, overwrite = T)
dbDisconnect(con)

