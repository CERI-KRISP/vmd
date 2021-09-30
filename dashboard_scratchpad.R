## import PHE/PHA variants file
#variants <- read_delim("variants_Africa.tsv", "\t",escape_double = FALSE, trim_ws = TRUE)

# dbFetch(dbSendQuery(con, 'SELECT country, count(variant) as "number of sequences", min(date) as "date of first sequence", max(date) as "date of last sequence", Cast ((JulianDay() - JulianDay(max(date))) As Integer) AS "No of days since last sampling" FROM tbl_metadata_Africa WHERE sadac == 1 AND variant == "Beta"  GROUP BY variant, country'))

# dbFetch(dbSendQuery(con, 'SELECT country, variant,  COUNT(variant) as "new sequences in the last 30 days" FROM tbl_metadata_Africa  WHERE date >= month AND sadac == 1 AND variant == "Beta"  GROUP BY variant, country'))
