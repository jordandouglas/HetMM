# Load required libraries
if (!require(jsonlite, quietly = T)) install.packages("jsonlite")
library(jsonlite)





# Define the input CSV file path and output JSON file path
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 1){
	stop(paste("Use: Rscript tsv2json.R <csvFile>"))
}
tsv_file <- args[1]
json_file <- gsub("[.]csv", ".json", tsv_file)



# Column numbers
column1_name <- 1
column2_name <- 2

# Read the TSV file into a data frame
data <- read.delim(tsv_file, sep = ",", header = TRUE)

# Create an empty list to store the XML data
xml_data = character(0)

# Loop through each row of the data frame
for (i in 1:nrow(data)) {
  # Extract the values from the specified columns
  variable1 <- data[i, column1_name]
  variable2 <- data[i, column2_name]


  # Ensure that rate is 0 when substrate is 0
  if (variable1 == 0){
  	variable2 = 0
  }
  
  # Create an XML snippet with the extracted values
  xml_snippet <- sprintf('<observation spec="hetmm.KineticDataPoint" reactant="X" a="%s" v="%s" />', variable1, variable2)
  
  # Append the XML snippet to the list
  xml_data = c(xml_data, xml_snippet)
}

xml_data = paste(xml_data, collapse="\n")
JSON = list()
JSON[["data"]] = xml_data
JSON[["name"]] = gsub(".+/", "", gsub("[.]tsv", "", tsv_file))

# Convert the list of XML snippets to a JSON array
json_array <- toJSON(JSON, auto_unbox = TRUE, pretty=TRUE)

# Save the JSON data to a file
cat(paste("Saving to", json_file, "\n"))
writeLines(json_array, json_file)

