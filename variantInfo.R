library(ghql)
library(jsonlite)

# create GraphQL client
cli <- GraphqlClient$new(
  url = "https://genetics-api.opentargets.io/graphql"
)

var <- "10_102877751_T_C"

qry <- Query$new()

qry$query("my_query", paste('query {variantInfo(variantId: "', var, '") {id, rsId, idB37, altAllele, gnomadNFE, gnomadFIN}}', sep=""))

res <- fromJSON(cli$exec(qry$queries$my_query), flatten = TRUE)

print(paste(res$data$variantInfo$id, res$data$variantInfo$rsId, res$data$variantInfo$idB37, res$data$variantInfo$altAllele, res$data$variantInfo$gnomadNFE, res$data$variantInfo$gnomadFIN, sep=" "))
