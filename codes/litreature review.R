source("https://bioconductor.org/biocLite.R")
biocLite("EBImage")
# then load metagear
library(metagear)
library(bib2df)

data(example_references_metagear)






# load a bibliographic dataset with the authors, titles, and abstracts of multiple study references
xyz <- bib2df("savedrecs.bib", separate_names = c(FALSE))
xyz = as.matrix(xyz)
# display the bibliographic variables in this dataset
xyzs["JOURNAL"]
theRefs <- effort_initialize(xyz)
names(theRefs)
theTeam <- c("ravi", "joerg")
theRefs_unscreened <- effort_distribute(theRefs, reviewers = theTeam)
theRefs_unscreened[c("STUDY_ID", "REVIEWERS")]
theRefs_unscreened <- effort_distribute(theRefs, reviewers = theTeam, effort = c(20, 80))
theRefs_unscreened[c("STUDY_ID", "REVIEWERS")]
theRefs_unscreened <- effort_distribute(theRefs, reviewers = theTeam, effort = c(20, 80), save_split = TRUE)
theRefs_unscreened[c("STUDY_ID", "REVIEWERS")]
list.files(pattern = "effort")
theRefs_screened <- effort_merge()
theRefs_screened[c("STUDY_ID", "REVIEWERS", "INCLUDE")]
theSummary <- effort_summary(theRefs_screened)
write.table(xyzs, file = "MyDatae.csv")
abstract_screener("effort_joerg.csv", aReviewer = "joerg")
