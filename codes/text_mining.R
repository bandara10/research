# read this page: eight2late.wordpress.com

# corpus means a collection of documents.
#myfiles <- list.files(path = "C:/Users/uqrdissa/Desktop/confirm docs", pattern = "txt",  full.names = TRUE)

#https://cran.r-project.org/web/packages/tm/vignettes/tm.pdf
DP <- DirSource("C:/Users/uqrdissa/Desktop/confirm docs")
ovid <-VCorpus(DP)
reuters <-ovid[43]
inspect(ovid[43])
meta(ovid[[43]], "id")
inspect(ovid[[43]])
lapply(ovid[43], as.character)


reuters <- tm_map(reuters, stripWhitespace)
reuters <- tm_map(reuters, removeWords, stopwords("english"))
tm_map(reuters, stemDocument)
dtm <- DocumentTermMatrix(reuters)
inspect(dtm[1:1, 1:1])
findFreqTerms(dtm, 10)
findAssocs(dtm, "subtype", 0.6)
inspect(removeSparseTerms(dtm, 0.4))


files <- list.files("C:\\Users\\uqrdissa\\Desktop\\confirm docs",pattern = "pdf$") # read all txt files
###
DP <- DirSource("C:\\Users\\uqrdissa\\Desktop\\confirm docs", pattern = "txt$") # has three txt files
ovid <-VCorpus(DP)
reuters <-ovid[1:3]
inspect(ovid[1:2])
meta(ovid[[1]], "id")
inspect(ovid[[1]])
DR <- lapply(ovid[1], as.character)
DR$`Confirmation Document_Suman_S43976371.txt`

###Text Mining Infrastructure in R | Journal of Statistical Software
#tdm <- TermDocumentMatrix(reuters, list(stemming = TRUE, stopwords = TRUE)) #create a term-document matrix, activate stemming and remove stopwords.
reuters <- tm_map(reuters, stripWhitespace)
reuters <- tm_map(reuters, removeWords, stopwords("english"))
tm_map(reuters, stemDocument)
dtm <- DocumentTermMatrix(reuters)
#plot(tdm, terms = findFreqTerms(tdm, lowfreq = 6)[1:25], corThreshold = 0.5) # need  package 'Rgraphviz. but not avilable.
inspect(dtm[1:3, 1:3])
findFreqTerms(dtm, 50) # find terms used 50 times in documents.

findAssocs(dtm, "risk", 0.9)
inspect(removeSparseTerms(dtm, 0.98))
tagPOS(reuters[[3]])

library(openNLP)
library(NLP)
tagPOS <-  function(x, ...) {
  s <- as.String(x)
  word_token_annotator <- Maxent_Word_Token_Annotator()
  a2 <- Annotation(1L, "sentence", 1L, nchar(s))
  a2 <- annotate(s, word_token_annotator, a2)
  a3 <- annotate(s, Maxent_POS_Tag_Annotator(), a2)
  a3w <- a3[a3$type == "word"]
  POStags <- unlist(lapply(a3w$features, `[[`, "POS"))
  POStagged <- paste(sprintf("%s/%s", s[a3w], POStags), collapse = " ")
  list(POStagged = POStagged, POStags = POStags)
}


##### https://eight2late.wordpress.com/2015/05/27/a-gentle-introduction-to-text-mining-using-r/
#Create Corpus
docs <- Corpus(DirSource("C:\\Users\\uqrdissa\\Desktop\\confirm docs", pattern = "txt$"))
docs <-ovid[3] # select 3rd document or go tot next step to see all 3 docs.
#Type in docs to see some information about the newly created corpus:
docs
#The summary() function gives more details, including a complete listing of files
summary(docs)
#inspect a particular document#prints the entire content of 3rd document in the corpus to the console
writeLines(as.character(docs[[1]])) # we chnage this to 1 document.
writeLines(as.character(docs[[3]]))
#The tm package offers a number of transformations that ease the tedium of cleaning data. To see the available transformations  type getTransformations() at the R prompt:
getTransformations()
#Here is the R code to build the content transformer, which  we will call toSpace:

#create the toSpace content transformer
toSpace <- content_transformer(function(x, pattern) {return (gsub(pattern, " ", x))})
#Now we can use  this content transformer to eliminate colons and hypens like so
docs <- tm_map(docs, toSpace, " ")
docs <- tm_map(docs, toSpace, " ")
#now apply the removePunctuation transformation.
#Remove punctuation - replace punctuation marks with " "
docs <- tm_map(docs, removePunctuation)
#Transform to lower case (need to wrap in content_transformer)
docs <- tm_map(docs,content_transformer(tolower))
#Strip digits (std transformation, so no need for content_transformer)
docs <- tm_map(docs, removeNumbers)
#remove stopwords using the standard list in tm
docs <- tm_map(docs, removeWords, stopwords("english"))
#remove all extraneous whitespaces using the stripWhitespace transformation:
#Strip whitespace (cosmetic?)
docs <- tm_map(docs, stripWhitespace)
#To see what stemming does, let's take a look at the  last few lines  of the corpus before and after stemming
#Here's what the last bit looks  like prior to stemming
writeLines(as.character(docs[[1]])) # for 1 document
writeLines(as.character(docs[[3]]))
#Now let's stem the corpus and reinspect it.
#load library
library(SnowballC)
#Stem document
docs <- tm_map(docs,stemDocument)
writeLines(as.character(docs[[1]])) # for 1 documents
writeLines(as.character(docs[[3]]))
#The next step in the process is the creation of the document term matrix  (DTM)- a matrix that lists all occurrences of words in the corpus, by document.
dtm <- DocumentTermMatrix(docs)
dtm
#This is a 3 x 3471 dimension matrix in which 57% of the rows are zero.it is clear that the large majority of words will appear only in a few documents. As a result a DTM is invariably sparse - that is, a large number of its entries are 0.
#To limit the information displayed in dtm, one can inspect a small section of it like so:
inspect(dtm[1,50:60]) # for 1 document
inspect(dtm[1:2,100:105])
#Mining the corpus
#Notice that in constructing the TDM, we have converted a corpus of text into a mathematical object that 
#can be analysed using quantitative techniques of matrix algebra.  It should be no surprise, 
#therefore, that the TDM (or DTM) is the starting point for quantitative text analysis.
#get the frequency of occurrence of each word in the corpus, we simply sum over all rows to give column sums:
freq <- colSums(as.matrix(dtm))
#first converted the TDM into a mathematical matrix using the as.matrix() function. 
#We have then summed over all rows to give us the totals for each column (term).
#length should be total number of terms
length(freq)
#create sort order (descending)
ord <- order(freq,decreasing=TRUE)
#Then list the most and least frequently occurring terms:
#inspect most frequently occurring terms
freq[head(ord)]
#inspect least frequently occurring terms
freq[tail(ord)]  
#Words like "can" and "one"  give us no information about the subject matter of the documents in which 
#they occur. They can therefore be eliminated without loss. Indeed, they ought to have been eliminated 
#by the stopword removal we did earlier. However, since such words occur very frequently - virtually in
#all documents - we can remove them by enforcing bounds when creating the DTM, like so:
dtmr <-DocumentTermMatrix(docs, control=list(wordLengths=c(4, 20),bounds = list(global = c(1,3))))
#Here we have told R to include only those words that occur in  1 to 3 documents. We have also enforced  
#lower and upper limit to length of the words included (between 4 and 20 characters).
#Inspecting the new DTM:
dtmr
#The dimension is reduced 3x3207
#Let's calculate the cumulative frequencies of words across documents and sort as before:
freqr <- colSums(as.matrix(dtmr))
#length should be total number of terms
length(freqr)
#create sort order (asc)
ordr <- order(freqr,decreasing=TRUE)
#inspect most frequently occurring terms
freqr[head(ordr)]
#inspect least frequently occurring terms
freqr[tail(ordr)]
#let's take get a list of terms that occur at least a  80 times in the entire corpus. This is easily done using the findFreqTerms() function as follows:
findFreqTerms(dtmr,lowfreq=80)
#check for correlations between some of these and other terms that occur in the corpus.  In this context, 
#correlation is a quantitative measure of the co-occurrence of words in multiple documents.
#if the correlation limit is 1, findAssocs() will return only  those words that always co-occur with the search term. A correlation limit
#of 0.5 will return terms that have a search term co-occurrence of at least  50% and so on.

findAssocs(dtmr,"group", .9)
#let's create a wordcloud for no other reason than everyone who can seems to be doing it.  The code for this is:

#wordcloud
library(wordcloud)
#setting the same seed each time ensures consistent look across clouds
set.seed(42)
#limit words by specifying min frequency
wordcloud(names(freqr),freqr, min.freq=5)
#Finally, one can make the wordcloud more visually appealing by adding colour as follows:

#.add color
wordcloud(names(freqr),freqr,min.freq=5,colors=brewer.pal(6,"Dark2"))
