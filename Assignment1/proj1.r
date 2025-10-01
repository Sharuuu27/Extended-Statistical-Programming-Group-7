# Team Members: 
# 1. Shivasharini Govindasamy (Sharuuu27) - S2766392
# 2. Yu-Hsuan Hung (YuHsuan07) - s2793274
# 3. Yasuhiro Hara (hiroh-git) - S2826059

# Team Contributions:
# Yasuhiro - Handled data pre-processing steps (steps 4a-4f). Implemented the weightage and probabilities in sampling the next word (step 7) and assisted with other steps.
# Yu-Hsuan - Implemented the tokenisation (step 5), matrix creation (step 6), provided useful insights and coding inputs on weight creation (step 7) and assisted with other steps.
# Shivasharini - Wrote the core part of `next.word` function (steps 7-9) and assisted with other steps.



## Read the file into R
#setwd("put/your/local/repo/location/here") ## comment out of submitted
a <- scan("shakespeare.txt", what="character", skip=83, nlines=196043-83,
          fileEncoding="UTF-8")



## Remove the stage directions by finding '[' and ']'  
i_ob <- grep("[", a, fixed=TRUE) ## locate all words in 'a' that contain '[' (note: ob is opening bracket)

i_uob <- c() ## create an empty list to store the locations of the unmatched ['
for (i in i_ob) {
  range<-min(i+100,length(a)) ## avoid errors
  i_uob_temp <- grep("]", a[i:range], fixed=TRUE) ## index of corresponding ']'
  
  if (length(i_uob_temp)==0) { ## if the corresponding ']' is not found within the next 100 words, 
    i_uob <- append(i_uob, i) ## store the location of '[' 
  }
}

a_s <- a
a_s[i_uob] <- gsub('\\[', '', a[i_uob]) ## delete the unmatched '['

i_obs <- grep("[", a_s, fixed=TRUE)
i_cbs <- grep("]", a_s, fixed=TRUE)
i_oc <- i_obs[i_obs %in% i_cbs] ## find index of words with ‘[’ and ‘]’

for (i in i_oc) {
  a_s[i] <- gsub("\\[.*?\\]", "", a_s[i]) ## delete stage directions within a word (e.g. "Drawer_.]—[_Singing_.]"->"Drawer_.]—" in line 234177)
}

i_obs2 <- grep("[", a_s, fixed=TRUE)
i_cbs2 <- grep("]", a_s, fixed=TRUE)

indices_to_remove <- unlist(mapply(seq, i_obs2, i_cbs2)) ## find remaining stage directions
a_s <- a_s[-indices_to_remove] ## delete remaining stage directions
## a_s --> "a" without stage directions


## Remove the words that are fully upper case
i_u <- which(
  (a_s == toupper(a_s)) &  ## locate words that are fully upper case, and numbers expressed as arabic numerals
    (a_s != 'I') &        ## except for 'I' and 'A'
    (a_s != 'A')
)

a_s2 <- a_s
a_s2 <- a_s[-i_u] ## delete character names and arabic numerals
## a_s2 --> "a" without stage directions, character names and arabic numerals


## Remove '-' and '_' from words by using gsub 
a_s2 <- gsub("([_-])", "", a_s2) ##remove “_” & "-"
## a_s2 --> "a" without stage directions, character names, arabic numerals, underscores ("_") and hyphens("-")


## Remove the punctuation marks from the word, and add the marks to punctuation marks as a new entry
split_punct <- function(x) {
  puncs <- c(",", ".", ";", "!", ":", "?") ## list of punctuations
  for (punc in puncs) {
    ii <- grep(punc, x, fixed=TRUE) ## which elements of text include punc.?
    if (length(ii)!=0) { ## if x includes corresponding punc
      xs <- rep("", length(ii)+length(x)) ## vector to store the words and puncs
      iis <- ii+1:length(ii) ## where should puncs go in xs?
      xs[iis] <- substr(x[ii], nchar(x[ii]), nchar(x[ii])) ## insert puncs
      xs[-iis] <- gsub(paste('\\', punc, sep=""), "", x) ## insert words without puncs
      x <- xs
    }
  }
  return(x)
} ## separates the punctuations from a word and makes the punctuations as a stand-alone word/token

x <- c("An", "omnishambles,", "in", "a", "headless", "chicken", "factory.")
split_punct(x) ## to check whether the function "split_punc()" works correctly


## Separate the punctuation marks from words
a_s3 <- split_punct(a_s2)
## a_s3 --> "a" without stage directions, character names, arabic numerals, underscores ("_"), hyphens("-")
##           and punctuations from the words



## Convert the cleaned word vector a ti lower case
a_s4 <- tolower(a_s3)
## a_s4 --> "a" without stage directions, character names, arabic numerals, underscores ("_"), hyphens("-"),
##           separated punctuations from the words and all words in lower case

a <- a_s4 ##rename 'a'


## Find the vector of unique words in a 
uq　<- unique(a) ## to obtain unique words from a


## 
uq_match <- match(a,uq) ##creating tokens for words in a


# Count up the occurrence of each unique word
ot<-data.frame(word=uq,count=tabulate(uq_match)) ## creating a frequency table for the unique words


# Create a vector b which contains the 1000 most 
b <- ot$word[order(ot$count, decreasing = TRUE)][1:1000] ## b = first-1000 common unique words


# Question 6a
M1 <- match(a,b,nomatch=NA) ## M1 = matching the 1000 common words with a


# Question 6b
n=length(a)
mlag=4
M <- matrix(0, nrow=n-mlag, ncol=mlag+1)
for (j in 0:mlag){
  M[,j+1] <- M1[(1 + j):(n - mlag + j)]
}
## M = 4-gram word Matrix based on the M1


# Question 7
next.word <- function(key, M, M1, w = rep(1, ncol(M) - 1)) {
  
  mlag <- ncol(M) - 1 ## define mlag from the shape of M
  M <- M[!is.na(M[, mlag+1]),] ## remove rows that last column is NA
  
  if (length(key)>mlag) {
    key <- tail(key,mlag) ## use only data from the end of key if it is too long, in order to be able to deal appropriately with any length of key
  }
  w <- w[1:length(key)] ## avoid errors with probs that is longer than the key length
  
  next.word.cands <- c() ##l ist of candidate words for each term = i length model
  p <- c() ## updated list of weights
  
  for (i in length(key):1) {
    keyword <- tail(key,i)  
    ## Ex: if the key = "besiege brow dig deep", the code will find for the 4-word match, then 3-word match, and so on until 1-word match
    ## The loop for i starts from max of 4 to 1 (if only 3 key are supplied, then the loop starts from 3 to 1)
    
    mc <- mlag - i + 1 
    ## the match will be found from the right side, Ex: using a 4-word match, it looks at column 1:4, using a 3-word match, it looks at columns 2:4, and so on
    
    ii <- colSums(!(t(M[, mc:mlag, drop = FALSE]) == keyword)) 
    ## this line finds the exact match from mc:mlag with the keyword, and returns FALSE if there is a perfect match. 
    
    matched_rows <- which(ii == 0 & is.finite(ii))
    ## rows with ALL FALSE (means = 0) is the perfect match with the keyword supplied. NOTE: There could be more than one rows with exact matches
    
    if(length(matched_rows) > 0) {
      
      next_word <- M[matched_rows, mlag + 1]
      ## it goes to the rows with exact matches, and fetches the word from the 5th column (mlag + 1)
      
      next.word.cands <- append(next.word.cands, next_word)
      ## there can be more than one word from the 5th column (m-lag +1) derived from the perfect match
      ## this collects all the potential candidates of words (from all the i-loop) to become the next word based on the weightage and probabilities
      
      p <- append(p, rep(w[length(key)-i+1]/length(next_word), length(next_word)))
      ## update the weight list
    }
  }
  
  if (length(next.word.cands)==0) {
    M1 <- M1[!is.na(M1)] ##remove NA
    next.word.cands <- append(next.word.cands, sample(M1,1)) 
  } ## when there is 0 potential candidate for the next word, randomly choose any words that is not NA from M1
  
  if (length(next.word.cands)==1) {
    next_word2 <- next.word.cands
    ## when there is only 1 potential candidate for the next word, it becomes the next word
    
  } else {
    next_word2 <- sample(next.word.cands, size=1, prob=p)
  } ## randomly choose a word to become the next word with involving weightage/biasness in the random sample
  
  return(next_word2) 
} ## returns the final output: the "next word"


# Question 8
## Randomly selecting the starting word to simulate the model.

puncs <- c(",", ".", ";", "!", ":", "?") ##list of punctuations
punctuation <- !is.na(match(b, puncs))
## using is.na = if there is a punctuation found in b, using Match, it returns FALSE 
##            = if it is word != punctuation, it returns TRUE
## ! = negates the statement, hence: False --> True, and True --> False
## Conclusion: punctuations are saved as TRUE 

word_only_token <- seq_along(b)[!punctuation] #contains only words, no punctuations


starting_token <- sample(seq_along(word_only_token), 1) #to randomly choose the starting index/word
print(starting_token)
print(b[starting_token])


# Question 8 Alternative: Using 'romeo' to start the model
get_word_token <- setNames(seq_along(b), b)

get_romeo <- which(names(get_word_token) == "romeo")
print(get_romeo)
print(b[get_romeo])

starting_token <- get_romeo ## o choose 'romeo' specifically as the starting word


# Question 9
set.seed(42) ## to enable reproducibility

generated <- starting_token

for (i in 1:100) { ## generates up to 100 words or until a full stop is generated (whichever comes first)
  token <- next.word(key = generated, M, M1)
  
  if (token == 2) {
    ## stop the loop if the generated token is a full stop
    generated <- c(generated, token)
    break
    
  } else {
    ## else keep generating the tokens until a full stop or reaching the i=100
    generated <- c(generated, token)
  }
}

full_sentence <- paste(names(get_word_token)[match(generated, get_word_token)], collapse = " ")
## converting the tokens generated from the "new.word" function into strings (clean sentence structure)

for (punc in puncs) {
  full_sentence <- gsub(paste(' \\', punc, sep=""), 
                        paste('\\', punc, sep=""), 
                        full_sentence)
} ## deletes space before punctuations to make the sentences more neat

full_sentence
