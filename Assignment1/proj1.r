
# 3. The following code will read the file into R. You will need to change the path in the setwd call to point to your local repo. Only use the given file name for the Shakespeare text file, to facilitate marking.

#setwd("put/your/local/repo/location/here") ## comment out of submitted
a <- scan("shakespeare.txt", what="character", skip=83, nlines=196043-83,
          fileEncoding="UTF-8")


# 4. Some pre-processing of a is needed.
## (a)

i_ob <- grep("[", a, fixed=TRUE) ##locate all words in 'a' that contain '[' (note: ob is opening bracket)

i_uob <- c() ##create an empty list to store the locations of the unmatched ['
for (i in i_ob) {
  range<-min(i+100,length(a)) ##avoid errors
  i_uob_temp <- grep("]", a[i:range], fixed=TRUE) ##index of corresponding ']'
  
  if (length(i_uob_temp)==0) { ##if the corresponding ']' is not found within the next 100 words, 
    i_uob <- append(i_uob, i) ##store the location of '[' 
  }
}

a_s <- a
a_s[i_uob] <- gsub('\\[', '', a[i_uob]) ##delete the unmatched '['

i_obs <- grep("[", a_s, fixed=TRUE)
i_cbs <- grep("]", a_s, fixed=TRUE)
i_oc <- i_obs[i_obs %in% i_cbs] ##find index of words with ‘[’ and ‘]’

for (i in i_oc) {
  a_s[i] <- gsub("\\[.*?\\]", "", a_s[i]) ##delete stage directions within a word (e.g. "Drawer_.]—[_Singing_.]"->"Drawer_.]—" in line 234177)
}

i_obs2 <- grep("[", a_s, fixed=TRUE)
i_cbs2 <- grep("]", a_s, fixed=TRUE)

indices_to_remove <- unlist(mapply(seq, i_obs2, i_cbs2)) ##find remaining stage directions
a_s <- a_s[-indices_to_remove] ##delete remaining stage directions


## (b)

i_u <- c()

for (i in 1:length(a_s)) {
  if ((a_s[i]==toupper(a_s[i])) &      ##locate words that are fully upper case, and numbers expressed as arabic numerals
      (a_s[i]!='I') & (a_s[i]!='A')) { ##except for 'I' and 'A'
    i_u <- append(i_u, i)
  }
}
#a_s[i_u[1:10]] ##check the first 10

a_s2 <- a_s
a_s2 <- a_s[-i_u] ##delete character names and arabic numerals


## (c)

a_s2 <- gsub("([_-])", "", a_s2) ##remove “_” & "-"


## (d)

split_punct <- function(x) {
  puncs <- c(",", ".", ";", "!", ":", "?") ##list of punctuations
  for (punc in puncs) {
    ii <- grep(punc, x, fixed=TRUE) ##which elements of text include punc.?
    if (length(ii)!=0) { ##if x includes corresponding punc
      xs <- rep("", length(ii)+length(x)) ##vector to store the words and puncs
      iis <- ii+1:length(ii) ##where should puncs go in xs?
      xs[iis] <- substr(x[ii], nchar(x[ii]), nchar(x[ii])) ##insert puncs
      xs[-iis] <- gsub(paste('\\', punc, sep=""), "", x) ##insert words without puncs
      x <- xs
    }
  }
  return(x)
}

#x <- c("An", "omnishambles,", "in", "a", "headless", "chicken", "factory.")
#split_punct(x)


## (e)

a_s3 <- split_punct(a_s2)
#a_s3[1:10]


## (f)

a_s4 <- tolower(a_s3)
#a_s4[1:10]

a <- a_s4 ##update 'a'

#5a
uq　<- unique(a)

#5b
uq_match <- match(a,uq) #--> token

#5c
ot<-data.frame(word=uq,count=tabulate(uq_match))

#5d
b <- ot$word[order(ot$count, decreasing = TRUE)][1:1000] # b= 1000 common unique words

#6a
M1 <- match(a,b,nomatch=NA) # M1= matching the 1000 common words with a


#6b
n=length(a)
mlag=4
M <- matrix(0, nrow=n-mlag, ncol=mlag+1)
for (j in 0:mlag){
  M[,j+1] <- M1[(1 + j):(n - mlag + j)]
}
# M= 4-gram word Matrix based on the M1

## Question 7

next.word_withNA <- function(key, M, M1, w = rep(1, ncol(M) - 1)) {
  
  mlag <- ncol(M) - 1 ##define mlag from the shape of M
  
  if (length(key)>mlag) {
    key <- key[(length(key)-mlag+1):length(key)] ##use only data from the end of key if it is too long, in order to be able to deal appropriately with any length of key
  }
  
  next.word.cands <- c() ##list of candidate words for each term length model
  
  for (i in length(key):1) {
    keyword <- key[(length(key) - i + 1): length(key)]  
    # Ex: if the key = besiege brow dig deep, the code will find for the 4-word match, if not found, then uses the 3-word, and so on
    # The loop for i starts from max of 4 to 1 (if only 3 key are supplied, then the loop starts from 3 to 1)
    
    mc <- mlag - i + 1 
    # the match will be found from the right side, Ex: using a 4-word match, it looks at column 1:4, using a 3-word match, it looks at columns 2:4, and so on
    
    ii <- colSums(!(t(M[, mc:mlag, drop = FALSE]) == keyword)) 
    # this line finds the exact match from mc:mlag with the keyword, and returns FALSE if there is a perfect match. 
    
    matched_rows <- which(ii == 0 & is.finite(ii))
    # rows with ALL FALSE (means = 0) is the perfect match with the keyword supplied. NOTE: There could be more than one rows with exact matches
    
    if(length(matched_rows) > 0) {
      
      next_word <- M[matched_rows, mlag + 1]
      # it goes to the rows with exact matches, and fetches the word from the 5th column (mlag + 1)
      
      next.word.cands <- append(next.word.cands, sample(next_word,1))
      # there can be more than one word from the 5th columns derived from the perfect match, hence from a collection of next-word, sample is used choose ONE word randomly
    } else {
      next.word.cands <- append(next.word.cands, NA) ##add NA to word candidate list if no matches are found
    }
  }
  if (length(key)==1) {
    next_word2 <- sample(next.word.cands,1)
    ##avoid errors with prob=1. prob must be a list.
  } else {
    w <- w[1:length(key)] ##avoid errors with probs that are longer than the key length
    next_word2 <- sample(next.word.cands, size=1, prob=w)
  }
  return(next_word2) ##randomly choose a term from each term length model
}

next.word <- function(key, M, M1, w = rep(1, ncol(M) - 1)) {
  result <- NA
  while(is.na(result)) { ##repeat generating while it generates NA
    result <- next.word_withNA(key, M, M1)
  }
  return(result)
}


# Question 8 --> Generate word from a or the 1000 common words?
## Randomly selecting the starting word to simulate the model.

puncs <- c(",", ".", ";", "!", ":", "?") ##list of punctuations
punctuation <- !is.na(match(b, puncs))
# using is.na = if there is a punctuation found in b, using Match, it returns FALSE 
#             = if it is word != punctuation, it returns TRUE
# ! = negates the statement, hence: False --> True, and True --> False
# Conclusion: punctuations are saved as TRUE 

word_only_token <- seq_along(b)[!punctuation] #contains only words, no punctuations
#word_only_token

starting_token <- sample(seq_along(word_only_token), 1)
print(starting_token)
print(b[starting_token])

## Alternative: Using 'romeo' to start the model

get_word_token <- setNames(seq_along(b), b)
#get_word_token

get_romeo <- which(names(get_word_token) == "romeo")
print(get_romeo)
print(b[get_romeo])

# Question 9

starting_token <- sample(seq_along(word_only_token), 1) #to randomly choose the starting index/word
starting_token <- get_romeo #to choose 'romeo' specifically as the starting word
generated <- starting_token


for (i in 1:100) {
  token <- next.word(key = generated, M, M1)
  #token  <- M1[next_w] # change from index to the token #<- maybe not needed?
  
  
  #if (is.na(token)) {
  #  # stop if the generated token is NA
  #  next
  #  
  #} else if (token == 2) { #I updated func. next.word not generate NA
  if (token == 2) {
    # stop the loop if the generated token is a full stop
    generated <- c(generated, token)
    break
    
  } else {
    # else keep generating the tokens until a full stop
    generated <- c(generated, token)
  }
}

full_sentence <- paste(names(get_word_token)[match(generated, get_word_token)], collapse = " ")
full_sentence

