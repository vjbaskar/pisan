getPrjName <- function(projectCSV){
    temp = strsplit(projectCSV, "[.]")[[1]]
    temp = temp[-length(temp)]
    temp
    pname = paste0(temp, sep="", collapse = ".")
    return(pname)
}

info <- function(text){
    text =  sapply(text, function(x) paste0(x, sep=" ", collapse = " "))
    d = date()
    cat(" [ Info: ", d, "]\n")
    cat(text,"\n")
}

errorMessage <- function(text){
    text =  sapply(text, function(x) paste0(x, sep=" ", collapse = " "))
    d = date()
    cat(" [ Error: ", d, "]\n")
    cat(text,"\n")
    cat("Quitting")
    message(" [ Error: ", d, "]\n")
    message(text,"\n")
    cat("Quitting")
    quit(save="yes", status=1)
}

createFileLink <- function(fileName){
    return(paste0(projectDir,"/", fileName))
}

stripNALines <- function(ss){
    # row wise. If Run
    temp = apply(ss, 2, unique)
    ss = ss[, - which(is.na(temp))]
    ss = ss [ !is.na(ss[,2]), ]
    ss = ss [ !is.na(ss[,1]), ]
    
    return(ss)
}

cleanSheet <- function(ss){
    
    ss = stripNALines(ss)
    
    for(j in 1:ncol(ss)){
        ss[,j] <- as.character(ss[,j])
    }
    
    head(ss)
    for(i in 1:nrow(ss)){
        for(j in 1:ncol(ss)){
            t <- as.character(ss[i,j])
            keepOnlyAlphaNumeric(t)
            
            ss[i,j] <- keepOnlyAlphaNumeric(t)
        }
    }
    
    cn <- colnames(ss)
    cn <- sapply(cn, keepOnlyAlphaNumeric)
    cn <- toupper(cn)
    colnames(ss) <- cn
    return(ss)
}



keepOnlyAlphaNumeric <- function(text){
    specialChars <- "a1~!@#$%^&*(){}_+:\"<>?./;'[]-= "
    # text = gsub("[[ ].$]", "",text)
    text = gsub("[ ]", "_",text)
    text = gsub("[^[:alnum:]^,^_^-]", "",text)
    text = gsub("_*$", "",text)
    return(text)
}