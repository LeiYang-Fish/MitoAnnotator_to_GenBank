library(tidyverse)
library(dplyr)
library(pgirmess)
library(seqinr)

# get taxon names from the mitoannotator result files in the folder
mitoannotator_file <- list.files(pattern=".txt")
mitoannotator_name <- gsub(".txt","",mitoannotator_file)

mito_read <- lapply(mitoannotator_file,read.delim)
db_seq <- seq(1,length(mito_read))

CDS_range <- lapply(lapply(mito_read,filter, TOPOLOGY=="CDS"),select,TOPOLOGY,X)
CDS_rows <- lapply(lapply(mito_read,filter, circular=="gene"),select, X.1)
CDS_gene <- mapply(bind_cols,CDS_range,CDS_rows,SIMPLIFY=F)

#nrows1 <- sapply(mito_read, nrow)
except_loc <- lapply(db_seq, function(x) which(mito_read[[x]]$circular=="transl_except"))
cds_loc <- lapply(db_seq, function(x) which(mito_read[[x]]$TOPOLOGY=="CDS"))
stop_loc <- lapply(db_seq, function(x) which(mito_read[[x]]$X.1=="Incomplete stop codon"))
stop_minus3 <- lapply(db_seq, function (x) stop_loc[[x]]-3)
except_minus4 <- lapply(db_seq, function (x) except_loc[[x]]-4)       # MitoAnnotator sometimes does not give transl_except information for incomplete stop codons

Full_stop_loc <- lapply(db_seq, function (x) which(cds_loc[[x]]%in%stop_minus3[[x]]==FALSE))           # which of the protein-coding genes has COMPLETE stop codons
Full_stop_loc1 <- lapply(db_seq, function (x) which(cds_loc[[x]]%in%except_minus4[[x]]==FALSE))         # same as above but consider those incomplete stop codons that have no transl_except info
incom_except_loc <- lapply(db_seq, function (x) which(cds_loc[[x]]%in%except_minus4[[x]]==TRUE))           # which of the protein-coding genes has trans_except info
incom_term <-lapply(db_seq, function (x) mito_read[[x]]$X.1[except_loc[[x]]])            # notes on incomplete stop codons, e.g. "(pos:7921,aa:TERM)"

Con_start <-lapply(db_seq, function (x) mito_read[[x]]$X.1[cds_loc[[x]]]) 

atypical_tRNA <- lapply(db_seq, function(x) which(str_detect(mito_read[[x]]$X.1, "Atypical", negate = FALSE)==TRUE))
atypical_note <- lapply(db_seq, function (x) mito_read[[x]]$X.1[atypical_tRNA[[x]]])            # notes on atypical codons for tRNAs, e.g. Atypical codon: UAA"

atypical_pos <- lapply(db_seq, function (x) atypical_tRNA[[x]]-1)
atypical_range <- lapply(db_seq, function (x) mito_read[[x]]$X[atypical_pos[[x]]])

atypical_start <- atypical_range%>%
  rapply(function(x) gsub("[a-z]", "", x), how = "replace")%>%
  rapply(function(x) gsub("[.]", "a", x), how = "replace")%>%
  rapply(function(x) gsub("aa.*", "", x), how = "replace")%>%
  rapply(function(x) gsub("[[:punct:]]","", x), how = "replace")


transTab_loc <- lapply(db_seq, function(x) which(str_detect(mito_read[[x]]$circular, "transl_table", negate = FALSE)==TRUE))          # location of "transl_table"
transTab <-lapply(db_seq, function (x) mito_read[[x]]$X.1[transTab_loc[[x]]])            # transl_table, e.g. transl_table=2 is for vertebrate mitochondrion

other_element <- lapply(lapply(mito_read,filter, TOPOLOGY %in% c("rRNA","tRNA","D-loop")),select,TOPOLOGY,X,X.1)
all_elements <- mapply(bind_rows,CDS_gene,other_element,SIMPLIFY=F)

all_elements_Final <- mapply(select,all_elements, "X", "TOPOLOGY", "X.1", SIMPLIFY=F) %>% 
  lapply(setNames, nm = c("Position","Topology","Element"))

Product_column0 <- mapply(select, all_elements, "TOPOLOGY", SIMPLIFY=F)
Product_column <- rapply(Product_column0, as.character, classes="factor", how="replace")%>%               # IMPORTANT! to change the "interger" type to "character" type
  lapply(function(x) replace(x, "PROD","product"))

Complem_loc <- lapply(db_seq, function (x) which(str_detect(all_elements_Final[[x]]$Position, "complement", negate = FALSE)==TRUE))
Complem_list <- lapply(db_seq, function (x) all_elements_Final[[x]]$Position[Complem_loc[[x]]])
Complem_range <- Complem_list%>%
  rapply(function(x) gsub("[a-z]", "", x), how = "replace")%>%
  rapply(function(x) gsub("[.]", "a", x), how = "replace")%>%
  rapply(function(x) gsub("[[:punct:]]","", x), how = "replace")%>%
  rapply(function(x) strsplit(x,"aa"), how = "replace")

Complem_range1 <- lapply(lapply(Complem_range,unlist), matrix, ncol=2, byrow=TRUE) %>% 
  lapply(as.data.frame) %>% 
  lapply(setNames, nm = c("End_ntC","Start_ntC"))


# split the location of each element into two numbers (start position, end position)
range <- lapply(all_elements_Final, "[","Position")%>%
  rapply(function(x) gsub("[a-z]", "", x), how = "replace")%>%
  rapply(function(x) gsub("[.]", "a", x), how = "replace")%>%
  rapply(function(x) gsub("[[:punct:]]","", x), how = "replace")%>%
  rapply(function(x) strsplit(x,"aa"), how = "replace")

range1 <- lapply(lapply(range,unlist), matrix, ncol=2, byrow=TRUE) %>% 
  lapply(as.data.frame) %>% 
  lapply(setNames, nm = c("Start_nt","End_nt"))

column3 <- mapply(select, range1,"Start_nt", SIMPLIFY=F) %>% 
  rapply(as.character,how = "replace") %>% 
  rapply(as.numeric,how = "replace")

column4 <- mapply(select, range1,"End_nt", SIMPLIFY=F) %>% 
  rapply(as.character,how = "replace") %>% 
  rapply(as.numeric,how = "replace")

Element_name <- lapply(all_elements_Final, "[","Element")
#Element_name <- mapply(select, all_elements_Final, "Element", SIMPLIFY=F)            # also works

Start_End <- mapply(cbind,column3,column4,SIMPLIFY=F)

mt_dt <- mapply(cbind,Start_End,Product_column, Element_name, SIMPLIFY=F)            # the 38-element annotations

# for protein-coding genes, need to annotate not only genes, but also proteins
CDS_rev <- lapply(mt_dt,filter, TOPOLOGY=="CDS") %>% 
  lapply(function(x) replace (x, "TOPOLOGY", "gene")) %>% 
  lapply(function(x) replace (x, "PROD", "gene"))

# change gene names to protein names for protein-coding genes in the 38-element annotations
Element_protein <- rapply(Element_name, function(x) gsub("ND", "NADH dehydrogenase subunit ", x), how = "replace") %>% 
  rapply(function(x) gsub("CO", "cytochrome c oxidase subunit ", x), how = "replace") %>% 
  rapply(function(x) gsub("ATPase ", "ATPase subunit ", x), how = "replace") %>% 
  rapply(function(x) gsub("Cyt b", "cytochrome b", x), how = "replace")

First_4cols <- lapply(db_seq, function(x) mt_dt[[x]][1:4])
#First_4cols <- mapply(select, mt_dt,"Start_nt", "End_nt","TOPOLOGY","PROD", SIMPLIFY=F)      # also works
New_rev <- mapply(cbind,First_4cols, Element_protein, SIMPLIFY=F)


# The following part focuses on the 13 protein-coding genes and their features
CDS_cds_only <- lapply(New_rev,filter, TOPOLOGY=="CDS")

nrows_cds <- sapply(CDS_cds_only, nrow)
CDS_cds_only1 <- lapply(db_seq, function (x) {CDS_cds_only[[x]][rep(1:nrows_cds[x], each = 5), ]})              # duplicate each line 5 times
nrows_cds2 <- sapply(CDS_cds_only1, nrow)

for (i in db_seq) {
  CDS_cds_only1[[i]][1:nrows_cds2[i] %% 5 == 0, 1:5]<- ""
  CDS_cds_only1[[i]][(1:nrows_cds2[i]+3) %% 5 == 0, 1:3]<- ""
  CDS_cds_only1[[i]][(1:nrows_cds2[i]+3) %% 5 == 0, 4]<- "transl_table"

  CDS_cds_only1[[i]][(1:nrows_cds2[i]+2) %% 5 == 0, 1:2]<- ""
  CDS_cds_only1[[i]][(1:nrows_cds2[i]+2) %% 5 == 0, 3]<- "codon_start=xx"
  CDS_cds_only1[[i]][(1:nrows_cds2[i]+2) %% 5 == 0, 4]<- "note"
  CDS_cds_only1[[i]][(1:nrows_cds2[i]+2) %% 5 == 0, 5]<- "Incomplete stop codon"
  
  CDS_cds_only1[[i]][(1:nrows_cds2[i]+1) %% 5 == 0, 1:3]<- ""
  CDS_cds_only1[[i]][(1:nrows_cds2[i]+1) %% 5 == 0, 4]<- "transl_except"
  CDS_cds_only1[[i]][Full_stop_loc1[[i]]*5-1, 4:5] <- ""
  CDS_cds_only1[[i]][Full_stop_loc[[i]]*5-2, 4:5] <- ""
}

cds_loc_dt <- lapply(db_seq, function(x) which(CDS_cds_only1[[x]]$TOPOLOGY=="CDS"))  # CDS gene location in the modified dataset CDS_cds_only1 

# adding incomplete stop codon information for each incomplete stop codon
for (i in db_seq) {
  for (j in 1:length(incom_term[[i]]))
    CDS_cds_only1[[i]][incom_except_loc[[i]][[j]]*5-1, 5]<- as.character(incom_term[[i]][[j]])
}

# adding tranlation table informtion for the 13 protein-coding genes
for (i in db_seq) {
  for (j in 1:length(transTab[[i]]))
    CDS_cds_only1[[i]][cds_loc_dt[[i]][[j]]+1, 5]<- as.character(transTab[[i]][[j]])
}

# adding condon start position informtion for the 13 protein-coding genes
for (i in db_seq) {
  for (j in 1:length(Con_start[[i]]))
    CDS_cds_only1[[i]][cds_loc_dt[[i]][[j]]+2, 3]<- paste0("codon_start=", as.character(Con_start[[i]][[j]]))
}

# separate columns 1-3 and columns 4&5
cds123 <- lapply(db_seq, function (x) CDS_cds_only1[[x]][1:3])
cds45 <- lapply(db_seq, function (x) CDS_cds_only1[[x]][4:5])

# add one blank row to the very front of columns 4&5
row0 <- lapply(db_seq, function (x) {cds45[[x]][rep(1, each = 1), ]})
for (i in db_seq) {
  row0[[i]][1:2]<- ""}
cds45_row0 <- mapply(bind_rows,row0,cds45,SIMPLIFY=F)

# add one blank row to the very end of columns 1-3
row0_2 <- lapply(db_seq, function (x) {cds123[[x]][rep(1, each = 1), ]})
for (i in db_seq) {
  row0_2[[i]][1:3]<- ""}

cds123_row0 <- mapply(bind_rows,cds123,row0_2, SIMPLIFY=F)

# bind all five columns together again
cds_final <- mapply(cbind,cds123_row0, cds45_row0, SIMPLIFY=F)

CDS_other <- lapply(New_rev,filter, TOPOLOGY!="CDS")

mt_dt1 <- mapply(bind_rows,CDS_rev,CDS_other,SIMPLIFY=F)          # the 51-element (38 genes +13 proteins) annotations
mt_dt2 <- lapply(db_seq,function(x) arrange(mt_dt1[[x]],mt_dt1[[x]]$Start_nt))          # sort according to start positions of genes



##### the following codes deal with blank rows, row names, and sample information row
nrows <- sapply(mt_dt2, nrow)
dt <- lapply(db_seq, function (x) {mt_dt2[[x]][rep(1:nrows[x], each = 3), ]})           # duplicate each line twice

nrows_dt <- sapply(dt, nrow)

for (i in db_seq) {
  dt[[i]][1:nrows_dt[i] %% 3 == 0, ]<- ""}            # make each of the 2nd duplicated line blank

for (i in db_seq) {
  dt[[i]][(1:nrows_dt[i]+1) %% 3 == 0, ]<- ""}              # make each of the 1st duplicated line blank

# separate columns 1-3 and columns 4&5
dt123 <- lapply(db_seq, function (x) dt[[x]][1:3])
dt45 <- lapply(db_seq, function (x) dt[[x]][4:5])

# add one blank row to the very front of columns 4&5
for (i in db_seq) {
  row0[[i]][1:2]<- ""}
dt45_row0 <- mapply(bind_rows,row0,dt45,SIMPLIFY=F)

# add one blank row to the very end of columns 1-3
for (i in db_seq) {
  row0_2[[i]][1:3]<- ""}
dt123_row0 <- mapply(bind_rows,dt123,row0_2, SIMPLIFY=F)

# bind all five columns together again
dt_final <- mapply(cbind,dt123_row0, dt45_row0, SIMPLIFY=F)

atypical_in_dt <- lapply(db_seq, function (y) lapply(1:length(atypical_start[[y]]), function(x) {which(dt_final[[y]]$Start_nt==atypical_start[[y]][x])}))
atypical_in_dt1 <- lapply(atypical_in_dt, unlist)

# some protein-coding and tRNA genes are encoded on the different strand, find them out and swap their start and end positions
Complem_len <- lapply(db_seq, function (x) nrow(Complem_range1[[x]]))
cds_comp_pos <- lapply(db_seq, function (y) lapply(1:Complem_len[[y]], function(x) {which(cds_final[[y]][,1]==Complem_range1[[y]][x,1])}))
dt_comp_pos <- lapply(db_seq, function (y) lapply(1:Complem_len[[y]], function(x) {which(dt_final[[y]][,1]==Complem_range1[[y]][x,1])}))
cds_comp_pos_unlist <- lapply(cds_comp_pos, unlist)
dt_comp_pos_unlist <- lapply(dt_comp_pos, unlist)

for (i in db_seq) {
  for (j in 1:length(cds_comp_pos_unlist[[i]]))
    swap(cds_final[[i]][cds_comp_pos_unlist[[i]][[j]],1],cds_final[[i]][cds_comp_pos_unlist[[i]][[j]],2])
    }

for (i in db_seq) {
  for (j in 1:length(dt_comp_pos_unlist[[i]]))
    swap(dt_final[[i]][dt_comp_pos_unlist[[i]][[j]],1],dt_final[[i]][dt_comp_pos_unlist[[i]][[j]],2])
}

# make one blank row to put sample information
Sample_info <- lapply(db_seq, function (x) {dt_final[[x]][rep(1, each = 1), ]})
for (i in db_seq) {
  Sample_info[[i]][1,1]<- paste0(">Feature ", mitoannotator_name[[i]])
  Sample_info[[i]][1,2]<-""
  Sample_info[[i]][1,3]<-""}

# add the sample information row to the very front
dt_final2 <- mapply(bind_rows,Sample_info, dt_final,cds_final, SIMPLIFY=F)

Pos_add_atypical <- lapply(db_seq, function (x) atypical_in_dt1[[x]]+3)


for (i in db_seq) {
  for (j in 1:length(atypical_note[[i]]))
    dt_final2[[i]][Pos_add_atypical[[i]][[j]], 5]<- as.character(atypical_note[[i]][[j]])
}

for (i in db_seq) {
  for (j in 1:length(atypical_note[[i]]))
  dt_final2[[i]][Pos_add_atypical[[i]][[j]], 4]<- "note"
}

# remove column names
for (i in db_seq) {
  colnames(dt_final2[[i]])<- NULL}


# write all transformed annotations into a single txt file
for (i in db_seq) {
  write.delim(dt_final2[[i]], "transformed_GenBank_output.txt", append= T)
}




