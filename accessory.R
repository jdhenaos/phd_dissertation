library(dplyr)

input_folder <- "/home/rstudio/project/input_data/"

amines <- readRDS(paste0(input_folder,"rds_data/amines.rds"))
acyl <- readRDS(paste0(input_folder,"rds_data/acyl.rds"))
non.pos.tgs <- readRDS(paste0(input_folder,"rds_data/non_pos_tgs.rds"))
pos.tgs <- readRDS(paste0(input_folder,"rds_data/pos_tgs.rds"))
proteomics <- readRDS(paste0(input_folder,"rds_data/proteomics.rds"))

global <- rbind(amines$annotation[,c("FullRunName","H2/Geburtsgewicht.[g]","BPDgrade","Gestational.age")], 
acyl$annotation[,c("FullRunName","H2/Geburtsgewicht.[g]","BPDgrade","Gestational.age")],
non.pos.tgs$annotation[,c("FullRunName","H2/Geburtsgewicht.[g]","BPDgrade","Gestational.age")],
pos.tgs$annotation[,c("FullRunName","H2/Geburtsgewicht.[g]","BPDgrade","Gestational.age")],
proteomics$annotation[,c("FullRunName","H2/Geburtsgewicht.[g]","BPDgrade","Gestational.age")]) |> 
  distinct()

###################################

e.n <- non.pos.tgs$annotation |> 
  filter(sample.times < 28) |> 
  select(FullRunName,`H2/Geburtsgewicht.[g]`,BPDgrade,Gestational.age, sample.times)

e.p <- proteomics$annotation |> 
  filter(sample.times <= 28) |> 
  select(FullRunName,`H2/Geburtsgewicht.[g]`,BPDgrade,Gestational.age, sample.times)

e.a <- acyl$annotation |> 
  filter(sample.times <= 28) |> 
  select(FullRunName,`H2/Geburtsgewicht.[g]`,BPDgrade,Gestational.age, sample.times)

early <- rbind(
non.pos.tgs$annotation |> 
  filter(sample.times < 28) |> 
  select(FullRunName,`H2/Geburtsgewicht.[g]`,BPDgrade,Gestational.age),

amines$annotation|>
  filter(sample.times < 28) |> 
  select(FullRunName,`H2/Geburtsgewicht.[g]`,BPDgrade,Gestational.age),

acyl$annotation |>
  filter(sample.times < 28) |> 
  select(FullRunName,`H2/Geburtsgewicht.[g]`,BPDgrade,Gestational.age),

pos.tgs$annotation |>
  filter(sample.times < 28) |> 
  select(FullRunName,`H2/Geburtsgewicht.[g]`,BPDgrade,Gestational.age),

proteomics$annotation |>
  filter(sample.times < 28) |> 
  select(FullRunName,`H2/Geburtsgewicht.[g]`,BPDgrade,Gestational.age)
) |> distinct()


