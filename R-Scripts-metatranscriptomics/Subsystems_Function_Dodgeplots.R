### Code by Eric Njiraini, E. Kariuki------
### Date: 22/11/21
### email/contact: ericmnjiraini@gmail.com, ericgathirwak@gmail.com


### Introduction------
### function to Load 'em packages-----

packeR  <-  function(pkg){
  new.pkg  <-  pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages_to_load <- c("optparse","tidyverse", "magrittr", "svDialogs", "ggsci") ## add list of packages depending on need

##install and/or loading packages
packeR(packages_to_load)

##Introduction message----
command<- dlg_message(message = "Welcome to this code. Your request is my Command. Hit 'OKAY' below", type="okcancel")$res
if(command=="cancel") {
  stop("You chose not to run the lines below")
}
###Defining parameters----
working_dir <- choose.dir(caption = "Choose the directory with the files in need:")
setwd(working_dir) # setting working directory, in my case ----- C:/Users/eric.njiraini/Downloads/Compressed/Sample_files

### Calling data-----

### Defining patterns to consider when pulling data 
control_files <- list.files(path = ".", pattern = "control_*")
experimental_files <- list.files(path = ".", pattern = "experimental_*")

### Control Data----
my_data_control <- map_df(control_files, 
                          ~read.table(.x,
                                      header = F, 
                                      quote = "", 
                                      sep = "\t", 
                                      fill = TRUE) %>% 
                            mutate(file = str_remove_all(str_remove_all(.x, ".reduced"), "control_"))) %>% 
  select_if(~ !any(is.na(.))) 

names(my_data_control) =  c("delete", "control", "Level4", "Level3", "Level2", "Level1", "file")

my_data_control %<>% select(-delete) 


### Experimental data ----
my_data_experimental <- map_df(experimental_files, 
                               ~read.table(.x, 
                                           header = F, 
                                           quote = "", 
                                           sep = "\t", 
                                           fill = TRUE) %>% 
                                 mutate(file = str_remove_all(str_remove_all(.x, "experimental_"), ".reduced"))) %>%
  select_if(~ !any(is.na(.)))


colnames(my_data_experimental) =   c("delete", "experimental", "Level4", "Level3", "Level2", "Level1", "file")

my_data_experimental %<>% select(-delete) 


#### Merging ctrl and experimental datasets based on levels------

col_to_use= paste0("Level",1) ## Which column to consider in the report

my_data_merged = my_data_control %>% 
  left_join(my_data_experimental, by=c("Level4", "Level3", "Level2", "Level1"), all.x= TRUE) %>% 
#my_data_merged %<>% 
  filter(nchar(.[,col_to_use])>0)  %>% 
  filter(!is.na(experimental) & !is.na(control)) %>% 
  rename(file=file.y)


#### Aggregating/averaging data----

my_data_plot<- my_data_merged %>% 
  select(col_to_use, control, experimental, file) %>% 
  mutate(file=str_remove_all(file, "[:digit:]")) %>% 
  group_by(Level=.[,col_to_use], file) %>% 
  summarize(experiment=mean(experimental),
            control=mean(control)) %>% 
  #filter(experimental_average!=control_average) %>% 
  pivot_longer(!Level & !file, names_to = "study", values_to = "Levels") %>% 
  arrange(desc(Level), study, desc(Levels)) %>% 
  ungroup() %>% 
  filter(str_detect(Level, fixed("RNA Metabolism", ignore_case = TRUE), negate = T)) %>%  
  head(300) 



### Plotting the data-----
my_data_plot %>% 
  ggplot(aes(x=Level, y=Levels, fill=study))+
  geom_bar(stat = "identity", position="dodge")+
  coord_flip() + xlab("Level 1 Functions") + ylab("Average Levels") +
  theme(axis.text= element_text(size=5, face="bold"),
        plot.title = element_text(color="black", size=14, face="bold"),
        plot.subtitle = element_text(color="black", size=12, face="italic"),
        legend.title = element_text(size=5.5, face="bold"), 
        axis.title.x = element_text(size=10, face="bold"),
        axis.title.y = element_text(size=10, face="bold"),
        legend.text = element_text(size=5))+
  labs(title= "SEED Subsystems Hierachical Functional Annotation",  ### Edit to desired header  
       subtitle= "Control vs. experimental subsystems hierarchy (Level 1) ", 
       fill = "Control/Experimental") + scale_fill_jco() +
  facet_wrap(~file, scales = "free_x")


dlg_message(message= "Graph output successful!Cheers")

### Export as .pdf for better resolution


