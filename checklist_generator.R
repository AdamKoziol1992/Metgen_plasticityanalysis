library(tidyverse)
library(data.table)

setwd('~/Documents/GitHub/Metgen_plasticityanalysis/')

my_table <- fread('./Checklist_GSC-MIxS host associated_1670950180471.csv')

metadata <- fread('./data/metadata.csv')

male <- metadata %>%
  filter(sex == 'MALE') %>%   
  select(sample_id) %>% 
  mutate(sample_id = str_to_lower(sample_id)) %>% 
  as.vector()

female <- metadata %>%
  filter(sex == 'FEMALE') %>%   
  select(sample_id) %>% 
  mutate(sample_id = str_to_lower(sample_id))

my_filenames <- fread('./filenames.txt', header = F) %>% 
  rename('Filename' = 'V1') %>%
  mutate(sample_title = str_extract(Filename, '[^.]+'),
         tax_id = if_else(str_detect(sample_title, '_s[1:9]'), 'Crocidura_russula', 'Apodemus_sylvaticus'),
         scientific_name = tax_id,
         sample_alias = str_remove(sample_title, '..$'),
         project_name = rep('ERA19381426'),
         geographic_location_country_andor_sea = rep('Spain'),
         geographic_location_latitude = rep('43.2848° N'),
         geographic_location_longitude = rep('2.1714° W'),
         environmental_medium = rep('Faeces'),
         read_number = str_extract(sample_title, '.$'),
         collection_date = case_when(str_detect(sample_alias, '(?i)pre') ~ '2019/10/09',
                                     str_detect(sample_alias, '(?i)post') ~ '2019/10/16',
                                     str_detect(sample_alias, '(?i)heat') ~ "2019/10/28",
                                     str_detect(sample_alias, 'cold') ~ "2019/11/11",
                                     str_detect(sample_alias, 'diet') ~ "2019/11/23",
                                     T ~ 'NA'),
         sample_description = case_when(str_detect(sample_alias, '(?i)pre') ~ 'Day1',
                                        str_detect(sample_alias, '(?i)post') ~ 'Acclimation',
                                        str_detect(sample_alias, '(?i)heat') ~ "Heat",
                                        str_detect(sample_alias, 'cold') ~ "Cold",
                                        str_detect(sample_alias, 'diet') ~ "Diet",
                                        T ~ 'NA')) %>%  
  filter(!str_detect(sample_alias, 'Antibiotic|fast|_anti|filenam|Mousechow|Shrew|Soil|reco|Diet2'))

fwrite(my_filenames, './my_filenames.csv')
