fig_source <- "C://Users//boennecd//Dropbox//skole_backup//phd//Moodys_DRD//figure"
rpres_file <- readChar("intro_seminar_phd.Rpres", file.info("intro_seminar_phd.Rpres")$size)

library(stringr)

all_pics <- str_extract_all(rpres_file, "(?<=\\]\\()[^\\)]+(?=\\))")
all_pics <- c(unlist(all_pics), unlist(str_extract_all(
  rpres_file, '(?<=src\\=\\")[^"]+\\.((jpg)|(jpeg)|(png))')))

all_pics <- str_replace(all_pics, "^figures/", "")
all_pics_regexp <- str_replace(all_pics, "\\.", "\\\\.")

for(s in seq_along(all_pics_regexp)){
  new_f <- list.files(path = fig_source, pattern = all_pics[s], recursive = TRUE, full.names = T)
  file.copy(new_f, paste0("figures/", all_pics[s]))
}
