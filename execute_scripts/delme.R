> access_list
function(list_to_open, condition = "a conditional statement", executable_code = "print('Insert string_code')", acr = c("tsv", "csv", "gz"), extra_args = NULL){
  
  access.dfs <- function(thislist, name.list = list()){
    
    lapply(names(thislist), function(thisname){
      
      #this item is a sublist or a df
      thisitem <- thislist[[thisname]]
      #name.list contains the name of the list you just entered
      
      
      if(eval(parse(text = condition))){
        # If you want to access a dataframe, make sure to use these
        # df <- thisitem[[1]]
        eval(parse(text = executable_code))
      }else{
        # store sublist's name
        name.list[[length(name.list) + 1]] <- thisname
        temp <- access.dfs(thislist = thisitem, name.list = name.list)
        names(temp) <- eval(parse(text = paste0("names(list_to_open$", paste(unlist(name.list), collapse = "$"), ")")))
        return(temp)
      }
    })
  }
  
  new.list <- access.dfs(thislist = list_to_open)
  names(new.list) <- names(list_to_open)
  return(new.list)
}