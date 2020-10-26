split_rat <- function(rat){
  split_list <- unlist(strsplit(rat,"[/]"))
  v1 <- first(split_list)
  v2 <- last(split_list)
  return(c(v1,v2))
}

pipeline_results <- function(which_response,plot_title="",iters = 0:110,best_iter=NULL){
  if(is.null(which_response)){stop("Pick a folder corresponding to a response variable!")}
  path1 <- paste("./scripts/",which_response,"/iter",sep="")
  grab_oob <- function(i){
    return(as.numeric(readLines(paste(path1,i,"_oob_results.txt",sep=""))[-1]))
  }
  length_rat <- function(i){
    length(readLines(paste(path1,i,"_ratios.txt",sep="")))
  }
  out_data <- tibble(auc_vec = lapply(iters,grab_oob),
                     num_rat = sapply(iters,length_rat),
                     iter = iters + 1)
  out_data <- out_data %>% 
    mutate(sims = lapply(iter,function(i){1:length(auc_vec[[i]])}))
  out_data <- out_data %>% 
    unnest(cols = c(auc_vec, sims))
  out_data <- out_data %>% 
    mutate(num_rat = as.factor(num_rat))
  mean_out_data <- out_data %>% 
    select(-sims) %>% 
    group_by(num_rat,iter) %>% 
    summarise_all(funs(mean,sd)) %>% 
    ungroup() %>% 
    mutate(iter_rev = 1:max(iter))
  if(is.null(best_iter)){print(best_iter <- with(mean_out_data,iter[which.max(mean)]))}
  p <- out_data %>%
    ggplot(aes(x=iter,y=auc_vec)) +
    geom_jitter(height=0,width=0) +
    theme_bw() +
    xlab("Iteration") +
    ylab("OOB AUROC") +
    ggtitle(plot_title) +
    theme(plot.title = element_text(hjust=0.5),
          axis.text.x = element_text(vjust=-.01,angle=45)) +
    geom_errorbar(data=mean_out_data,aes(ymin=mean-2*sd,
                                         ymax=mean+2*sd,
                                         x=iter),
                  width=0.1,inherit.aes = FALSE) +
    geom_point(data=mean_out_data,mapping=aes(x=iter,y=mean),colour="red",size=1) +
    geom_line(data=mean_out_data,mapping=aes(x=iter,y=mean)) +
    geom_vline(xintercept=best_iter)
  return(p)
}
