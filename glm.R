
library(lme4)


fitGLM_fixed <- function(res, condition,   subject_effect = TRUE, pairwise = TRUE, fixed_only = FALSE, verbose = TRUE){
    
    fit_random <- list()
    fit_fixed <- list()
    for(i in 1:ncol(res$nstar)){
        # idx <- indexes_list[[i]]
        if(verbose){
            if(i%%10==0){
                print(paste("fitting GLM...", i))
            }
        }
        
        
        glm_df <-  cbind(res$info[,1:2], res$nstar[,i])
        
        # glm_df <- melt(glm_df)
        colnames(glm_df) <- c("cellTypes", "subject", "cell_count")
        glm_df$cond <- condition
        
        # get degree of freedom
        mod <- stats::glm(cell_count ~ cellTypes   + cond +  cellTypes:cond + subject,
                          data = glm_df,
                          family = poisson(link=log))
        
        # get degree of freedom, need this for mice::pool() otherwise weird error
        degree_freedom = mod$df.null
        
        if(subject_effect){
            if(pairwise){
                
                fit_fixed[[i]] <- stats::glm(cell_count ~ cellTypes   + cond +  cellTypes:cond + subject,
                                             data = glm_df,
                                             family = poisson(link=log))
                if(!fixed_only){
                    fit_random[[i]] <- lme4::glmer(cell_count ~ cellTypes   + cond +  
                                                       cellTypes:cond + (1 | subject ),
                                                   data = glm_df, family = poisson(link=log),
                                                   control = glmerControl(nAGQ = 0L))
                }
            }else{
                fit_fixed[[i]] <- stats::glm(cell_count ~ cellTypes  + cond +  cellTypes:cond + subject,
                                             data = glm_df, family = poisson(link=log))
                if(!fixed_only){
                    
                    fit_random[[i]] <- lme4::glmer(cell_count ~ cellTypes   + cond +
                                                       cellTypes*cond + (1 | subject ), data = glm_df,
                                                   family = poisson(link=log),
                                                   control = glmerControl(nAGQ = 0L))
                }
            }
        }else{
            fit_fixed[[i]] <- stats::glm(cell_count ~ cellTypes   + cond +  cellTypes:cond, data = glm_df,
                                         family = poisson(link=log))
        }
        
    }
    
    if(!subject_effect){
        fixed_only = TRUE
    }
    
    if(!fixed_only){
        pool_res_random = mice::pool(fit_random  , dfcom = degree_freedom)
        pool_res_fixed = mice::pool(fit_fixed , dfcom = degree_freedom)
        return(list(pool_res_random = pool_res_random,
                    pool_res_fixed = pool_res_fixed,
                    fit_random = fit_random,
                    fit_fixed = fit_fixed))
    }else{
        pool_res_fixed = mice::pool(fit_fixed , dfcom = degree_freedom)
        return(list(pool_res_fixed = pool_res_fixed,
                    fit_fixed = fit_fixed))
    }
    
}
