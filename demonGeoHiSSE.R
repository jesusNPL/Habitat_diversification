##### Function to perform GeoHiSSE analysis on a multiPhylo object #####

# demonGeoHiSSE is a function to perform analysis of geographical-state diversification
# based on the GeoHiSSE model.
# see details in Caetano et al. 2018 Evolution (https://onlinelibrary.wiley.com/doi/abs/10.1111/evo.13602) 
# Also, this function is a modification of the function "evaluate.models" of Caetano et al. 2018
# available at https://figshare.com/collections/Data_for_Hidden_state_models_improve_state-dependent_diversification_approaches_including_biogeographical_models/4069580/2

demonGeoHiSSE <- function(model, phy, ranges, f = c(1, 1, 1), number_of_trees, outfile){
  
  posteriors <- sample(phy, number_of_trees)
  
  for(i in 1:length(posteriors)){
    svMisc::progress(i, max.value = length(posteriors))
    
    if (length(posteriors) == 1){phy <- phy} else {phy <- posteriors[[i]]}
    
    ## Model 1 - Dispersal parameters vary only, no range-dependent diversification.
    if (model == 1){
      speciation <- c(1, 1, 1)
      extirpation <- c(1, 1)
      trans.rate <- TransMatMakerGeoHiSSE(hidden.areas = 0, include.jumps = FALSE, separate.extirpation = FALSE)
      model1 <- try(GeoHiSSE(phy, ranges, f = f, speciation = speciation, extirpation = extirpation, 
                           hidden.areas = FALSE, trans.rate = trans.rate, assume.cladogenetic = TRUE))
      saveRDS(model1, file = paste0(outfile, i , "_model1.rds"))
      
      capture.output(print(paste0("Analysis - ", outfile, i, " - model 1 done.")), file = paste0(outfile, i, ".log"), append = TRUE)
      
    }
    ## Model 2. Canonical GeoSSE model, range effect on diversification
    if(model == 2){
      speciation <- c(1, 2, 3)
      extirpation <- c(1, 2)
      trans.rate <- TransMatMakerGeoHiSSE(hidden.areas = 0, include.jumps = FALSE, separate.extirpation = FALSE)
      model2 <- try(GeoHiSSE(phy, ranges, f = f, speciation = speciation, extirpation = extirpation, 
                           hidden.areas = FALSE, trans.rate = trans.rate, assume.cladogenetic = TRUE) )
      saveRDS(model2, file = paste0(outfile, i, "_model2.rds"))
      
      capture.output(print(paste0("Analysis - ", outfile, i,  " - model 2 done.")), file = paste0(outfile, i, ".log"), append = TRUE)
    }
    ## Model 3. Heterogeneous diversification, not tied to range evolution.
    ## Assumes three distinct diversification rates.
    if(model == 3){
      speciation <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
      extirpation <- c(1, 1, 2, 2, 3, 3)
      trans.rate <- TransMatMakerGeoHiSSE(hidden.areas = 2, make.null = TRUE, include.jumps = FALSE, separate.extirpation = FALSE)
      model3 <- try(GeoHiSSE(phy, ranges, f = f, speciation = speciation, extirpation = extirpation, 
                             hidden.areas = TRUE, trans.rate = trans.rate, assume.cladogenetic = TRUE) )
      saveRDS(model3, file = paste0(outfile, i, "_model3.rds"))
      
      capture.output(print(paste0("Analysis - ", outfile, i, " - model 3 done.")), file = paste0(outfile, i, ".log"), append = TRUE)
    }
    ## Model 4. Heterogeneous diversification, tied to range evolution.
    ## Assumes 6 distinct diversification rates.
    if(model == 4){
      print(4)
      speciation <- c(1, 2, 3, 4, 5, 6)
      extirpation <- c(1, 2, 3, 4)
      trans.rate <- TransMatMakerGeoHiSSE(hidden.areas = 1, include.jumps = FALSE, separate.extirpation = FALSE)
      model4 <- try(GeoHiSSE(phy, ranges, f = f, speciation = speciation, extirpation = extirpation, 
                           hidden.areas = TRUE, trans.rate = trans.rate, assume.cladogenetic = TRUE))
      saveRDS(model4, file = paste0(outfile, i, "_model4.rds"))
      
      capture.output(print(paste0("Analysis - ", outfile, i, " - model 4 done.")), file = paste0(outfile, i, ".log"), append = TRUE)
    }
    ## Model 5. Heterogeneous diversification, not tied to range evolution. Assumes 5 distinct diversification rates.
    if(model == 5){
      speciation <- c(rep(1 ,3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3))
      extirpation <- c(rep(1, 2), rep(2, 2), rep(3, 2), rep(4, 2), rep(5, 2))
      trans.rate <- TransMatMakerGeoHiSSE(hidden.areas = 4, make.null = TRUE, include.jumps = FALSE, separate.extirpation = FALSE)
      model5 <- try(GeoHiSSE(phy, ranges, f = f, speciation = speciation, extirpation = extirpation, 
                              hidden.areas = TRUE, trans.rate = trans.rate, assume.cladogenetic = TRUE))
      saveRDS(model5, file = paste0(outfile, i, "_model5.rds"))
      capture.output(print(paste0("Analysis - ", outfile, i, " - model 5 done.")), file = paste0(outfile, i, ".log"), append = TRUE)
    }
    ## Model 6. Heterogeneous diversification, not tied to range evolution. Assumes two distinct diversification rates.
    if(model == 6){
      speciation <- c(1, 1, 1, 2, 2, 2)
      extirpation <- c(1, 1, 2, 2)
      trans.rate <- TransMatMakerGeoHiSSE(hidden.areas = 1, make.null = TRUE, include.jumps = FALSE, separate.extirpation = FALSE)
      model6 <- try(GeoHiSSE(phy, ranges, f = f, speciation = speciation, extirpation = extirpation, 
                              hidden.areas = TRUE, trans.rate = trans.rate, assume.cladogenetic = TRUE))
      saveRDS(model6, file=paste0(outfile, i, "_model6.rds"))
      capture.output(print(paste0("Analysis - ", outfile, i, " - model 6 done.")), file = paste0(outfile, i, ".log"), append = TRUE)
    }
    ###############################################################################
    ## Block of cladogenetic models not GeoSSE-like.
    ## Here extirpation is NOT linked to range reduction.
    ## So range reduction is different from the extinction of an endemic lineage.
    ## Jumps between endemic areas are not allowed.
    ###############################################################################
    
    ## Model 1 - Dispersal parameters vary only, no range-dependent diversification.
    if(model == 7){
      speciation <- c(1, 1, 1)
      extirpation <- c(1, 1)
      trans.rate <- TransMatMakerGeoHiSSE(hidden.areas = 0, include.jumps = FALSE, separate.extirpation = TRUE)
      model7 <- try(GeoHiSSE(phy, ranges, f = f, speciation = speciation, extirpation = extirpation, 
                              hidden.areas = FALSE, trans.rate = trans.rate, assume.cladogenetic = TRUE))
      saveRDS(model7, file = paste0(outfile, i, "_model7.rds"))
      capture.output(print(paste0("Analysis - ", outfile, i, " - model 7 done.")), file = paste0(outfile, i, ".log"), append = TRUE)
    }
    ## Model 2. Canonical GeoSSE model, range effect on diversification
    if(model == 8){
      print(8)
      speciation <- c(1, 2, 3)
      extirpation <- c(1, 2)
      trans.rate <- TransMatMakerGeoHiSSE(hidden.areas = 0, include.jumps = FALSE, separate.extirpation = TRUE)
      model8 <- try(GeoHiSSE(phy, ranges, f = f, speciation = speciation, extirpation = extirpation, 
                              hidden.areas = FALSE, trans.rate = trans.rate, assume.cladogenetic = TRUE))
      saveRDS(model8, file = paste0(outfile, i, "_model8.rds"))
      capture.output(print(paste0("Analysis - ", outfile, i, " - model 8 done.")), file = paste0(outfile, i, ".log"), append = TRUE)
    }
    ## Model 3. Heterogeneous diversification, not tied to range evolution.
    ## Assumes three distinct diversification rates.
    if(model == 9){
      speciation <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
      extirpation <- c(1, 1, 2, 2, 3, 3)
      trans.rate <- TransMatMakerGeoHiSSE(hidden.areas = 2, make.null = TRUE, include.jumps = FALSE, separate.extirpation = TRUE)
      model9 <- try(GeoHiSSE(phy, ranges, f = f, speciation = speciation, extirpation = extirpation, 
                             hidden.areas = TRUE, trans.rate = trans.rate, assume.cladogenetic = TRUE))
      saveRDS(model9, file = paste0(outfile, i, "_model9.rds"))
      capture.output( print( paste0("Analysis - ", outfile, i, " - model 9 done.")), file = paste0(outfile, i, ".log"), append = TRUE)
    }
    ## Model 4. Heterogeneous diversification, tied to range evolution.
    ## Assumes 6 distinct diversification rates.
    if(model == 10){
      speciation <- c(1, 2, 3, 4, 5, 6)
      extirpation <- c(1, 2, 3, 4)
      trans.rate <- TransMatMakerGeoHiSSE(hidden.areas = 1, include.jumps = FALSE, separate.extirpation = TRUE)
      model10 <- try(GeoHiSSE(phy, ranges, f = f, speciation = speciation, extirpation = extirpation, 
                               hidden.areas = TRUE, trans.rate = trans.rate, assume.cladogenetic = TRUE))
      saveRDS(model10, file = paste0(outfile, i, "_model10.rds"))
      capture.output(print(paste0("Analysis - ", outfile, i, " - model 10 done.")), file = paste0(outfile, i, ".log"), append = TRUE)
    }
    ## Model 5. Heterogeneous diversification, not tied to range evolution. Assumes 5 distinct diversification rates.
    if(model == 11){
      speciation <- c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3))
      extirpation <- c(rep(1, 2), rep(2, 2), rep(3, 2), rep(4, 2), rep(5, 2))
      trans.rate <- TransMatMakerGeoHiSSE(hidden.areas = 4, make.null = TRUE, include.jumps = FALSE, separate.extirpation = TRUE)
      model11 <- try(GeoHiSSE(phy, ranges, f = f, speciation = speciation, extirpation = extirpation, 
                              hidden.areas = TRUE, trans.rate = trans.rate, assume.cladogenetic = TRUE))
      saveRDS(model11, file = paste0(outfile, i, "_model11.rds"))
      capture.output(print(paste0("Analysis - ", outfile, i, " - model 11 done.")), file = paste0(outfile, i, ".log"), append = TRUE)
    }
    ## Model 6. Heterogeneous diversification, not tied to range evolution. Assumes two distinct diversification rates.
    if(model == 12){
      speciation <- c(1, 1, 1, 2, 2, 2)
      extirpation <- c(1, 1, 2, 2)
      trans.rate <- TransMatMakerGeoHiSSE(hidden.areas = 1, make.null = TRUE, include.jumps = FALSE, separate.extirpation = TRUE)
      model12 <- try(GeoHiSSE(phy, ranges, f = f, speciation = speciation, extirpation = extirpation, 
                              hidden.areas = TRUE, trans.rate = trans.rate, assume.cladogenetic = TRUE) )
      saveRDS(model12, file = paste0(outfile, i, "_model12.rds"))
      capture.output(print(paste0("Analysis - ", outfile, i, " - model 12 done.")), file = paste0(outfile, i, ".log"), append = TRUE)
    }
    
    ###############################################################################
    ## Second block of anagenetic models MuSSE-like.
    ## These are very liberal models. They are really GeoSSE-like without cladogenetic events
    ###############################################################################
    
    ## Model 1. Transitions only. No character effect on diversification
    if(model == 13){
      speciation <- c(1, 1, 1)
      extirpation <- c(1, 1, 1)
      trans.rate <- TransMatMakerGeoHiSSE(hidden.areas = 0, include.jumps = FALSE, separate.extirpation = TRUE)
      model13 <- try(GeoHiSSE(phy, ranges, f = f, speciation = speciation, extirpation = extirpation, 
                          hidden.areas = FALSE, trans.rate=trans.rate, assume.cladogenetic = FALSE))
      saveRDS(model13, file = paste0(outfile, i, "_model13.rds"))
      
      capture.output(print(paste0("Analysis - ", outfile, i, " - model 13 done.")), file = paste0(outfile, i, ".log"), append = TRUE)
    }
    ## Model 2. Character effect on diversification.
    if(model == 14){
      speciation <- c(1, 2, 3)
      extirpation <- c(1, 2, 3)
      trans.rate <- TransMatMakerGeoHiSSE(hidden.areas = 0, include.jumps = FALSE, separate.extirpation = TRUE)
      model14 <- try(GeoHiSSE(phy, ranges, f = f, speciation = speciation, extirpation = extirpation, 
                              hidden.areas = FALSE, trans.rate = trans.rate, assume.cladogenetic = FALSE))
      saveRDS(model14, file = paste0(outfile, i, "_model14.rds"))
      
      capture.output(print(paste0("Analysis - ", outfile, i, " - model 14 done.")), file = paste0(outfile, i, ".log"), append = TRUE)
    }
    ## Model 3. No character effect on diversification.
    if(model == 15){
      speciation <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
      extirpation <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
      trans.rate <- TransMatMakerGeoHiSSE(hidden.areas = 2, include.jumps = FALSE, separate.extirpation = TRUE, make.null = TRUE)
      model15 <- try(GeoHiSSE(phy, ranges, f = f, speciation = speciation, extirpation = extirpation, 
                          hidden.areas = TRUE, trans.rate = trans.rate, assume.cladogenetic = FALSE))
      saveRDS(model15, file = paste0(outfile, i, "_model15.rds"))
      
      capture.output(print(paste0("Analysis - ", outfile, i, " - model 15 done.")), file = paste0(outfile, i, ".log"), append = TRUE)
    }
    ## Model 4. Character effect on diversification, with a hidden state
    if(model == 16){
      speciation <- c(1, 2, 3, 4, 5, 6)
      extirpation <- c(1, 2, 3, 4, 5, 6)
      trans.rate <- TransMatMakerGeoHiSSE(hidden.areas = 1, include.jumps = FALSE, separate.extirpation = TRUE)
      model16 <- try(GeoHiSSE(phy, ranges, f = f, speciation = speciation, extirpation = extirpation, 
                          hidden.areas = TRUE, trans.rate = trans.rate, assume.cladogenetic = FALSE))
      saveRDS(model16, file = paste0(outfile, i, "_model16.rds"))
      
      capture.output(print(paste0("Analysis - ", outfile, i, " - model 16 done.")), file = paste0(outfile, i, ".log"), append = TRUE)
    }
    ## Model 5. No character effect on diversification, multiple shifts
    if(model == 17){
      speciation <- c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3))
      extirpation <- c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3))
      trans.rate <- TransMatMakerGeoHiSSE(hidden.areas = 4, include.jumps = FALSE, separate.extirpation = TRUE, make.null = TRUE)
      model17 <- try(GeoHiSSE(phy, ranges, f = f, speciation = speciation, extirpation = extirpation, 
                              hidden.areas = TRUE, trans.rate = trans.rate, assume.cladogenetic = FALSE)
      saveRDS(model17, file = paste0(outfile, i, "_model17.rds"))
      
      capture.output(print(paste0("Analysis - ", outfile, i, " - model 17 done.")), file = paste0(outfile, i, ".log"), append = TRUE)
    }
    ## Model 6*. No character effect on diversification, multiple shifts.
    if(model == 18){
      speciation <- c(rep(1, 3), rep(2, 3))
      extirpation <- c(rep(1, 3), rep(2, 3))
      trans.rate <- TransMatMakerGeoHiSSE(hidden.areas = 1, include.jumps = FALSE, separate.extirpation = TRUE, make.null = TRUE)
      model18 <- try(GeoHiSSE(phy, ranges, f = f, speciation = speciation, extirpation = extirpation, 
                          hidden.areas = TRUE, trans.rate = trans.rate, assume.cladogenetic = FALSE))
      saveRDS(model18, file = paste0(outfile, i, "_model18.rds"))
      
      capture.output(print(paste0("Analysis - ", outfile, i, " - model 18 done.")), file = paste0(outfile, i, ".log"), append = TRUE)
    }
  }
  print("Analyses with GeoHiSSE model completed ! Check your results in the working directory.")
}

