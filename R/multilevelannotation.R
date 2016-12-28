multilevelannotation <- function(dataA,
           max.mz.diff = 10,
           max.rt.diff = 10,
           cormethod = "pearson",
           num_nodes = 4,
           queryadductlist = c("all"),
           mode = "pos",
           outloc,
           db_name = "HMDB",
           adduct_weights = NA,
           num_sets = 30,
           allsteps = TRUE,
           corthresh = 0.7,
           NOPS_check = TRUE,
           customIDs = NA,
           missing.value = NA,
           hclustmethod = "complete",
           deepsplit = 2,
           networktype = "unsigned",
           minclustsize = 10,
           module.merge.dissimilarity = 0.2,
           filter.by = c("M+H"),
           redundancy_check = TRUE,
           min_ions_perchem = 1,
           biofluid.location = NA,
           origin = NA,
           status = "Detected",
           boostIDs = NA,
           max_isp = 5,
           MplusH.abundance.ratio.check = FALSE,
           customDB = NA,
           HMDBselect = "union",
           mass_defect_window = 0.01,
           mass_defect_mode = "pos",
           pathwaycheckmode = "pm")
  {
    options(warn = -1)
    browser()
    allowWGCNAThreads(nThreads = num_nodes)

    data(adduct_table)

    suppressWarnings(if (is.na(customIDs) == FALSE) {
      customIDs <- as.data.frame(customIDs)

    })

    suppressWarnings(dir.create(outloc))
    setwd(outloc)

    suppressWarnings(suppressWarnings(sink(file = NULL)))


    targeted_feat_raw <- {
    }
    num_nodes <- min(num_nodes, detectCores())

    x <- date()
    x <- strsplit(x, split = " ")
    x1 <- unlist(x)


    x1 <- gsub(x1, pattern = ":", replacement = "_")
    fname <- paste(x1[2:4], collapse = "")
    fname <- gsub(fname, pattern = ":", replacement = "_")
    #fname<-paste(fname,x1[5],sep="")
    x1[5] <- gsub(x1[5], pattern = ":", replacement = "_")
    fname <- paste(fname, x1[5], sep = "_")

    fname <- paste(outloc, "/Log_", fname, ".txt", sep = "")
    print(
      paste(
        "Program is running. Please check the logfile for runtime status: ",
        fname,
        sep = ""
      )
    )

    #不输出结果到txt文件
    # sink(fname)

    dbAllinf = NA
    ##max_diff_rt是同位素峰和加合物峰的RT的最大差值
    max_diff_rt <- max.rt.diff
    cutheight = module.merge.dissimilarity
    time_step = 1
    step1log2scale = FALSE

    adduct_table$Adduct <-
      gsub(adduct_table$Adduct,
           pattern = "2M\\+2H",
           replacement = "M+H")
    adduct_table$Adduct <-
      gsub(adduct_table$Adduct,
           pattern = "3M\\+3H",
           replacement = "M+H")

    adduct_table$Adduct <-
      gsub(adduct_table$Adduct,
           pattern = "2M\\+2Na",
           replacement = "M+Na")
    adduct_table$Adduct <-
      gsub(adduct_table$Adduct,
           pattern = "3M\\+3Na",
           replacement = "M+Na")

    adduct_table$Adduct <-
      gsub(adduct_table$Adduct,
           pattern = "2M\\+2K",
           replacement = "M+K")
    adduct_table$Adduct <-
      gsub(adduct_table$Adduct,
           pattern = "3M\\+3K",
           replacement = "M+K")

    adduct_table <- unique(adduct_table)



    if (queryadductlist == "all" & mode == "pos") {
      print(
        "Warning: You are using all positive adduct forms. Some adducts might not be applicable based on the gradient."
      )
      adduct_names <-
        adduct_table$Adduct[(adduct_table$Mode == "positive")]

      adduct_table <-
        adduct_table[which(adduct_table$Adduct %in% adduct_names),]

    } else{
      if (queryadductlist == "all" & mode == "neg") {
        print(
          "Warning: You are using all negative adduct forms. Some adducts might not be applicable based on the gradient."
        )
        adduct_names <-
          adduct_table$Adduct[(adduct_table$Mode == "negative")]

        #-adduct_table$Adduct[(adduct_table$Type=="S" & adduct_table$Mode=="negative") | (adduct_table$Type==gradienttype & adduct_table$Mode=="negative")]
        adduct_table <-
          adduct_table[which(adduct_table$Adduct %in% adduct_names),]
      } else{
        adduct_names <-
          adduct_table$Adduct[which(adduct_table$Adduct %in% queryadductlist)]

        adduct_table <-
          adduct_table[which(adduct_table$Adduct %in% adduct_names),]

      }
    }

    adduct_names <- unique(adduct_names)

    #print("Using the following adducts for annotation:")
    #print(adduct_names)

    if (db_name == "HMDB") {
      #print("Using HMDB database")
      #data(hmdbAllinfv3.5)


      data(hmdbAllinf)
      dbAllinf <- hmdbAllinf
      rm(hmdbAllinf)
      hmdbAllinfv3.6 = dbAllinf


      #load("/home/kuppal2/Documents/Projects/xMSannotator/hmdbAllinfv3.6.Rda")

      #load("/Users/karanuppal/Documents/Emory/JonesLab/Projects/xMSannotator/hmdbAllinfv3.6.Rda")

      suppressWarnings(if (is.na(customIDs) == TRUE) {
        if (is.na(biofluid.location) == FALSE & is.na(origin) == TRUE) {
          gres <-
            gregexpr(
              hmdbAllinfv3.6$BioFluidLocation,
              pattern = biofluid.location,
              ignore.case = TRUE
            )

          if (is.na(status) == FALSE) {
            gres3 <-
              gregexpr(hmdbAllinfv3.6$HMDBStatus,
                       pattern = status,
                       ignore.case = TRUE)

            gres3_good <- which(gres3 == 1)
            gres_good <- which(gres == 1)
            if (HMDBselect == "intersect") {
              gres <- intersect(gres_good, gres3_good)
            } else{
              gres <- c(gres_good, gres3_good)
              gres <- unique(gres)

            }


          } else{
            gres <- which(gres == 1)

          }


          customIDs <- hmdbAllinfv3.6[gres, c(1, 20)]




        } else{
          if (is.na(biofluid.location) == FALSE & is.na(origin) == FALSE) {
            gres <-
              gregexpr(
                hmdbAllinfv3.6$BioFluidLocation,
                pattern = biofluid.location,
                ignore.case = TRUE
              )

            gres2 <-
              gregexpr(hmdbAllinfv3.6$Origin,
                       pattern = origin,
                       ignore.case = TRUE)

            gres_good <- which(gres == 1)
            gres2_good <- which(gres2 == 1)




            if (is.na(status) == FALSE) {
              gres3 <-
                gregexpr(
                  hmdbAllinfv3.6$HMDBStatus,
                  pattern = status,
                  ignore.case = TRUE
                )

              gres3_good <- which(gres3 == 1)


              if (HMDBselect == "intersect") {
                gres <- intersect(gres_good, gres2_good, gres3_good)
              } else{
                gres <- c(gres_good, gres2_good, gres3_good)
                gres <- unique(gres)

              }
            } else{
              if (HMDBselect == "intersect") {
                gres <- intersect(gres_good, gres2_good)
              } else{
                gres <- c(gres_good, gres2_good)
                gres <- unique(gres)

              }
            }





            customIDs <- hmdbAllinfv3.6[gres, c(1, 20)]



          } else{
            if (is.na(biofluid.location) == TRUE & is.na(origin) == FALSE) {
              gres <-
                gregexpr(hmdbAllinfv3.6$Origin,
                         pattern = origin,
                         ignore.case = TRUE)


              if (is.na(status) == FALSE) {
                gres3 <-
                  gregexpr(
                    hmdbAllinfv3.6$HMDBStatus,
                    pattern = status,
                    ignore.case = TRUE
                  )

                gres3_good <- which(gres3 == 1)
                gres_good <- which(gres == 1)
                #gres<-intersect(gres_good,gres3_good)
                if (HMDBselect == "intersect") {
                  gres <- intersect(gres_good, gres3_good)
                } else{
                  gres <- c(gres_good, gres3_good)
                  gres <- unique(gres)

                }
              } else{
                gres <- which(gres == 1)

              }


              #customIDs<-hmdbAllinfv3.6[gres,c(1)]
              customIDs <- hmdbAllinfv3.6[gres, c(1, 20)]

            } else{
              if (is.na(status) == FALSE) {
                gres3 <-
                  gregexpr(
                    hmdbAllinfv3.6$HMDBStatus,
                    pattern = status,
                    ignore.case = TRUE
                  )

                gres3_good <- which(gres3 == 1)

                gres <- gres3_good




                customIDs <- hmdbAllinfv3.6[gres, c(1, 20)]

              }


            }


          }
        }

      })


      data(hmdbCompMZ)


      suppressWarnings(if (is.na(customIDs) == FALSE) {
        customIDs <- unique(customIDs)


        hmdbCompMZ <-
          hmdbCompMZ[which(hmdbCompMZ$HMDBID %in% customIDs[, 1]),]

        hmdbCompMZ <-
          hmdbCompMZ[which(hmdbCompMZ$Adduct %in% adduct_names),]
      })

      hmdbCompMZ <-
        hmdbCompMZ[which(hmdbCompMZ$Adduct %in% adduct_names),]

      chemCompMZ <- hmdbCompMZ


      print("Dimension of precomputed HMDB m/z database")
      print(dim(chemCompMZ))

      rm(hmdbCompMZ)

    } else{
      if (db_name == "KEGG") {
        #print("Using KEGG database")

        data(keggCompMZ)

        #data(keggotherinf)
        #dbAllinf<-keggotherinf
        #rm(keggotherinf) #<-dbAllinf



        suppressWarnings(if (is.na(customIDs) == FALSE) {
          customIDs <- unique(customIDs)


          keggCompMZ <-
            keggCompMZ[which(keggCompMZ$KEGGID %in% customIDs[, 1]),]

          keggCompMZ <-
            keggCompMZ[which(keggCompMZ$Adduct %in% adduct_names),]
        })

        keggCompMZ <-
          keggCompMZ[which(keggCompMZ$Adduct %in% adduct_names),]
        chemCompMZ <- keggCompMZ

        print("Dimension of precomputed KEGG m/z database")
        print(dim(chemCompMZ))

        rm(keggCompMZ)
      } else{
        if (db_name == "LipidMaps") {
          data(lipidmapsCompMZ)
          lipidmapsCompMZ <-
            lipidmapsCompMZ[which(lipidmapsCompMZ$Adduct %in% adduct_names),]
          chemCompMZ <- lipidmapsCompMZ

          print("Dimension of precomputed LipidMaps m/z database")
          print(dim(chemCompMZ))

          rm(lipidmapsCompMZ)

        } else{
          if (db_name == "T3DB") {
            #load("/Users/karanuppal/Documents/Emory/JonesLab/Projects/xMSannotator/t3dbotherinf.Rda")
            #load("/Users/karanuppal/Documents/Emory/JonesLab/Projects/xMSannotator/t3dbCompMZ.Rda")

            #t3dbotherinf<-dbAllinf

            data(t3dbCompMZ)
            t3dbCompMZ <-
              t3dbCompMZ[which(t3dbCompMZ$Adduct %in% adduct_names),]
            chemCompMZ <- t3dbCompMZ

            print("Dimension of precomputed T3DB m/z database")
            print(dim(chemCompMZ))

            rm(t3dbCompMZ)

          } else{
            if (db_name == "Custom") {
              inputmassmat <- customDB

              if (length(which(duplicated(inputmassmat[, 1]) == TRUE)) > 0) {
                inputmassmat <-
                  inputmassmat[-which(duplicated(inputmassmat[, 1]) == TRUE),]
              }
              mz_search_list <- {

              }
              mz_search_list <-
                lapply(1:dim(inputmassmat)[1], function(m) {
                  mz_search_list <-
                    get_mz_by_monoisotopicmass(
                      monoisotopicmass = as.numeric(as.character(inputmassmat[m, 4])),
                      dbid = inputmassmat[m, 1],
                      name = as.character(inputmassmat[m, 2]),
                      formula = as.character(inputmassmat[m, 3]),
                      queryadductlist = adduct_names
                    )

                  return(mz_search_list)
                })

              mz_search_list <- ldply(mz_search_list, rbind)

              chemCompMZ <- mz_search_list
              rm(mz_search_list)
              rm(customDB)


            } else{
              stop("db_name should be: KEGG, HMDB, T3DB, LipidMaps, or Custom")

            }


          }


        }
      }
    }
    data(adduct_table)

    if (is.na(adduct_weights) == TRUE) {
      data(adduct_weights)
      adduct_weights1 <- matrix(nrow = 2, ncol = 2, 0)
      adduct_weights1[1,] <- c("M+H", 1)
      adduct_weights1[2,] <- c("M-H", 1)
      adduct_weights <- as.data.frame(adduct_weights1)
      colnames(adduct_weights) <- c("Adduct", "Weight")
      adduct_weights <- as.data.frame(adduct_weights)

    }


    cnames <- colnames(chemCompMZ)

    cnames[2] <- "chemical_ID"
    colnames(chemCompMZ) <- cnames


    l1 <- list.files(outloc)

    check_step1 <- which(l1 == "step1_results.Rda")


    if (length(check_step1) < 1) {
      print("Status 1: Computing modules using WGCNA")

      check_levelA <- which(l1 == "xMSannotator_levelA_modules.Rda")

      if (is.na(missing.value) == FALSE) {
        dataA <-
          replace(as.matrix(dataA), which(dataA == missing.value), NA)
        dataA <- as.data.frame(dataA)
      }

      system.time(global_cor <-
                    WGCNA::cor(
                      t(dataA[, -c(1:2)]),
                      nThreads = num_nodes,
                      method = cormethod,
                      use = 'p'
                    ))

      setwd(outloc)

      mzid <- paste(dataA$mz, dataA$time, sep = "_")
      colnames(global_cor) <- mzid
      rownames(global_cor) <- mzid

      global_cor <- round(global_cor, 3)

      dataA <- unique(dataA)

      clustmethod = "WGCNA"
      if (length(check_levelA) < 1) {
        if (clustmethod == "WGCNA") {
          levelA_res <-
            get_peak_blocks_modulesvhclust(
              dataA,
              simmat = global_cor,
              adjacencyfromsimilarity = TRUE,
              time_step = time_step,
              max.rt.diff = max_diff_rt,
              outloc,
              column.rm.index = NA,
              cor.thresh = NA,
              deepsplit = deepsplit,
              minclustsize = minclustsize,
              cutheight = cutheight,
              cormethod = cormethod,
              networktype = networktype,
              num_nodes = num_nodes,
              step1log2scale = step1log2scale,
              hclustmethod = hclustmethod
            )

browser()
          alldegrees <- levelA_res[, c(1:4)]


          levelA_res <- levelA_res[, -c(1:4)]

        } else{
          if (clustmethod == "graph")
          {
            levelA_res <-
              get_peak_blocks_graph(
                dataA,
                simmat = global_cor,
                adjacencyfromsimilarity = TRUE,
                time_step = 3,
                max.rt.diff = max_diff_rt,
                outloc,
                column.rm.index = NA,
                cor.thresh = NA,
                cormethod = cormethod,
                networktype = networktype,
                num_nodes = num_nodes
              )

          } else{
            if (clustmethod == "hclust") {
              levelA_res <-
                get_peak_blocks_hclust(
                  dataA,
                  time_step = time_step,
                  max.rt.diff = max_diff_rt,
                  outloc,
                  column.rm.index = NA,
                  cor.thresh = NA,
                  deepsplit = deepsplit,
                  minclustsize = minclustsize,
                  cutheight = cutheight,
                  cormethod = cormethod
                )
            }

          }

        }

        setwd(outloc)
        save(levelA_res, file = "xMSannotator_levelA_modules.Rda")
        save(alldegrees, file = "alldegrees.Rda")
        write.table(levelA_res,
                    file = "Stage1.txt",
                    sep = "\t",
                    row.names = FALSE)

      } else{
        setwd(outloc)
        load("xMSannotator_levelA_modules.Rda")
        if (clustmethod == "WGCNA")
        {
          load("alldegrees.Rda")
        }
      }

      setwd(outloc)


      randsetsize <- 1000
      if (dim(dataA)[1] < 1000) {
        randsetsize <- dim(dataA)[1]
      }

      randmzset <- sample(x = 1:dim(dataA)[1], size = randsetsize)
      quant1 <-
        apply(global_cor[randmzset,], 1, function(x) {
          return(quantile(x, 0.99, na.rm = TRUE))
        })

      #corthresh<-max(0.6,mean(quant1,na.rm=TRUE))

      chemCompMZ <- as.data.frame(chemCompMZ)

      setwd(outloc)
      l1 <- list.files(outloc)
      check_levelB <- which(l1 == "xMSannotator_levelB.Rda")

      if (length(check_levelB) < 1) {
        #mz.annot.res<-parLapply(cl,s1,Annotationbychemical_IDschild,dataA,queryadductlist,adduct_table,max.mz.diff,outloc, otherdbs,otherinfo,keggCompMZ)

        cl <- makeSOCKcluster(num_nodes)

        clusterEvalQ(cl, library(XML))
        clusterEvalQ(cl, library(R2HTML))
        clusterEvalQ(cl, library(RCurl))
        clusterEvalQ(cl, library(SSOAP))
        clusterEvalQ(cl, library(limma))

        clusterEvalQ(cl, library(plyr))

        clusterEvalQ(cl, "processWSDL")
        clusterEvalQ(cl, library(png))
        clusterExport(cl, "Annotationbychemical_IDschild")
        #		clusterExport(cl, "chemCompMZ")
        clusterExport(cl, "find.Overlapping.mzs")
        #clusterEvalQ(cl, "chemCompMZ")
        clusterExport(cl, "getVenn")
        if (FALSE) {
          l2 <- mclapply(1:length(adduct_names), function(j) {
            annot_res <-
              Annotationbychemical_IDschild(
                dataA = dataA,
                queryadductlist = c(adduct_names[j]),
                adduct_table = adduct_table,
                max.mz.diff = max.mz.diff,
                outloc = outloc,
                keggCompMZ = chemCompMZ,
                otherdbs = FALSE,
                otherinfo = FALSE
              )

            Sys.sleep(4)
            return(annot_res)
          }, mc.cores = num_nodes, mc.preschedule = FALSE)
        }

        s1 <- seq(1, length(adduct_names))
        print("Status 2: Mapping m/z to metabolites:")


        l2 <-
          parLapply(
            cl,
            s1,
            Annotationbychemical_IDschild,
            dataA = dataA,
            queryadductlist = c(adduct_names),
            adduct_table = adduct_table,
            max.mz.diff = max.mz.diff,
            outloc = outloc,
            keggCompMZ = chemCompMZ,
            otherdbs = FALSE,
            otherinfo = FALSE
          )

        stopCluster(cl)


        rm(chemCompMZ)
        levelB_res <- {

        }

        for (j in 1:length(l2)) {
          if (length(l2[[j]]) > 1) {
            levelB_res <- rbind(levelB_res, l2[[j]])
          }
        }


        rm(l2)

        #m2$mz<-as.numeric(as.character(m2$mz))
        #levelB_res[which(levelB_res$mz>175.117 & levelB_res$mz<175.119),1:8]

        levelB_res$mz <- as.numeric(as.character(levelB_res$mz))

        levelB_res$time <- as.numeric(as.character(levelB_res$time))

        levelB_res <- as.data.frame(levelB_res)


        #save(levelB_res,file="xMSannotator_levelB_1.Rda")

        uniq_formula <- as.character(unique(levelB_res$Formula))

        bad_formula <- which(is.na(uniq_formula) == TRUE)
        if (length(bad_formula) > 0) {
          uniq_formula <- uniq_formula[-c(bad_formula)]
        }

        cl <- makeSOCKcluster(num_nodes)


        clusterExport(cl, "check_golden_rules")
        clusterExport(cl, "check_element")
        #clusterExport(cl, "uniq_formula")
        #clusterExport(cl, "NOPS_check")

        levelB_res_check <-
          parLapply(cl, 1:length(uniq_formula), function(j, uniq_formula, NOPS_check) {
            curformula <- as.character(uniq_formula[j])
            return(check_golden_rules(curformula, NOPS_check = NOPS_check))

          }, uniq_formula = uniq_formula, NOPS_check = NOPS_check)
        stopCluster(cl)

        #save(levelB_res_check,file="xMSannotator_levelB_check.Rda")
        levelB_res_check2 <- ldply(levelB_res_check, rbind)

        levelB_res_check3 <-
          levelB_res_check2[which(levelB_res_check2[, 2] == 1),]


        levelB_res <-
          levelB_res[which(levelB_res$Formula %in% as.character(levelB_res_check3[, 1])),]

        water_adducts <- c("M+H-H2O", "M+H-2H2O", "M-H2O-H")

        water_adduct_ind <-
          which(levelB_res$Adduct %in% water_adducts)

        cl <- makeSOCKcluster(num_nodes)


        clusterExport(cl, "check_element")



        if (length(water_adduct_ind) > 0) {
          levelB_res2 <- levelB_res[c(water_adduct_ind),]

          levelB_res <- levelB_res[-c(water_adduct_ind),]

          sind1 <- seq(1:dim(levelB_res2)[1])

          levelB_res_check3 <- parLapply(cl, sind1, function(j) {
            adduct <- as.character(levelB_res2$Adduct[j])
            curformula <- as.character(levelB_res2$Formula[j])

            numoxygens <- check_element(curformula, "O")

            if (numoxygens > 0) {
              bool_check <- 1
            } else{
              bool_check <- 0
            }

            res <- cbind(curformula, bool_check)
            res <- as.data.frame(res)
            return(res)


          })

          levelB_res_check4 <- ldply(levelB_res_check3, rbind)

          valid_form <- {

          }

          if (length(which(levelB_res_check4[, 2] == 1)) > 0) {
            levelB_res_check4 <-
              levelB_res_check4[which(levelB_res_check4[, 2] == 1),]


            valid_form <-
              which(levelB_res2$Formula %in% as.character(levelB_res_check4[, 1]))
          }
          if (length(valid_form) > 0) {
            levelB_res2 <- levelB_res2[valid_form,]
            levelB_res <- rbind(levelB_res, levelB_res2)
          }

        }

        #save(levelB_res,file="xMSannotator_levelB.Rda")



      } else{
        load("xMSannotator_levelB.Rda")
      }


      no_match <- dataA[-which(dataA$mz %in% levelB_res$mz), 1:2]

      no_match_annot <-
        matrix(nrow = dim(no_match)[1], ncol = 8, "-")

      no_match_1 <-
        cbind(no_match[, 1], no_match_annot, no_match[, -c(1)])

      annot_chid_res_1 <- levelB_res[, 1:10]

      colnames(no_match_1) <- colnames(annot_chid_res_1)

      annot_chid_res_1 <- rbind(annot_chid_res_1, no_match_1)

      last_col_ind <- dim(annot_chid_res_1)[2]

      annot_chid_res_1 <- annot_chid_res_1[, c(1, 10, 2:8)]


      #diffmatA<-as.data.frame(levelA_res)




      levels_A <- levels(levelA_res$Module_RTclust)

      levels_A <- table(levelA_res$Module_RTclust)
      levels_A <- names(levels_A)


      #levelA_res_orig<-levelA_res
      dataA_orig <- dataA

      levelA_res <- levelA_res[order(levelA_res$mz),]

      if (length(which(duplicated(mzid) == TRUE)) > 0) {
        levelA_res <- levelA_res[-which(duplicated(mzid) == TRUE),]

        dataA <- dataA[-which(duplicated(mzid) == TRUE),]
      }

      levelA_res <- levelA_res[order(levelA_res$Module_RTclust),]
      levels_A <- levels(unique(levelA_res$Module_RT))



      module_num <-
        gsub(levelA_res$Module_RTclust,
             pattern = "_[0-9]{1,}",
             replacement = "")


      levelA_res_all <- levelA_res[, c(1:2)]

      orig_module_labels <- levelA_res$Module_RTclust

      levelA_res_all$Module_RTclust <- module_num



      #mzdefect<-1*((levelA_res$mz-round(levelA_res$mz)))

      mzdefect <- 1 * ((levelA_res$mz - floor(levelA_res$mz)))
      #print(summary(mzdefect))


      d1 <-
        density(
          mzdefect,
          bw = "nrd",
          from = min(mzdefect),
          to = (0.01 + max(mzdefect))
        )




      t1 <- d1$x


      gid <- {

      }

      gid <- paste("ISgroup", dim(dataA)[1], sep = "")

      levelA_res2 <- cbind(mzdefect, levelA_res)



      #massdefect_cor_groups<-sapply(list(myData1=levelA_res2),function(x)  split(x,cut(levelA_res2$mzdefect,breaks=seq(-0.5,0.5,0.01))))

      massdefect_cor_groups <-
        sapply(list(myData1 = levelA_res2), function(x)
          split(x, cut(
            levelA_res2$mzdefect, breaks = seq(0, 1, 0.01)
          )))


      diffmatB <- {

      }
      for (gnum in 1:length(massdefect_cor_groups)) {
        cur_group <- {

        }

        #print(dim(massdefect_cor_groups[[gnum]]))
        if (dim(massdefect_cor_groups[[gnum]])[1] > 0) {
          ISgroup <-
            paste("ISgroup",
                  massdefect_cor_groups[[gnum]]$Module_RTclust,
                  gnum,
                  sep = "_")

          #print(massdefect_cor_groups[[gnum]])
          cur_group <- as.data.frame(massdefect_cor_groups[[gnum]])
          cur_group <- cbind(ISgroup, cur_group)


          if (length(cur_group) > 0) {
            cur_group <- cur_group[order(cur_group$mz),]


            if (length(diffmatB) < 1) {
              diffmatB <- cur_group
            } else{
              if (dim(cur_group)[1] > 0) {
                diffmatB <- rbind(diffmatB, cur_group)
              }
            }

          }
        }
      }



      diffmatB <- as.data.frame(diffmatB)



      if (dim(diffmatB)[1] > 0) {
        cnames <- colnames(diffmatB)
        cnames[1] <- "ISGroup"

        colnames(diffmatB) <- cnames
      }


      if (dim(diffmatB)[1] < dim(dataA)[1]) {
        diffmatC <- levelA_res2[-which(levelA_res$mz %in% diffmatB$mz),]

        t1 <- table(diffmatB[, 1])
        isop_last <-
          paste("ISgroup_", diffmatC$Module_RTclust, "_", 1, sep = "")

        diffmatC <- cbind(isop_last, diffmatC)
        colnames(diffmatC) <- colnames(diffmatB)


        diffmatD <- rbind(diffmatB[, c(1:10)], diffmatC[, c(1:10)])
        rm(levelA_res2)
        rm(diffmatC)
        rm(diffmatB)
      } else{
        diffmatD <- diffmatB#[,c(1:10)]
        rm(diffmatB)
      }

      diffmatD <- as.data.frame(diffmatD)
      diffmatD$mz <- as.numeric(diffmatD$mz)

      diffmatD <- diffmatD[order(diffmatD$mz),]

      levelA_res <- levelA_res[order(levelA_res$mz),]

      levelA_res1 <- cbind(diffmatD[, 1], levelA_res)

      mean_int_vec <-
        apply(levelA_res1[, -c(1:4)], 1, function(x) {
          mean(x, na.rm = TRUE)
        })

      isop_res_md <-
        cbind(diffmatD[, c(4, 5, 1, 3)], mean_int_vec, diffmatD[, c(2)])

      colnames(isop_res_md) <-
        c("mz",
          "time",
          "ISgroup",
          "Module_RTclust",
          "AvgIntensity",
          "MD")

      isop_res_md[which(isop_res_md$mz > 238.09 &
                          isop_res_md$mz < 241),]


      #save(isop_res_md,file="isop_res.Rda")




      MD <- diffmatD[, c(2)]
      levelA_res1 <- cbind(levelA_res1[, c(1:4)], mean_int_vec, MD)

      cnames <- colnames(levelA_res1)
      cnames[1] <- "ISgroup"

      colnames(levelA_res1) <- cnames

      #decomposeIsotopes(masses=c(174.1116,175.112),intensities=c(100,5),ppm=5)


      m1 <- merge(annot_chid_res_1, levelA_res1, by = "mz")

      rm(annot_chid_res_1)

      m2 <- m1[order(m1$Module_RTclust),]

      rm(m1)

      #m3<-merge(m2[,c(1:12)],diffmatB,by="Module_RTclust")

      t2 <- table(m2$mz)

      same1 <- which(t2 == 1)

      uniquemz <- names(same1)

      m2$MatchCategory = rep("Multiple", dim(m2)[1])

      m2$MatchCategory[which(m2$mz %in% uniquemz)] <- "Unique"

      munique <- m2[which(m2$MatchCategory == "Unique"),]

      t3chem <- table(munique$chemical_ID)


      setwd(outloc)
      #save(m2,file="xMSannotator_levelABC.Rda")


      tmultimatches <- table(m2$chemical_ID, m2$Adduct)
      smutimatches <- apply(tmultimatches, 1, sum)


      multiresmat <- m2[order(m2$Module_RTclust, m2$chemical_ID),]



      dupmz <-
        multiresmat$mz[which(duplicated(multiresmat$mz) == TRUE)]

      tablemz <- table(multiresmat$Module_RTclust, multiresmat$mz)

      multiresmat$MatchCategory <-
        rep("Multiple", dim(multiresmat)[1])

      multiresmat$MatchCategory[-which(multiresmat$mz %in% dupmz)] <-
        "Unique"


      multiresmat2 <-
        multiresmat[-which(multiresmat$chemical_ID == "-"),]


      multiresmat <- multiresmat2

      tall <- table(multiresmat$chemical_ID, multiresmat$mz)

      tall_mznames <- colnames(tall)



      tall_checkmultichems <- apply(tall, 2, sum)

      tall_multichems <- tall[, which(tall_checkmultichems > 1)]

      names_tall_multichems <- colnames(tall_multichems)

      tall_checkmultimz <- apply(tall, 1, sum)

      tall_multimzperchem <- tall[which(tall_checkmultimz > 1),]
      tall_unimzperchem <- tall[which(tall_checkmultimz == 1),]

      names_tall_multimzperchem <- rownames(tall_multimzperchem)

      chemids <-
        names_tall_multimzperchem #chemids[which(t1chemids>0)]

      single_mz_chemids <- rownames(tall_unimzperchem)


      chemids <- chemids
      chem_score <- new("list")
      i <- 1


      del_mz <- {

      }


      multiresmat_filt <- {

      }

      dataA <- dataA_orig #m2[,-c(3:12)]
      cnames <- colnames(dataA)
      cnames[2] <- "time"
      colnames(dataA) <- as.character(cnames)
      dataA <- unique(dataA)

      cnames <- colnames(multiresmat)
      cnames[10] <- "ISgroup"
      colnames(multiresmat) <- cnames




      level_module_isop_annot <-
        levelA_res1 #merge(isop_res_md,levelB_res,by="mz")


      chem_ids <- table(multiresmat$chemical_ID)
      chem_ids_1 <- chem_ids[which(chem_ids >= 1)]
      chem_ids_1 <- names(chem_ids_1)

      chem_ids_2 <- chem_ids_1

      #mz=((nmol*M)+adductMass)/charge

      #chemids<-c("C00018","C00068","C00086","C00114","C00120","C00153","C00187","C00194","C00213","C00250","C00253","C00255","C00262","C00300","C00314","C00327","C00328","C00366","C00378","C00415","C00440","C00504","C00534","C00544")

      #chemids<-c("C00544")


      chemids <- chem_ids_2


      cnames <- colnames(multiresmat)
      cnames[2] <- "time"
      colnames(multiresmat) <- cnames


      rm(levelA_res)
      rm(levelB_res)
      rm(m2)

      mchemdata <-
        multiresmat   #[which(multiresmat$chemical_ID%in%chemids),]

      rm(multiresmat)


      chemids <- unique(mchemdata$chemical_ID)

      #chemids<-c("HMDB15352","HMDB14017","HMDB60549","HMDB00641", "HMDB02352", "HMDB00472", "HMDB14426", "HMDB00517",
      #"HMDB00162","HMDB00883","HMDB00214","HMDB03405","HMDB01847","HMDB00168", "HMDB00043", "HMDB00123","HMDB00148",
      #"HMDB00161","HMDB00159","HMDB00167","HMDB00168","HMDB00177","HMDB00172","HMDB00181","HMDB00182","HMDB00187",
      #"HMDB00191","HMDB00192","HMDB00214","HMDB00450", "HMDB00504","HMDB00748")


      chemids <- unique(chemids)


      chemscoremat <- {

      }

      if (length(chemids) > 20000) {
        num_sets <- 2000

      }

      list_winsize <- num_sets

      list_size <- round(length(chemids) / list_winsize)


      #list_size<-20 #round(length(chemids)/list_winsize)


      #list_winsize<-round(length(chemids)/list_size)

      if (length(chemids) > list_winsize) {
        g <- seq(1, length(chemids), list_size)
        g <- factor(g)
        chemids_split <- split(1:length(chemids), f = g)
        split_size <- 1:list_winsize

      } else{
        chemids_split <- split(1:length(chemids), f = length(chemids))
        split_size <- c(1)
      }
      #print(length(chemids))

      num_sets = length(chemids_split)

      write.table(mchemdata,
                  file = "Stage2.txt",
                  sep = "\t",
                  row.names = FALSE)

      save(
        list = c(
          "mchemdata",
          "chemids",
          "adduct_table",
          "global_cor",
          "mzid",
          "max_diff_rt",
          "isop_res_md",
          "corthresh",
          "level_module_isop_annot",
          "chemids_split",
          "dataA",
          "corthresh",
          "max.mz.diff",
          "outloc",
          "num_sets",
          "db_name",
          "num_nodes",
          "num_sets",
          "adduct_weights",
          "alldegrees",
          "filter.by",
          "max.rt.diff",
          "max_isp",
          "MplusH.abundance.ratio.check",
          "mass_defect_window",
          "mass_defect_mode"
        ),
        file = "step1_results.Rda"
      )

    } else{
      print("Status 1: Skipping step 1.")
      print("Status 2: Using existing step1_results.Rda file.")

      load("step1_results.Rda")

    }

    rm(chemids)
    rm(global_cor)
    rm(isop_res_md)
    rm(level_module_isop_annot)
    rm(dataA)
    rm(alldegrees)

    if (allsteps == TRUE) {
      print("Status 3: Calculating scores for individual chemicals/metabolites")


      #if(FALSE)

      if (num_sets > 1)
      {
        cl <- makeSOCKcluster(num_nodes)

        #clusterEvalQ(cl, "multilevelannotationstep2")

        clusterExport(cl, "multilevelannotationstep2", envir = environment())

        clusterEvalQ(cl, "library(plyr)")

        browser()
        clusterExport(cl, "get_chemscorev1.6.71", envir = environment())

        clusterExport(cl, "ldply", envir = environment())
        clusterExport(cl, "get_confidence_stage2", envir = environment())
        clusterExport(cl, "find.Overlapping.mzs", envir = environment())
        clusterExport(cl, "check_element", envir = environment())
        clusterExport(cl, "group_by_rt_hist", envir = environment())
        clusterExport(cl, "adduct_table", envir = environment())
        clusterExport(cl, "adduct_weights", envir = environment())


        parLapply(cl, 1:num_sets, function(arg1) {
          suppressWarnings(multilevelannotationstep2(outloc1 = outloc, list_number =
                                                       arg1))
        })

        #max.time.diff=max.rt.diff,filter.by=filter.by,max_isp=max_isp,numnodes=1,MplusH.abundance.ratio.check=MplusH.abundance.ratio.check,mass_defect_window=mass_defect_window,mass_defect_mode=mass_defect_mode

        stopCluster(cl)
      } else{
        #lapply(1:num_sets,function(arg1){

        #suppressWarnings(multilevelannotationstep2(outloc1=outloc,list_number=arg1))})

        cl <- makeSOCKcluster(num_nodes)

        clusterExport(cl, "multilevelannotationstep2", envir = environment())

        clusterEvalQ(cl, "library(plyr)")

        clusterExport(cl, "get_chemscorev1.6.71", envir = environment())

        clusterExport(cl, "ldply", envir = environment())
        clusterExport(cl, "get_confidence_stage2", envir = environment())
        clusterExport(cl, "find.Overlapping.mzs", envir = environment())
        clusterExport(cl, "check_element", envir = environment())
        clusterExport(cl, "group_by_rt_hist", envir = environment())
        clusterExport(cl, "adduct_table", envir = environment())
        clusterExport(cl, "adduct_weights", envir = environment())


        #  parLapply(cl,1,function(arg1){


        #    suppressWarnings(multilevelannotationstep2(outloc1=outloc,list_number=arg1))})

        #max.time.diff=max.rt.diff,filter.by=filter.by,max_isp=max_isp,numnodes=1,MplusH.abundance.ratio.check=MplusH.abundance.ratio.check,mass_defect_window=mass_defect_window,mass_defect_mode=mass_defect_mode

        lapply(1:num_sets, function(arg1) {
          suppressWarnings(multilevelannotationstep2(outloc1 = outloc, list_number =
                                                       arg1))
        })

        stopCluster(cl)

      }
browser()
      print("Status 4: Pathway evaluation")
      suppressWarnings(
        annotres <-
          multilevelannotationstep3(
            outloc,
            adduct_weights = adduct_weights,
            boostIDs = boostIDs,
            pathwaycheckmode = pathwaycheckmode
          )
      )

      print("Status 5: Assigning confidence levels")
      suppressWarnings(
        annotresstage4 <-
          multilevelannotationstep4(
            outloc,
            max_diff_rt = max_diff_rt,
            filter.by = filter.by,
            adduct_weights = adduct_weights,
            min_ions_perchem = min_ions_perchem,
            boostIDs = boostIDs,
            max_isp = max_isp
          )
      )


      if (redundancy_check == TRUE) {
        print("Status 6: Redundancy filtering")


        suppressWarnings(
          annotresstage5 <-
            multilevelannotationstep5(
              outloc,
              max_diff_rt = max_diff_rt,
              filter.by = filter.by,
              adduct_weights = adduct_weights,
              min_ions_perchem = min_ions_perchem,
              boostIDs = boostIDs,
              max_isp = max_isp
            )
        )


        print("Processing complete")
        print("Stage 4 confidence level distribution for unique chemical/metabolite IDs")
        print(table(annotresstage4$Confidence[-which(duplicated(annotresstage4$chemical_ID) ==
                                                       TRUE)]))


        suppressWarnings(if (annotresstage5 != 0) {
          print("Stage 5 confidence level distribution for unique chemical/metabolite IDs")
          print(table(annotresstage5$Confidence[-which(duplicated(annotresstage5$chemical_ID) ==
                                                         TRUE)]))
        })


      } else{
        annotresstage5 <- {

        }

        print("Processing complete")
        print("Stage 4 confidence level distribution for unique chemical/metabolite IDs")
        print(table(annotresstage4$Confidence[-which(duplicated(annotresstage4$chemical_ID) ==
                                                       TRUE)]))


      }


    } else{
      annotres <- mchemdata
    }
    rm(mchemdata)


    outloc2 <- paste(outloc, "/stage2/", sep = "")

    unlink(outloc2, force = TRUE, recursive = TRUE)

    file.remove(dir(outloc2, full.names = TRUE))

    outloc2 <- paste(outloc, "/stage3/", sep = "")
    unlink(outloc2, force = TRUE, recursive = TRUE)
    file.remove(dir(outloc2, full.names = TRUE))
    outloc2 <- paste(outloc, "/stage4/", sep = "")
    file.remove(dir(outloc2, full.names = TRUE))
    unlink(outloc2, force = TRUE, recursive = TRUE)
    outloc2 <- paste(outloc, "/stage5/", sep = "")
    file.remove(dir(outloc2, full.names = TRUE))
    unlink(outloc2, force = TRUE, recursive = TRUE)

    unlink("xMSannotator_levelA_modules.Rda")
    unlink("alldegrees.Rda")

    unlink("step1_results.Rda")
    unlink("plot.pdf")

    suppressWarnings(suppressWarnings(sink(file = NULL)))

    # fname <- paste(outloc, "/README.txt", sep = "")

    # suppressWarnings(suppressWarnings(sink(file = fname)))
    print("Output files description:")
    print(
      "Stage 1 includes clustering of features based on intensity and retention time without annotations"
    )
    print("Stage 2 includes clustering results along with simple m/z based database matching")
    print("Stage 3 includes scores for annotations assigned in stage 2 based on multiple criteria")
    print(
      "Stages 4 and 5 include the confidence levels before and after redundancy (multiple matches) filtering, respectively"
    )
    suppressWarnings(suppressWarnings(sink(file = NULL)))

    print("################")
    print("Processing complete")
    print(paste("README.txt and results are located under ", outloc, sep = ""))

    if (allsteps == TRUE) {
      return(list("stage4" = annotresstage4, "stage5" = annotresstage5))
    } else{
      return(annotres)
    }

  }
