multilevelannotationstep3 <-
  function(outloc1,
           adduct_weights = NA,
           num_sets = NA,
           boostIDs = NA,
           pathwaycheckmode = "p",
           dbAllinf = NA) {
    setwd(outloc1)
    browser()
    if (is.na(num_sets) == FALSE) {
      num_sets_1 = num_sets
    } else{
      num_sets_1 = NA
    }
    load("step1_results.Rda")
    rm(dataA)
    setwd(outloc1)

    num_sets <- length(chemids_split)
    #load("xMSannotator_levelABC.Rda")
    if (is.na(num_sets_1) == FALSE) {
      num_sets = num_sets_1
    }


    rm(dataA)
    rm(levelA_res1)

    #print(ls())

    if (num_sets >= length(chemids_split)) {
      num_sets <- length(chemids_split)


    }



    outloc <- outloc1


    if (is.na(adduct_weights) == TRUE) {
      data(adduct_weights)
      adduct_weights <- as.data.frame(adduct_weights)
      #print(dim(adduct_weights))

      adduct_weights1 <- matrix(nrow = 2, ncol = 2, 0)
      adduct_weights1[1, ] <- c("M+H", 1)
      adduct_weights1[2, ] <- c("M-H", 1)
      adduct_weights <- as.data.frame(adduct_weights1)
      colnames(adduct_weights) <- c("Adduct", "Weight")

    }
    outloc1 <- paste(outloc, "/stage2/", sep = "")
    dir.create(outloc1)
    setwd(outloc1)

    chemscoremat <- {
    }

    #print("num_sets")
    #print(num_sets)
    #print(length(chemids_split))

    for (sind in seq(1, num_sets))
    {
      cur_fname <- paste("chem_score", sind, ".Rda", sep = "")

      load(cur_fname)


      #curchemscoremat<-curchemscoremat[,1:14]
      #print(dim(curchemscoremat))
      #print(colnames(curchemscoremat))
      chemscoremat <- rbind(chemscoremat, curchemscoremat)

    }

    chemscoremat <- as.data.frame(chemscoremat)

    write.table(chemscoremat,
                file = "xMSannotator_scoremat_stage2.txt",
                sep = "\t",
                row.names = FALSE)

    #print("Final chemscoremat")
    #print(dim(chemscoremat))
    otherchems <- mchemdata[which(chemscoremat$mz %in% chemscoremat)]

    cnames <- colnames(chemscoremat)
    cnames[1] <- "score"
    colnames(chemscoremat) <- cnames

    #chemscoremat$score[which(chemscoremat$Adduct%in%adduct_weights[,1])]<-chemscoremat$score[which(chemscoremat$Adduct%in%adduct_weights[,1])]*10^(adduct_weights[which(adduct_weights[,1]%in%chemscoremat$Adduct),2])
    #

    #chemscoremat$score<-(as.numeric(chemscoremat$score)/1000000)
    #chemscoremat$score[which(chemscoremat$Adduct%in%adduct_weights[,1])]<-(chemscoremat$score[which(chemscoremat$Adduct%in%adduct_weights[,1])]*(10^as.numeric(as.character(adduct_weights[which(adduct_weights[,1]%in%chemscoremat$Adduct),2]))))

    #chemscoremat$score<-round(chemscoremat$score,1)
    #write.table(chemscoremat,file="xMSannotator_scoremat_stage1.txt",sep="\t",row.names=FALSE)



    outloc2 <- paste(outloc, "/stage3/", sep = "")
    dir.create(outloc2)
    setwd(outloc2)

    scorethresh <- 0 #corthresh*(1/(max_diff_rt))*2

    if (db_name == "KEGG")
    {
      data(keggotherinf)


      #keggotherinf<-dbAllinf

      m1 <- apply(keggotherinf, 1, function(x) {
        chemid <- x[1]

        g1 <- gregexpr(x[4], pattern = "map")
        regexp_check <- attr(g1[[1]], "match.length")
        if (regexp_check[1] < 0) {
          pathid = "-"

          return(cbind(chemid, pathid))
        } else{
          pathid <- strsplit(x = x[4], split = ";")

          pathid <- unlist(pathid)
          return(cbind(chemid, pathid))
        }
      })

      m2 <- ldply(m1, rbind)

      chemscoremat <-
        merge(chemscoremat,
              keggotherinf,
              by.x = "chemical_ID",
              by.y = "KEGGID")


      #write.table(chemscoremat,file=paste("xMSannotator_",db_name,"scorematotherinf_stage1.txt",sep=""),sep="\t",row.names=FALSE)
      #chemscoremat<-read.table("xMSannotator_KEGGscorematotherinf_stage1.txt",sep="\t",header=TRUE)



      chemids <- as.character(chemscoremat$chemical_ID)

      pathway_ids <- as.character(m2[, 2])

      pathway_ids <- unique(pathway_ids)

      module_num <-
        gsub(chemscoremat$Module_RTclust,
             pattern = "_[0-9]*",
             replacement = "")

      chemscoremat <- cbind(chemscoremat, module_num)

      chemscoremat_orig <- chemscoremat

      chemscoremat <- chemscoremat_orig




      for (path_id in pathway_ids)
      {
        if (path_id != "-" & path_id != "map01100") {
          pathway_chemicals <-
            m2[which(m2[, 2] %in% path_id), 1] #m2[which(m2[,2]%in%pathwayscurchemical),1]


          curmchemicaldata1 <-
            chemscoremat[which(
              chemscoremat$chemical_ID %in% pathway_chemicals &
                chemscoremat$score >= 0.1 &
                chemscoremat$Adduct %in% as.character(adduct_weights[, 1])
            ), ]

          # hmdbAllinfv3.5[which(as.character(hmdbAllinfv3.5$hmdbid)==as.character(curmchemicaldata$chemical_ID[4])),]

          num_chems <- length(unique(curmchemicaldata1$chemical_ID)) ^ 5

          t1 <- table(curmchemicaldata1$module_num)

          module_counts <- t1[which(t1 > 0)]

          module_names <- names(module_counts)

          pathway_chemicals_1 <- curmchemicaldata1$chemical_ID

          for (c in pathway_chemicals_1) {
            curmchemicaldata <-
              curmchemicaldata1[which(as.character(curmchemicaldata1$chemical_ID) == c), ]

            #chemscoremat[which(as.character(chemscoremat$chemical_ID)==c),] #& chemscoremat$Adduct%in%as.character(adduct_weights[,1])),]
            #cur_module<-max(curmchemicaldata$module_num)[1]

            t2 <- table(curmchemicaldata$module_num)
            cur_module <- names(t2[which(t2 == max(t2)[1])])

            mzid_cur <-
              paste(curmchemicaldata$mz, curmchemicaldata$time, sep = "_")

            dweights <- alldegrees[which(mzid %in% mzid_cur), 1]

            if (nrow(curmchemicaldata) > 0) {
              pathwayscurchemical <- m2[which(m2[, 1] == c), 2]
              #for(cur_module in cur_modules)
              {
                if (pathwaycheckmode == "pm") {
                  num_chems <- t1[as.character(cur_module)]
                }
                #if(num_chems>1 && length(which(dweights>min(alldegrees[,1])))>=1){
                #	num_chems=10
                #}

                if (is.na(curmchemicaldata$score[1]) == TRUE) {
                  diff_rt <- max(curmchemicaldata$time) - min(curmchemicaldata$time)

                  if (diff_rt > max_diff_rt) {
                    if (length(which(t2 > 1)) == 1) {
                      curmchemicaldata$score <- rep(0.1, length(curmchemicaldata$score))
                    } else{
                      curmchemicaldata$score <- rep(0, length(curmchemicaldata$score))
                    }
                  } else{
                    curmchemicaldata$score <- rep(0, length(curmchemicaldata$score))
                  }
                }

                if (curmchemicaldata$score[1] < 0.01) {
                  #chemscoremat$score[which(as.character(chemscoremat$chemical_ID)==c & chemscoremat$module_num==cur_module & chemscoremat$Adduct%in%as.character(adduct_weights[,1]))]=as.numeric(chemscoremat$score[which(as.character(chemscoremat$chemical_ID)==c & chemscoremat$Adduct%in%as.character(adduct_weights[,1]))][1])+num_chems

                  chemscoremat$score[which(
                    as.character(chemscoremat$chemical_ID) == c &
                      chemscoremat$Adduct %in% as.character(adduct_weights[, 1])
                  )] = as.numeric(chemscoremat$score[which(
                    as.character(chemscoremat$chemical_ID) == c &
                      chemscoremat$Adduct %in% as.character(adduct_weights[, 1])
                  )][1]) + num_chems

                } else{
                  #chemscoremat$score[which(as.character(chemscoremat$chemical_ID)==c & chemscoremat$module_num==cur_module)]=chemscoremat$score[which(as.character(chemscoremat$chemical_ID)==c)][1]+num_chems

                  chemscoremat$score[which(as.character(chemscoremat$chemical_ID) ==
                                             c)] = chemscoremat$score[which(as.character(chemscoremat$chemical_ID) ==
                                                                              c)][1] + num_chems


                }
              }
              #chemscoremat[which(as.character(chemscoremat$chemical_ID)==c),]
            }
          }

        }
      }
      cnames <- colnames(chemscoremat)

      cnames <- gsub(cnames, pattern = ".x", replacement = "")

      colnames(chemscoremat) <- cnames

      #write.table(chemscoremat,file=paste("xMSannotator_",db_name,"scorematotherinf_stage2.txt",sep=""),sep="\t",row.names=FALSE)
    }
    else{
      if (db_name == "HMDB")
      {
        data(hmdbAllinf)

        hmdbAllinfv3.5 <- hmdbAllinf[, -c(26:27)] #dbAllinf[,-c(26:27)]
        rm(hmdbAllinf)

        #hmdbAllinfv3.5<-gsub(x=hmdbAllinfv3.5,pattern="\n[\\s]+",replacement="")

        #replace(as.matrix(hmdbAllinfv3.5),which(hmdbAllinfv3.5=="\n"),"")

        m1 <- apply(hmdbAllinfv3.5, 1, function(x) {
          chemid <- x[1]

          g1 <- gregexpr(x[14], pattern = "SM")
          regexp_check <- attr(g1[[1]], "match.length")
          if (regexp_check[1] < 0) {
            pathid = "-"

            return(cbind(chemid, pathid))
          } else{
            pathid <- strsplit(x = x[14], split = ";")

            pathid <- unlist(pathid)
            return(cbind(chemid, pathid))
          }
        })

        m2 <- ldply(m1, rbind)

        chemscoremat <-
          merge(chemscoremat,
                hmdbAllinfv3.5,
                by.x = "chemical_ID",
                by.y = "HMDBID")

        #write.table(chemscoremat,file=paste("xMSannotator_",db_name,"scoremat_stage1.txt",sep=""),sep="\t",row.names=FALSE)


        chemids <- as.character(chemscoremat$chemical_ID)

        pathway_ids <- as.character(m2[, 2])

        pathway_ids <- unique(pathway_ids)

        module_num <-
          gsub(chemscoremat$Module_RTclust,
               pattern = "_[0-9]*",
               replacement = "")

        chemscoremat <- cbind(chemscoremat, module_num)

        chemscoremat_orig <- chemscoremat

        chemscoremat <- chemscoremat_orig

        for (path_id in pathway_ids)
        {
          if (path_id != "-") {
            pathway_chemicals <-
              m2[which(m2[, 2] %in% path_id), 1] #m2[which(m2[,2]%in%pathwayscurchemical),1]


            #curmchemicaldata<-chemscoremat[which(chemscoremat$chemical_ID%in%pathway_chemicals & chemscoremat$Adduct%in%as.character(adduct_weights[,1])),]

            curmchemicaldata <-
              chemscoremat[which(
                chemscoremat$chemical_ID %in% pathway_chemicals &
                  chemscoremat$score >= 0.1 &
                  chemscoremat$Adduct %in% as.character(adduct_weights[, 1])
              ), ]

            # hmdbAllinfv3.5[which(as.character(hmdbAllinfv3.5$hmdbid)==as.character(curmchemicaldata$chemical_ID[4])),]

            num_chems <- length(unique(curmchemicaldata$chemical_ID)) ^ 5

            t1 <- table(curmchemicaldata$module_num)

            module_counts <- t1[which(t1 > 0)]

            module_names <- names(module_counts)

            for (c in pathway_chemicals) {
              curmchemicaldata <-
                chemscoremat[which(as.character(chemscoremat$chemical_ID) == c), ] #& chemscoremat$Adduct%in%as.character(adduct_weights[,1])),]

              #cur_module<-max(curmchemicaldata$module_num)[1]
              t2 <- table(curmchemicaldata$module_num)

              cur_module <- names(t2[which(t2 == max(t2)[1])])

              mzid_cur <-
                paste(curmchemicaldata$mz, curmchemicaldata$time, sep = "_")

              dweights <-
                alldegrees[which(mzid %in% mzid_cur), 1]

              if (nrow(curmchemicaldata) > 0) {
                pathwayscurchemical <- m2[which(m2[, 1] == c), 2]
                #for(cur_module in cur_modules)
                {
                  if (pathwaycheckmode == "pm") {
                    num_chems <- t1[as.character(cur_module)]
                  }
                  # if(num_chems>1 && length(which(dweights>min(alldegrees[,1])))>=1){
                  #        num_chems=10
                  #}

                  if (is.na(curmchemicaldata$score[1]) ==
                      TRUE) {
                    diff_rt <- max(curmchemicaldata$time) - min(curmchemicaldata$time)

                    if (diff_rt > max_diff_rt) {
                      if (length(which(t2 > 1)) == 1) {
                        curmchemicaldata$score <- rep(0.1, length(curmchemicaldata$score))
                      } else{
                        curmchemicaldata$score <- rep(0, length(curmchemicaldata$score))
                      }
                    } else{
                      curmchemicaldata$score <- rep(0, length(curmchemicaldata$score))
                    }
                  }

                  if (curmchemicaldata$score[1] <
                      0.01) {
                    #chemscoremat$score[which(as.character(chemscoremat$chemical_ID)==c & chemscoremat$module_num==cur_module & chemscoremat$Adduct%in%as.character(adduct_weights[,1]))]=as.numeric(chemscoremat$score[which(as.character(chemscoremat$chemical_ID)==c & chemscoremat$Adduct%in%as.character(adduct_weights[,1]))][1])+num_chems

                    chemscoremat$score[which(
                      as.character(chemscoremat$chemical_ID) == c &
                        chemscoremat$Adduct %in% as.character(adduct_weights[, 1])
                    )] = as.numeric(chemscoremat$score[which(
                      as.character(chemscoremat$chemical_ID) == c &
                        chemscoremat$Adduct %in% as.character(adduct_weights[, 1])
                    )][1]) + num_chems

                  } else{
                    #chemscoremat$score[which(as.character(chemscoremat$chemical_ID)==c & chemscoremat$module_num==cur_module)]=chemscoremat$score[which(as.character(chemscoremat$chemical_ID)==c)][1]+num_chems

                    chemscoremat$score[which(as.character(chemscoremat$chemical_ID) ==
                                               c)] = chemscoremat$score[which(as.character(chemscoremat$chemical_ID) ==
                                                                                c)][1] + num_chems


                  }
                }


              }
            }

          }
        }
        cnames <- colnames(chemscoremat)

        cnames <- gsub(cnames, pattern = ".x", replacement = "")

        colnames(chemscoremat) <- cnames

        #write.table(chemscoremat,file=paste("xMSannotator_",db_name,"scoremat_stage2.txt",sep=""),sep="\t",row.names=FALSE)
      } else{
        # if(db_name=="T3DB"){

        #load("/Users/karanuppal/Documents/Emory/JonesLab/Projects/xMSannotator/t3dbotherinf.Rda")
        #chemscoremat_otherinf<-merge(chemscoremat,t3dbotherinf,by.x="chemical_ID",by.y="T3DB.ID")

        #write.table(chemscoremat_otherinf,file=paste("xMSannotator_",db_name,"scoremat_stage1.txt",sep=""),sep="\t",row.names=FALSE)
        #rm(chemscoremat_otherinf)


        #}

      }
    }


    cnames <- colnames(chemscoremat)

    cnames <- gsub(cnames, pattern = ".x", replacement = "")

    colnames(chemscoremat) <- cnames
    multiresmat <- chemscoremat #mchemdata_orig

    cnames <- colnames(chemscoremat)

    chemscoremat$score[which(chemscoremat$chemical_ID %in% boostIDs)] <-
      max(chemscoremat$score, na.rm = TRUE) * 100



    good_ind <- which(chemscoremat$score >= scorethresh)

    chemscoremat_highconf <- {
    }
    if (length(good_ind) > 0) {
      chemscoremat_highconf <-
        chemscoremat #[which(chemscoremat$score>=scorethresh),]
      rm(chemscoremat)
    }


    chemscoremat_highconf <- chemscoremat_highconf[, c(1:12, 14:15)]

    #save(list=c("chemscoremat","chemscoremat_highconf"),file="chemscore_res.Rda")

    write.table(
      chemscoremat_highconf,
      file = "../Stage3.txt",
      sep = "\t",
      row.names = FALSE
    )



    rm(
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
      "filter.by"
    )

    return(chemscoremat_highconf)
  }
