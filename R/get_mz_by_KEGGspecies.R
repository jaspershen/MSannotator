get_mz_by_KEGGspecies <-
function(kegg_species_code="hsa",queryadductlist=c("M+H"),syssleep=0.01){
    
    
    data(adduct_table)
    adduct_table<-as.data.frame(adduct_table)
    
    
    path_list<-keggList(kegg_species_code,database="pathway")
    
    path_ids<-names(path_list)
    
    path_ids<-gsub(path_ids,pattern="path:",replacement="")
    
    map_res<-get_mz_by_KEGGpathwayIDs(kegg_pathway_ids=path_ids,outloc=outloc,queryadductlist=queryadductlist,syssleep=syssleep)
 
     map_res<-map_res[,c(1:6)]
     colnames(map_res)<-c("mz","KEGG.Compound.ID","Name","Formula","MonoisotopicMass","Adduct")
     return(map_res)
}
