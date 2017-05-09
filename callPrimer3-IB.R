#' call primer3 for a given set of DNAstringSet object
#' Forked from https://gist.github.com/al2na/8540391
#' TODO: Add support for target amplicon region (Maybe as [] in the fasta input)
#' @param seq DNAstring object, one DNA string for the given amplicon
#' @param sequence_target: A space separated list of target pairs: 'starting position,target length starting position2,target length2'
#' @param SSR_exclude_region: Usually same as target, make sure the SSR is not considered as an internal oligo
#' @param size_range default: '151-500'
#' @param Tm melting temprature parameters default:c(55,57,60)
#' @param GC min and max %GC content parameters default: c(30,65)
#' @param name name of the amplicon in chr_start_end format
#' @param primer3 primer3 location
#' @param thermo.param thermodynamic parameters folder
#' @param settings text file for parameters
#' @author Ido Bar forked from Altuna Akalin which  modified Arnaud Krebs' original function
#' @example
#'
#' primers=mapply(.callP3NreadOrg,seq=transcripts$sequence[match(filtered_PSR$Seq_ID, transcripts$transcript_id)],
#'              name = filtered_PSR$Seq_ID,
#'      report = paste0(filtered_PSR$Seq_ID, sprintf("_primer3_report_%s.txt", format(Sys.Date(), "%d_%m_%y"))),
#'               sequence_target = paste(filtered_PSR$start, (filtered_PSR$stop-filtered_PSR$start), sep=","),
#'               MoreArgs=list(size_range='100-350',drive = "I") ,SIMPLIFY = FALSE,USE.NAMES = FALSE)
#'
#'
.callP3NreadOrg<-function(seq,size_range='151-500',Tm=c(55,57,60),GC=c(30,65),name,sequence_target=NULL,
                          SSR_exclude_region=TRUE,report=NULL, primer3_path="primer3",primer_num=3, liberal_bases=TRUE,
        primer3=file.path(primer3_path,"bin","primer3_core.exe"),
        thermo.param=file.path(primer3_path, "bin", "primer3_config\\"),
        settings=file.path(primer3_path,"primer3_v1_1_4_Ido_Bar_settings.txt")){

  sequence_exclude_region <- ifelse(SSR_exclude_region,sprintf("SEQUENCE_INTERNAL_EXCLUDED_REGION=%s\n", sequence_target), NULL)
  #print(excluded.regions)
  # make primer 3 input file
  p3.input=tempfile()
  p3.output=tempfile()
  write(
    paste0( sprintf("SEQUENCE_ID=%s\n",name  ),
            sprintf("SEQUENCE_TEMPLATE=%s\n",as.character(seq)),
            sprintf("SEQUENCE_TARGET=%s\n",sequence_target),
            sequence_exclude_region,
            sprintf("PRIMER_LIBERAL_BASE=%d\n", as.numeric(liberal_bases)),
            sprintf("PRIMER_NUM_RETURN=%d\n", primer_num),
            "PRIMER_TASK=pick_detection_primers\n",
            "PRIMER_PICK_LEFT_PRIMER=1\n" ,
            "PRIMER_PICK_INTERNAL_OLIGO=0\n",
            "PRIMER_PICK_RIGHT_PRIMER=1\n"  ,
            "PRIMER_EXPLAIN_FLAG=1\n"  ,
            "PRIMER_PAIR_MAX_DIFF_TM=3\n",
            sprintf("PRIMER_MIN_GC=%s\n" ,GC[1]),
            sprintf("PRIMER_MAX_GC=%s\n" ,GC[2]),
            sprintf("PRIMER_MIN_TM=%s\n" ,Tm[1]),
            sprintf("PRIMER_OPT_TM=%s\n" ,Tm[2]),
            sprintf("PRIMER_MAX_TM=%s\n" ,Tm[3]),
            sprintf("PRIMER_PRODUCT_SIZE_RANGE=%s\n" ,size_range),
            sprintf("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s\n" ,thermo.param),
            "=")
    ,  p3.input
  )
  #call primer 3 and store the output in a temporary file

#   try(system2(primer3 ,args=c(sprintf("-p3_settings_file=%s",settings), p3.input,
#                               ">", p3.output))
#   )

  try(system2(primer3 ,args=c(sprintf("-p3_settings_file=%s -output=%s",settings, p3.output), p3.input)))

  # sprintf("-p3_settings_file=%s", settings),
  #import and parse the output into a dataframe named designed.primers
  out=read.delim(p3.output, sep='=', header=FALSE)


  returned.primers=as.numeric(as.vector(out[out[,1]=='PRIMER_PAIR_NUM_RETURNED',][,2]))
  if (length(returned.primers)==0){ warning('primers not detected for ',name,call. = FALSE);return(NA)}
  if ((returned.primers)==0){ warning('primers not detected for ',name,call. = FALSE);return(NA)}
  if (returned.primers>0){
    designed.primers=data.frame()
    for (i in 0:(returned.primers-1)){
      #IMPORT SEQUENCES
      id=sprintf(  'PRIMER_LEFT_%i_SEQUENCE',i)
      PRIMER_LEFT_SEQUENCE=as.character(out[out[,1]==id,][,2])
      id=sprintf(  'PRIMER_RIGHT_%i_SEQUENCE',i)
      PRIMER_RIGHT_SEQUENCE=as.character(out[out[,1]==id,][,2])

      #IMPORT PRIMING POSITIONS
      id=sprintf(  'PRIMER_LEFT_%i',i)
      PRIMER_LEFT=as.numeric(unlist(strsplit(as.vector((out[out[,1]==id,][,2])),',')))
      #PRIMER_LEFT_LEN=as.numeric(unlist(strsplit(as.vector((out[out[,1]==id,][,2])),',')))
      id=sprintf(  'PRIMER_RIGHT_%i',i)
      PRIMER_RIGHT=as.numeric(unlist(strsplit(as.vector((out[out[,1]==id,][,2])),',')))
      #IMPORT Tm
      id=sprintf(  'PRIMER_LEFT_%i_TM',i)
      PRIMER_LEFT_TM=as.numeric(as.vector((out[out[,1]==id,][,2])),',')
      id=sprintf(  'PRIMER_RIGHT_%i_TM',i)
      PRIMER_RIGHT_TM=as.numeric(as.vector((out[out[,1]==id,][,2])),',')

      res=out[grep(i,out[,1]),]
      extra.inf=t(res)[2,,drop=FALSE]
      colnames(extra.inf)=sub( paste("_",i,sep=""),"",res[,1])
      extra.inf=extra.inf[,-c(4:9),drop=FALSE] # remove redundant columns
      extra.inf=apply(extra.inf,2,as.numeric)
      #Aggegate in a dataframe
      primer.info=data.frame(Seq_ID=name,i=i+1,
                             PRIMER_LEFT_SEQUENCE,PRIMER_RIGHT_SEQUENCE,
                             PRIMER_LEFT_TM, PRIMER_RIGHT_TM,
                             PRIMER_LEFT_pos=PRIMER_LEFT[1],
                             PRIMER_LEFT_len=PRIMER_LEFT[2],
                             PRIMER_RIGHT_pos=PRIMER_RIGHT[1],
                             PRIMER_RIGHT_len=PRIMER_RIGHT[2],
                             t(data.frame(extra.inf))

      )
      rownames(primer.info)=NULL
      designed.primers=rbind(designed.primers, primer.info)
      #print(primer.info)
    }



    #colnames(designed.primers)=c('PrimerID',
    #                             'Fwseq','Rvseq',
    #                             'FwTm','RvTm',
    #                             'FwPos','Fwlen',
    #                             'RvPos','Rvlen',
    #                             'fragLen' )
    if (!is.null(report)) {
      try(system2(primer3 ,args=c(sprintf("-p3_settings_file=%s -format_output -output=%s",settings, report), p3.input)))
    }
  }
  unlink(c(p3.input,p3.output) ) # delete temp files
  return(designed.primers)
}


#' convert primer list to genomic intervals
#'
#' function returns a GRanges or data frame object from a list of primers designed
#' by \code{designPrimers} function after calculating genomic location of the amplicon
#' targeted by the primers.
#'
#' @param primers a list of primers returned by \code{designPrimers} function
#' @param as.data.frame logical indicating if a data frame should be returned
#'        instead of \code{GRanges} object.
#'
#' @examples
#'  data(bisPrimers)
#'  # remove data with no primers found
#'  bisPrimers=bisPrimers[!is.na(bisPrimers)]
#'  gr.pr=primers2ranges(bisPrimers) # convert primers to GRanges
#'
#' @seealso \code{\link{filterPrimers}}, \code{\link{designPrimers}}
#'
#' @export
#' @docType methods
primers2ranges<-function(primers,as.data.frame=FALSE){

  if(any(is.na(primers))){
    warning( "There are targets without primers\nfiltering those before conversion")
    primers=primers[ !is.na(primers) ]
  }
  df=do.call("rbind",primers) # get primers to a df
  locs=gsub("\\|\\.+","",rownames(df)) # get the coordinates from list ids
  temp=do.call("rbind",strsplit(locs,"_")) #

  start=as.numeric(temp[,2])
  chr=as.character(temp[,1])



  amp.start= start + as.numeric(df$PRIMER_LEFT_pos)
  amp.end  = start + as.numeric(df$PRIMER_RIGHT_pos)
  res=data.frame(chr=chr,start=amp.start,end=amp.end,df)
  #saveRDS(res,file="/work2/gschub/altuna/projects/DMR_alignments/all.designed.primers.to.amps.rds")
  if(as.data.frame)
  {
    return(res)
  }
  gr=GRanges(seqnames=res[,1],ranges=IRanges(res[,2],res[,3]) )
  values(gr)=res[,-c(1,2,3)]
  gr
}

