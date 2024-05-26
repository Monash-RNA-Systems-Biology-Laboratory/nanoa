
rev_comp <- function(str) as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(str)))

rep_char <- function(c,n) c |> rep(n) |> paste(collapse="")

read_fastq <- function(filename) {
    reads <- ShortRead::readFastq(filename)
    tibble(
        name=as.character(reads@id), 
        seq=as.character(reads@sread), 
        qual=as.character(reads@quality@quality))
}

score_cigar <- function(cigar, score_match=1, score_mismatch=-2, gap_open=-1, gap_extend=-2) { 
    #Actually FASTA's modified cigar
    score <- 0
    i <- 0
    for(char in strsplit(cigar,"")[[1]]) {
        if (char == "M") {
            score <- score + i*score_match
            i <- 0
        } else if (char == "X") {
            score <- score + i*score_mismatch
            i <- 0
        } else if (char == "I" || char == "D") {
            score <- score + gap_open + i*gap_extend
            i <- 0
        } else {
            i <- i*10 + as.numeric(char)
        }
    }

    score
}

# Returns one hit for each query x subject
# Defaults to local-local alignment
# Can also use program="glsearch36" for global-local alignment
run_fasta <- function(query, subject, program="ssearch36") {
    int_query <- is.null(names(query))
    if (int_query)
        names(query) <- seq_along(query)
    
    int_subject <- is.null(names(subject))
    if (int_subject)
        names(subject) <- seq_along(subject)

    query_file <- tempfile()
    subject_file <- tempfile()
    hits_file <- tempfile()

    ShortRead::writeFasta(Biostrings::DNAStringSet(query), query_file)
    #writeFasta(DNAStringSet(subject), subject_file)

    n <- length(subject)

    # -n nucleotide
    # -3 forward strand only (I hope)
    # -r match/mismatch
    # -f gap open
    # -g gap extend
    ## -Z search as if there are this many subject sequences
    ## -b output results from this many subject sequences
    # -E expectation value cutoff
    ## -z -1 seems to make it output 1 hit per subject
    # -m output format
    # number ktup
    
    chunk_size <- 10000
    chunks <- split(seq_len(n), ceiling(seq_len(n)/chunk_size))
    hits <- map_df(chunks, function(indexes) {
        ShortRead::writeFasta(Biostrings::DNAStringSet(subject[indexes]), subject_file)
        system(paste0(
            program," ",
            #"-n -3 -r +1/-2 -f 0 -g -2 -E ",E_per_subject*n," -m 8CD ",
            #For  more consistent scores "-k 1000000 ",
            "-n -3 -r +1/-2 -f -1 -g -2 -z -1 -b =",length(indexes)," -m 8CD ",
            query_file, " ", subject_file, " >", hits_file))
        
        readr::read_tsv(hits_file, comment="#", col_names=FALSE, col_types="ccniiiiiiinnc")
    })

    unlink(query_file)
    unlink(subject_file)
    unlink(hits_file)

    colnames(hits) <- c("query","subject","pct_identity","length","mismatches","gaps","qstart","qend","sstart","send","evalue","bitscore", "cigar")
    
    hits$score <- purrr::map_dbl(hits$cigar, score_cigar)

    if (int_query) hits$query <- as.integer(hits$query)
    if (int_subject) hits$subject <- as.integer(hits$subject)

    hits
}


view_read <- function(seq, qual, hits) {
    nmax <- nchar(seq)
    step <- 90
    i <- 1
    while(i <= nmax) {
        j <- i+step-1
        cat(paste0(substring(qual,i,j),"\n"))
        cat(paste0(substring(seq,i,j),"\n"))

        for(k in seq_len(nrow(hits)))
        if (hits$sstart[k] <= j && hits$send[k] >= i) {
            a <- max(i,hits$sstart[k])
            b <- min(j,hits$send[k])
            cat(paste0(
                rep_char(" ",a-i),
                rep_char("=",b-a+1),
                rep_char(" ",j-b),
                " ",hits$query[k],
                " ",hits$score[k],
                " ",hits$cigar[k],
                "\n"
            ))
        }

        cat("\n")
        i <- j+1
    }
}


view_reads <- function(reads, hits) {
    for(i in seq_len(nrow(reads))) {
        cat(paste0("\n\n\n\nRead: ",reads$name[i],"\n\nLength: ",nchar(reads$seq[i]),"\n\n\n\n"))
        view_read(reads$seq[i], reads$qual[i], filter(hits, subject == i))
    }
}


second_best <- function(vec, minimum=0) {
    if (length(vec) >= 2)
        sort(vec,decreasing=TRUE)[2]
    else
        minimum
}


# For better_by, note each mismatch decreases score by 3.
demux <- function(reads, samples, starts, ends, min_score_start=10, min_score_end=10, better_by=6, min_length=20) {
    names(starts) <- samples
    names(ends) <- samples
    
    reads_both_strands <- tibble(
        name = c(reads$name, reads$name),
        reverse = rep(c(FALSE, TRUE), each=nrow(reads)),
        seq = c(reads$seq, rev_comp(reads$seq)),
        quality = c(reads$qual, stringi::stri_reverse(reads$qual)))
    
    hits_start <- run_fasta(starts, reads_both_strands$seq) |>
        transmute(
            read_row=subject, 
            start_start = sstart, 
            start_end = send, 
            start_cigar=cigar, 
            start_score=score, 
            sample=query)
    
    hits_end <- run_fasta(ends, reads_both_strands$seq) |>
        transmute(
            read_row=subject, 
            end_start = sstart, 
            end_end = send, 
            end_cigar=cigar, 
            end_score=score,
            sample=query)
    
    result <- inner_join(hits_start, hits_end, by=c("read_row","sample"), relationship="many-to-many") |>
        filter(start_end < end_start) |>
        mutate(score = start_score + end_score) |>
        mutate(
            name=reads_both_strands$name[read_row], 
            reverse=reads_both_strands$reverse[read_row], 
            seq=reads_both_strands$seq[read_row],
            quality=reads_both_strands$quality[read_row]) |>
        group_by(name) |>
        filter(score >= second_best(score)+better_by) |>
        filter(n() == 1) |>
        ungroup() |>
        mutate(start=start_end+1, end=end_start-1, clipped_length=end-start+1) |>
        filter(
            clipped_length >= min_length,
            start_score >= min_score_start,
            end_score >= min_score_end)
    
    result
}

mkdirs <- function(fp) {
    if(!file.exists(fp)) {
        mkdirs(dirname(fp))
        dir.create(fp)
    }
}

write_out <- function(out_dir, hits, samples, append=FALSE) {
    mkdirs(out_dir)
    
    csv_filename = paste0(out_dir,"/","reads.csv")
    if (!append) {
        con <- file(csv_filename, "w")
        writeLines("sample,reverse,length,clipped_length,name", con=con)
    } else {
        con <- file(csv_filename, "a")
    }
    writeLines(paste(
        hits$sample,
        hits$reverse,
        nchar(hits$seq),
        hits$clipped_length,
        hits$name,
        sep=","), con=con)
    close(con)
    
    for(this_sample in samples) {
        this_hits <- filter(hits, sample == this_sample)
        
        if (this_sample == "") 
            this_sample <- "reads"
        filename <- paste0(out_dir,"/",this_sample,".fastq.gz")
    
        if (!append) 
            writeLines(character(0), filename)
        
        from <- this_hits$start
        to <- this_hits$end
        seqs <- substring(this_hits$seq, from, to) |> Biostrings::DNAStringSet()
        quals <- substring(this_hits$quality, from, to) |> Biostrings::BStringSet()
        
        ShortRead::ShortReadQ(
                sread=seqs,
                quality=quals,
                id=Biostrings::BStringSet(this_hits$name)) |>
            ShortRead::writeFastq(
                filename, 
                mode="a")
    }
}


#' Demultiplex nanopore reads
#'
#' Separate a set of nanopore reads into separate samples, trimmed of adaptor sequence and correctly stranded.
#'
#' For alignment scores: Each matching base scores 1. Each mismatch scores -2. Each indel scores -1-2*length_of_indel.
#'
#' @param samples A data frame with three columns. The first column is the sample name (will be used in the output filename). The second column is the left adaptor sequence. The third column is the right adaptor sequence.
#' @param filenames A character vector of fastq read filenames.
#' @param out_dir Directory to put output files in.
#' @param min_score_left Minimum alignment score for the left primer.
#' @param min_score_right Minimum alignment score for the left primer.
#' @param better_by Total score of the best sample must be this much better than the runner-up. Note each mismatch reduces the score by three, so the default of 6 means the best sample has to have two less mismatches than the runner-up.
#' @param min_length Only output reads containing this many bases after adaptor trimming.
#'
#' @export
demultiplex <- function(
        samples, filenames, out_dir, 
        min_score_left=10, min_score_right=10, better_by=6, min_length=20) {
    samples <- tibble::as_tibble(samples)
    assertthat::assert_that(ncol(samples) == 3)
    sample_names <- samples[[1]]
    starts <- samples[[2]]
    ends <- rev_comp( samples[[3]] )
    
    first <- TRUE
    for(i in seq_along(filenames)) {
        filename <- filenames[i]
        message("Processing file ",i," of ",length(filenames)," ",filename)
        
        reads <- read_fastq(filename)
        hits <- demux(
            reads, sample_names, starts, ends,
            min_score_start=min_score_left, min_score_end=min_score_right, 
            better_by=better_by, min_length=min_length)
        
        write_out(out_dir, hits, sample_names, append=!first)
        first <- FALSE
    }
}

