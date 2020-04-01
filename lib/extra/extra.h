//
// Created by borja on 31/3/20.
//

#ifndef VIADBG_TRIAL_EXTRA_H
#define VIADBG_TRIAL_EXTRA_H

/*
 * Parameters
 */
struct Parameters
{
    static Parameters& get(){
        static Parameters param;
        return param;
    }
    void set(size_t param, std::string param_read)
    {
        if (param_read == "KmerSize")
            Parameters::get().kmerSize = param;
        if (param_read == "Accumulative")
            Parameters::get().accumulative_h = param;
        if (param_read == "missmatches")
            Parameters::get().missmatches = param;
        if (param_read == "numKmers")
            Parameters::get().num_unique_kmers = param;
        if (param_read == "genome_size")
            Parameters::get().genome_size = param;
        if (param_read == "full_info")
            Parameters::get().full_info = param;
        if (param_read == "metagenomic")
            Parameters::get().metagenomic = param;
        if (param_read == "remove_duplicates")
            Parameters::get().remove_duplicates = param;
        if (param_read == "numThreads")
            Parameters::get().numThreads = param;
        if (param_read == "polish")
            Parameters::get().polish = param;
        if (param_read == "postProcess")
            Parameters::get().postProcess = param;
        if (param_read == "gfa")
            Parameters::get().gfa = param;
    }
    double missmatches;
    size_t genome_size = 0;
    double num_unique_kmers = 0;
    size_t accumulative_h = 0;
    bool full_info = false, metagenomic = false,
            postProcess = false, remove_duplicates = true, polish = true, gfa = false;
    size_t kmerSize;
    size_t numThreads;
    bool show = false;
};

#endif //VIADBG_TRIAL_EXTRA_H
